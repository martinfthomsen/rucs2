# -*- coding: utf-8 -*-
"""
    PCR library
    ~~~~~~~~~~~~~~~~~~~~~

    Contain PCR handling functionalities.

    * AdvancedDictionary (dict-based class) - This class expands on the dictionary class by making it hierachical
    * multirange (class) - Class for keeping track of marked ranges. This makes it possible to keep track of which positions have been marked without creating a huge list for each position.
    * find_validated_primer_pairs - Design and identify primer pairs for the provided template, which are validated against the possitive and negative references.
    * pair_sort - Method for deciding sorting rank of primers.
    * compute_binding_sites - Compute binding sites by calculating heterodimer thermodynamics for all alignments and filtering those with a melting temperature below 0 degrees
    * grade_binding_sites - Grade the binding sites for each primer sequence based on the thermodynamics according to the predefined threshold criteria
    * estimate_primer_rank - Calculate the penalty for the alignments given the scheme (positive or negative)
    * validate_primer_pairs - Validate and score primer pairs
    * find_pcr_products - Find primer pair products
    * present_pairs_full - Presents essential aspects of the primer pairs
    * predict_pcr_results - Predict the results of a PCR
    * print_pcr_results - Print PCR results
    * get_primer_and_probe_bindsites - Extract the sequences and the primer and probe location found in the
    sequences.
    * print_primer_and_probe_bindsites - print primer and probe binding sites
    * color_seq - fct for applying coloring to print statements from print_primer_and_probe_bindsites

"""
# (c) 2023 Martin Thomsen

import sys, os, re, types, bisect, pickle
import numpy as np

from rucs import settings, num2binarray, binarray2num, round_sig, get_directories
from rucs.log import log
from rucs.file import seqs_from_file, save_as_fasta, parse_tsv, open_
from rucs.seq import reverse_complement, align_seqs
from rucs.primer3 import compute_primer_pairs, configure_p3_thermoanalysis
from rucs.blast import blast_to_ref, get_blast_annotations


# CLASSES
class AdvancedDictionary(dict):
    """ This class expands on the dictionary class by making it hierachical

    Hierachical Dictionary Methods:
        get_tree - Extract nested key combinations
        sort     - Go through the key value pairs, and sort the value arrays
        update   - Add k,v pairs where identical keys just adds multiple values
                   to the value list
        invert   - returns the invert of the dictionary, so that all values
                   become keys and vise versa

    USAGE
        >>> import types
        >>> d = AdvancedDictionary([('A',1),('B',5),('B',2)], sort=True)
        >>> d
        {'B': [2, 5], 'A': [1]}
        >>> d.update([('C', 1)])
        >>> d.invert()
        {1: ['A', 'C'], 2: ['B'], 5: ['B']}
        >>> d = AdvancedDictionary({'A':{'B':5}})
        >>> d.get_tree(['A','B'])
        5
    """
    def __init__(self, data, sort=False, key=None, reverse=False):
        if isinstance(data, (list, types.GeneratorType)):
            self.update(data)
        elif isinstance(data, dict):
            for k in data:
                self[k] = data[k]
        else:
            raise ValueError('Error: Non supported input format! '
                             'AdvancedDictionary only supports three types of '
                             'input; <dict>s, <list>s and <generator>s.')
        if sort: self.sort(key=key, reverse=reverse)
    def update(self, list_of_kv_pairs):
        ''' Add list of key-value pairs to the hierachical dictionary

        USAGE
            >>> mydict = AdvancedDictionary([('A',1),('B',5),('B',2)])
            >>> mydict.update([('C',3),('C',1),('A',2),('A',4),('A',3),('B',3)])
            >>> mydict
            {'A': [1, 2, 4, 3], 'C': [3, 1], 'B': [5, 2, 3]}
        '''
        for k, v in list_of_kv_pairs:
            if not k in self:
                self[k] = [v]
            elif isinstance(self[k], (str, int, float)):
                self[k] = [self[k], v]
            elif not isinstance(self[k], list):
                raise ValueError('Error: Non supported value format! Values may '
                                 'only be numerical, strings, or lists of numbers '
                                 'and strings.')
            else:
                self[k].append(v)
    def get_tree(self, list_of_keys):
        """ gettree will extract the value from a nested tree

        INPUT
            list_of_keys: a list of keys ie. ['key1', 'key2']
        USAGE
            >>> # Access the value for key2 within the nested dictionary
            >>> ad = AdvancedDictionary({'key1': {'key2': 'value'}})
            >>> ad.gettree(['key1', 'key2'])
            'value'
        """
        cur_obj = self
        for key in list_of_keys:
            cur_obj = cur_obj.get(key)
            if not cur_obj: break
        return cur_obj
    def invert(self):
        ''' Return inverse mapping of dictionary with sorted values.

        USAGE
            >>> # Switch the keys and values
            >>> AdvancedDictionary({
            ...     'A': [1, 2, 3],
            ...     'B': [4, 2],
            ...     'C': [1, 4],
            ... }).invert()
            {1: ['A', 'C'], 2: ['A', 'B'], 3: ['A'], 4: ['B', 'C']}
        '''
        inv_map = {}
        for k, v in self.items():
            if isinstance(v, (str, int, float)): v = [v]
            elif not isinstance(v, list):
                raise ValueError('Error: Non supported value format! Values may '
                                 'only be numerical, strings, or lists of numbers '
                                 'and strings.')
            for val in v:
                inv_map[val] = inv_map.get(val, [])
                inv_map[val].append(k)
                inv_map[val].sort()
        return inv_map
    def sort(self, key=None, reverse=False):
        ''' Sort the value-arrays in place using the built-in method "sorted".

        USAGE
            >>> a = AdvancedDictionary({
            ...     'A': [1, 3, 2],
            ...     'B': [4, 2],
            ...     'C': [1, 4],
            ... })
            >>> a.sort()
            >>> a
            {'B': [2, 4], 'A': [1, 2, 3], 'C': [1, 4]}
        '''
        for k in self.keys():
            self[k] = sorted(self[k], key=key, reverse=reverse)


class multirange():
    ''' Class for keeping track of marked ranges

    This makes it possible to keep track of which positions have been marked
    without creating a huge list for each position.

    AUTHOR: Martin CF Thomsen

    USAGE:
        >>> import bisect
        >>> # Create a multirange object and print or iterate over it
        >>> marked_positions = multirange()
        >>> marked_positions.add(1,2)
        >>> marked_positions.add(3,4)
        >>> marked_positions.add(7,8)
        >>> print(marked_positions)
        [(1, 4), (7, 8)]
        >>> for start, end in marked_positions:
        ...     print(start, end)
        ...
        1 4
        7 8
        >>> # Add two or more multiranges together, and see where the ranges overlap etc.
        >>> first = multirange([(0,11),(21,28),(41,65)])
        >>> second = multirange([(6,23),(35,42)])
        >>> third = multirange([(1,2),(68,70)])
        >>> print(first + second + third)
        [(0, 0, 4), (1, 2, 5), (3, 5, 4), (6, 11, 6), (12, 20, 2), (21, 23, 6),
         (24, 28, 4), (35, 40, 2), (41, 42, 6), (43, 65, 4), (68, 70, 1)]
    TEST:
        >>> import bisect
        >>> multirange().__test__()
        95 out of 95 successful tests!
    '''
    def __init__(self, ranges=None):
        self.index = 0
        self.tags = None
        try:
            self.starts, self.ends, self.tags = map(list,zip(*ranges))
        except:
            try:
                self.starts, self.ends = map(list,zip(*ranges))
            except:
                if ranges is None:
                    self.starts, self.ends = [], []
                else:
                    print(ranges)
                    raise
    def __iter__(self):
        for i in range(len(self.starts)):
            if self.tags is not None:
                yield self.starts[i], self.ends[i], self.tags[i]
            else:
                yield self.starts[i], self.ends[i]
    def __next__(self):
        self.index += 1
        if self.index <= len(self.starts):
            if self.tags is not None:
                return self.starts[self.index-1], self.ends[self.index-1], self.tags[self.index-1]
            else:
                return self.starts[self.index-1], self.ends[self.index-1]
        else:
            raise StopIteration
    def __reset__(self):
        self.index = 0
    def __eq__(self, other):
        return self._get() == other
    def __str__(self):
        return str(self._get())
    def _get(self):
        if self.tags is not None:
            return list(zip(self.starts, self.ends, self.tags))
        else:
            return list(zip(self.starts, self.ends))
    def __len__(self):
        return sum((e-s+1 for s, e in self))
    def __sub__(self, other):
        ''' Return ranges which are not covered by other '''
        return [(s, e) for s, e, t in self.__add__(other) if t == 1]
    def add(self, start, end):
        ''' add a new start and end position to the stored ranges '''
        # Identify placement of new range
        idx_start = bisect.bisect(self.starts, start) -1
        idx_end   = bisect.bisect_left(self.ends, end)
        # Nose dipping?
        if not idx_start < 0 and start-1 <= self.ends[idx_start]:
            if start > self.starts[idx_start]:
               start = self.starts[idx_start]
            idx_start -= 1

        # Tail dipping?
        if not idx_end >= len(self.ends) and end+1 >= self.starts[idx_end]:
            if end < self.ends[idx_end]:
               end = self.ends[idx_end]
            idx_end += 1

        # Remove covered ranges
        del self.starts[idx_start+1:idx_end], self.ends[idx_start+1:idx_end]
        # Insert new range
        self.starts.insert(idx_start+1, start)
        self.ends.insert(idx_start+1, end)
    def __add__(self, other):
        ''' add two multiranges together

        The output is a multirange object with tags in the third column.
        The tag indicate which multirange object overlap the range:

        EXAMPLE
            >>> print(multirange([(1,2),(6,8)]) +     # 3. bit in tag
                      multirange([(5,7),(10,12)]) +   # 2. bit in tag
                      multirange([(12,12)]))          # 1. bit in tag
            [(1, 2, 4), (5, 5, 2), (6, 7, 6), (8, 8, 4), (10, 11, 2), (12, 12, 3)]
        '''
        segments = [] # [(start, end, tag), ...]
        try:
            tmp = next(self)
            s_start, s_end = tmp[:2]
            if len(tmp) == 3:
                # Binary Left Shift Operation of previous tag
                tag_size = tmp[2] << 1
            else:
                tag_size = 2
        except: s_start = None

        try: o_start, o_end = next(other)
        except: o_start = None

        incision_pos = 0
        while s_start is not None or o_start is not None:
            # get next ranges if previous is passed
            # Find next incision point
            tag = 0
            if s_start is not None:
                if incision_pos > s_end:
                    try:
                        tmp = next(self)
                        s_start, s_end = tmp[:2]
                        if len(tmp) == 3:
                            # Binary Left Shift Operation of previous tag
                            tag_size = tmp[2] << 1
                        else:
                            tag_size = 2
                    except: s_start = None
                if s_start is not None:
                    if s_start > incision_pos:
                        next_incision_dist = s_start - incision_pos -1
                    elif s_end >= incision_pos:
                        next_incision_dist = s_end - incision_pos
                        tag += tag_size
                elif o_start is not None:
                    next_incision_dist = o_end
            else:
                next_incision_dist = o_end

            if o_start is not None:
                if incision_pos > o_end:
                    try: o_start, o_end = next(other)
                    except: o_start = None
                if o_start is not None:
                    if o_start > incision_pos:
                        dist = o_start - incision_pos -1
                        if dist < next_incision_dist:
                            next_incision_dist = dist
                    elif o_end >= incision_pos:
                        dist = o_end - incision_pos
                        if dist < next_incision_dist:
                            next_incision_dist = dist
                        tag += 1

            # Store segment
            if tag > 0:
                next_incision_pos = incision_pos + next_incision_dist
                segments.append((incision_pos, next_incision_pos, tag))
            # Set next incision point
            incision_pos += next_incision_dist +1

        return multirange(segments)
    def __test__(self):
        '''  Tests '''
        test_start_values = self.starts, self.ends, self.tags
        # Set answers
        a = {}
        a[(0,0)] = [(0, 0),(2,2),(4,4),(8,10)]
        a[(0,1)] = a[(0,2)] = [(0,2),(4,4),(8,10)]
        a[(0,3)] = a[(0,4)] = [(0,4),(8,10)]
        a[(0,5)] = [(0,5),(8,10)]
        a[(0,6)] = [(0,6),(8,10)]
        a[(0,7)] = a[(0,10)] = a[(0,8)] = a[(0,9)] = [(0,10)]
        a[(0,11)] = [(0,11)]
        a[(0,12)] = [(0,12)]
        a[(0,13)] = [(0,13)]
        a[(1,1)] = a[(1,2)] = [(1,2),(4,4),(8,10)]
        a[(1,3)] = a[(1,4)] = [(1,4),(8,10)]
        a[(1,5)] = [(1,5),(8,10)]
        a[(1,6)] = [(1,6),(8,10)]
        a[(1,7)] = a[(1,8)] = a[(1,9)] = a[(1,10)] = [(1,10)]
        a[(1,11)] = [(1,11)]
        a[(1,12)] = [(1,12)]
        a[(2,2)] = a[(4,4)] = [(2,2),(4,4),(8,10)]
        a[(2,3)] = a[(2,4)] = a[(3,3)] = a[(3,4)] = [(2,4),(8,10)]
        a[(2,5)] = a[(3,5)] = [(2,5),(8,10)]
        a[(2,6)] = a[(3,6)] = [(2,6),(8,10)]
        a[(2,7)] = a[(2,8)] = a[(2,9)] = a[(2,10)] = [(2,10)]
        a[(2,11)] = a[(3,11)] = [(2,11)]
        a[(2,12)] = a[(3,12)] = [(2,12)]
        a[(3,7)] = a[(3,8)] = a[(3,9)] = a[(3,10)] = [(2,10)]
        a[(4,5)] = a[(5,5)] = [(2,2),(4,5),(8,10)]
        a[(4,6)] = a[(5,6)] = [(2,2),(4,6),(8,10)]
        a[(4,7)] = a[(4,8)] = a[(4,9)] = a[(4,10)] = [(2,2),(4,10)]
        a[(4,11)] = a[(5,11)] = [(2,2),(4,11)]
        a[(4,12)] = a[(5,12)] = [(2,2),(4,12)]
        a[(5,7)] = a[(5,8)] = a[(5,9)] = a[(5,10)] = [(2,2),(4,10)]
        a[(6,6)] = [(2, 2), (4, 4), (6, 6), (8, 10)]
        a[(6,7)] = a[(6,8)] = a[(6,9)] = a[(6,10)] = [(2,2),(4,4),(6,10)]
        a[(6,11)] = [(2,2),(4,4),(6,11)]
        a[(6,12)] = [(2,2),(4,4),(6,12)]
        a[(7,7)] = a[(7,8)] = a[(7,9)] = a[(7,10)] = [(2,2),(4,4),(7,10)]
        a[(7,11)] = [(2,2),(4,4),(7,11)]
        a[(7,12)] = [(2,2),(4,4),(7,12)]
        a[(8,8)] = a[(8,9)] = a[(8,10)] = [(2,2),(4,4),(8,10)]
        a[(9,9)] = a[(9,10)] = a[(10,10)] = [(2,2),(4,4),(8,10)]
        a[(8,11)] = a[(9,11)] = a[(10,11)] = a[(11,11)] = [(2,2),(4,4),(8,11)]
        a[(8,12)] = a[(9,12)] = a[(10,12)] = a[(11,12)] = [(2,2),(4,4),(8,12)]
        a[(12,12)] = [(2, 2), (4, 4), (8, 10), (12, 12)]
        # Set starts and ends for each test
        starts, ends = [2, 4, 8], [2, 4, 10]
        # Test all
        successes = 0
        failures = 0
        for start, end in ((i,j+1) for i in range(0,12) for j in range(i-1,12)):
            test = multirange(list(zip(starts, ends)))
            test.add(start, end)
            if test == a.get((start, end)):
                successes += 1
            else:
                failures += 1
                print('%s %s -> Failure: %s --> %s'%(start, end, test, a.get((start, end))))

        # Single tests
        test = multirange([(1,2),(6,7),(8,9)])
        start, end = 4, 4
        test.add(start, end)
        answer = [(1,2),(4,4),(6,7),(8,9)]
        if test == answer:
            successes += 1
        else:
            failures += 1
            print('%s %s -> Failure: %s --> %s'%(start, end, test, answer))

        test = multirange([(1,2),(6,7),(8,9)])
        start, end = 4, 4
        test.add(start, end)
        answer = [(1,2),(4,4),(6,7),(8,9)]
        if test == answer:
            successes += 1
        else:
            failures += 1
            print('%s %s -> Failure: %s --> %s'%(start, end, test, answer))

        test = multirange([(5256,5276),(7541,7561),(7863,7883)])
        start, end = 5883, 5903
        test.add(start, end)
        answer = [(5256,5276),(5883,5903),(7541,7561),(7863,7883)]
        if test == answer:
            successes += 1
        else:
            failures += 1
            print('%s %s -> Failure: %s --> %s'%(start, end, test, answer))

        test = multirange([(5255,5280),(5883,5903),(7221,7246),(7541,7562),
                           (7856,7884),(7906,7935)])
        start, end = 5877, 5898
        test.add(start, end)
        answer = [(5255,5280),(5877,5903),(7221,7246),(7541,7562),(7856,7884),
                  (7906,7935)]
        if test == answer:
           successes += 1
        else:
            failures += 1
            print('%s %s -> Failure: %s --> %s'%(start, end, test, answer))

        r1, r2, r3 = (multirange([(1,2),(6,8)]),
                      multirange([(5,7),(10,12)]),
                      multirange([(12,12)]))

        test = r1 + r2 + r3
        answer = [(1, 2, 4), (5, 5, 2), (6, 7, 6), (8, 8, 4), (10, 11, 2),
                  (12, 12, 3)]
        if test == answer:
            successes += 1
        else:
            failures += 1
            print('%s + %s + %s -> Failure: %s --> %s'%(r1, r2, r3, test, answer))

        print("%s out of %s successful tests!"%(successes, successes+failures))
        # Return to original values
        self.starts, self.ends, self.tags = test_start_values


def find_validated_primer_pairs(contig_file, p_refs, n_refs,
                                contig_names=None, annotate=True,
                                skip_start=0, max_contigs=0,
                                max_primer_pairs=10000, reuse=False):
    ''' if you do not have scaffolds, the contigs can be passed here

    USAGE
    >>> contig_file = 'contigs_tem.fa'
    >>> p_refs = ['genome_pos_1.fna']
    >>> n_refs = ['genome_neg_1.fna']
    >>> bam_name = 'myseq'
    >>> bwa_settings = None
    >>> threshold_tm = 45
    >>> pairs = find_validated_primer_pairs(
    ...    contig_file, p_refs, n_refs, bam_name,
    ...    bwa_settings, threshold_tm)
    >>> print(present_pairs(pairs[:3]))
    # product_size	tm	forward primer	reverse primer	test score
    515	60.0	CAACATTTTCGTGTCGCCCTT	TCGTCGTTTGGTATGGCTTCA	0
    636	60.0	CAACATTTTCGTGTCGCCCTT	TGCAACTTTATCCGCCTCCAT	0
    339	60.0	TGAAGCCATACCAAACGACGA	GAGGCACCTATCTCAGCGATC	0
    '''
    # Get directories
    global work_dir, ref_dir, result_dir
    work_dir, ref_dir, result_dir = get_directories(get_ref_dir=True)

    p3_args = settings['pcr']['priming']['primer3']
    if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'], list):
        if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0], list):
            product_length_limits = [lim for lims in p3_args['PRIMER_PRODUCT_SIZE_RANGE'] for lim in lims]
            min_seq_len = min(product_length_limits)
        else:
            min_seq_len = p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0]
    else:
        exit('Settings error: invalid input format for PRIMER_PRODUCT_SIZE_RANGE, expected "list", got: ' + p3_args['PRIMER_PRODUCT_SIZE_RANGE'])

    qpcr = p3_args['PRIMER_PICK_INTERNAL_OLIGO'] == 1
    to_upper = settings['input']['to_upper']
    filter_grade = settings['pcr']['priming']['threshold_grade']
    tm_thresholds = settings['pcr']['priming']['tm_thresholds']
    buffer = settings['input']['use_ram_buffer']
    # Input validation
    assert (contig_names is None or isinstance(contig_names, list)), \
           ('Invalid value (contig_names)! Only list allowed.')

    if (contig_names is None and 'seq_selection' in settings['pcr']
        and settings['pcr']['seq_selection'] is not None):
        contig_names = list(map(str, settings['pcr']['seq_selection']))

    if log is not None:
        log.progress.add('pcr', 'Find valid primer pairs for PCR', 'main')
    else:
        print('Find valid primer pairs for PCR')
    # WHILE UCS and need more good PP
    plen = len(p_refs)
    nlen = len(n_refs)
    too_short = 0
    skipped = 0
    ignored = 0
    no_pairs = 0
    max_pp_reached = False
    primer_pairs = []
    good_pp = []
    if reuse:
        # Try fetching previously computed values
        try:
            with open(f'{work_dir}good_primer_pairs.pkl', 'rb') as f:
                reuse_i, too_short, skipped, ignored, no_pairs, good_pp = pickle.load(f)
            with open(f'{work_dir}primer_pairs.pkl', 'rb') as f:
                max_pp_reached, primer_pairs = pickle.load(f)
        except Exception as e:
            reuse_i = -1
            log.progress.add('reuse', f'Could not load primer pair dumps from previous run. {e}', 'pcr')
        else:
            log.progress.add('reuse', f'Reuse successful. Primer pair data from the first {reuse_i+1} contigs was loaded.', 'pcr')

    for i, (seq, name, desc) in enumerate(seqs_from_file(contig_file,
                                                         to_upper=to_upper,
                                                         use_ram_buffer=buffer)):
        if reuse and i <= reuse_i:
            # Skip previously processed contigs
            continue
        if i < skip_start:
            # Ignore beginning contig, that are requested omitted
            skipped += 1
            continue
        if max_contigs > 0 and i >= max_contigs:
            # Ignore contig, since the desired number of contigs have been analysed
            break
        if contig_names is not None and name not in contig_names and str(i) not in contig_names:
            # Ignore unspecified contigs
            ignored += 1
            continue
        seqlen = len(seq)
        if seqlen < min_seq_len:
            # Ignore sequences which are too short for Primer3 settings
            too_short += 1
            continue
        if seqlen - seq.count('n') < 50:
            # Ignore scarce sequences content
            too_short += 1
            continue
        if len(good_pp) >= settings['pcr']['no_good_pp']:
            # Ignore contig, since the desired number of good pairs is reached
            skipped += 1
            continue

        if log is not None:
            # Prepare logs for primer pairs statistics
            log.stats.add_table('prims%s'%i,
                                'P3 primer scan - %s'%name,
                                ['Right', 'Left', 'Probe', 'Category'])
            log.stats.add_table('pair%s'%i,
                                'P3 pair scan - %s'%name,
                                ['Pairs', 'Category'])
            log.progress.add('scan%s'%i,
                             f'Scanning contig #{i} {name} for primer pairs (size:{seqlen})', 'pcr')
        else:
            print(f' Scanning contig #{i} {name} for primer pairs (size:{seqlen})')

        # COMPUTE primer pairs for sequence
        pairs, notes = compute_primer_pairs(seq)

        if log is not None:
            # Log primer pairs statistics for sequence
            log.progress['scan%s'%i].log_time()
            # Process primer3 notes
            right = dict([x.rsplit(' ', 1) for x in notes['right'].split(', ')])
            left = dict([x.rsplit(' ', 1) for x in notes['left'].split(', ')])
            probe = dict([x.rsplit(' ', 1) for x in notes.get('internal', '').split(', ') if x != ''])
            pair = dict([x.rsplit(' ', 1) for x in notes['pair'].split(', ')])
            headers = ['considered', 'too many Ns', 'GC content failed',
                       'long poly-x seq', 'low tm', 'high tm',
                       'high any compl', 'high hairpin stability', 'ok']
            headers_pair = ['considered', 'unacceptable product size',
                            'high any compl', 'high end compl',
                            'tm diff too large', 'no internal oligo', 'ok']
            # Warn about headers unaccounted for
            if miss_head := ', '.join(h for h in left if not h in headers):
                print(('Warning: One or more forward primer headers were not '
                       'reported: %s')%miss_head)
            if miss_head := ', '.join(h for h in right if not h in headers):
                print(('Warning: One or more reverse primer headers were not '
                       'reported: %s')%miss_head)
            if miss_head := ', '.join(h for h in probe if not h in headers):
                print(('Warning: One or more probe headers were not '
                       'reported: %s')%miss_head)
            if miss_head := ', '.join(h for h in pair if not h in headers_pair):
                print(('Warning: One or more pair headers were not '
                       'reported: %s')%miss_head)
            # Log Primer3 notes as statistics
            for h in headers:
                log.stats.add_row('prims%s'%i, [right.get(h, 0), left.get(h, 0), probe.get(h, 0), h])

            for h in headers_pair:
                log.stats.add_row('pair%s'%i, [pair.get(h, 0), h])

        if len(pairs) == 0:
            # Skip contigs with no primer pair candidates
            no_pairs += 1
            continue

        # EXTRACT primer pair sequences in prep for Virtual PCR
        p_fw = []
        p_rv = []
        probes = []
        for p in pairs:
            p_fw.append(p['left']['sequence'])
            p_rv.append(p['right']['sequence'])
            if qpcr: probes.append(p['internal']['sequence'])
            # Set sequence name which the primer was based on
            p['sequence_id'] = name

        # Extract all primers from primer3 results
        primers = dict([(p, {'pos':[],'neg':[],'penalty':None})
                        for p in p_fw + list(map(reverse_complement, p_rv)) + probes])

        # Store primers as fastq in prep for alignment algorithm
        primers_fa = f'{work_dir}primers.fa'
        save_as_fasta(primers, primers_fa)

        if log is not None:
            # Prepare logs for vitual PCR based on alignment (positives)
            log.stats.add_table('val%s'%i,
                                'Primer Prediction Overview (%s) for Contig %s'%(
                                   len(primers), name),
                                ['Aligned','Grade 1', 'Grade 2', 'Grade 3',
                                 'Grade 4', 'Grade 5', 'Reference'])
            log.progress.add('pos%s'%i, 'Aligning primers to positive references',
                             'pcr')
        else:
            print(' Aligning primers to positive references')

        # Create primer-reference mapping matrix for positive references
        for j, ref in enumerate(p_refs):
            if log is not None:
                log.progress.add("%s_pos_%s"%(i, j),
                                 'Aligning %s to %s'%(
                                     os.path.basename(primers_fa),
                                     os.path.basename(ref)),
                                 'pos%s'%i)
            else:
                print('  Aligning %s to %s'%(os.path.basename(primers_fa),
                                             os.path.basename(ref)))

            # Align sequences to reference
            alignments = blast_to_ref(ref, primers_fa,
                                      settings['pcr']['priming']['blastn_settings'],
                                      buffer)
            prim_log = [sum((len(m) > 0 for m in alignments.values())),
                        0 if filter_grade <= 1 else '-',
                        0 if filter_grade <= 2 else '-',
                        0 if filter_grade <= 3 else '-',
                        0 if filter_grade <= 4 else '-',
                        0, os.path.basename(ref)]

            # Filter primer matches which fail the thermodynamic test
            binding_sites = compute_binding_sites(ref, alignments, probes, filter=True)
            # print(f"Stats: {sum(map(len, binding_sites.values()))} binding sites found for {ref}")
            binding_sites = grade_binding_sites(binding_sites, probes, tm_thresholds, filter_grade)
            # print(f"Stats: {sum(map(len, binding_sites.values()))} graded binding sites found for {ref}")

            # Update primer-reference matrix
            for p in primers:
                primers[p]['pos'].append(binding_sites[p] if p in binding_sites else [])

            if log is not None:
                log.progress["%s_pos_%s"%(i, j)].log_time()
                # Update prim_log with grade data
                for v in binding_sites.values():
                    max_g = 0
                    for a in v:
                        #(contig_name, strand, pos, tm, grade)
                        if a[-1] > max_g:
                            max_g = a[-1]

                    if max_g >= filter_grade:
                        prim_log[max_g] += 1

                log.stats.add_row('val%s'%i, prim_log)

        if log is not None:
            log.progress['pos%s'%i].log_time()

            # Prepare logs for vitual PCR based on alignment (negatives)
            log.progress.add('neg%s'%i, 'Aligning primers to negative references',
                             'pcr')
        else:
            print(' Aligning primers to negative references')

        # Create primer-reference mapping matrix for negative references
        for j, ref in enumerate(n_refs):
            if log is not None:
                log.progress.add("%s_neg_%s"%(i, j),
                                 'Aligning %s to %s'%(
                                     os.path.basename(primers_fa),
                                     os.path.basename(ref)),
                                 'neg%s'%i)
            else:
                print('  Aligning %s to %s'%(os.path.basename(primers_fa),
                                           os.path.basename(ref)))
            # Align sequences to reference
            alignments = blast_to_ref(ref, primers_fa,
                                      settings['pcr']['priming']['blastn_settings'],
                                      buffer)
            prim_log = [sum((len(m) > 0 for m in alignments.values())),
                        0 if filter_grade <= 1 else '-',
                        0 if filter_grade <= 2 else '-',
                        0 if filter_grade <= 3 else '-',
                        0 if filter_grade <= 4 else '-',
                        0, os.path.basename(ref)]

            # Filter primer matches which fail the thermodynamic test
            binding_sites = compute_binding_sites(ref, alignments, probes, filter=True)
            # print(f"Stats: {sum(map(len, binding_sites.values()))} binding sites found for {ref}")
            binding_sites = grade_binding_sites(binding_sites, probes, tm_thresholds, filter_grade)
            # print(f"Stats: {sum(map(len, binding_sites.values()))} graded binding sites found for {ref}")

            # Update primer-reference matrix
            for p in primers:
                primers[p]['neg'].append(binding_sites[p] if p in binding_sites else [])

            if log is not None:
                log.progress["%s_neg_%s"%(i, j)].log_time()
                # Update prim_log with grade data
                for v in binding_sites.values():
                    max_g = 0
                    for a in v:
                        #(contig_name, strand, pos, tm, grade)
                        if a[-1] > max_g:
                           max_g = a[-1]

                    if max_g >= filter_grade:
                        prim_log[max_g] += 1

                log.stats.add_row('val%s'%i, prim_log)

        if log is not None:
            log.progress['neg%s'%i].log_time()
            log.progress.add('rank%s'%i, 'Ranking primer pairs', 'pcr')
        else:
            print(' Ranking primer pairs')

        # Compute primer rank based on the number of hits to the negative references
        for p in primers:
            # Rank primers based on their performance
            # Compute the sum of the primer's penalties to each reference
            penalty_pos = sum([estimate_primer_rank(aln, scheme='positive')
                            for aln in primers[p]['pos']])
            penalty_neg = sum([estimate_primer_rank(aln, scheme='negative')
                            for aln in primers[p]['neg']])
            primers[p]['penalty'] = (penalty_pos, penalty_neg)

        # Compute primer pair penalty
        for p in pairs:
            p_fw = p['left']['sequence']
            p_rv = reverse_complement(p['right']['sequence'])
            # Store ranks for forward and reverse primer
            p['left']['pen_p'] = primers[p_fw]['penalty'][0]
            p['left']['pen_n'] = primers[p_fw]['penalty'][1]
            p['right']['pen_p'] = primers[p_rv]['penalty'][0]
            p['right']['pen_n'] = primers[p_rv]['penalty'][1]
            # Compute ranks for forward and reverse primer
            if not 'test' in p: p['test'] = {}
            penalty_pos = (primers[p_fw]['penalty'][0] + 1) * (primers[p_rv]['penalty'][0] + 1) - 1
            penalty_neg = (primers[p_fw]['penalty'][1] + 1) * (primers[p_rv]['penalty'][1] + 1) - 1
            if qpcr:
                probe = p['internal']['sequence']
                # Store ranks for probe
                p['internal']['pen_p'] = primers[probe]['penalty'][0]
                p['internal']['pen_n'] = primers[probe]['penalty'][1]
                # Compute ranks for probe
                penalty_pos += primers[probe]['penalty'][0]
                penalty_neg += primers[probe]['penalty'][1]
            # Store pair rank
            p['test']['penalty'] = penalty_pos + penalty_neg

        # Sort primer pairs according to their RUCS and Primer3 penalty scores
        if log is not None:
            log.progress['rank%s'%i].log_time()
            log.progress.add('sort%s'%i, 'Sorting primer pairs', 'pcr')
        else:
            print(' Sorting primer pairs')
        pairs = sorted(pairs, key=lambda p: (p['test']['penalty'],
                                             p['pair']['penalty']))

        # Validate and score primer pairs
        if log is not None:
            log.progress['sort%s'%i].log_time()

        good_pp.extend(validate_primer_pairs(pairs,p_refs,n_refs,primers,i))

        if not max_pp_reached:
            # Add evaluated pairs to the final list of primer pairs
            for p in pairs:
                if ('test' in p and
                    'penalty' in p['test']):
                    primer_pairs.append(p)

            # Check if max_pp is reached
            if len(primer_pairs) > max_primer_pairs:
                max_pp_reached = True
                log.progress.add('max_pp_reached',
                                 'Maximum number of primer_pairs stored. Now only good primer pairs will be kept!', 'pcr')
            with open(f'{work_dir}primer_pairs.pkl', 'wb') as f:
                pickle.dump((max_pp_reached, primer_pairs), f)

        if log is not None:
            log.progress.add('status_%s'%(name),
                             'Current good pp: %s'%(len(good_pp)), 'pcr')
        else:
            print(' Current good pp: %s'%(len(good_pp)))

        # Dump good primer pairs as backup for reuse purposes, if an error happens
        with open(f'{work_dir}good_primer_pairs.pkl', 'wb') as f:
            pickle.dump((i, too_short, skipped, ignored, no_pairs, good_pp), f)

    # Annotate PCR products with flanks (regions)
    if annotate:
        if log is not None:
            log.progress.add('annotate', 'Annotating PCR product environment', 'pcr')
        else:
            print(' Annotating PCR product environment')
        flank_size = settings['pcr']['annotation']['flank_size']
        regions = {}
        for p in good_pp:
            sid = p['sequence_id']
            start = p['left']['position'] - flank_size
            end   = p['right']['position'] + flank_size + 1
            if sid in regions:
                new_regions = []
                for region in regions[sid]:
                    # Merge overlapping regions into current region
                    if region[0] >= start and region[0] <= end:
                        if region[1] > end:
                            end = region[1]
                    elif region[1] >= start and region[1] <= end:
                        if region[0] < start:
                            start = region[0]
                    elif region[0] <= start and region[1] >= end:
                        start = region[0]
                        end = region[1]
                    else:
                        # Add non-overlapping region to the list
                        new_regions.append(region)
                # Add the current region to the list
                if start < 0: start = 0
                new_regions.append((start, end))
                # Update regions list
                regions[sid] = new_regions
            else:
                if start < 0: start = 0
                regions[sid] = [(start, end)]

        # Extract sequences
        sequences = []
        contigs = dict((n, seq) for seq, n, d in seqs_from_file(contig_file, use_ram_buffer=buffer))
        sid_registry = {}
        for sid, regs in regions.items():
            for i, (start, end) in enumerate(sorted(regs)):
                seq = contigs[sid][start:end]
                name = "%s_%s"%(sid, i)
                desc = ''
                sequences.append((seq, name, desc))
                # keep track of what names exist for each sequence ID
                if not sid in sid_registry:
                    sid_registry[sid] = []
                sid_registry[sid].append(name)

        # Store sequences as fasta
        env_fasta = f'{work_dir}pcr_products_with_skirts.fa'
        save_as_fasta(sequences, env_fasta, add_unique_id=False)

        # Analyse PCR product environments
        annotations = get_blast_annotations(env_fasta,
            settings['pcr']['annotation']['blastx_settings'],
            settings['pcr']['annotation']['blast_db_path'])
        for p in good_pp:
            p_notes = []
            sid     = p['sequence_id']
            start   = p['left']['position']
            end     = p['right']['position'] + 1
            for hit_region_name in sid_registry[sid]:
                for h_start, h_end, notes in annotations.get(hit_region_name, []):
                    if (  (h_start >= start and h_start <= end)
                        or (h_end >= start and h_end <= end)
                        or (h_start <= start and h_end >= end)  ):
                        p_notes.extend(notes)

            # Store the mapped annotations with the pair
            p['annotations'] = list(set(p_notes))

        if log is not None:
            log.progress['annotate'].log_time()
        else:
            print('Done')

    if too_short:
        print("%s sequence were too short to be analysed (<%s)!"%(
              too_short, min_seq_len))
    if no_pairs:
        print(("%s sequence had no identifiable primer pairs!")%(no_pairs))
    if skipped:
        print(("%s sequence were skipped since the requested number of pairs "
               "is reached!")%(skipped))
    if ignored:
        print(("%s sequence were ignored since they were not in the sequence "
               "selection!")%(ignored))
    # SORT and return pairs according to validation and penalty score
    if log is not None:
        log.progress['pcr'].log_time()

    return (sorted(primer_pairs, key=pair_sort, reverse=True),
            sorted(good_pp,      key=pair_sort, reverse=True))


def pair_sort(p):
    ''' Method for deciding sorting rank of primers. '''
    if 'test' in p:
        t = p['test']
        u = sum(map(int, list("{0:b}".format(t['unique_flags']).zfill(3))))
        return (t['sensitivity'], t['specificity'], u, 1/(1+t['noise']), 1/((1+t['penalty'])*(1+p['pair']['penalty'])))
    else:
        return (0,0,0,0,0)


def compute_binding_sites(ref, alignments, probes=[], filter=True):
    ''' Compute binding sites by calculating heterodimer thermodynamics for all
    alignments and filtering those with a melting temperature below 0 degrees

    INPUTS:
        * ref is a path to the reference file.
        * alignments is a dictionary of aligned sequences as keys and a list of
          alignments as values consisting of the contig name, the cigar value of
          the alignment, the strand which was aligned to and the position of the
          alignment.
          EG. {'TGATAAGGCGATGGCAGATGA': [(contig_name, cigar_value, strand, position),...]}
        * probes is a list of aligned sequences which are probes sequences and
          should be handled as such. All sequences not in this list are
          considered primer sequences.

    OUTPUTS:
        * binding_sites is a dictionary of primer sequences as keys and a
          list of alignments as values consisting of the contig name, the strand
          which was aligned to, the position of the alignment and the melting
          tempetrature of the alignment.
          EG. {'TGATAAGGCGATGGCAGATGA': [(contig_name, strand, position, Tm),...]}

    EXAMPLE:
        >>> ref = 'old/reference_tem.fa'
        >>> alignments = {'CAACATTTTCGTGTCGCCCTT': [('reference', '', '+', 1043),
        ...                                         ('reference', '', '+', 300)]}
        >>> compute_binding_sites(ref, alignments)
        {'CAACATTTTCGTGTCGCCCTT': [('reference', '+', 1043, 50.132)]}
    '''
    buffer = settings['input']['use_ram_buffer']

    # Initiate the thermodynamic analyses (defining: p3_primer, p3_probe)
    p3_primer, p3_probe = configure_p3_thermoanalysis()

    # EXTRACT contigs from reference
    contigs = dict((name, seq) for seq, name, desc in seqs_from_file(ref, use_ram_buffer=buffer))

    # Identify binding sites
    binding_sites = {}
    for primer in alignments:
        binding_sites[primer] = []
        # COMPUTE primer details
        primer_rc = reverse_complement(primer)
        primer_len = len(primer)
        for (contig_name, cigar, strand, position) in alignments[primer]:
            tm = None
            # EXTRACT sequence
            seq = contigs[contig_name][position:position+primer_len].upper()
            # CALCULATE the alignment hetero-dimer annealing temperature
            seq2 = primer
            if strand == '+':
                seq2 = primer_rc
            if len(re.findall('[^M0-9]', cigar)) > 0:
                # verify that primer3 is able to analyse the alignment
                s_aln, p_aln, sim = align_seqs(seq, reverse_complement(seq2))
                if sim < settings['pcr']['priming']['alignment_similarity']:
                    tm = 0
            if tm is None:
                if primer in probes:
                    tm = p3_probe.calc_heterodimer(seq, seq2).tm
                else:
                    tm = p3_primer.calc_heterodimer(seq, seq2).tm

            # Filter binding_sites with Tm of 0 or below
            if not filter or tm > 0:
                binding_sites[primer].append((contig_name, strand, position, tm))

    # dump computed binding sites to temporary pickle_file
    name = os.path.basename(ref).rsplit('.',1)[0]
    with open_(f'{work_dir}binding_sites_{name}.pkl', 'wb') as f:
        pickle.dump(binding_sites, f)

    return binding_sites


def grade_binding_sites(binding_sites, probes=[], thresholds=None, min_grade=None):
    ''' Grade the binding sites for each primer sequence based on the
    thermodynamics according to the predefined threshold criteria

    +-------+--------------------+--------------------------+
    | Grade | Description        | Thresholds/boundaries    |
    +-------+--------------------+--------------------------+
    |   0   | Negligible binding |      Tm =  0             |
    |   1   | Unlikely binding   |  0 < Tm < 20             |
    |   2   | Poor binding       | 20 < Tm < 40             |
    |   3   | Detectable binding | 40 < Tm < 55 and Tm > 65 |
    |   4   | Good binding       | 55 < Tm < 65             |
    |   5   | Optimal binding    | 59 < Tm < 61             |
    +-------+--------------------+--------------------------+

    * Note: thresholds above can be modified. Defaults are set around an anneal
            temperature of 60^C.

    INPUTS:
        * binding_sites is a dictionary of primer sequences as keys and a
          list of alignments as values consisting of the contig name, the strand
          which was aligned to, the position of the alignment and the melting
          tempetrature of the alignment.
          EG. {'TGATAAGGCGATGGCAGATGA': [(contig_name, strand, position, Tm), ...]}
        * thresholds is a tuple or list of thresholds for the temperature, which
          if matched increases the grade by one. Each threshold must be either a
          number, which the temperature must be above or a list/tuple of
          two numbers which the temperature must be above and below accordingly
          to be matched. Default is the following 5 thresholds:
          (0, 20, 40, [55,65], [59,61]).
        * min_grade is the threshold grade for the filtering. Any binding site
          with a lower grade than this minimum is not returned.

    OUTPUTS:
        * graded_binding_sites
          EG. {'TGATAAGGCGATGGCAGATGA': [(contig_name, strand, position, Tm, grade), ...]}

    USAGE:
        >>> binding_sites = {'CAACATTTTCGTGTCGCCCTT': [('reference', '+', 1043, 50.132),
        ...                                             'reference', '+', 300, 19.132)]}
        >>> grade_binding_sites(binding_sites)
        {'CAACATTTTCGTGTCGCCCTT': [('reference', '+', 1043, 50.132, 3)]}

    '''
    p3_set = settings["pcr"]["priming"]["primer3"]
    thresholds = (0, 20, 40, [55,65], [59,61]) if thresholds is None else thresholds
    thresholds_probe = (0, 20, 40, p3_set["PRIMER_MIN_TM"], [p3_set["PRIMER_INTERNAL_MIN_TM"], p3_set["PRIMER_INTERNAL_MAX_TM"]])

    min_grade  = 0 if min_grade is None else min_grade
    # Grade and filter binding sites
    graded_binding_sites = {}
    for primer in binding_sites:
        graded_binding_sites[primer] = []
        for (contig_name, strand, position, tm) in binding_sites[primer]:
            # Compute binding grade based on the thermodynamics and defined thresholds
            grade = 0
            ths = thresholds_probe if primer in probes else thresholds
            for threshold in ths:
                try:
                    if isinstance(threshold, list):
                        if tm > threshold[0] and tm < threshold[1]:
                            grade += 1
                    else:
                        if tm > threshold:
                            grade += 1
                except:
                    exit(f"Configuration error: Threshold for grade_binding_sites contained unsupported format! ({threshold})")

            if grade >= min_grade:
                graded_binding_sites[primer].append((contig_name, strand,
                                                     position, tm, grade))
            # else:
            #     print('BAD:', tm, primer, alignment, ref)

    return graded_binding_sites


def estimate_primer_rank(alignments, scheme='positive'):
    ''' Calculate the penalty for the alignments given the scheme (positive or negative)

    For an implicit primer, summarize the alignments of the primer to the
    implicit reference. Based on the alignment summary and the given scheme,
    calculate the penalty based on the penalty rules and settings.

    Penalty rules:
        * Penalize primers not binding well to a positive reference
        * Penalize primers with multiple binding sites to a positive reference
        * Penalize primers that binds to a negative reference
    '''
    penalty = settings["pcr"]["priming"]["penalties"][scheme]
    # Get grade summary
    grade_sum = {0:0,1:0,2:0,3:0,4:0,5:0}
    for aln in alignments:
        # contig_name, strand, pos, tm, grade
        grade_sum[aln[-1]] += 1

    # Compute primer penalty
    penalty_score = 0
    if scheme == 'positive':
        if grade_sum[5] == 0: # No optimal binding site
            if grade_sum[4] == 0: # No good binding site
                if grade_sum[3] == 0: # No expected binding site
                    penalty_score += penalty['no_grade_3']
                else:
                    penalty_score += penalty['no_grade_4']
            else:
                penalty_score += penalty['no_grade_5']
        if grade_sum[5] > 1: # Multiple optimal binding sites
            penalty_score += penalty['multi_grade_5']
        elif grade_sum[4] > 1: # Multiple good binding sites
            penalty_score += penalty['multi_grade_4']
        elif grade_sum[3] > 1: # Multiple expected binding sites
            penalty_score += penalty['multi_grade_3']
        elif grade_sum[2] > 1: # Multiple non-negligible binding sites
            penalty_score += penalty['multi_grade_2']
        elif grade_sum[1] > 1: # Multiple gonegligibleod binding sites
            penalty_score += penalty['multi_grade_1']
    elif scheme == 'negative':
        if grade_sum[5] > 0: # one or more optimal binding sites
            penalty_score += penalty['grade_5']
        elif grade_sum[4] > 0: # one or more good binding sites
            penalty_score += penalty['grade_4']
        elif grade_sum[3] > 0: # one or more expected binding sites
            penalty_score += penalty['grade_3']
        elif grade_sum[2] > 0: # one or more non-negligible binding sites
            penalty_score += penalty['grade_2']
        elif grade_sum[1] > 0: # one or more negligible binding sites
            penalty_score += penalty['grade_1']
    else:
        raise ValueError('Error: Unknown scheme %s'%scheme)

    return penalty_score


def validate_primer_pairs(pairs, p_refs, n_refs, primers, seq_id=None):
    ''' Validate and score primer pairs '''
    qpcr = settings['pcr']['priming']['primer3']['PRIMER_PICK_INTERNAL_OLIGO'] == 1
    mg = settings['pcr']['min_pcr_grade']
    if log is not None:
        log.progress.add('val%s'%seq_id, 'Validating primer pairs', 'pcr')
    plen = len(p_refs)
    nlen = len(n_refs)
    # For each primer pair
    good_pp = []
    pos_counts = [0,0,0,0,0]
    neg_counts = [0,0,0,0,0]
    for p in pairs:
        # COMPUTE validation score
        p['products'] = {'pos': [], 'neg': []}

        # FIND possible matches which fit the product_size criteria
        p_fw = p['left']['sequence']
        p_rv = p['right']['sequence']
        probe = p['internal']['sequence'] if qpcr else ''
        target = p['pair']['product_size']
        hb = target + settings['pcr']['product_deviation']
        lb = target - settings['pcr']['product_deviation']
        sensi = 0
        speci = 0
        noise = 0
        # Set all unique flags true for the pair
        p['test']['unique_flags'] = 7 if qpcr else 3

        # EVALUATE pcr for positive references
        total_count = 0
        small_count = 0
        large_count = 0
        therm_count = 0
        match_count = 0
        for i, ref in enumerate(p_refs):
            # GET primer match data
            try:
                rc_p_rv = reverse_complement(p_rv)
                fw_locs = AdvancedDictionary((
                    (contig_name, (strand, (grade, pos, tm)))
                    for contig_name, strand, pos, tm, grade in primers[p_fw]['pos'][i]))
                rv_locs = AdvancedDictionary((
                    (contig_name, (strand, (grade, pos, tm)))
                    for contig_name, strand, pos, tm, grade in primers[rc_p_rv]['pos'][i]))
                if qpcr:
                    probe_locs = AdvancedDictionary((
                       (contig_name, (strand, (grade, pos, tm)))
                       for contig_name, strand, pos, tm, grade in primers[probe]['pos'][i]))
                else:
                    probe_locs = AdvancedDictionary({})
            except:
                sys.stderr.write("error: Primer matching failed for positive reference %s!\n"%(i))
                raise

            # FIND potential pcr products (OBS: Not filtered on grade yet!)
            products, counts = find_pcr_products(fw_locs, rv_locs, probe_locs,
                                                 len(p_fw), len(p_rv), len(probe),
                                                 filter_on_grade=False)
            total_count += counts[0]
            small_count += counts[1]
            large_count += counts[2]

            # Post process products
            filtered_products = []
            # Find best target match
            best_target = [False, False, False] # [forward, reverse, probe]
            for product in products:
                # product format: [fw_contig_name, fw_strand, fw_position, rv_contig_name, rv_strand, rv_position, [fw_grade, rv_grade], [fw_tm, rv_tm], product_size]
                if product[-1] > lb and product[-1] < hb:
                    target_match = [product[-3][0] >= mg, product[-3][1] >= mg,
                                    qpcr and product[-3][2] >= mg]
                    if sum(target_match) > sum(best_target):
                        best_target = target_match

                # Filter products with bad grades
                if any(g < mg for g in product[-2]):
                    therm_count += 1
                else:
                    filtered_products.append(product)

            # Update primer/probe uniqueness flags
            current_flags = np.array(num2binarray(p['test']['unique_flags'], 3))
            p['test']['unique_flags'] = binarray2num(current_flags * best_target)

            # Update products
            products = filtered_products
            match_count += len(products)

            # UPDATE pair with PCR products
            p['products']['pos'].append(products)

            # Compute sensitivity
            sensi += any(True for x in products if x[-1] > lb and x[-1] < hb)

            # Compute noise
            noise += sum(True for x in products if x[-1] <= lb or x[-1] >= hb)

        pos_counts[0] += total_count
        pos_counts[1] += small_count
        pos_counts[2] += large_count
        pos_counts[3] += therm_count
        pos_counts[4] += match_count

        # Store sensitivity
        p['test']['sensitivity'] = round_sig(sensi / plen, 3)

        # EVALUATE pcr for negative references
        total_count = 0
        small_count = 0
        large_count = 0
        therm_count = 0
        match_count = 0
        for i, ref in enumerate(n_refs):
            # GET primer match data
            try:
                rc_p_rv = reverse_complement(p_rv)
                fw_locs = AdvancedDictionary((
                    (contig_name, (strand, (grade, pos, tm)))
                    for contig_name, strand, pos, tm, grade in primers[p_fw]['neg'][i]))
                rv_locs = AdvancedDictionary((
                    (contig_name, (strand, (grade, pos, tm)))
                    for contig_name, strand, pos, tm, grade in primers[rc_p_rv]['neg'][i]))
                if qpcr:
                    probe_locs = AdvancedDictionary((
                       (contig_name, (strand, (grade, pos, tm)))
                       for contig_name, strand, pos, tm, grade in primers[probe]['neg'][i]))
                else:
                    probe_locs = AdvancedDictionary({})
            except:
                print("error: Primer matching failed for negative reference %s!"%(i))
                raise

            # FIND potential pcr products (OBS: Not filtered on grade yet!)
            products, counts = find_pcr_products(fw_locs, rv_locs, probe_locs,
                                                 len(p_fw), len(p_rv), len(probe),
                                                 filter_on_grade=False)
            total_count += counts[0]
            small_count += counts[1]
            large_count += counts[2]

            # Post process products
            filtered_products = []
            for product in products:
                # product format: [fw_contig_name, fw_strand, fw_position, rv_contig_name, rv_strand, rv_position, [fw_grade, rv_grade], [fw_tm, rv_tm], product_size]
                # Identify binding primers/probe
                good_binders = np.zeros(3, dtype=bool)
                for j, g in enumerate(product[-3]):
                    if g >= mg: good_binders[j] = True

                if good_binders.any():
                    # Update primer/probe uniqueness flags
                    remove_flags = np.logical_not(good_binders)
                    current_flags = np.array(num2binarray(p['test']['unique_flags'], 3))
                    p['test']['unique_flags'] = binarray2num(current_flags * remove_flags)

                # Filter products with bad grades
                if any(g < mg for g in product[-3]):
                    therm_count += 1
                else:
                    filtered_products.append(product)

            # Update products
            products = filtered_products
            match_count += len(products)

            # UPDATE pair with PCR products
            p['products']['neg'].append(products)

            # Compute specificity
            # Use product size deviation limit, to determine when the product size is indistiguishable from target size.
            hb = target + settings['pcr']['product_deviation']
            lb = target - settings['pcr']['product_deviation']
            speci += not any(True for x in products if x[-1] > lb and x[-1] < hb)

            # Compute noise
            noise += sum(True for x in products if x[-1] <= lb or x[-1] >= hb)

        neg_counts[0] += total_count
        neg_counts[1] += small_count
        neg_counts[2] += large_count
        neg_counts[3] += therm_count
        neg_counts[4] += match_count

        # Store specificity
        p['test']['specificity'] = round_sig(speci / nlen, 3) if nlen > 0 else 1.0

        # Store noise
        p['test']['noise'] = round_sig(noise / (plen + nlen), 3)

        # CHECK pair performance
        if p['test']['sensitivity'] == 1 and p['test']['specificity'] == 1:
            too_close_to_other_pairs = False
            # Check distance to previous pairs
            for gp in good_pp:
                dist = (abs(p['left']['position']  - gp['left']['position']) +
                        abs(p['right']['position'] - gp['right']['position']))
                if dist < settings['pcr']['min_pair_distance']:
                    too_close_to_other_pairs = True
                    break

            if not too_close_to_other_pairs:
                # Add to list of good primer pair candidates
                good_pp.append(p)

    if log is not None:
        log.progress['val%s'%seq_id].log_time()
        log.stats.add_row('pair%s'%seq_id, [len(good_pp), 'Non-overlapping pairs'])
        # Create PCR product statistics table
        headers = ['considered', 'too small', 'too big', 'low tm', 'possible']
        log.stats.add_table('pcr%s'%seq_id,
                            'PCR product analysis %s'%seq_id,
                            ['Attribute', 'Positive', 'Negative'])
        for k, h in enumerate(headers):
            log.stats.add_row('pcr%s'%seq_id, [h, pos_counts[k], neg_counts[k]])

    return good_pp


def find_pcr_products(fw_locs, rv_locs, probe_locs, fw_len, rv_len, probe_len,
                       filter_on_grade=True):
    ''' Find primer pair products

    PRIMER/PCR GRADES
        0 = no priming/probe  (alignment_similarity < 0.6)
        1 = unlikely priming  (Tm < threshold_tm)
        2 = potential priming (Tm > threshold_tm)
        3 = probable priming  (Tm = target_tm +- 5*C)
        4 = definite priming  (Tm = target_tm +- 1*C)

    # If quantitative PCR (qPCR), also known as real-time PCR is chosen, the
    # probe is included in the grading process

    products returned are in the following format:
    [forward_contig_name, forward_strand_indicator, forward_position, reverse_contig_name, reverse_strand_indicator, reverse_position, pcr_grades, pcr_tms, product_size]

    contig_name = name of the sequences in the reference fasta file
    strand_indicator = + or -
    position = the location of the match in the reference
    pcr_grades = [grade_of_forward_primer, grade_of_reverse_primer]
    pcr_tms = [Tm_of_forward_primer, Tm_of_reverse_primer]
    product_size = The size of the predicted PCR product/amplicon

    '''
    p3_args = settings['pcr']['priming']['primer3']
    mg = settings['pcr']['min_pcr_grade']
    qpcr = p3_args['PRIMER_PICK_INTERNAL_OLIGO'] == 1
    if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'], list):
        if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0], list):
            product_length_limits = [lim for lims in p3_args['PRIMER_PRODUCT_SIZE_RANGE'] for lim in lims]
            min_len, max_len = min(product_length_limits), max(product_length_limits)
        else:
            min_len, max_len = p3_args['PRIMER_PRODUCT_SIZE_RANGE']
    else:
        exit('Settings error: invalid input format for PRIMER_PRODUCT_SIZE_RANGE, expected "list", got: ' + p3_args['PRIMER_PRODUCT_SIZE_RANGE'])
    products = []
    total_count = 0
    small_count = 0
    therm_count = 0
    # print(f"min grade: {mg}\nqpcr: {qpcr}\nmin_len: {min_len}\nmax_len: {max_len}")
    # For each contig
    for contig_name in fw_locs:
        if contig_name not in rv_locs: continue # Skip if reverse is unmapped
        total_count += len(fw_locs[contig_name]) * len(rv_locs[contig_name])
        # Devide into + and - strand and sort the positions
        p_fw_list = AdvancedDictionary(fw_locs[contig_name], sort=True, key=lambda x: x[1])
        p_rv_list = AdvancedDictionary(rv_locs[contig_name], sort=True, key=lambda x: x[1])
        if qpcr:
            if contig_name not in probe_locs: continue # Skip if probe is unmapped
            probe_list = AdvancedDictionary(probe_locs[contig_name])
            probe_plus = np.array(probe_list['+']) if '+' in probe_list else None
            probe_minus = np.array(probe_list['-']) if '-' in probe_list else None
        # Find pairs on the + strand
        if '+' in p_fw_list and '+' in p_rv_list:
            for g_fw, p_fw, t_fw in p_fw_list['+']: # grade, pos, tm
                for g_rv, p_rv, t_rv in p_rv_list['+']: # grade, pos, tm
                    # Get PCR grade and product size
                    product_size = p_rv + rv_len - p_fw
                    pcr_grades = [g_fw, g_rv]
                    pcr_tms = [t_fw, t_rv]
                    if qpcr:
                        if probe_plus is not None:
                            # Select probes that are within the forward and reverse primer bind site
                            probes = probe_plus[np.logical_and(probe_plus[:,1]>(p_fw+fw_len), probe_plus[:,1]<(p_rv-probe_len))]
                        else:
                            probes = np.array([])
                        # Append the best PCR grade for the probe
                        pcr_grades.append(0 if probes.shape[0] == 0 else int(probes[:,0].max()))
                        # Append the best PCR Tm for the probe
                        pcr_tms.append(0 if probes.shape[0] == 0 else probes[:,2].max())
                    # Filter bad products
                    if filter_on_grade and any(g < mg for g in pcr_grades):
                        therm_count += 1
                        continue
                    # Filter small products
                    if product_size < min_len:
                        small_count += 1
                        continue
                    # Filter big products
                    if product_size > max_len:
                        break
                    # Add match to alignment_list
                    products.append([contig_name, '+', p_fw,
                                     contig_name, '+', p_rv,
                                     pcr_grades, pcr_tms, product_size])

        # Find pairs on the - strand
        if '-' in p_fw_list and '-' in p_rv_list:
            p_fw_list.sort(key=lambda x: x[1], reverse=True)
            p_rv_list.sort(key=lambda x: x[1], reverse=True)
            for g_fw, p_fw, t_fw in p_fw_list['-']: # grade, pos, tm
                for g_rv, p_rv, t_rv in p_rv_list['-']: # grade, pos, tm
                    product_size = p_fw + fw_len - p_rv
                    pcr_grades = [g_fw, g_rv]
                    pcr_tms = [t_fw, t_rv]
                    if qpcr:
                        if probe_minus is not None:
                            # Select probes that are within the forward and reverse primer bind site
                            probes = probe_minus[np.logical_and(probe_minus[:,1]>(p_rv+rv_len), probe_minus[:,1]<(p_fw-probe_len))]
                        else:
                            probes = np.array([])
                        # Append the best PCR grade for the probe
                        pcr_grades.append(0 if probes.shape[0] == 0 else int(probes[:,0].max()))
                        # Append the best PCR Tm for the probe
                        pcr_tms.append(0 if probes.shape[0] == 0 else probes[:,2].max())
                    # Filter bad products
                    if filter_on_grade and any(g < mg for g in pcr_grades):
                        therm_count += 1
                        continue
                    # Filter small products
                    if product_size < min_len:
                        small_count += 1
                        continue
                    # Filter big products
                    if product_size > max_len:
                        break
                    # Add match to alignment_list
                    products.append([contig_name, '-', p_fw,
                                     contig_name, '-', p_rv,
                                     pcr_grades, pcr_tms, product_size])

    large_count = total_count - small_count - therm_count - len(products)
    return products, (total_count, small_count, large_count, therm_count)


def present_pairs_full(pairs):
    ''' Presents essential aspects of the primer pairs
    USAGE
        >>> print(present_pairs_full(pairs[:3]))
        # product_size	tm	forward primer	reverse primer	test score
        515	60.0	CAACATTTTCGTGTCGCCCTT	TCGTCGTTTGGTATGGCTTCA	0
        605	60.0	CAACATTTTCGTGTCGCCCTT	AATTGTTGCCGGGAAGCTAGA	0
        636	60.0	CAACATTTTCGTGTCGCCCTT	TGCAACTTTATCCGCCTCCAT	0
    '''
    def _repr(i, p, h):
        return '\t'.join(["%s"]*len(h))%(
            i,
            p['sequence_id'],
            p['pair']['product_size'],
            p['test']['unique_flags'],
            p['test']['sensitivity'],
            p['test']['specificity'],
            p['test']['noise'],
            p['test']['penalty'],
            p['pair']['penalty'],
            p['left']['sequence'],
            "%.1f"%p['left']['tm'],
            p['left']['length'],
            "%.1f"%p['left']['gc_percent'],
            p['left']['position'],
            p['right']['sequence'],
            "%.1f"%p['right']['tm'],
            p['right']['length'],
            "%.1f"%p['right']['gc_percent'],
            p['right']['position'],
            p['internal']['sequence'] if p['internal'] else '',
            "%.1f"%p['internal']['tm'] if p['internal'] else '',
            p['internal']['length'] if p['internal'] else '',
            "%.1f"%p['internal']['gc_percent'] if p['internal'] else '',
            p['internal']['position'] if p['internal'] else '',
            '; '.join(p['annotations']) if 'annotations' in p else ''
            )
    h = ['pair', 'sequence_id', 'product_size', 'unique_flags', 'sensitivity',
         'specificity', 'noise', 'penalty', 'p3_penalty', 'forward_primer',
         'forward_tm', 'forward_length', 'forward_gc%', 'forward_position',
         'reverse_primer', 'reverse_tm', 'reverse_length', 'reverse_gc%',
         'reverse_position', 'probe', 'probe_tm', 'probe_length',
         'probe_gc%', 'probe_position', 'annotation']
    return "#%s\n%s"%('\t'.join(h), '\n'.join((_repr(i, p, h) for i, p in enumerate(pairs))))


def predict_pcr_results(refs, pairs, fail_on_non_match=False, tm_thresholds=None, min_grade=None):
    ''' Predict the results of a PCR

    USAGE
        >>> predict_pcr_results(refs, pairs)
        Forward primer has %s binding locations hereof %s are good binders.
        Reverse primer has %s binding locations hereof %s are good binders.
        Potential products sizes: %s
    '''
    # Get directories
    global work_dir, ref_dir, result_dir
    work_dir, ref_dir, result_dir = get_directories(get_ref_dir=True)

    buffer = settings['input']['use_ram_buffer']
    # Prepare primers and probes
    primers = []
    probes = []
    for pair in pairs:
        primers.append((pair[0], {}))
        primers.append((reverse_complement(pair[1]), {}))
        if len(pair) == 3 and pair[2] is not None:
            primers.append((pair[2], {}))
            probes.append(pair[2])

    # Store primer and probe sequences as fasta (required for BLASTing)
    primers_fa = f'{work_dir}primers.fa'
    save_as_fasta(dict(primers), primers_fa)

    # Align primers to references (BLASTing)
    pcr_results = [['NA' for r in refs] for p in pairs]
    for r, ref in enumerate(refs):
        alignments = blast_to_ref(ref, primers_fa,
                                  settings['pcr']['priming']['blastn_settings'],
                                  buffer)
        print(f"Stats: Of the {len(primers)} blasted primer/probe sequences, {len(alignments)} sequences had matches with a total of     {sum(map(len, alignments.values()))} alignments for {ref}")
        # Filter primer matches which fail the thermodynamic test
        binding_sites = compute_binding_sites(ref, alignments, probes, filter=False)
        print(f"Stats: {sum(map(len, binding_sites.values()))} binding sites found for {ref}")
        binding_sites = grade_binding_sites(binding_sites, probes, tm_thresholds, min_grade)
        print(f"Stats: {sum(map(len, binding_sites.values()))} graded binding sites found for {ref}")
        for p, pair in enumerate(pairs):
            fw = pair[0]
            rv = pair[1]
            if len(pair) == 3 and pair[2] is not None:
                probe = pair[2]
            else:
                probe = ''

            rv_rc = reverse_complement(rv)
            # Validate and score primer pair
            try:
                fw_locs = AdvancedDictionary((
                    (contig_name, (strand, (grade, pos, tm)))
                    for contig_name, strand, pos, tm, grade in binding_sites[fw]))
                rv_locs = AdvancedDictionary((
                    (contig_name, (strand, (grade, pos, tm)))
                    for contig_name, strand, pos, tm, grade in binding_sites[rv_rc]))
                if probe:
                    probe_locs = AdvancedDictionary((
                        (contig_name, (strand, (grade, pos, tm)))
                        for contig_name, strand, pos, tm, grade in binding_sites[probe]))
                else:
                    probe_locs = AdvancedDictionary({})
            except:
                if fail_on_non_match:
                    print("error: Primer matching failed for PCR prediction (%s)!"%(ref))
                    raise
                else:
                    # Skip pairs, where one of the primers do not match anything
                    pcr_results[p][r] = ''
                    continue

            # FIND potential pcr products
            products, counts = find_pcr_products(fw_locs, rv_locs, probe_locs,
                                                 len(fw), len(rv), len(probe))
            total_count = counts[0]
            small_count = counts[1]
            large_count = counts[2]
            therm_count = counts[3]
            match_count = len(products)
            print(f"Stats: PCR counts total={total_count}, small={small_count}, large={large_count}, therm={therm_count},     match={match_count} for {ref}")

            pcr_results[p][r] = products

    return pcr_results


def print_pcr_results(refs, pcr_results, output=None):
    """ Print PCR results

    Fetch pair index (i), product size (p[-1]) and calculate binding grade for the product based on the binding grade for fwd and rvs primer (p[6][0] and p[6][1])
    * ppp = products per pair (rows of output)
    * ppr = products per reference (columns of output)
    * p   = product
    product format: [fw_contig_name, fw_strand, fw_position, rv_contig_name, rv_strand, rv_position, [fw_grade, rv_grade (, pr_grade)], [fw_tm, rv_tm (, pr_tm)], product_size]

    """
    # Check if alternative IDs for the input files exists
    try:
        input_details = {os.path.basename(x['filepath']):x['id'] for x in parse_tsv('inputs/inputs.tsv') if 'filepath' in x}
    except:
        input_details = {}
    if output is not None:
        if '/' not in output or os.path.exists(os.path.dirname(output)):
            with open_(output, 'w') as f:
                f.write(f"Pair\\ref\t" + '\t'.join((input_details.get(os.path.basename(r), os.path.basename(r)) for r in refs)))
                f.write('\n')
                f.write('\n'.join([f"{i}:\t" + '\t'.join((','.join([f"{p[-1]}({round_sig(p[6][0]*p[6][1]/25*p[6][2]/5) if len(p[6]) == 3 else round_sig(p[6][0]*p[6][1]/25)})" for p in ppr]) for ppr in ppp)) for i, ppp in enumerate(pcr_results)]))
                f.write('\n')
        else:
            f.write(f"Pair\\ref\t" + '\t'.join((input_details.get(os.path.basename(r), os.path.basename(r)) for r in refs)))
            print('\n'.join([f"{i}:\t" + '\t'.join((','.join([f"{p[-1]}({round_sig(p[6][0]*p[6][1]/25*p[6][2]/5) if len(p[6]) == 3 else round_sig(p[6][0]*p[6][1]/25)})" for p in ppr]) for ppr in ppp)) for i, ppp in enumerate(pcr_results)]))
            raise OSError('Error: Output path directory do not exist!')
    else:
        f.write(f"Pair\\ref\t" + '\t'.join((input_details.get(os.path.basename(r), os.path.basename(r)) for r in refs)))
        print('\n'.join([f"{i}:\t" + '\t'.join((','.join([f"{p[-1]}({round_sig(p[6][0]*p[6][1]/25*p[6][2]/5) if len(p[6]) == 3 else round_sig(p[6][0]*p[6][1]/25)})" for p in ppr]) for ppr in ppp)) for i, ppp in enumerate(pcr_results)]))


def get_primer_and_probe_bindsites(results_file, unique_only=False):
    ''' Extract the sequences and the primer and probe location found in the
    sequences.
    '''
    count = 0
    pribind_locs = {}
    probind_locs = {}
    with open_(results_file) as f:
        for l in f:
            if l.isspace(): continue
            if l.startswith('#'): continue
            tmp = list(map(lambda x: x.strip(), l.split('\t')))
            sequence_id = tmp[0]
            target_size = int(tmp[1])
            unique_flags = num2binarray(int(tmp[2]),3)
            sensitivity = float(tmp[3])
            specificity = float(tmp[4])
            noise = float(tmp[5])
            if not sequence_id in pribind_locs:
                pribind_locs[sequence_id] = multirange()
                probind_locs[sequence_id] = multirange()

            if not unique_only or unique_flags[0]:
                fw_len = int(tmp[10])
                fw_pos = int(tmp[12])
                pribind_locs[sequence_id].add(fw_pos,fw_pos+fw_len-1)
            # else:
            #     print('Ignored Forward: %s %s %s'%(tmp[10], tmp[12], unique_flags))

            if not unique_only or unique_flags[1]:
                rv_len = int(tmp[15].strip())
                rv_pos = int(tmp[17].strip())
                pribind_locs[sequence_id].add(rv_pos-rv_len+1,rv_pos)
            # else:
            #     print('Ignored Reverse: %s %s %s'%(tmp[10], tmp[12], unique_flags))

            if not unique_only or unique_flags[2]:
                try:
                    pb_len = int(tmp[20].strip())
                    pb_pos = int(tmp[22].strip())
                except: pass
                else:
                    probind_locs[sequence_id].add(pb_pos,pb_pos+pb_len-1)

    return pribind_locs, probind_locs


def print_primer_and_probe_bindsites(seq_file, pribind_locs, probind_locs={},
                                     tag1='31', tag2='32', contigs=None):
    '''  '''
    buffer = settings['input']['use_ram_buffer']
    # Fetch Sequence data
    for seq, name, desc in seqs_from_file(seq_file, use_ram_buffer=buffer):
        if contigs is not None:
            if not name in contigs:
                continue
        elif not name in pribind_locs: continue

        print('SEQUENCE NAME: %s'%name)
        check1 = name in pribind_locs
        check2 = name in probind_locs
        if check1 and check2:
            print(color_seq(seq, pribind_locs[name] + probind_locs[name],
                            tag_mods=[['0'], ['0','1',tag1], ['0','1',tag2],
                                      ['0','1',tag1,tag2]]))
        elif check1:
            print(color_seq(seq, pribind_locs[name],
                            tag_mods=[['0'], ['0','1',tag1], ['0','1',tag2],
                                      ['0','1',tag1,tag2]]))
        elif check2:
            print(color_seq(seq, probind_locs[name],
                            tag_mods=[['0'], ['0','1',tag1], ['0','1',tag2],
                                      ['0','1',tag1,tag2]]))
        else:
            print(seq)
        print('\n')


def color_seq(seq, ranges,
              tag_mods=[['0'], ['0','1','31'], ['0','1','43'], ['0','1','31','43']]):
    ''' default tag color is red '''
    default_text_mods = '\x1B[%sm'%';'.join(tag_mods[0])
    seq2 = []
    incision_pos = 0
    for start, end, tag in ranges:
        if start > incision_pos:
            seq2.append(default_text_mods + seq[incision_pos:start])

        text_mods = '\x1B[%sm'%';'.join(tag_mods[tag])
        seq2.append(text_mods + seq[start:end+1])
        incision_pos = end +1

    # Reset coloring
    seq2.append(default_text_mods + seq[incision_pos:] +'\x1B[0m')
    # print Sequence
    return ''.join(seq2)
