# -*- coding: utf-8 -*-
"""
    k-mer library
    ~~~~~~~~~~~~~~~~~~~~~

    Contains functionalities relating to k-mer operations.


"""
# (c) 2023 Martin Thomsen

import os, pickle
import numpy as np

from rucs import settings, get_directories
from rucs.log import log
from rucs.file import open_, check_file_type, save_as_fastq, save_as_fasta, seqs_from_file
from rucs.bwa import align_to_ref
from rucs.seq import analyse_genome, reverse_complement, split_scaffolds_into_contigs


def find_unique_core_sequences(positives, negatives, reference, kmer_size=20):
    ''' Find unique core sequences which are sequences only found in the positive genomes

    >>> cs_files, ucs_files = find_unique_core_sequences(positives, negatives, reference, kmer_size):
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

    if log is not None:
        # Initialize k-mer and sequences statistics logging
        log.progress.add('ucs', 'Finding unique core sequences', 'main')
        log.stats.add_table('kmer', 'k-mer analyses', ['k-mers', 'Process'])
        log.stats.add_table('seqs', 'Sequence Analyses',
                            ['Fasta file', 'Sequences', 'Size in bases',
                             'Seqs >%s'%min_seq_len, 'Size >%s'%min_seq_len])
        # Analyse reference sequences
        counts = analyse_genome(reference, min_seq_len)
        log.stats.add_row('seqs', [os.path.basename(reference)] + counts)

    # Identify core sequences through computation of intersecting k-mers from positive references
    core_kmers, cs_files = compute_core_sequences(positives, reference, kmer_size, min_seq_len)

    if log is not None:
        # Analyse contigs of core sequences
        counts = analyse_genome(cs_files[1], min_seq_len)
        log.stats.add_row('seqs', [os.path.basename(cs_files[1])] + counts)

    # Identify unique sequences through computation of complementing k-mers to negative references
    unique_core_kmers, ucs_files = remove_common_sequences(core_kmers, negatives, reference, kmer_size, min_seq_len)

    if log is not None:
        # Analyse contigs of unique core sequences
        counts = analyse_genome(ucs_files[1], min_seq_len)
        log.stats.add_row('seqs', [os.path.basename(ucs_files[1])] + counts)
        log.progress['ucs'].log_time()

    return cs_files, ucs_files


def compute_core_sequences(positives, reference, kmer_size=20, min_seq_len=100):
    ''' Compute the intersecting k-mers from positive references, and
    generate the core sequences from the common kmers

    >>> compute_core_sequences(positives, negatives, reference, kmer_size, min_seq_len)
    '''
    # Compute intersecting k-mers from positive references
    core_kmers = find_intersecting_kmers(positives, kmer_size, min_seq_len)
    if len(core_kmers) == 0:
        raise UserWarning('No intersecting k-mers could be identified!')

    # Store k-mers as fastq
    c_kmers_fq = 'core_kmers.fq'
    save_as_fastq(core_kmers, c_kmers_fq)

    # Align k-mers to reference (Note: Unmapped k-mers are lost in this process)
    core_kmers = align_to_ref(reference, c_kmers_fq, settings['ucs']['bwa_settings'],
                              'core_kmers', settings['ucs']['sam_flags_ignore'])

    # Compute core sequences
    if log is not None:
        log.progress.add('make_cs', 'Computing k-mer contigs and scaffolds','ucs')

    cs_files = compute_consensus_sequences(core_kmers, reference, kmer_size,
                                           settings['ucs']['charspace'], 'core_sequences')
    if log is not None:
        log.progress['make_cs'].log_time()

    return core_kmers, cs_files


def remove_common_sequences(core_kmers, negatives, reference, kmer_size=20, min_seq_len=100):
    ''' Compute and remove complementing k-mers to negative references, and
    generate the unique core sequences from the remaining kmers

    >>> unique_core_kmers, ucs_files = remove_common_sequences(core_kmers, negatives, reference, kmer_size, min_seq_len)
    '''
    # Compute complementing k-mers to negative references
    unique_core_kmers = find_complementing_kmers(core_kmers, negatives, kmer_size, min_seq_len)
    if len(unique_core_kmers) == 0:
        raise UserWarning('No complementing k-mers could be identified!')

    # Compute unique core sequences
    if log is not None:
        log.progress.add('make_ucs', 'Computing k-mer contigs and scaffolds',
                         'ucs')
    ucs_files = compute_consensus_sequences(unique_core_kmers, reference,
                                            kmer_size, settings['ucs']['charspace'],
                                            'unique_core_sequences')
    if log is not None:
        log.progress['make_ucs'].log_time()

    return unique_core_kmers, ucs_files


def find_intersecting_kmers(files, kmer_size=20, min_seq_len=100):
    ''' Compute the intersection of kmer sequences across the set of files

    * get kmers from file (RAM)
    * store kmer object as tmp file (cpickle)
    * clean memory
    * retrieve and merge tmp file objects (RAM)
      - Filter/ignore those kmers with low counts
    * store object as file (cpickle)
    * clean up memory
    * delete tmp files
    '''
    if log is not None:
        log.progress.add('core', 'Computing k-mer intersection', 'ucs')
    #Extract kmers
    kmers = {}
    for i, genome in enumerate(files):
        # ADD GENOME TO SET
        gname = genome.split('/')[-1]
        # EXTRACT K-MERS FROM SEQEUNCE DATA FROM INPUT FILES
        to_upper = settings['input']['to_upper']
        buffer = settings['input']['use_ram_buffer']
        kmers_i = extract_kmers_from_file(genome, i, kmer_size, '',
                                          settings['ucs']['min_kmer_count'],
                                          settings['ucs']['rev_comp'],
                                          min_seq_len, to_upper, buffer, 'core')
        # COMPUTE INTERSECTION
        if i > 0: kmers = set_op('intersection', kmers, kmers_i, method='list')
        else: kmers = kmers_i

    if log is not None:
        log.progress['core'].log_time()
        log.stats.add_row('kmer',
                          [len(kmers), 'Computation of k-mer intersection'])

    return kmers

def extract_kmers_from_file(filename, genome_prefix='', kmer_size=20,
                            kmer_prefix='', fastq_kmer_count_threshold=20,
                            revcom=True, min_seq_len=100, to_upper=False,
                            buffer=False, log_entry=None, reuse=False):
    '''
    NAME      Extract K-mers from sequence
    AUTHOR    Martin CF Thomsen
    DATE      18 Jul 2013
    DESCRIPTION
        Extract K-mers from a file with sequence data, and return a list of unique
        k-mers.
    ARGUMENTS
        filename: string which contain a path to the input file
        genome_prefix: A prefix which is added to the position note of the
           K-mers.
        kmer_size: The size of the generated kmers. (Sequences shorter than the
           K-mer size will not generate K-mers.)
        kmer_prefix: A prefix used to filter which K-mers are stored. (Only
           K-mers starting with this prefix will be added.)
    SUPPORTED FILE FORMATS
        fasta
        fastq
    USAGE
        >>> kmers = extract_kmers_from_file('test.fa', '1', 13, '')
        # Extracting Kmers from inputfile test.fa...
        # Found 4 sequences and 4 Kmers in 0.000451 seconds
        >>> kmers
        {'this_is_seq_4': 1, 'this_is_seq_3': 1,
         'this_is_seq_2': 1, 'this_is_seq_1': 1}
    '''
    if reuse and os.path.exists(f'{work_dir}ext_kmer_reuse_{os.path.basename(filename)}.pkl'):
        # Extract K-mers from previous result file
        if log is not None and log_entry is not None:
            log.progress.add(f'ext_{os.path.basename(filename)}',
                             f'Reloading previous extracted k-mers for {os.path.basename(filename)}', log_entry)

        with open_(f'{work_dir}ext_kmer_reuse_{os.path.basename(filename)}.pkl', 'rb') as f:
            kmers = pickle.load(f)

        if log is not None and log_entry is not None:
            kmersum = sum(kmers.values())
            log.stats.add_row('kmer', [len(kmers), f'Extracted k-mers from {os.path.basename(filename)}'])
            log.stats.add_row('kmer', [kmersum - len(kmers), f'Doublicate k-mers from {os.path.basename(filename)}'])
            log.progress[f'ext_{os.path.basename(filename)}'].log_time()
    else:
        if log is not None and log_entry is not None:
            log.progress.add(f'ext_{os.path.basename(filename)}',
                             f'Extracting k-mers from {os.path.basename(filename)}', log_entry)
        # Extract K-mers from the sequences
        kmers = {}
        seqcount = 0
        kmerstot = 0
        file_type = check_file_type([filename])
        for i, (seq, name, desc) in enumerate(seqs_from_file(filename,
                                                             to_upper=to_upper,
                                                             use_ram_buffer=buffer)):
            if len(seq) < min_seq_len: continue # Skip small sequences

            # Extract kmers from sequence (GenomePrefix_SequencePrefix_KmerPosition) to 'kmers'
            extract_kmers(kmers, seq, kmer_size, kmer_prefix)
            kmerstot += len(seq) - kmer_size + 1

            # Extract kmers from the reverse complement sequence
            if revcom:
                extract_kmers(kmers, reverse_complement(seq), kmer_size, kmer_prefix)
                kmerstot += len(seq) - kmer_size + 1

            seqcount += 1

        if log is not None and log_entry is not None:
            kmersum = sum(kmers.values())
            log.stats.add_row('kmer', [len(kmers), f'Extracted k-mers from {os.path.basename(filename)}'])
            log.stats.add_row('kmer', [kmerstot - kmersum, f"Skipped k-mers from {os.path.basename(filename)}"])
            log.stats.add_row('kmer', [kmersum - len(kmers), f'Doublicate k-mers from {os.path.basename(filename)}'])
            log.progress[f'ext_{os.path.basename(filename)}'].log_time()

        # File type dependend filtering
        if file_type == 'fastq':
            # Filter K-mers with low occurences
            if log is not None and log_entry is not None:
                log.progress.add(f'flt_{os.path.basename(filename)}',
                                 f'Filtering k-mers with low coverage', log_entry)
            sum_ = 0
            kmer_count = len(kmers)
            for kmer, count in kmers.items():
                if count < fastq_kmer_count_threshold:
                    del kmers[kmer]
                    sum_ += count
            kmer_count -= len(kmers)
            if log is not None and log_entry is not None:
                log.stats.add_row('kmer', [kmer_count, f'Filtered k-mers from {os.path.basename(filename)} with an average       coverage of {sum_/kmer_count:.4f}'])
                log.progress[f'flt_{os.path.basename(filename)}'].log_time()

        if reuse:
            # Store k-mer extraction results for later reuse
            with open_(f'{work_dir}ext_kmer_reuse_{os.path.basename(filename)}.pkl', 'wb') as f:
                pickle.dump(kmers, f)

    return kmers

def extract_kmers(kmers, seq, size=20, kmer_prefix='', charspace=list('ATGC')):
    '''
    NAME      Extract K-mers from sequence
    AUTHOR    Martin CF Thomsen
    DATE      18 Jul 2013
    DESCRIPTION
        Extract K-mers from seqeunce data (seq), and add it to the K-mer
        dictionary (kmers).
    ARGUMENTS
        kmers: The dictionary where the found K-mers are stored.
        seq: The sequence from which the kmers are generated.
        size: The size of the generated kmers. (Sequences shorter than the K-mer
            size will not generate K-mers.)
        kmer_prefix: A prefix used to filter which K-mers are stored. (Only
            K-mers starting with this prefix will be added.)
    USAGE
        >>> let = {}; extract_kmers(let, "abcdefghijklmnopqrstuvwxyz", 5, '', 'efghijkl'); let
        {'efghi': 1, 'fghij': 1, 'ghijk': 1, 'hijkl': 1}
        >>> kmers = {}; extract_kmers(kmers, 'ATGTCSTAGAGGNGCTA', 5); kmers
        {'ATGTC': 1, 'TAGAG': 1, 'AGAGG': 1}
        >>> kmers = {}; extract_kmers(kmers, 'AAAGGNGCTA', 5); kmers
        {'AAAGG': 1}
        >>> kmers = {}; extract_kmers(kmers, 'ACGNAATTGCACTGATACCGCCGGCNTGNGAGAGGCAAGCGATGACGAG'); kmers
        {'AATTGCACTGATACCGCCGG': 1, 'ATTGCACTGATACCGCCGGC': 1, 'GAGAGGCAAGCGATGACGAG': 1}
        >>> kmers = {}; extract_kmers(kmers, 'NACGAATTGCACTGCCGCCGGCNTGNGAGAGGCAAGCGATGACGAG'); kmers
        {'ACGAATTGCACTGCCGCCGG': 1, 'CGAATTGCACTGCCGCCGGC': 1, 'GAGAGGCAAGCGATGACGAG': 1}
        >>> kmers = {}; extract_kmers(kmers, 'AAAGTNNAAAA',5); kmers
        {'AAAGT': 1}
        >>> kmers = {}; extract_kmers(kmers, 'AAAGTNAAAAN',5); kmers
        {'AAAGT': 1}
    '''
    # ITERATE THE SEQUENCE TO FIND K-MERS
    kmer_prefix_len = len(kmer_prefix)
    # Find a valid starting sequence
    j = 0
    h = j + size - 2
    while h >= j:
        if seq[h] in charspace:
            h -= 1
        else:
            # Reset and search again
            j = h + 1
            h = j + size - 2
            if h > len(seq):
                j = len(seq)
                break
    # Find k-mers
    while j <= len(seq) - size:
        if seq[j+size-1] not in charspace:
            # Skip size ahead to get past the invalid character
            j += size
            # Find the first valid sequence
            if j < len(seq)-size+1:
                eol = False
                h = j + size - 1
                while h >= j:
                    if seq[h] in charspace:
                        h -= 1
                    else:
                        # Reset and search again
                        j = h + 1
                        h = j + size - 1
                        if h > len(seq) - 1:
                            eol = True
                            j = len(seq)
                            break
            else:
                # No more possible kmers - skip the rest
                j = len(seq)
                eol = True
            if eol:
                break
        kmer = seq[j:j+size]
        # CHECK IF K-MER STARTS WITH SPECIFIED PREFIX
        if kmer_prefix == kmer[:kmer_prefix_len]:
            # ADD K-MER
            if kmer in kmers: kmers[kmer] += 1
            else: kmers[kmer] = 1
        j += 1


def set_op(op, a, b, method='number'):
    ''' Calculate the intersection for sets or dictionaries.

    Set Operations:
                     union (|): {1,2,3} | {3,4,5} = {1, 2, 3, 4, 5}
              intersection (&): {1,2,3} & {3,4,5} = {3}
      symmetric difference (^): {1,2,3} ^ {3,4,5} = {1, 2, 4, 5}
                complement (-): {1,2,3} - {3,4,5} = {1, 2}

    USAGE
        >>> set_op('union', {1,2,3}, {3,4,5})
        {1, 2, 3, 4, 5}
        >>> set_op('intersection', {1,2,3}, {3,4,5})
        {3}
        >>> set_op('symmetric difference', {1,2,3}, {3,4,5})
        {1, 2, 4, 5}
        >>> set_op('complement', {1,2,3}, {3,4,5})
        {1, 2}
        >>> setA = {1:1,2:2,3:3,4:4}
        >>> setB = {3:3,4:4,5:5,6:6}
        >>> set_op('i', setA, setB)
        {3: 6, 4: 8}
        >>> setA = {1:[1],2:[2],3:[3],4:[4]}
        >>> setB = {3:[3],4:[4],5:[5],6:[6]}
        >>> set_op('i', setA, setB, 'list')
        {3: [3, 3], 4: [4, 4]}
    '''
    if op in ['union', 'uni', 'u', 'U', 'join', 'or', '|']:
        op = 'union'
        collection = set(a) | set(b)
    elif op in ['intersection', 'int', 'i', 'I', 'and', '&']:
        op = 'intersection'
        collection = set(a) & set(b)
    elif op in ['symmetric difference', 'sym', 's', 'xor', '^']:
        op = 'symmetric difference'
        collection = set(a) ^ set(b)
    elif op in ['complement', 'difference', 'diff', 'd', 'subtract', 'not', '-', '!']:
        op = 'complement'
        collection = set(a) - set(b)
    else:
        raise ValueError("operation '%s' is not a recognised option!"%op)
    if isinstance(a, dict) and isinstance(b, dict):
        dict_collection = {}
        # Set default value by method
        if method == 'number':
            default = 0
        elif method == 'list':
            default = []
        else:
            raise ValueError("method %s is not a recognised option!"%method)
        # Compute value sum for the sets intersection
        if op in ['union', 'intersection']:
            for key in collection:
                dict_collection[key] = a.get(key, default) + b.get(key, default)
        elif op == 'complement':
            for key in collection:
                dict_collection[key] = a.get(key, default)
        elif op == 'symmetric difference':
            for key in collection:
                if key in a:
                    dict_collection[key] = a.get(key, default)
                elif key in b:
                    dict_collection[key] = b.get(key, default)
                else:
                    raise ValueError('Impossible case: key not in a nor b!')
        return dict_collection
    else:
        return collection


def compute_consensus_sequences(kmers, reference, kmer_size=20,
                                charspace='nATGC', name_prefix='consensus',
                                buffer=False):
    ''' Method for computing scaffolds and contigs from a kmers object using a
    reference sequence as template
    '''
    if not 'work_dir' in globals():
        # Get directories
        global work_dir, ref_dir, result_dir
        work_dir, ref_dir, result_dir = get_directories(get_ref_dir=True)

    scaffolds_file = f'{result_dir}{name_prefix}.scaffolds.fa'
    contigs_file = f'{result_dir}{name_prefix}.contigs.fa'
    auxiliary_file= f'{result_dir}{name_prefix}.aux.tsv'
    dissected_scafs_file = f'{result_dir}{name_prefix}.disscafs.fa'
    clen = len(charspace)
    max_ns = settings['ucs']['max_ns']

    # Compute scaffold length
    scaffold_lengths = dict((name, (len(seq), i))
                            for i, (seq, name, desc) in enumerate(seqs_from_file(reference, use_ram_buffer=buffer)))

    # Divide k-mers in scaffolds
    kmer_dict = {}
    for kmer in kmers:
        for contig, cigar, strand, position in kmers[kmer]:
            if not contig in kmer_dict: kmer_dict[contig] = {}
            if not kmer in kmer_dict[contig]: kmer_dict[contig][kmer] = []
            kmer_dict[contig][kmer].append(position)

    # Compute scaffold consensuses
    with open_(auxiliary_file, 'w') as f:
        _ = f.write('# Auxiliary details on the k-mer-based scaffolds\n')
        _ = f.write('# Base\tContig\tPosition\tDepth\tConfidence\n')
        scaffolds = []
        DNAbins = np.asarray(list(charspace))[np.newaxis].T
        for scaffold in sorted(kmer_dict.keys(), key=lambda x: scaffold_lengths[x][1]):
            # Build scaffold matrix
            mat = np.zeros([scaffold_lengths[scaffold][0],clen], dtype=int)
            for kmer in kmer_dict[scaffold]:
                for p in kmer_dict[scaffold][kmer]:
                    dlen = scaffold_lengths[scaffold][0] - p
                    if dlen > 0:
                        if dlen >= kmer_size: dlen = kmer_size
                        mat[p:p+dlen,:] += (DNAbins==list(kmer[:dlen])).T
            # Extract scaffold consensus
            scaffold_consensus = ''.join([charspace[j] for j in np.argmax(mat, 1)])

            # Update scaffold object
            scaffolds.append((scaffold_consensus, scaffold, ""))

            # Compute Depth and Confidence
            depths = mat.sum(1)
            confidences = mat.max(1) / kmer_size

            # Store K-mer depths and confidence
            for j, char in enumerate(scaffold_consensus):
                _ = f.write("%s\t%s\t%s\t%s\t%s\n"%(char, scaffold, j, depths[j],
                                                    confidences[j]))

    # Store scaffolds as fasta (position annotation is 0-indexed)
    save_as_fasta(scaffolds, scaffolds_file)

    # Split scaffolds into contigs
    contigs = split_scaffolds_into_contigs(scaffolds, charspace[0])

    # Store contigs as fasta (position annotation is 0-indexed)
    save_as_fasta(contigs, contigs_file)

    # Split scaffolds into dissected scaffolds
    dissected_scaffolds = sorted((
        (z,"%s_%s"%(se[1], se[0].index(z)),'') for z, se in
        ((y.strip(charspace[0]), se) for x, se in
         map(lambda se: (se[0].split(charspace[0]*max_ns),se), scaffolds) for y in x if y)
        ), key=lambda x: len(x[0])-x[0].count(charspace[0]), reverse=True)

    # Store contigs as fasta (position annotation is 0-indexed)
    save_as_fasta(dissected_scaffolds, dissected_scafs_file)

    return (scaffolds_file, contigs_file, auxiliary_file, dissected_scafs_file)


def find_complementing_kmers(kmers, files, kmer_size=20, min_seq_len=100):
    ''' Compute the complementing kmers of a kmers object to a set of files
    containing sequences

    '''
    if log is not None:
        log.progress.add('pan',
                         'Computing k-mer complement to the negative genomes',
                         'ucs')

    # Loop through set B files
    for i, genome in enumerate(files):
        # Extract K-mers from the B file
        to_upper = settings['input']['to_upper']
        buffer = settings['input']['use_ram_buffer']
        negative_kmers = extract_kmers_from_file(genome, i, kmer_size, '',
                                                 settings['ucs']['min_kmer_count'],
                                                 settings['ucs']['rev_comp'],
                                                 min_seq_len, to_upper, buffer, 'pan')
        # TODO: Filter low counts?

        # Remove K-mers from A which are found in B
        if log is not None:
            log.progress.add(f'comp_{os.path.basename(genome)}', 'Filtering the negative k-mers found in {os.path.basename(genome)}','pan')
        kmers = set_op('complement', kmers, negative_kmers, method='list')
        if log is not None:
            log.progress[f'comp_{os.path.basename(genome)}'].log_time()

    if log is not None:
        log.progress['pan'].log_time()
        log.stats.add_row('kmer', [len(kmers), 'Computation of k-mer complement'])

    return kmers
