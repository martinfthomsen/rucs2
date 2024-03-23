# -*- coding: utf-8 -*-
"""
    Method library
    ~~~~~~~~~~~~~~~~~~~~~

    Contains the primary methods.

    * main - RUCS - a tool for designing PCR primer pairs suitable for distinguishing closely related strains
    * find_ucs - This script computes the core sequences of the positive genomes and removes sequences covered by any of the negative genomes, thus creating unique core sequences.
    * find_primer_pairs - This methods allows you to run Primer Identification for a one or more fasta entries
    * annotate - BLAST the fasta sequences to the annotation database
    * virtual_pcr - Run a virtual PCR / simulation of a PCR in silico
    * pcr_stats - Show PCR Statistics
    * show_pcr_stats - PCR Stat Analysis
    * explore - This script computes the over- and underrepresentation of k-mers in the positive genomes versus the negative genomes. This feature will create two output fasta files; the first containing the overrepresented sequences, and the second containing the underrepresented sequences. The fasta description header will contain the genome names where the given sequence is found.
    * show_primer_probe_locs - Show the location of primer bindingsites (red) and probebinding sites
    (green) using shell text coloring.
    * compute_tm - Compute the Tm of each of the sequences, and the hetero dimer complex.
    * generate_sequence_atlas_data_file - Generate Sequence Atlas data file (json)
    * print_multiranges - Markup the template with data from a multirange object with up to 4 targets.


"""
# (c) 2023 Martin Thomsen

import sys, os
import json, pickle              # refactor away

from rucs import settings, log, load_global_settings, setup_directories, text_table, round_sig
from rucs.log import log
from rucs.file import seqs_from_file, save_as_fasta, create_symbolic_files, open_, clean_up
from rucs.explore import explore_representation
from rucs.blast import get_blast_annotations
from rucs.seq import find_unique_core_sequences, analyse_genome, reverse_complement
from rucs.pcr import print_pcr_results, present_pairs_full, predict_pcr_results, configure_p3_thermoanalysis, get_primer_and_probe_bindsites, print_primer_and_probe_bindsites, color_seq, multirange, find_validated_primer_pairs

def main(positives, negatives, ref_input=None, kmer_size=None, quiet=False,
         clean_files=settings["clean_files"], annotate=True, settings_file=None, name=None,
         reuse=False):
    ''' RUCS - a tool for designing PCR primer pairs suitable for
    distinguishing closely related strains

    Part 1: This script computes the core sequences of the positive genomes and
    removes sequences covered by any of the negative genomes, thus creating
    unique core sequences.
    Part 2: The unique core sequences are then processed for potential primer
    pairs
    Part 3: The primer pair candidates are then validated by mapping the pairs
    against, first the positive genomes, to ensure that the target PCR is found,
    and that only that target is found. Secondly the candidates are mapped to
    the negative genomes to ensure that there are no false positives.

    The final primer pair candidates are sorted according to their suitability
    for PCR.
    '''
    # Set Globals
    global log
    log.set_quiet(quiet)
    if settings_file is not None:
       global settings
       settings = load_global_settings(settings_file)

    if kmer_size is None: kmer_size = settings['ucs']['kmer_size']

    # Validate Input
    if positives is None or not positives:
       raise UserWarning('No Positive genomes provided!')
    if negatives is None:
       negatives = []

    # Create run subdirectories to store temporary, result files, etc.
    work_dir, ref_dir, result_dir = setup_directories(get_ref_dir=True)

    # Define result file paths
    name = ''  if name is None else f'{name}_'
    stats_file = f'{result_dir}{name}stats.log'
    pairs_json_file = f'{result_dir}{name}pairs.json'
    products_file = f'{result_dir}{name}products.tsv'
    full_products_file = f'{result_dir}{name}products_full.tsv'
    results_file = f'{result_dir}{name}results.tsv'
    results_file_best = f'{result_dir}{name}results_best.tsv'

    try:
        log.progress.add('main', 'Running RUCS', None)

        # Save indexed Reference
        reference = f"{work_dir}reference.fa"
        if ref_input is None: ref_input = positives[0]

        log.progress.add('ref',
                         'Prepare reference: %s'%os.path.basename(ref_input),
                         'main')
        to_upper = settings['input']['to_upper']
        buffer = settings['input']['use_ram_buffer']
        save_as_fasta([seq for seq, n, d in seqs_from_file(ref_input,
                                                           to_upper=to_upper,
                                                           use_ram_buffer=buffer)],
                     reference)
        log.progress['ref'].log_time()

        # Create symlinks for all reference
        log.progress.add('input',
                        'Prepare inputs: %s positive and %s negative genomes'%(
                        len(positives), len(negatives)),
                        'main')
        reference = create_symbolic_files([reference], ref_dir, reuse=True)[0]
        positives = create_symbolic_files(positives, ref_dir, reuse=True)
        negatives = create_symbolic_files(negatives, ref_dir, reuse=True)
        log.progress['input'].log_time()

        # Find unique core sequences
        cs_f, ucs_files = find_unique_core_sequences(positives, negatives,
                                                     reference, kmer_size)

        # Identify Primer Pairs
        contig_file = ucs_files[3]

        pairs, good_pp = find_validated_primer_pairs(contig_file, positives,
                                                     negatives,
                                                     annotate=annotate,
                                                     reuse=reuse)

        # Check how many primer pairs have been tested (only tested pp in results)
        no_tested_pp = 0
        for p in pairs:
            if 'test' in p:
                no_tested_pp += 1

        # Check if any pairs were found
        if len(pairs) == 0:
            raise UserWarning('No valid PCR primer pairs could be found!')

        # Create json file containing all pairs
        with open_(pairs_json_file, 'w') as f:
            json.dump(pairs[:no_tested_pp], f)

        # Create test summary file
        with open_(full_products_file, 'w') as f:
            rp = range(len(positives))
            rn = range(len(negatives))
            f.write('#%s\t%s\n'%('\t'.join(('P_%s'%(x+1) for x in rp)),
                                 '\t'.join(('N_%s'%(x+1) for x in rn))))
            f.write('\n'.join(
                ('\t'.join(  # old version, kept for later debugging
                    [','.join(str(y[-1]) for y in x) for x in p['products']['pos']]
                    + [','.join(str(y[-1]) for y in x) for x in p['products']['neg']]
                    if 'test' in p else []
                    ) for p in pairs[:no_tested_pp])
                ))

        # Prepare reference and pcr_products list
        refs = positives + negatives
        lp = len(positives)
        pcr_products = []
        for p, pair in enumerate(good_pp):
            pcr_products.append([[] for r in refs])
            for r, ref in enumerate(positives):
                pcr_products[p][r] = pair['products']['pos'][r]
            for r, ref in enumerate(negatives):
                pcr_products[p][lp+r] = pair['products']['neg'][r]

        # Pickle dump of products
        with open_(f'{work_dir}pcr_products.pkl', 'wb') as f:
            pickle.dump(pcr_products, f)

        # Print PCR results to products file
        print_pcr_results(refs, pcr_products, products_file)

        # Create tab separated summary file of good pairs
        with open_(results_file_best, 'w') as f:
            f.write(present_pairs_full(good_pp))
        # Create tab separated summary file of results
        with open_(results_file, 'w') as f:
            f.write(present_pairs_full(pairs[:no_tested_pp]))
    except UserWarning as msg:
        # Clean up (to reduce space usage)
        clean_up(ref_dir, clean_files)
        sys.stderr.write("%s\n"%msg)
    except:
        # Clean up (to reduce space usage)
        log.progress['main'].log_time()
        log.progress.summary()
        log.stats.summary()
        clean_up(ref_dir, clean_files)
        raise
    else:
        # Clean up (to reduce space usage)
        clean_up(ref_dir, clean_files)
    finally:
        log.progress['main'].log_time()
        # Store time and stats in debug.log
        with open_(stats_file, 'w') as f:
            log.progress.summary(f)
            log.stats.summary(f)


def find_ucs(positives, negatives, ref_input=None, kmer_size=None, quiet=False,
             clean_files=settings["clean_files"], settings_file=None, name=None):
    ''' This script computes the core sequences of the positive genomes and
    removes sequences covered by any of the negative genomes, thus creating
    unique core sequences.
    '''
    # Set Globals
    global log
    log.set_quiet(quiet)
    if settings_file is not None:
        global settings
        settings = load_global_settings(settings_file)

    # Validate Input
    if positives is None or not positives:
        raise UserWarning('No Positive genomes provided!')
    if negatives is None:
        negatives = []

    # Create run subdirectories to store temporary, result files, etc.
    work_dir, ref_dir, result_dir = setup_directories(get_ref_dir=True)

    stats_file = '%sstats.log'%("%s_"%name if name is not None else '')

    if kmer_size is None:
        kmer_size = settings['ucs']['kmer_size']

    p3_args = settings['pcr']['priming']['primer3']
    if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'], list):
        if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0], list):
            product_length_limits = [lim for lims in p3_args['PRIMER_PRODUCT_SIZE_RANGE'] for lim in lims]
            min_seq_len = min(product_length_limits)
        else:
            min_seq_len = p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0]
    else:
        exit('Settings error: invalid input format for PRIMER_PRODUCT_SIZE_RANGE, expected "list", got: ' + p3_args['PRIMER_PRODUCT_SIZE_RANGE'])

    try:
        log.progress.add('main', 'Running find_ucs', None)

        # Save Sorted Reference
        reference = "reference.fa"
        if ref_input is None: ref_input = positives[0]

        log.progress.add('ref',
                         'Prepare reference: %s'%os.path.basename(ref_input),
                         'main')
        to_upper = settings['input']['to_upper']
        buffer = settings['input']['use_ram_buffer']
        save_as_fasta([seq for seq, n, d in seqs_from_file(ref_input,
                                                           to_upper=to_upper,
                                                           use_ram_buffer=buffer)],
                      reference)
        log.progress['ref'].log_time()

        # Create symlinks for all reference
        log.progress.add('input',
                         'Prepare inputs: %s positive and %s negative genomes'%(
                         len(positives), len(negatives)),
                         'main')
        reference = create_symbolic_files([reference], ref_dir, reuse=True)[0]
        positives = create_symbolic_files(positives, ref_dir, reuse=True)
        negatives = create_symbolic_files(negatives, ref_dir, reuse=True)
        log.progress['input'].log_time()

        # Find unique core sequences
        cs_files, ucs_files = find_unique_core_sequences(positives, negatives,
                                                         reference, kmer_size)

        # Print Sequence Analysis Table
        title = 'Sequence Analysis'
        headers = ['', 'Sequences', 'Size in bases', 'Seqs >%s'%min_seq_len,
                   'Size >%s'%min_seq_len]
        rows = [['Reference']+analyse_genome(reference, min_seq_len),
                ['Core Sequences']+analyse_genome(cs_files[1], min_seq_len),
                ['Unique Core Sequences']+analyse_genome(ucs_files[1], min_seq_len)]
        print(text_table(title, headers, rows))

    except UserWarning as msg:
        # Clean up (to reduce space usage)
        clean_up(ref_dir, clean_files)
        sys.stderr.write("%s\n"%msg)
    except:
        # Clean up (to reduce space usage)
        log.progress['main'].log_time()
        log.progress.summary()
        log.stats.summary()
        clean_up(ref_dir, clean_files)
        raise
    else:
        # Clean up (to reduce space usage)
        clean_up(ref_dir, clean_files)
    finally:
        log.progress['main'].log_time()
        # Store time and stats in stats.log
        with open_(stats_file, 'w') as f:
            log.progress.summary(f)
            log.stats.summary(f)


def find_primer_pairs(contig_file, positives, negatives, contig_names=None,
                      quiet=False, clean_files=settings["clean_files"], annotate=True,
                      settings_file=None, name=None, reuse=False):
    ''' This methods allows you to run Primer Identification for a one or more
        fasta entries '''
    # Set Globals
    global log
    log.set_quiet(quiet)
    if settings_file is not None:
        global settings
        settings = load_global_settings(settings_file)

    # Create run subdirectories to store temporary, result files, etc.
    work_dir, ref_dir, result_dir = setup_directories(get_ref_dir=True)

    # Define result file paths
    name = ''  if name is None else f'{name}_'
    stats_file = f'{result_dir}{name}stats.log'
    pairs_json_file = f'{result_dir}{name}pairs.json'
    products_file = f'{result_dir}{name}products.tsv'
    full_products_file = f'{result_dir}{name}products_full.tsv'
    results_file = f'{result_dir}{name}results.tsv'
    results_file_best = f'{result_dir}{name}results_best.tsv'

    try:
        log.progress.add('main', 'Running RUCS', None)

        # Create symlinks for all reference
        log.progress.add('input',
                         'Prepare inputs: %s positive and %s negative genomes'%(
                         len(positives), len(negatives)),
                         'main')
        positives = create_symbolic_files(positives, ref_dir, reuse=True)
        negatives = create_symbolic_files(negatives, ref_dir, reuse=True)
        log.progress['input'].log_time()

        if contig_names is not None and isinstance(contig_names, str):
            contig_names = [contig_names]
        pairs, good_pp = find_validated_primer_pairs(contig_file, positives,
                                                     negatives,
                                                     contig_names=contig_names,
                                                     annotate=annotate,
                                                     reuse=reuse)

        # Check how many primer pairs have been tested (only tested pp in results)
        no_tested_pp = 0
        for p in pairs:
            if 'test' in p:
                no_tested_pp += 1

        # Check if any pairs were found
        if len(pairs) == 0:
            raise UserWarning('No valid PCR primer pairs could be found!')

        # Create json file containing all pairs
        with open_(pairs_json_file, 'w') as f:
            json.dump(pairs[:no_tested_pp], f)

        # Create test summary file
        with open_(full_products_file, 'w') as f:
            rp = range(len(positives))
            rn = range(len(negatives))
            f.write('#%s\t%s\n'%('\t'.join(('P_%s'%(x+1) for x in rp)),
                                 '\t'.join(('N_%s'%(x+1) for x in rn))))
            f.write('\n'.join(
                ('\t'.join(  # old version, kept for later debugging
                    [','.join(str(y[-1]) for y in x) for x in p['products']['pos']]
                    + [','.join(str(y[-1]) for y in x) for x in p['products']['neg']]
                    if 'test' in p else []
                    ) for p in pairs[:no_tested_pp])
                ))

        # Prepare reference and pcr_products list
        refs = positives + negatives
        lp = len(positives)
        pcr_products = []
        for p, pair in enumerate(good_pp):
            pcr_products.append([[] for r in refs])
            for r, ref in enumerate(positives):
                pcr_products[p][r] = pair['products']['pos'][r]
            for r, ref in enumerate(negatives):
                pcr_products[p][lp+r] = pair['products']['neg'][r]

        # Pickle dump of products
        with open_(f'{work_dir}pcr_products.pkl', 'wb') as f:
            pickle.dump(pcr_products, f)

        # Print PCR results to products file
        print_pcr_results(refs, pcr_products, products_file)

        # Create tab separated summary file of good pairs
        with open_(results_file_best, 'w') as f:
            f.write(present_pairs_full(good_pp))
        # Create tab separated summary file of results
        with open_(results_file, 'w') as f:
            f.write(present_pairs_full(pairs[:no_tested_pp]))
    except UserWarning as msg:
        # Clean up (to reduce space usage)
        clean_up(ref_dir, clean_files)
        sys.stderr.write("%s\n"%msg)
    except:
        # Clean up (to reduce space usage)
        log.progress['main'].log_time()
        log.progress.summary()
        log.stats.summary()
        clean_up(ref_dir, clean_files)
        raise
    else:
        # Clean up (to reduce space usage)
        clean_up(ref_dir, clean_files)
    finally:
        log.progress['main'].log_time()
        # Store time and stats in debug.log
        with open_(stats_file, 'w') as f:
            log.progress.summary(f)
            log.stats.summary(f)


def annotate(fasta, blast_settings=None, dbpath=''):
    ''' BLAST the fasta sequences to the annotation database (default db is the swissprot protein database)

    USAGE
        >>> from rucs.blast import get_blast_annotations
        >>> annotations = annotate('unique_core_sequences.contigs.fa', dbpath='/path/to/blastdb/')
    '''
    # Create run subdirectories to store temporary, result files, etc.
    work_dir, ref_dir, result_dir = setup_directories(get_ref_dir=True)

    return get_blast_annotations(fasta, blast_settings, dbpath)


def virtual_pcr(references, pairs, output=None, name=None, tm_thresholds=None, min_grade=None):
    ''' Virtual PCR

    USAGE
        >>> refs = ['/full/path/to/ref1.fa', ...]
        >>> pairs = [['forward_primer1_seq', 'reverse_primer1_seq'], ...]
        >>> virtual_pcr(refs, pairs, output='my_pcr_results.tsv')
    '''
    # Create run subdirectories to store temporary, result files, etc.
    work_dir, ref_dir, result_dir = setup_directories(get_ref_dir=True)

    # Define
    if output is None or not isinstance(output, str):
        name = '' if name is None else f'{name}_'
        output = f'{result_dir}{name}products.tsv'

    # Create symlinks for all reference
    refs = create_symbolic_files(references, ref_dir, reuse=True)

    fail_on_non_match=False
    tm_thresholds = settings['pcr']['priming']['tm_thresholds'] if tm_thresholds is None else tm_thresholds
    min_grade = settings['pcr']['priming']['threshold_grade'] if min_grade is None else min_grade
    pcr_products = predict_pcr_results(refs, pairs, fail_on_non_match, tm_thresholds, min_grade)

    # Pickle dump of products
    with open_(f'{work_dir}pcr_products.pkl', 'wb') as products_file:
        pickle.dump(pcr_products, products_file)

    # Print PCR results
    print_pcr_results(refs, pcr_products, output)


def pcr_stats(pairs, template=None, to_upper=True, buffer=False):
    ''' Show PCR Statistics

    Analyses scope:
     - Product/amplicon size
     - Primer and probe
       - Length
       - GC-content
       - Melting temperature
       - Homo dimer temperature
       - Hairpin temperature
     - Hetero dimer temperature between primers and probe
     - 5' ends on G or C check
     - Probe distance to primer

    Template is required for:
     - estimating product size
     - check probe distance to primers
    '''
    # Get first template entry (TODO: Check all template sequences for primer matches)
    if template is not None:
        try:
            for seq, n, d in seqs_from_file(template, to_upper=to_upper,
                                            use_ram_buffer=buffer):
                template = seq
                break
        except:
            sys.stderr.write('Template: %s\n'%template)
            sys.stderr.write('Warning: Template not valid!\n')
            template = None

    # Show PCR stats for each pair
    for i, seqs in enumerate(pairs):
        forward = seqs[0]
        reverse = seqs[1]
        probe   = seqs[2] if len(seqs) == 3 else None
        show_pcr_stats(forward, reverse, probe, template,
                       title='PCR Stat Analysis for pair %s'%(i+1))


def show_pcr_stats(forward, reverse, probe=None, template=None,
                   settings_file=None, title='PCR Stat Analysis'):
    '''
    USAGE:
        >>> settings_file = 'settings.default.cjson'
        >>> forward = 'TGAACATACACGGCACAGA'
        >>> reverse = 'CGTAAGCCACGCCTAGT'
        >>> probe = 'TGTTCGTCGTCGGTGAGACGGC'
        >>> show_pcr_stats(forward, reverse, probe, settings_file=settings_file)

PCR Stat Analysis

PCR product size: N/A

         Sequence
Forward: TGAACATACACGGCACAGA
Reverse: CGTAAGCCACGCCTAGT
Probe:   TGTTCGTCGTCGGTGAGACGGC

             Forward  Reverse  Probe
Length:         19       17     22
GC% content:    47.4     58.8   63.6
Melting Tm:     56.7     56.2   67.1
Homo dimer:      -        -      -
Hairpin:         -        -     62.5

              Fw-Rv  Fw-Pr  Rv-Pr
Hetero dimer:   -      -      -

Probe 5' end not G or C?  PASSED
Probe distance to primer: N/A
    '''
    def neg2mask(num):
        ''' Change all negative values to "-  " '''
        return num if num > 0 else '-  '

    if settings_file is not None:
        global settings
        settings = load_global_settings(settings_file)

    # Initiate the thermodynamic analyses (defining: p3_primer, p3_probe)
    p3_primer, p3_probe = configure_p3_thermoanalysis()

    prim_sets = settings['pcr']['priming']
    threshold_tm = prim_sets['threshold_tm'] if 'threshold_tm' in prim_sets else 47
    max_probe_dist = prim_sets['max_probe_dist'] if 'max_probe_dist' in prim_sets else 30
    min_primer_tm = prim_sets['primer3']['PRIMER_MIN_TM'] if 'PRIMER_MIN_TM' in prim_sets['primer3'] else 57
    max_primer_tm = prim_sets['primer3']['PRIMER_MAX_TM'] if 'PRIMER_MAX_TM' in prim_sets['primer3'] else 62
    min_probe_tm = prim_sets['primer3']['PRIMER_INTERNAL_MIN_TM'] if 'PRIMER_INTERNAL_MIN_TM' in prim_sets['primer3'] else 67
    max_probe_tm = prim_sets['primer3']['PRIMER_INTERNAL_MAX_TM'] if 'PRIMER_INTERNAL_MAX_TM' in prim_sets['primer3'] else 72

    n = "\x1B[0m"
    h = "\x1B[1;4m"
    r = "\x1B[31m"
    g = "\x1B[32m"
    print("\n%s%s%s\n"%(h, title, n))
    # Check PCR product size
    if template is None:
        size = 'N/A'
    else:
        try:
            start = template.index(forward)
            end = template.index(reverse_complement(reverse)) + len(reverse)
            size = "%sbp"%(end - start)
        except ValueError as e:
            try:
                tmp = forward
                forward = reverse
                reverse = tmp
                start = template.index(forward)
                end = template.index(reverse_complement(reverse)) + len(reverse)
                size = "%sbp"%(end - start)
            except ValueError as e:
                exit('Error: Unable to map forward or reverse primer to template!')
    print("PCR product size: %s\n"%(size))

    # Print sequences
    print("         Sequence")
    print("Forward: %s\nReverse: %s\nProbe:   %s\n"%(forward, reverse, probe))

    print("             Forward  Reverse  Probe")
    # Print length
    l1, l2 = len(forward), len(reverse)
    l3 = len(probe) if probe is not None else '  N/A'
    print("Length:      %5s    %5s    %3s  "%(l1, l2, l3))

    # Count GC%
    gc1 = round((forward.count('G') + forward.count('C')) / l1 * 100, 1)
    gc2 = round((reverse.count('G') + reverse.count('C')) / l2 * 100, 1)
    if probe is not None:
        gc3 = round((probe.count('G') + probe.count('C')) / l3 * 100, 1)
    else:
        gc3 = 'N/A'
    print("GC%% content: %7s  %7s  %5s"%(gc1, gc2, gc3))

    # Melting temperature
    tm1 = neg2mask(round(p3_primer.calc_tm(forward),1))
    tm2 = neg2mask(round(p3_primer.calc_tm(reverse),1))
    tm3 = neg2mask(round(p3_primer.calc_tm(probe),1)) if probe is not None else 'N/A'
    m1, m2, m3 = n, n, n
    if isinstance(tm1, float): m1 = g if tm1 > min_primer_tm and tm1 < max_primer_tm else r
    if isinstance(tm2, float): m2 = g if tm2 > min_primer_tm and tm2 < max_primer_tm else r
    if isinstance(tm3, float): m3 = g if tm3 > min_probe_tm and tm3 < max_probe_tm else r
    print("Melting Tm:  %s%7s  %s%7s  %s%5s%s"%(m1, tm1, m2, tm2, m3, tm3, n))

    # Homo dimer temperature
    tm1 = neg2mask(round(p3_primer.calc_homodimer(forward).tm,1))
    tm2 = neg2mask(round(p3_primer.calc_homodimer(reverse).tm,1))
    tm3 = neg2mask(round(p3_primer.calc_homodimer(probe).tm,1)) if probe is not None else 'N/A'
    m1, m2, m3 = n, n, n
    if isinstance(tm1, float): m1 = r if tm1 > threshold_tm else g
    if isinstance(tm2, float): m2 = r if tm2 > threshold_tm else g
    if isinstance(tm3, float): m3 = r if tm3 > threshold_tm else g
    print("Homo dimer:  %s%7s  %s%7s  %s%5s%s"%(m1, tm1, m2, tm2, m3, tm3, n))

    # Hairpin temperature
    tm1 = neg2mask(round(p3_primer.calc_hairpin(forward).tm,1))
    tm2 = neg2mask(round(p3_primer.calc_hairpin(reverse).tm,1))
    tm3 = neg2mask(round(p3_primer.calc_hairpin(probe).tm,1)) if probe is not None else 'N/A'
    m1, m2, m3 = n, n, n
    if isinstance(tm1, float): m1 = r if tm1 > threshold_tm else g
    if isinstance(tm2, float): m2 = r if tm2 > threshold_tm else g
    if isinstance(tm3, float): m3 = r if tm3 > threshold_tm else g
    print("Hairpin:     %s%7s  %s%7s  %s%5s%s\n"%(m1, tm1, m2, tm2, m3, tm3, n))

    print("              Fw-Rv  Fw-Pr  Rv-Pr")
    # Hetero dimer temperatur
    tm1 = neg2mask(round(p3_primer.calc_heterodimer(forward, reverse).tm,1))
    tm2 = neg2mask(round(p3_primer.calc_heterodimer(forward, probe).tm,1)) if probe is not None else 'N/A'
    tm3 = neg2mask(round(p3_primer.calc_heterodimer(reverse, probe).tm,1)) if probe is not None else 'N/A'
    m1, m2, m3 = n, n, n
    if isinstance(tm1, float): m1 = r if tm1 > threshold_tm else g
    if isinstance(tm2, float): m2 = r if tm2 > threshold_tm else g
    if isinstance(tm3, float): m3 = r if tm3 > threshold_tm else g
    print("Hetero dimer: %s%5s  %s%5s  %s%5s%s\n"%(m1, tm1, m2, tm2, m3, tm3, n))

    # Check Probe not starting with GC
    if probe is not None:
        state = "%sFAILED"%r if probe[0] in 'GC' else "%sPASSED"%g
        print("Probe 5' end not G or C?  %s%s"%(state, n))

        # Check Probe distance to primer
        if template is None:
            dist = 'N/A'
            m = n
        else:
            if probe in template:
                start = template.index(forward) + len(forward)
                end = template.index(probe)
            else:
                rctem = reverse_complement(template)
                start = rctem.index(reverse) + len(reverse)
                end = rctem.index(probe)
            dist = "%sbp"%(end - start)
            m = r if end - start > max_probe_dist else g
        print("Probe distance to primer: %s%s%s\n"%(m, dist, n))


def explore(positives, negatives, kmer_size=None, quiet=False, clean_files=settings["clean_files"],
            settings_file=None, name=None, reuse=False):
    ''' This script computes the over- and underrepresentation of k-mers in the
    positive genomes versus the negative genomes. This feature will create two
    output fasta files; the first containing the overrepresented sequences, and
    the second containing the underrepresented sequences. The fasta description
    header will contain the genome names where the given sequence is found.
    '''
    # Set Globals
    global log
    log.set_quiet(quiet)
    if settings_file is not None:
        global settings
        settings = load_global_settings(settings_file)

    # Validate Input
    if positives is None or not positives:
        raise UserWarning('No Positive genomes provided!')
    if negatives is None:
        negatives = []

    stats_file = '%sstats.log'%("%s_"%name if name is not None else '')

    if kmer_size is None:
        kmer_size = settings['ucs']['kmer_size']

    p3_args = settings['pcr']['priming']['primer3']
    if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'], list):
        if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0], list):
            product_length_limits = [lim for lims in p3_args['PRIMER_PRODUCT_SIZE_RANGE'] for lim in lims]
            min_seq_len = min(product_length_limits)
        else:
            min_seq_len = p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0]
    else:
        exit('Settings error: invalid input format for PRIMER_PRODUCT_SIZE_RANGE, expected "list", got: ' + p3_args['PRIMER_PRODUCT_SIZE_RANGE'])

    # Create run subdirectories to store temporary, result files, etc.
    work_dir, ref_dir, result_dir = setup_directories(get_ref_dir=True)
    try:
        log.progress.add('main', 'Running exploration', None)

        # Create symlinks for all reference
        log.progress.add('input',
                         'Prepare inputs: %s positive and %s negative genomes'%(
                         len(positives), len(negatives)),
                         'main')
        positives = create_symbolic_files(positives, ref_dir, reuse=True)
        negatives = create_symbolic_files(negatives, ref_dir, reuse=True)
        log.progress['input'].log_time()

        # Find over- and under-represented sequences
        ors_files, urs_files = explore_representation(positives, negatives, kmer_size, reuse=reuse)

        # Print Sequence Analysis Table
        title = 'Sequence Analysis'
        headers = ['', 'Sequences', 'Size in bases', 'Seqs >%s'%min_seq_len,
                   'Size >%s'%min_seq_len]
        rows = [['Over Represented Sequences']+analyse_genome(ors_files[1], min_seq_len),
                ['Under Represented Sequences']+analyse_genome(urs_files[1], min_seq_len)]
        print(text_table(title, headers, rows))

    except UserWarning as msg:
        # Clean up (to reduce space usage)
        clean_up(ref_dir, clean_files)
        sys.stderr.write("%s\n"%msg)
    except:
        # Clean up (to reduce space usage)
        log.progress['main'].log_time()
        log.progress.summary()
        log.stats.summary()
        clean_up(ref_dir, clean_files)
        raise
    else:
        # Clean up (to reduce space usage)
        clean_up(ref_dir, clean_files)
    finally:
        log.progress['main'].log_time()
        # Store time and stats in stats.log
        with open_(stats_file, 'w') as f:
            log.progress.summary(f)
            log.stats.summary(f)


def show_primer_probe_locs(contigs=None, name=None):
    ''' Show the location of primer bindingsites (red) and probebinding sites
    (green) using shell text coloring.
    '''
    # Create run subdirectories to store temporary, result files, etc.
    work_dir, ref_dir, result_dir = setup_directories(get_ref_dir=True)

    # Define result file paths
    name = '' if name is None else f'{name}_'
    # stats_file = f'{result_dir}{name}stats.log'
    # pairs_json_file = f'{result_dir}{name}pairs.json'
    products_file = f'{result_dir}{name}products.tsv'
    results_file = f'{result_dir}{name}results.tsv'
    # results_file_best = f'{result_dir}{name}results_best.tsv'
    seq_file    = f'{work_dir}{name}unique_core_sequences.disscafs.fa'
    pribind_locs, probind_locs = get_primer_and_probe_bindsites(results_file)
    print_primer_and_probe_bindsites(seq_file, pribind_locs, probind_locs, contigs=contigs, tag1='41', tag2='32')


def compute_tm(seq1, seq2=None):
    ''' Compute the Tm of each of the sequences, and the hetero dimer complex.

    USAGE:
        >>> compute_tm('CAACATTTTCGTGTCGCCCTT')
        (61.9, 61.9, 61.8)
        >>> compute_tm('CAACATTTTCGTGTCGCCCTT', 'AAGGGCGACACGAGAATGTTG')
        (61.9, 63.1, 51.6)
    '''
    # Initiate the thermodynamic analyses (defining: p3_primer, p3_probe)
    p3_primer, p3_probe = configure_p3_thermoanalysis()

    tm_seq1 = round_sig(p3_primer.calc_tm(seq1))
    if seq2 is not None:
        tm_seq2 = round_sig(p3_primer.calc_tm(seq2))
        tm_dimer = round_sig(p3_primer.calc_heterodimer(seq1, seq2).tm)
    else:
        tm_seq2 = tm_seq1
        tm_dimer = round_sig(p3_primer.calc_heterodimer(seq1,
                                                       reverse_complement(seq1)).tm)

    return tm_seq1, tm_seq2, tm_dimer


def generate_sequence_atlas_data_file(reference_file, core_sequence_file, unique_core_sequence_file, aux_file, file_name='results.json'):
    ''' Generate Sequence Atlas data file (json)

    >>> generate_sequence_atlas_data_file(
    ...     'reference.fa',
    ...     'core_sequences.contigs.fa',
    ...     'unique_core_sequences.contigs.fa',
    ...     'unique_core_sequences.aux.tsv')
    '''
    buffer = settings['input']['use_ram_buffer']
    # read reference   -> contigs: name(id), length
    data = {'lanes': ['reference', 'cs', 'ucs', 'depth', 'confidence'],
            'attributes': ['name', 'length'], 'data': []}
    for i, (seq, name, desc) in enumerate(seqs_from_file(reference_file, use_ram_buffer=buffer)):
        data['data'].append({
            'attributes': {
                'name': str(i),
                'length': len(seq),
                },
            'lanes': {
                'reference': [{'name': str(i),
                               'position': 0,
                               'length': len(seq)}],
                'depth': [],
                'confidence': []
                }
            })

    # read cs          -> contigs: length, position, name
    css = {}
    for seq, name, desc in seqs_from_file(core_sequence_file, use_ram_buffer=buffer):
        id = int(name.split('_')[0])
        if not id in css: css[id] = []
        att = dict(map(lambda x: (y.strip() for y in x.split('=')), desc.split(',')))
        try: pos = int(att.get('position'))
        except: pos = None
        css[id].append({'name': len(css[id]),
                        'position': pos,
                        'length': len(seq)})

    for id, cs in css.items():
        data['data'][id]['lanes']['cs'] = cs

    # read ucs         -> contigs: length, position, name
    ucss = {}
    for seq, name, desc in seqs_from_file(unique_core_sequence_file, use_ram_buffer=buffer):
        id = int(name.split('_')[0])
        if not id in ucss: ucss[id] = []
        att = dict(map(lambda x: (y.strip() for y in x.split('=')), desc.split(',')))
        try: pos = int(att.get('position'))
        except: pos = None
        ucss[id].append({'name': len(ucss[id]),
                         'position': pos,
                         'length': len(seq)})

    for id, ucs in ucss.items():
        data['data'][id]['lanes']['ucs'] = ucs

    # read tsv         -> depth, confidence
    with open_(aux_file) as f:
        for l in f:
            l = l.strip()
            if l == '': continue
            if l[0] == '#': continue
            # Base	Contig	Position	Depth	Confidence
            tmp = l.split('\t')
            id = int(tmp[1])
            ldata = data['data'][id]['lanes']
            ldata['depth'].append(int(tmp[3]))
            ldata['confidence'].append(float(tmp[4]))

    # Filter empty depth/confidence
    for id in range(len(data['data'])):
        adata = data['data'][id]['attributes']
        ldata = data['data'][id]['lanes']
        if adata['length'] > 0:
            if len(ldata['depth']) == 0: del ldata['depth']
            elif len(ldata['depth']) != adata['length']: print('ALERT!')

            if len(ldata['confidence']) == 0: del ldata['confidence']
            elif len(ldata['confidence']) != adata['length']: print('ALERT!')

    # Dump data to file
    with open_(file_name, 'w') as f: json.dump(data, f)


def print_multiranges(template, multiranges, entries=None, buffer=False):
    ''' Markup the template with data from a multirange object with up to 4
    targets.

    Text Modifications:
        The last target: Red
        The second last: Green
        The third last : Blue
        The fourth last: Underlined (_)
        Red   + Green  = Yellow
        Red   + Blue   = Magenta
        Green + Blue   = Cyan
        R + G + B      = Greyish/White

    Any modification above is also followed up by making the text bold.

    USAGE
        * Template must be a fasta file
        * Multiranges must be a dictionary, where the keys are the names of the
          fasta entries, and where each value is a multirange object combining up
          to 4 multirange objects through addition.
    '''
    tag_mods = [['0'],              # No Modifications
                ['0','1','31'],     # Bold Red
                ['0','1','32'],     # Bold Green
                ['0','1','33'],     # Bold Yellow
                ['0','1','34'],     # Bold Blue
                ['0','1','35'],     # Bold Magenta
                ['0','1','36'],     # Bold Cyan
                ['0','1','37'],     # Bold Whiteish/Grey
                ['0','1','4'],      # Bold Underlined
                ['0','1','4','31'], # Bold Underlined Red
                ['0','1','4','32'], # Bold Underlined Green
                ['0','1','4','33'], # Bold Underlined Yellow
                ['0','1','4','34'], # Bold Underlined Blue
                ['0','1','4','35'], # Bold Underlined Magenta
                ['0','1','4','36'], # Bold Underlined Cyan
                ['0','1','4','37']] # Bold Underlined Whiteish/Grey

    # Fetch Template data
    for seq, name, desc in seqs_from_file(template, use_ram_buffer=buffer):
        if entries is not None:
            if not name in entries:
                continue
        elif not name in multiranges: continue

        print('SEQUENCE NAME: %s'%name)
        if name in multiranges:
            try:
                print(color_seq(seq, multiranges[name], tag_mods))
            except ValueError:
                print(color_seq(seq, multirange() + multiranges[name], tag_mods))
        else:
            print(seq)
        print('\n')
