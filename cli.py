#!/usr/local/bin/python3 -u
# -*- coding: utf-8 -*-
"""
    Commandline Interface
    ~~~~~~~~~~~~~~~~~~~~~

    Main entry point for command line invocation.

"""
# (c) 2023 Martin Thomsen

import sys, os, json

from rucs import settings, load_global_settings, get_directories
from rucs.methods import main, find_ucs, find_primer_pairs, annotate, virtual_pcr, pcr_stats, explore
from rucs.file import get_pairs, get_fasta_files

# Entry points for methods
def full(args):
    ''' Run full diagnostic: fucs + fppp '''
    main(args.positives, args.negatives, args.reference, quiet=args.quiet,
         clean_files=settings["clean_files"], annotate=True, reuse=args.reuse)

def fucs(args):
    ''' Find Unique Core Sequences '''
    find_ucs(args.positives, args.negatives, args.reference, quiet=args.quiet,
             clean_files=settings["clean_files"])

def fppp(args):
    ''' Find PCR Primer Pairs '''
    settings['pcr']['seq_selection'] = None
    find_primer_pairs(args.template, args.positives, args.negatives,
                      quiet=args.quiet, clean_files=settings["clean_files"], annotate=True,
                      reuse=args.reuse)

def anno(args):
    ''' Annotate sequences using a protein BLAST DB '''
    annotations = annotate(args.template,
        settings['pcr']['annotation']['blastx_settings'],
        settings['pcr']['annotation']['blast_db_path'])
    print(annotations)

def vpcr(args):
    ''' Virtual PCR - Simulate PCR and predict PCR product for the provided
    primer pairs against the provided references.

    MINIMUM CMDLINE ARGUMENT INPUTS:
       --pairs pair_file.tsv        - This is a tab-separated file listing the forward and reverse primer sequences (and probes)
       --references references/*    - This is a list of all the reference samples to run the virtual PCR against
    '''
    pairs = get_pairs(args.pairs)
    min_grade = settings['pcr']['priming']['threshold_grade']
    virtual_pcr(args.references, pairs, output=None, name=None, tm_thresholds=args.tm_thresholds, min_grade=min_grade)

def pcrs(args):
    ''' Show PCR Statistics - Analyse the temperature and more for the provided
    primer pairs.

    MINIMUM CMDLINE ARGUMENT INPUTS:
       --pairs pair_file.tsv        - This is a tab-separated file listing the forward and reverse primer sequences (and probes)
       --template template.fa       - This is a fasta file containing the template where the primers and probes align to
    '''
    # Get pairs
    pairs = get_pairs(args.pairs)
    to_upper = settings['input']['to_upper']
    buffer = settings['input']['use_ram_buffer']
    pcr_stats(pairs, args.template, to_upper, buffer)

def expl(args):
    ''' Explore the positive and negative genomes for overrepresented k-mers'''
    explore(args.positives, args.negatives, quiet=args.quiet, clean_files=settings["clean_files"], reuse=args.reuse)

def test(args):
    ''' Test run - Run the standard test set through the main (full run)
    function to test that all the different dependencies operate as expected.
    '''
    rucs_dir = os.path.dirname(os.path.realpath(__file__))
    pos = ["%s/testdata/bla.fa"%(rucs_dir)]
    neg = ["%s/testdata/sul.fa"%(rucs_dir)]
    main(pos, neg, None, quiet=args.quiet, clean_files=None, annotate=True,
         reuse=args.reuse)

    # Get directories
    work_dir, ref_dir, result_dir = get_directories(get_ref_dir=True)
    if os.path.exists(f'{result_dir}results.tsv'):
        print('Test completed successfully!')
    else:
        print('Test completed with errors!')


if __name__ == '__main__':
    from argparse import ArgumentParser

    # PARSE COMMAND LINE ARGUMENTS
    parser = ArgumentParser()
    # Entry point
    parser.add_argument("entry_point", nargs=1)
    # Main arguments
    parser.add_argument("--positives", nargs='+', default=None,
                        help=("File paths for the positive dataset"))
    parser.add_argument("--negatives", nargs='+', default=None,
                        help=("File paths for the nagetive dataset"))
    parser.add_argument("--references", nargs='+', default=None,
                        help=("File paths for the references to be tested"))
    parser.add_argument("--reference", default=None,
                        help=("The reference file to which the k-mers should be "
                              "mapped."))
    parser.add_argument("--template", default=None,
                        help=("File path for the template file"))
    parser.add_argument("--pairs", default=None,
                        help=("File containing pcr primer pair sets (1 pair per "
                              "line, tab-separated sequences, format: forward, "
                              "reverse, probe*) *optional"))
    parser.add_argument("--settings_file", default="settings.default.cjson",
                        help=("File containing the default settings for this "
                              "program"))
    # Settings-modification options
    # The following options if passed will overwrite the default settings in the
    # settings_file
    parser.add_argument("--kmer_size", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--read_length", default=None,
                        help=("This option will modify min_seq_len_pos and "
                              "min_seq_len_neg in the settings file"))
    parser.add_argument("--tm_threshold", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--primer_pair_max_diff_tm", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--dna_conc", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--dna_conc_probe", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--salt_monovalent_conc", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--salt_divalent_conc", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--dntp_conc", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--anneal_tm", default=None,
                        help=("This will modify the appropriate values in the settings"))
    parser.add_argument("--tm_thresholds", default=None,
                        help=(f"Overwrite the settings default value of '{settings['pcr']['priming']['tm_thresholds']}'. List of Tm     thresholds for grading primer binding to binding site. each item result in one grade increase. For single     numbers the Tm must be higher than the threshold to pass, for list items, the Tm must be between the two     numbers to pass."))
    parser.add_argument("--max_3end_gc", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--product_size_min", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--product_size_max", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--pick_probe", action='store_true',
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--annotation_evalue", default=None,
                        help=("This will overwrite the set value in the settings"))
    parser.add_argument("--kmer_count_threshold", default=None,
                        help=("kmers found below this limit within each file are"
                              " ignored. Setting this argument will overwrite the"
                              " set value in the settings"))
    parser.add_argument("--sensitivity_threshold", default=None,
                        help=("The sensitivity threshold defines how often a kmer"
                              " must be found in the positive set to be included "
                              "in the results. If provided as an integer the "
                              "number is considered as a minimum count. Setting "
                              "this argument will overwrite the set value in the "
                              "settings"))
    parser.add_argument("--fallout_threshold", default=None,
                        help=("The fall-out threshold defines how often k-mers "
                              "may be found in the negative set and still be "
                              "included in the results. If provided as an integer"
                              " the number is considered as a maximum count. "
                              "Setting this argument will overwrite the set value"
                              " in the settings"))
    parser.add_argument("--align_percent_threshold", default=None,
                        help=("The alignment percent threshold defines the "
                              "acceptable amount of kmers to not be aligned to a "
                              "contig. These k-mers are lost from further "
                              "analysis to speed up the process. Set to 0, if you"
                              " want as much data as possible. Setting this "
                              "argument will overwrite the set value in the "
                              "settings"))
    # Standard arguments
    parser.add_argument("-r", "--reuse", default=False, action='store_true',
                        help=("This option allows the reuse of some work files"
                              " making subsequent rerun of the tool faster"))
    parser.add_argument("-v", "--verbose", default=False, action='store_true',
                        help=("This option write live information of the "
                              "progress to the screen"))
    args = parser.parse_args()

    entry_points = ['full', 'fucs', 'fppp', 'vpcr', 'anno', 'pcrs', 'expl', 'test']
    args.entry_point = args.entry_point[0]

    if args.entry_point not in entry_points:
        sys.stderr.write('Unknown entry point provided, for help use the -h '
                         'option!\n')
        sys.exit(1)

    # Load Default Settings
    settings = load_global_settings(args.settings_file)

    # Update Settings
    if args.kmer_size is not None:
        settings['ucs']['kmer_size'] = int(args.kmer_size)
    if args.read_length is not None:
        settings['ucs']['min_seq_len_pos'] = int(args.read_length) * 2
        settings['ucs']['min_seq_len_neg'] = int(args.read_length) * 2
    if args.tm_threshold is not None:
        tm = int(args.tm_threshold)
        settings['pcr']['priming']['threshold_tm'] = tm
        p3_settings = settings['pcr']['priming']['primer3']
        p3_settings['PRIMER_MAX_HAIRPIN_TH'] = tm
        p3_settings['PRIMER_PAIR_MAX_COMPL_ANY_TH'] = tm
        p3_settings['PRIMER_PAIR_MAX_COMPL_END_TH'] = tm
        p3_settings['PRIMER_MAX_SELF_ANY_TH'] = tm
        p3_settings['PRIMER_MAX_SELF_END_TH'] = tm
        p3_settings['PRIMER_INTERNAL_MAX_HAIRPIN_TH'] = tm
        p3_settings['PRIMER_INTERNAL_MAX_SELF_ANY_TH'] = tm
        p3_settings['PRIMER_INTERNAL_MAX_SELF_END_TH'] = tm
    if args.primer_pair_max_diff_tm is not None:
        settings['pcr']['priming']['primer3']['PRIMER_PAIR_MAX_DIFF_TM'] = float(args.primer_pair_max_diff_tm)
    if args.dna_conc is not None:
        settings['pcr']['priming']['primer3']['PRIMER_DNA_CONC'] = float(args.dna_conc)
    if args.dna_conc_probe is not None:
        settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_DNA_CONC'] = float(args.dna_conc_probe)
    if args.salt_monovalent_conc is not None:
        settings['pcr']['priming']['primer3']['PRIMER_SALT_MONOVALENT'] = float(args.salt_monovalent_conc)
        settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_SALT_MONOVALENT'] = float(args.salt_monovalent_conc)
    if args.salt_divalent_conc is not None:
        settings['pcr']['priming']['primer3']['PRIMER_SALT_DIVALENT'] = float(args.salt_divalent_conc)
        settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_SALT_DIVALENT'] = float(args.salt_divalent_conc)
    if args.dntp_conc is not None:
        settings['pcr']['priming']['primer3']['PRIMER_DNTP_CONC'] = float(args.dntp_conc)
        settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_DNTP_CONC'] = float(args.dntp_conc)
    if args.tm_thresholds is not None:
        # Convert input string to list
        try:
            args.tm_thresholds = json.loads(args.tm_thresholds)
        except:
            exit(f"input parameter value of --tm_thresholds invalid! '{args.tm_thresholds}'\nPlease use a valid JSON input, fx '[0,     20, 40, [55,65], [59,61]]'")
    if args.anneal_tm is not None:
        ta = float(args.anneal_tm)
        # Modify primer3 temperature settings according to the user specified anneal Temperature (Ta)
        settings['pcr']['priming']['primer3']['PRIMER_MIN_TM'] = ta - 3
        settings['pcr']['priming']['primer3']['PRIMER_OPT_TM'] = ta
        settings['pcr']['priming']['primer3']['PRIMER_MAX_TM'] = ta + 3
        settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_MIN_TM'] = ta + 8
        settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_OPT_TM'] = ta + 10
        settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_MAX_TM'] = ta + 12
        # Set grading thresholds for melting temperature of primer binding sites
        if args.tm_thresholds is None:
            args.tm_thresholds = (0, 20, 40, [ta-5,ta+5],[ta-1,ta+1])
    if args.max_3end_gc is not None:
        settings['pcr']['priming']['primer3']['PRIMER_MAX_END_GC'] = int(args.max_3end_gc)
    if args.product_size_min is not None and args.product_size_max is not None:
        pr = [int(args.product_size_min), int(args.product_size_max)]
        settings['pcr']['priming']['primer3']['PRIMER_PRODUCT_SIZE_RANGE'] = pr
    settings['pcr']['priming']['primer3']['PRIMER_PICK_INTERNAL_OLIGO'] = 1 if args.pick_probe else 0
    if args.annotation_evalue is not None:
        settings['pcr']['annotation']['blastx_settings']['evalue'] = float(args.annotation_evalue)
    if args.kmer_count_threshold is not None:
        settings['explore']['kmer_count_threshold'] = float(args.kmer_count_threshold)
    if args.sensitivity_threshold is not None:
        settings['explore']['sensitivity_threshold'] = float(args.sensitivity_threshold)
    if args.fallout_threshold is not None:
        settings['explore']['fall-out_threshold'] = float(args.fallout_threshold)
    if args.align_percent_threshold is not None:
        settings['explore']['align_percent_threshold'] = float(args.align_percent_threshold)

    # Handle wildcards in positives, negatives and references
    if args.positives is not None:
        args.positives = get_fasta_files(args.positives, subdir='positives')
    if args.negatives is not None:
        args.negatives = get_fasta_files(args.negatives, subdir='negatives')
    if args.references is not None:
        args.references = get_fasta_files(args.references, subdir='references')

    args.quiet = False if args.verbose else True

    # Run Service
    locals().get(args.entry_point)(args)
