# -*- coding: utf-8 -*-
"""
    Explore feature lib
    ~~~~~~~~~~~~~~~~~~~~~

    Contains relevant functionalities for the explore method.

"""
# (c) 2023 Martin Thomsen

import os, pickle, gc

from rucs import settings, get_directories
from rucs.log import log
from rucs.file import open_, seqs_from_file, save_as_fastq
from rucs.seq import analyse_genome, reverse_complement
from rucs.bwa import align_to_ref
from rucs.kmer import compute_consensus_sequences, extract_kmers_from_file


def explore_representation(positives, negatives, kmer_size=None, reuse=False):
    ''' Explore Over- and underrepresentation of k-mer in the positive genomes
    versus the negative genomes.
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

    kmer_count_threshold = settings['explore']['kmer_count_threshold']
    sensitivity_threshold = settings['explore']['sensitivity_threshold']
    fall_out_threshold = settings['explore']['fall-out_threshold']
    align_percent_threshold = settings['explore']['align_percent_threshold']
    kmer_seps = ['A', 'T', 'G', 'C']
    p_count = len(positives)
    n_count = len(negatives)
    ors_sens_threshold = sensitivity_threshold * p_count if type(sensitivity_threshold) is float else sensitivity_threshold
    ors_fo_threshold = fall_out_threshold * n_count if type(fall_out_threshold) is float else fall_out_threshold
    urs_sens_threshold = sensitivity_threshold * n_count if type(sensitivity_threshold) is float else sensitivity_threshold
    urs_fo_threshold = fall_out_threshold * p_count if type(fall_out_threshold) is float else fall_out_threshold

    if log is not None:
        # Initialize k-mer and sequences statistics logging
        log.progress.add('explore', 'Exploring over- and underrepresentation', 'main')
        log.stats.add_table('kmer', 'k-mer analyses', ['k-mers', 'Process'])
        log.stats.add_table('seqs', 'Sequence Analyses',
                            ['Fasta file', 'Sequences', 'Size in bases',
                             'Seqs >%s'%min_seq_len, 'Size >%s'%min_seq_len])

    # Count k-mers from positive references
    if log is not None:
        log.progress.add('pos', 'Counting Positive k-mers', 'explore')

    kmer_counts_pos = count_kmers(positives, kmer_size, kmer_count_threshold,
                                  min_seq_len=settings['ucs']['min_seq_len_pos'],
                                  sep_on_first=1, log_entry='pos', reuse=reuse)

    pos_kmer_count = sum(map(len, kmer_counts_pos.values()))
    if pos_kmer_count == 0:
        raise UserWarning('No positive k-mers were found!')

    # Pickle the kmer-keys and kmer_counts_pos dictionary and remove object to save space
    if log is not None:
        log.stats.add_row('kmer',[pos_kmer_count, 'Number of positive k-mers'])
        log.progress.add('store_pos_kmercounts', 'store k-mers counts', 'pos')

    for sep in kmer_seps:
        with open_(f'{work_dir}kmer_counts_pos_{sep}.pkl', 'wb') as f:
            pickle.dump(kmer_counts_pos[sep], f)

    del kmer_counts_pos
    gc.collect()

    # Count k-mers from negative references
    if log is not None:
        log.progress['store_pos_kmercounts'].log_time()
        log.progress['pos'].log_time()
        log.progress.add('neg', 'Counting Negative k-mers', 'explore')

    kmer_counts_neg = count_kmers(negatives, kmer_size, kmer_count_threshold,
                                  min_seq_len=settings['ucs']['min_seq_len_pos'],
                                  sep_on_first=1, log_entry='neg', reuse=reuse)

    neg_kmer_count = sum(map(len, kmer_counts_neg.values()))
    if neg_kmer_count == 0:
        raise UserWarning('No negative k-mers were found!')

    # Pickle the kmer-keys and kmer_counts_neg dictionary and remove object to save space
    if log is not None:
        log.stats.add_row('kmer', [neg_kmer_count, 'Number of negative k-mers'])
        log.progress.add('store_neg_kmercounts', 'store k-mers counts', 'neg')

    for sep in kmer_seps:
        with open_(f'{work_dir}kmer_counts_neg_{sep}.pkl', 'wb') as f:
            pickle.dump(kmer_counts_neg[sep], f)

    del kmer_counts_neg
    gc.collect()

    # Filter insignificant kmers
    if log is not None:
        log.progress['store_neg_kmercounts'].log_time()
        log.progress['neg'].log_time()
        log.progress.add('filter', 'Filtering k-mers', 'explore')

    pos_kmer_count = 0
    neg_kmer_count = 0
    # Process the kmer-sets independently
    for sep in kmer_seps:
        if log is not None:
            log.progress.add(f'filter_{sep}', f'Processing k-mers starting with "{sep}"', 'filter')

        # Load kmer counts
        with open_(f'{work_dir}kmer_counts_pos_{sep}.pkl', 'rb') as f:
            kmer_counts_pos = pickle.load(f)

        with open_(f'{work_dir}kmer_counts_neg_{sep}.pkl', 'rb') as f:
            kmer_counts_neg = pickle.load(f)

        # Identify significantly over- or under-represented k-mers
        for kmer in kmer_counts_pos.keys() | kmer_counts_neg.keys():
            pkc = kmer_counts_pos[kmer] if kmer in kmer_counts_pos else 0
            nkc = kmer_counts_neg[kmer] if kmer in kmer_counts_neg else 0

            # Filter k-mers not passing sensitivity and specificity threshold
            if pkc > 0 and (pkc < ors_sens_threshold or nkc > ors_fo_threshold):
                del kmer_counts_pos[kmer]

            if nkc > 0 and (nkc < urs_sens_threshold or pkc > urs_fo_threshold):
                del kmer_counts_neg[kmer]

        pos_kmer_count += len(kmer_counts_pos)
        neg_kmer_count += len(kmer_counts_neg)

        gc.collect()
        # Pickle the significant k-mer counts dictionary
        with open_(f'{work_dir}kmer_counts_pos_sig_{sep}.pkl', 'wb') as f:
            pickle.dump(set(kmer_counts_pos.keys()), f)

        with open_(f'{work_dir}kmer_counts_neg_sig_{sep}.pkl', 'wb') as f:
            pickle.dump(set(kmer_counts_neg.keys()), f)

        del kmer_counts_pos, kmer_counts_neg
        gc.collect()

        if log is not None:
            log.progress[f'filter_{sep}'].log_time()

    if log is not None:
        log.progress['filter'].log_time()
        log.stats.add_row('kmer', [pos_kmer_count, 'Number of significant positive k-mers'])
        log.stats.add_row('kmer', [neg_kmer_count, 'Number of significant negative k-mers'])

    # Create consensus sequences for the over represented k-mers
    if log is not None:
        log.progress.add('make_ors', 'Computing k-mer contigs and scaffolds (ORS)','explore')

    # Load all significant positive kmer counts
    sig_pos_kmers = set()
    for sep in kmer_seps:
        with open_(f'{work_dir}kmer_counts_pos_sig_{sep}.pkl', 'rb') as f:
            sig_pos_kmers.update(pickle.load(f))

    # Align to pos refs until enough k-mers has been aligned to a ref
    max_pos = len(sig_pos_kmers)
    combined_contigs_ors = 'ors.contigs.fa'
    combined_disscafs_ors = 'ors.disscafs.fa'
    with open_(combined_contigs_ors, 'w') as f_cont, open_(combined_disscafs_ors, 'w') as f_diss:
        for ref in positives:
            align_percent = len(sig_pos_kmers) / max_pos
            if max_pos and align_percent > align_percent_threshold:
                prefix = 'ors_%s'%(os.path.basename(ref).rsplit('.', 1)[0])
                # Store k-mers as fastq
                if log is not None:
                    log.progress.add(f'make_fq_{ref}', f'Preparing k-mers for alignment (make fastq)', 'make_ors')
                pos_kmers_fq = 'pos_kmers.fq'
                save_as_fastq(sig_pos_kmers, pos_kmers_fq)
                if log is not None:
                    log.progress[f'make_fq_{ref}'].log_time()
                # Align k-mers to reference (Note: Unmapped k-mers are lost in this process)
                pos_kmers = align_to_ref(ref, pos_kmers_fq, settings['ucs']['bwa_settings'],
                                         prefix, settings['ucs']['sam_flags_ignore'],
                                         log_entry='make_ors')
                # Compute sequences
                if log is not None:
                    log.progress.add(f'make_sequences_{ref}', 'Building contigs and scaffolds', 'make_ors')
                ors_files = compute_consensus_sequences(pos_kmers, ref, kmer_size,
                                                        settings['ucs']['charspace'],
                                                        prefix+'cons')
                # Add contigs from ref for ref to combined contigs
                try: f_cont.write('\n'.join(">%s_%s\n%s"%(prefix,n,s) for s,n,d in seqs_from_file(ors_files[1]) if len(s) >= kmer_size))
                except IOError: pass
                else: f_cont.write('\n')

                # Add dissected scaffolds for ref to combined dissected scaffolds
                try: f_diss.write('\n'.join(">%s_%s\n%s"%(prefix,n,s) for s,n,d in seqs_from_file(ors_files[3]) if len(s) >= kmer_size))
                except IOError: pass
                else: f_diss.write('\n')

                # Remove aligned k-mers and thier reverse complement
                if log is not None:
                    log.progress[f'make_sequences_{ref}'].log_time()
                    log.progress.add(f'remove_kmers_{ref}', f'Removing the k-mers that was aligned sucessfully', 'make_ors')

                sig_pos_kmers -= pos_kmers.keys()
                sig_pos_kmers -= set(reverse_complement(k) for k in pos_kmers.keys())

                if log is not None:
                    log.progress[f'remove_kmers_{ref}'].log_time()
                    log.progress.add(f'status_{ref}', f'{len(sig_pos_kmers)} k-mers out of {max_pos} ({int(len(sig_pos_kmers) / max_pos  *     100)}%) remains to be aligned.', 'make_ors')

    if log is not None:
        log.progress['make_ors'].log_time()
        counts = analyse_genome(combined_contigs_ors, min_seq_len)
        log.stats.add_row('seqs', [os.path.basename(combined_contigs_ors)] + counts)

    # Delete sig_pos_kmers list to save space
    del sig_pos_kmers
    gc.collect()

    # Create consensus sequences for the under represented k-mers
    if log is not None:
        log.progress.add('make_urs', 'Computing k-mer contigs and scaffolds (URS)','explore')

    # Load all significant negative kmer counts
    sig_neg_kmers = set()
    for sep in kmer_seps:
        with open_(f'{work_dir}kmer_counts_neg_sig_{sep}.pkl', 'rb') as f:
            sig_neg_kmers.update(pickle.load(f))

    # WHILE loop neg refs until enough k-mers has been aligned to a ref
    max_neg = len(sig_neg_kmers)
    combined_contigs_urs = 'urs.contigs.fa'
    combined_disscafs_urs = 'urs.disscafs.fa'
    with open_(combined_contigs_urs, 'w') as f_cont, open_(combined_disscafs_urs, 'w') as f_diss:
        for ref in negatives:
            align_percent = len(sig_neg_kmers) / max_neg
            if max_neg and align_percent > align_percent_threshold:
                prefix = 'urs_%s'%(os.path.basename(ref).rsplit('.', 1)[0])
                # Store k-mers as fastq
                neg_kmers_fq = 'neg_kmers.fq'
                save_as_fastq(sig_neg_kmers, neg_kmers_fq)
                # Align k-mers to reference (Note: Unmapped k-mers are lost in this process)
                neg_kmers = align_to_ref(ref, neg_kmers_fq, settings['ucs']['bwa_settings'],
                                         prefix, settings['ucs']['sam_flags_ignore'],
                                         log_entry='make_urs')
                # Compute sequences
                urs_files = compute_consensus_sequences(neg_kmers, ref, kmer_size,
                                                        settings['ucs']['charspace'],
                                                        prefix+'cons')
                # Add contigs from ref for ref to combined contigs
                try: f_cont.write('\n'.join(">%s_%s\n%s"%(prefix,n,s) for s,n,d in seqs_from_file(urs_files[1]) if len(s) >= kmer_size))
                except IOError: pass
                else: f_cont.write('\n')

                # Add dissected scaffolds for ref to combined dissected scaffolds
                try: f_diss.write('\n'.join(">%s_%s\n%s"%(prefix,n,s) for s,n,d in seqs_from_file(urs_files[3]) if len(s) >= kmer_size))
                except IOError: pass
                else: f_diss.write('\n')

                # Remove aligned k-mers and thier reverse complement
                if log is not None:
                    log.progress.add(f'remove_kmers_{ref}', f'Removing the k-mers that was aligned sucessfully', 'make_urs')

                sig_neg_kmers -= neg_kmers.keys()
                sig_neg_kmers -= set(reverse_complement(k) for k in neg_kmers.keys())

                if log is not None:
                    log.progress.add(f'status_{ref}', f'{len(sig_neg_kmers)} k-mers out of {max_neg} ({int(len(sig_neg_kmers) / max_neg    *100)}%)  remains to be aligned.', 'make_urs')

    if log is not None:
        log.progress['make_urs'].log_time()
        counts = analyse_genome(combined_contigs_urs, min_seq_len)
        log.stats.add_row('seqs', [os.path.basename(combined_contigs_urs)] + counts)
        log.progress['explore'].log_time()

    return ((combined_contigs_ors, combined_disscafs_ors),
           (combined_contigs_urs, combined_disscafs_urs))


def count_kmers(files, kmer_size=20, kmer_count_threshold=1, min_seq_len=500,
                sep_on_first=0, log_entry=None, reuse=False):
    ''' Count in how many files, a given k-mer is found. '''
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
                                          min_seq_len, to_upper, buffer, log_entry, reuse=reuse)
        # COMPUTE K-MER COUNT
        if log is not None and log_entry is not None:
            log.progress.add(f'kmercount_{os.path.basename(genome)}', f'Updating k-mer counts', log_entry)

        for kmer in kmers_i:
            if sep_on_first > 0:
                sep = kmer[0:sep_on_first].upper()
                if not sep in kmers:
                    kmers[sep] = {}
                if not kmer in kmers[sep]:
                    kmers[sep][kmer] = 1
                else:
                    kmers[sep][kmer] += 1
            else:
                if not kmer in kmers:
                    kmers[kmer] = 1
                else:
                    kmers[kmer] += 1
            if log is not None and log_entry is not None:
                log.progress[f'kmercount_{os.path.basename(genome)}'].log_time()

    # Filter low occuring k-mers
    if kmer_count_threshold > 1:
        if log is not None and log_entry is not None:
            log.progress.add(f'filter_kmers_{log_entry}', f'Filter low occuring k-mers', log_entry)
        if sep_on_first > 0:
            kc_before = sum(len(kmers[sep]) for sep in kmers)
            for sep in kmers:
                for kmer in kmers[sep]:
                    if kmers[kmer] < kmer_count_threshold:
                        del kmers[kmer]
            kc_after = sum(len(kmers[sep]) for sep in kmers)
        else:
            kc_before = len(kmers)
            for kmer in kmers:
                if kmers[kmer] < kmer_count_threshold:
                    del kmers[kmer]
            kc_after = len(kmers)
        if log is not None and log_entry is not None:
            log.stats.add_row('kmer',[kc_before - kc_after, 'Number of low occuring k-mers filtered'])
            log.progress[f'filter_kmers_{log_entry}'].log_time()

    return kmers
