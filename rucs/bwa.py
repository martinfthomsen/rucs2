# -*- coding: utf-8 -*-
"""
    BWA library
    ~~~~~~~~~~~~~~~~~~~~~

    Contains functionalities relating to BWA.

"""
# (c) 2023 Martin Thomsen

import sys, os
from subprocess import Popen, PIPE

from rucs import settings
from rucs.log import log
from rucs.file import which, open_
from rucs.seq import reverse_complement

# CLASSES
class DependencyError(EnvironmentError):
    '''raise this when there's an unsolved program dependency error'''
    pass


# Functions
def align_to_ref(reference, fastq, bwa_settings=None, output_prefix='aln',
                 ignore_flags=4, log_entry='ucs'):
    ''' Align the sequences in the fastq to the reference

    output is a dictionary of the sequences containing all the alignment hits as
    contigname_position pairs.

    USAGE
        >>> sequences = ['CAACATTTTCGTGTCGCCCTT', 'TGAAGCCATACCAAACGACGA']
        >>> save_as_fastq(sequences, 'sequences.fastq')
        >>> align_to_ref('tem_reference.fasta', 'sequences.fastq')
        {'CAACATTTTCGTGTCGCCCTT': ['TEM-1b_10'], 'TGAAGCCATACCAAACGACGA': ['TEM-1b_504']}
    '''
    if log is not None:
        log.progress.add(output_prefix,
                         "Aligning %s to %s"%(os.path.basename(fastq),
                                              os.path.basename(reference)),
                         log_entry)
    # Align sequences to reference
    bam = bwa(reference, fastq, None, bwa_settings, output_prefix)
    # Extract alignment details
    alignments = extract_alignment_details(bam, ignore_flags)

    if log is not None:
        log.progress[output_prefix].log_time()
        log.stats.add_row('kmer', [len(alignments),
                                   "Aligning %s to %s"%(
                                    os.path.basename(fastq),
                                    os.path.basename(reference))
                                   ])
    return alignments

def bwa(ref, fastq, paired_fastq=None, bwa_settings=None, output_prefix='aln'):
    ''' Wrapper for running BWA commandline '''
    # Create path class and object and method to extract them
    path = type('path_object', (object,), {})()
    all_elements = lambda obj: [p for p in dir(obj) if p[:2] != '__']

    # Check Bowtie can be found in the system environment
    path.bwa = which('bwa')
    path.samtools = which('samtools')
    for p in all_elements(path):
        if path.__getattribute__(p) is None:
            raise DependencyError(1, '%s could not be found'%p)

    # Validate inputs -M? -R? -N?
    # http://bio-bwa.sourceforge.net/bwa.shtml
    defaults = { # Name: (type, arg, default_value)
        'max_mismatch':       (int, '-n',   6), # edit distance
        'max_mismatch_score': (int, '-M',  10), # edit distance
        'max_gaps':           (int, '-o',   1), # max number of gaps
        'max_gap_ext':        (int, '-e',   0), # max gap extensions
        'gap_penalty_open':   (int, '-O',  10), # gap penalty open
        'gap_penalty_ext':    (int, '-E',  10), # gap penalty extension
        'seed_length':        (int, '-l',  35), # Length of the seed
        'max_mismatch_seed':  (int, '-k',   0), # edit distance in seed
        'threads':            (int, '-t',   4)  # paralelisation using multi threads
    }
    if bwa_settings is None: bwa_settings = {}
    for s, v in bwa_settings.items():
        assert s in defaults, 'unknown setting encountered (%s)!'%s
        assert isinstance(v, defaults[s][0]), ('invalid type for %s, expected %s,'
                                               ' found %s!')%(s, defaults[s][0],
                                                              type(v))

    # Create BWA index $ bwa index <REFERENCE>
    if not os.path.exists("%s.sa"%(ref)):
        cmd = [path.bwa, 'index', ref]
        # _ = sys.stderr.write("#CMD=%s\n"%(' '.join(cmd)))
        so_path = '%s_index.out.txt'%output_prefix
        se_path = '%s_index.err.txt'%output_prefix
        with open_(so_path, 'w') as so, open_(se_path, 'w') as se:
            ec = Popen(cmd, stdout=so, stderr=se).wait()
            if ec != 0: raise RuntimeError('bwa indexing failed during execution')

    # Align seqs to reference $ bwa aln <REFERENCE> <FASTQ>
    cmd = [path.bwa, 'aln']
    for s, (t, a, v) in defaults.items():
        if s in bwa_settings: cmd.extend([a, str(bwa_settings[s])])
        elif v is not None: cmd.extend([a, str(v)])

    cmd.extend([ref, fastq])
    # _ = sys.stderr.write("#CMD=%s\n"%(' '.join(cmd)))
    sai_file = '%s_1.sai'%output_prefix
    se_path = '%s_sai.err.txt'%output_prefix
    with open_(sai_file, 'w') as so, open_(se_path, 'w') as se:
        ec = Popen(cmd, stdout=so, stderr=se).wait()
        if ec != 0: raise RuntimeError('bwa aln failed during execution')

    # Get BAM file
    path.sam = '%s.sam'%(output_prefix)
    path.bam = '%s.bam'%(output_prefix)
    path.bwalog = '%s_bwa.out.txt'%(output_prefix)
    if paired_fastq is None:
        # Align using BWA (single end) $ bwa samse <REFERENCE> <SAI> <FASTQ>
        cmd1 = [path.bwa, 'samse', '-n', '99', ref, sai_file, fastq]
        cmd2 = [path.samtools, 'view', '-Sb', path.sam]
        with open_(path.sam, 'w') as sam, open_(path.bwalog, 'w') as se:
            # Run BWA
            p1 = Popen(cmd1, stderr=se, stdout=sam)
            ec1 = p1.wait()
            if ec1 != 0: raise RuntimeError('bwa failed during execution (return code: %s)\n%s'%(ec1, ' '.join(cmd1)))
        with open_(path.bam, 'w') as bam, open_(path.bwalog, 'a') as se:
            # Convert SAM to BAM output for disk space efficiency
            p2 = Popen(cmd2, stdout=bam, stderr=se)
            # Wait for the programs to finish
            ec2 = p2.wait()
            if ec2 != 0: raise RuntimeError('samtools failed during execution (return code: %s)\n%s'%(ec2,
            ' '.join(cmd2)))
        os.unlink(path.sam)
    else:
        # Align paired seqs to reference $ bwa aln <REFERENCE> <FASTQ>
        cmd = [path.bwa, 'aln']
        for s, (t, a, v) in defaults.items():
            if s in bwa_settings: cmd.extend([a, str(bwa_settings[s])])
            elif v is not None: cmd.extend([a, str(v)])

        cmd.extend([ref, paired_fastq])
        # _ = sys.stderr.write("#CMD=%s\n"%(' '.join(cmd)))
        sai_file_2 = '%s_2.sai'%output_prefix
        se_path_2 = '%s_2_sai.err.txt'%output_prefix
        with open_(sai_file_2, 'w') as so, open_(se_path_2, 'w') as se:
            ec = Popen(cmd, stdout=so, stderr=se).wait()
            if ec != 0: raise RuntimeError('bwa aln failed during execution')

        # Align using BWA (paired end) $ bwa sampe <REFERENCE> <SAI_1> <SAI_2> <FASTQ_1> <FASTQ_2>      -P?
        cmd = [path.bwa, 'sampe', '-n', '99', '-o', '2000', ref, sai_file, sai_file_2, fastq, paired_fastq]
        # sys.stderr.write("#CMD=%s\n"%(' '.join(cmd)))
        with open_(path.bam, 'w') as so, open_(path.bwalog, 'w') as se:
            # Run bowtie
            p1 = Popen(cmd, stderr=se, stdout=PIPE)
            # Convert SAM to BAM output for disk space efficiency
            p2 = Popen([path.samtools, 'view', '-Sb'], stdin=p1.stdout, stdout=so)
            # Wait for the programs to finish
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
            ec1 = p1.wait()
            ec2 = p2.wait()
            if ec1 != 0: raise RuntimeError('bwa failed during execution')
            if ec2 != 0: raise RuntimeError('samtools failed during execution')

    # Clean up
    os.unlink(sai_file)
    return path.bam


def extract_alignment_details(bam, ignore_flags=4):
    ''' Adds contig and position details from the sequence mapping (sam) to the
    sequence dictionary.

    FLAGS   DESCRIPTION  (NOTE: the flag is a binary combination: 20 = 4 + 16)
       1   the read is paired in sequencing
       2   the read is mapped in a proper pair
       4   the query sequence itself is unmapped
       8   the mate is unmapped
      16   strand of the query (1 for reverse)
      32   strand of the mate
      64   the read is the first read in a pair
     128   the read is the second read in a pair
     256   the alignment is not primary
     512   QC failure
    1024   optical or PCR duplicate
    2048   supplementary alignment

    SAM COLUMNS:
    QNAME Query template/pair NAME
    FLAG  bitwise FLAG
    RNAME Reference sequence NAME
    POS   1-based leftmost POSition/coordinate of clipped sequence
    MAPQ  MAPping Quality (Phred-scaled)
    CIGAR extended CIGAR string
    MRNM  Mate Reference sequence NaMe ('=' if same as RNAME)
    MPOS  1-based Mate POSistion
    LEN   inferred Template LENgth (insert size)
    SEQ   query SEQuence on the same strand as the reference
    QUAL  query QUALity (ASCII-33 gives the Phred base quality)
    OPT   variable OPTional fields in the format TAG:VTYPE:VALUE
         XA:  Alternative hits
         NM:  Edit distance

    XA FORMAT:
    (RNAME,POS,CIGAR,NM;)
    Note: The sign in front of POS (+ or -) determines the strand
    '''
    # Create path class and object and method to extract them
    path = type('path_object', (object,), {})()
    all_elements = lambda obj: [p for p in dir(obj) if p[:2] != '__']

    # Check Bowtie can be found in the system environment
    path.samtools = which('samtools')
    for p in all_elements(path):
        if path.__getattribute__(p) is None:
            raise DependencyError(1, '%s could not be found'%p)

    # Parse BAM output -> create tab-list of matches
    # with pysam.AlignmentFile("ex1.bam", "rb") as si:
    #    for alignment in fin.fetch():
    #http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
    #pysam.view(path.bam).split('\n')
    # Extract all matches with flag: 99 = first mate of proper paired reads
    cmd = [path.samtools, 'view', '-F', str(ignore_flags), bam]
    p1 = Popen(cmd, stderr=PIPE, stdout=PIPE) #stdin=si,
    stdout, stderr = p1.communicate()
    # Extract possible PCR products for the queries
    ignore_flags = list(x == '1' for x in reversed("{0:b}".format(ignore_flags).zfill(12)))
    alignments = {}
    for l in stdout.decode("utf-8").split('\n'):
        if l.strip() == '': continue
        try:
            tmp = l.split('\t')
            flags = list(x == '1' for x in reversed("{0:b}".format(int(tmp[1])).zfill(12)))
            if flags[2] and ignore_flags[2]: # Ignore unmatched
                continue
            if flags[4]:
                strand = '-'
            else:
                strand = '+'
            seq = tmp[9] if strand == '+' else reverse_complement(tmp[9])
            contig_name = tmp[2]
            position = int(tmp[3]) - 1
            cigar = tmp[5]
            # Update Sequence Position(s)
            if not (flags[4] and ignore_flags[4]): # Ignore negative strand
                if not seq in alignments: alignments[seq] = []
                alignments[seq].append((contig_name, cigar, strand, position))
            for i in range(11, len(tmp)):
                # Skip all options which are not alternative positions
                if tmp[i][:2] != 'XA': continue
                # Alternative sequences found
                XA = tmp[i].split(':')[-1].split(';')
                for alt in XA:
                    if alt == '': continue # skip empty alternatives
                    a = alt.split(',')
                    contig_name = a[0]
                    cigar = a[2]
                    strand = a[1][0]
                    pos_alt = int(a[1][1:])-1
                    if strand == '-' and ignore_flags[4]: # Ignore negative strand
                        continue
                    # Add alternative position to hit_list
                    if not seq in alignments: alignments[seq] = []
                    alignments[seq].append((contig_name, cigar, strand, pos_alt))
        except:
            print(l)
            raise

    return alignments
