# -*- coding: utf-8 -*-
"""
    Sequence library
    ~~~~~~~~~~~~~~~~~~~~~

    Contains sequence handling functionalities.

"""
# (c) 2023 Martin Thomsen

from difflib import SequenceMatcher

from rucs import settings
from rucs.file import seqs_from_file


def analyse_genome(genome, min_seq_len=300):
    ''' Analyse the genome, and annotate sequence and base data. '''
    buffer = settings['input']['use_ram_buffer']
    # counts (seqs, bases, seqs >threshold, bases >threshold)
    counts = [0, 0, 0, 0]
    try:
        for seq, n, d in seqs_from_file(genome, use_ram_buffer=buffer):
            seqlen = len(seq)
            counts[0] += 1
            counts[1] += seqlen
            if seqlen >= min_seq_len:
                counts[2] += 1
                counts[3] += seqlen
    except IOError as e:
        pass

    return counts


def reverse_complement(seq):
    ''' Compute the reverse complementary DNA string.
    >>> print(reverse_complement('ATGCTTGCGATCGTGTCAGTGCGTATGCTAGCTTCTAGGTTCGA'))
    TCGAACCTAGAAGCTAGCATACGCACTGACACGATCGCAAGCAT
    '''
    return seq.translate(str.maketrans("ATGC","TACG"))[::-1]


def split_scaffolds_into_contigs(scaffolds, gap='n'):
    ''' Expands scaffolds into contigs, by splitting the sequences by gap

    USAGE
        >>> scaffolds = [('nnnnGTTAnnnTAGCGTACAnn', 'name', '')]
        >>> split_scaffolds_into_contigs(scaffolds)
        [('TAGCGTACA', 'name_1', 'position=11'), ('GTTA', 'name_0', 'position=4')]
    '''
    contigs = []
    contig_count = 0
    for seq, name, desc in scaffolds:
        if not gap in seq:
            # Scaffold is a contig
            contigs.append((seq, name, 'position=0'))
            contig_count += 1
        else:
            # Extract sub contigs
            for contig, position in extract_contigs(seq, gap):
                desc = 'position=%s'%(position)
                contigs.append((contig, "%s_%s"%(name, contig_count), desc))
                contig_count += 1

    # Sort and return the processed sequences according to sequence length
    contigs.sort(key=lambda x: len(x[0]), reverse=True)
    return contigs


def extract_contigs(seq, gap='n'):
    ''' Extract sub-sequences (and their position) which are divided by gaps

    USAGE
        >>> seq = 'nnnTACGTACGTACGGCATnnnGTTCATGCTAnnn'
        >>> extract_contigs(seq)
        [('TACGTACGTACGGCAT', 3), ('GTTCATGCTA', 22)]
    '''
    contigs = []
    start = 0
    end = 0
    in_gap = True
    for i, c in enumerate(seq):
        if c == gap:
            in_gap = True
            end = i - 1
            if end > start:
                contig = seq[start:end+1]
                contigs.append((contig, start))
            start = i
        elif in_gap:
            in_gap = False
            start = i

    if not in_gap and len(seq) >= start:
        contig = seq[start:]
        contigs.append((contig, start))

    return contigs


def align_seqs(a, b=None):
    ''' Align 2 sequences, '-' are used for gaps.

    USAGE
        >>> from difflib import SequenceMatcher
        >>> seq1 = 'TGGCCGCAAATGCATTCCCCA'
        >>> seq2 = 'GGCAATGCATTATGTGGGCAT'
        >>> print("%s\n%s\nSimilarity score: %s"%align_seqs(seq1, seq2))
        TGGCCGCAAATGCATTCCC----CA-
        -GGC----AATGCATTATGTGGGCAT
        Similarity score: 0.619
    '''
    if b is None: a, b = a
    s = [[], []]
    p = [0, 0]
    match = SequenceMatcher(a=a, b=b)
    similarity = round(match.ratio(), 3)
    for op, a0, a1, b0, b1 in match.get_opcodes():
        diff = (b1-b0)-(a1-a0)
        s[0].append('-'*(a0-p[0]))
        s[1].append('-'*(b0-p[1]))
        s[0].append(a[a0:a1])
        s[1].append(b[b0:b1])
        if diff > 0:
            s[0].append('-'*(diff))
        elif diff < 0:
            s[1].append('-'*(-diff))
        p[0] = a1
        p[1] = b1

    return (''.join(s[0]), ''.join(s[1]), similarity)
