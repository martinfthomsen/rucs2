# -*- coding: utf-8 -*-
"""
    Blast library
    ~~~~~~~~~~~~~~~~~~~~~

    Contains BLAST functionalities.

"""
# (c) 2023 Martin Thomsen

import sys, os, re, pickle
import numpy as np
from subprocess import Popen, PIPE, DEVNULL
from difflib import SequenceMatcher

from rucs import settings, get_directories
from rucs.file import open_, seqs_from_file


def blast_to_ref(reference, fasta, blast_settings=None, buffer=False):
    ''' BLAST the fasta sequences to the reference

    -num_alignments: Number of database sequences to show alignments for
    -max_hsps      : Set maximum number of HSPs per subject sequence to save for
                     each query
    -query_loc     : Location on the query sequence in 1-based offsets
                     (Format: start-stop)
    -strand        : Strand: both, minus or plus
    -evalue        : Expectation value (E) threshold for saving hits
    -word_size     : Word size for wordfinder algorithm
                     (length of best perfect match)
    -gapopen       : Cost to open a gap
    -gapextend     : Cost to extend a gap
    -penalty       : Penalty for a nucleotide mismatch
    -reward        : Reward for a nucleotide match
    -dbtype        : 'nucl' for DNA and 'prot' for protein
    -dust          : Filter query sequence with DUST
                     (Format: 'yes', 'level window linker', or 'no' to disable)
    -num_threads   : paralelisation using multiple threads

    USAGE
        >>> alignments = blast_to_ref(reference, fasta)
    '''
    # Validate inputs
    assert ' ' not in reference, 'BLAST cannot handle spaces in reference path!'
    defaults = { # Name: (type, arg, default_value)
        'num_alignments': (int,   '-num_alignments',   1000),
        'max_hsps':       (int,   '-max_hsps',         None),
        'query_loc':      (str,   '-query_loc',        None),
        'strand':         (str,   '-strand',         'both'),
        'evalue':         (float, '-evalue',              8.0),
        'word_size':      (int,   '-word_size',          10),
        'gapopen':        (int,   '-gapopen',             0),
        'gapextend':      (int,   '-gapextend',           2),
        'penalty':        (int,   '-penalty',            -2),
        'reward':         (int,   '-reward',              1),
        'dbtype':         (str,   '-dbtype',         'nucl'),
        'dust':           (str,   '-dust',             'no'),
        'num_threads':    (int,   '-num_threads',         4)
    }
    if blast_settings is None: blast_settings = {}
    # Filter None Items
    blast_settings = {k:v for k,v in blast_settings.items() if v is not None}

    for s, v in blast_settings.items():
        assert s in defaults, 'unknown setting encountered (%s)!'%s
        assert isinstance(v, defaults[s][0]), ('invalid type for %s, expected %s,'
                                               ' found %s!')%(s, defaults[s][0],
                                                              type(v))

    # Get directories
    global work_dir, ref_dir, result_dir
    work_dir, ref_dir, result_dir = get_directories(get_ref_dir=True)

    # Create BLAST DB
    name = os.path.basename(reference).rsplit('.',1)[0]
    if not os.path.exists("%s.nin"%(reference)):
        cmd = ['makeblastdb', '-in', '-', '-title', name, '-out', reference, '-parse_seqids']
        if 'dbtype' in blast_settings:
            cmd.extend(['-dbtype', blast_settings['dbtype']])
        else:
            cmd.extend(['-dbtype', defaults['dbtype'][-1]])
        so_path = f'{work_dir}blast_db.out.txt'
        se_path = f'{work_dir}blast_db.err.txt'
        with open_(reference, 'r') as six, open_(so_path, 'w') as so, open_(se_path, 'w') as se:
            so.write(f"#CMD ")
            if reference[-3:] == '.gz':
                # gunzip input
                so.write(f"gzip -cd {reference} | ")
                pgz = Popen(['gzip', '-cd', reference], stdout=PIPE, stderr=DEVNULL)
                si = pgz.stdout
            else:
                so.write(f"cat {reference} | ")
                si = six
            so.write(f"{' '.join(cmd)}\n")
            so.flush()
            ec = Popen(cmd, stdin=si, stdout=so, stderr=se).wait()
            if ec != 0: raise RuntimeError(f'BLAST formatdb failed during execution ({ec})')
    del defaults['dbtype']

    # BLAST fasta to DB
    # blastn output format options (https://www.metagenomics.wiki/tools/blast/blastn-output-format-6)
    # qseqid    query or source (gene) sequence id
    # sseqid    subject or target (reference genome) sequence id
    # pident    percentage of identical positions
    # length    alignment length (sequence overlap)
    # mismatch  number of mismatches
    # gapopen   number of gap openings
    # qstart    start of alignment in query
    # qend      end of alignment in query
    # sstart    start of alignment in subject
    # send      end of alignment in subject
    # evalue    expect value
    # bitscore  bit score
    # sstrand   Subject Strand
    # qseq      Aligned part of query sequence
    # sseq      Aligned part of subject sequence
    # The E-value (expectation value) is a corrected bit-score adjusted to the sequence database size. The E-value therefore depends on the size of the used sequence database. Since large databases increase the chance of false positive hits,  the E-value corrects for the higher chance. It's a correction for multiple comparisons. This means that a sequence hit would get a better E-value when present in a smaller database.
    cmd = ['blastn', '-db', reference, '-query', fasta,
           '-outfmt', '6 qseqid sseqid sstrand qstart qend sstart send qseq sseq'] # evalue
    for s, (t, a, v) in defaults.items():
        if s in blast_settings: cmd.extend([a, str(blast_settings[s])])
        elif v is not None: cmd.extend([a, str(v)])

    # print(' '.join(cmd))
    so_path = f'{work_dir}blastn.{name}.tsv'
    se_path = f'{work_dir}blastn.err.txt'
    with open_(so_path, 'w') as so, open_(se_path, 'w') as se:
        ec = Popen(cmd, stdout=so, stderr=se).wait()
        if ec != 0:
            sys.stderr.write('Failing command:\n%s\n\n'%(' '.join(cmd)))
            raise RuntimeError('BLASTn failed during execution')

    # Extract alignments
    primers = dict((i, seq) for i,(seq, n, d) in enumerate(seqs_from_file(fasta, use_ram_buffer=buffer)))
    alignments = {}
    with open_(so_path) as f:
        for l in f:
            l = l.strip()
            if l == '': continue
            if l[0] == '#': continue
            try:
                # Parse BLAST line
                flds = l.split('\t')
                query = int(flds[0])
                contig_name = flds[1].split('|')[1] if '|' in flds[1] else flds[1]
                strand = '-' if flds[2] == 'minus' else '+'
                qstart, qend, sstart, send = map(int, flds[3:7])
                qseq, sseq = flds[7:]
            except:
                print('Error while parsing BLAST output:')
                print(l)
                raise
            else:
                primer = primers[query]
                if strand == '+':
                    position = sstart-qstart
                else:
                    position = send-len(primer)+qend-1
                # Correct position for gap offset
                gaps = [qseq.count('-'), sseq.count('-')]
                splts = ([len(x) for x in qseq.split('-')],
                         [len(x) for x in sseq.split('-')])
                tmp = (np.arange(len(splts[0])), np.arange(len(splts[1])))
                splsq = sum((tmp[0]-tmp[0][::-1])*splts[0])
                splss = sum((tmp[1]-tmp[1][::-1])*splts[1])
                if (splsq < 0 and strand == '-') or (splsq >= 0 and strand == '+'):
                    position += qseq.count('-')
                elif splsq == 0 and len(splts[0]) > 1:
                    position -= qseq[:round(len(qseq)/2)].count('-')

                if (splss >= 0 and strand == '+') or (splss < 0 and strand == '-'):
                    position -= sseq.count('-')
                elif splss == 0 and len(splts[1]) > 1:
                    position -= sseq[:round(len(sseq)/2)].count('-')

                # Add alignment
                if not primer in alignments: alignments[primer] = []
                alignments[primer].append((contig_name, 'NA', strand, position))

    # dump computed alignments to temporary pickle_file
    with open_(f'{work_dir}alignments_{name}.pkl', 'wb') as f:
        pickle.dump(alignments, f)

    return alignments


def get_blast_annotations(fasta, blast_settings=None, dbpath=''):
    ''' BLAST the fasta sequences to the reference

    FAST MODE: -task blastx-fast
    PARALELISM: -num_threads 4
    REDUCE REDUNDANT RESULTS:
        -culling_limit 5
        -best_hit_overhang 0.4
        -best_hit_score_edge 0.4

    GENETIC CODE TRRANSLATION:
       -query_gencode <integer>  # DEFAULT 1
    #    https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes
    #     1. The Standard Code
    #     2. The Vertebrate Mitochondrial Code
    #     3. The Yeast Mitochondrial Code
    #     4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
    #     5. The Invertebrate Mitochondrial Code
    #     6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
    #     9. The Echinoderm and Flatworm Mitochondrial Code
    #    10. The Euplotid Nuclear Code
    #    11. The Bacterial, Archaeal and Plant Plastid Code
    #    12. The Alternative Yeast Nuclear Code
    #    13. The Ascidian Mitochondrial Code
    #    14. The Alternative Flatworm Mitochondrial Code
    #    16. Chlorophycean Mitochondrial Code
    #    21. Trematode Mitochondrial Code
    #    22. Scenedesmus obliquus Mitochondrial Code
    #    23. Thraustochytrium Mitochondrial Code
    #    24. Pterobranchia Mitochondrial Code
    #    25. Candidate Division SR1 and Gracilibacteria Code
    #    26. Pachysolen tannophilus Nuclear Code

    USAGE
        >>> import re
        >>> annotations = get_blast_annotations(fasta)
        get_blast_annotations('unique_core_sequences.contigs.fa', dbpath='/path/to/blastdb/')
    '''
    class RegexTermReplacementObj():
        r''' For a given string, replace all keys from a given dictionary to their
        corresponding value.

        USAGE
            >>> import re
            >>> text = 'MULTI: hello--world!, partial [la la land]'
            >>> replace_terms = {'\[.+\]':'','--':' ','MULTI:':'',', partial':''}
            >>> regex = RegexTermReplacementObj(replace_terms)
            >>> print(regex.replace_terms(text).strip())
            hello world!
            '''
        def __init__(self, term_dict):
            self.term_dict = term_dict
            self.regex = re.compile(r"(%s)" % r"|".join(term_dict))
        def replace_terms(self, term):
            return self.regex.sub(self.mo2term, term)
        def mo2term(self, mo):
            match = mo.string[mo.start():mo.end()]
            if match in self.term_dict:
                return self.term_dict[match]
            else:
                for term in self.term_dict:
                    if re.match(term, match):
                        return self.term_dict[term]
                raise KeyError('could not find %s amongst %s'%(match, self.term_dict.keys()))

    def extract_annotations(so_path, ignore_terms, replace_terms, max_cov):
        ''' Extract BLAST annotations from tsv file '''
        annotations = {}

        # Prepare a regular expression from the replace terms
        regex = RegexTermReplacementObj(replace_terms)

        with open_(so_path) as f:
            for l in f:
                l = l.strip()
                if l == '': continue
                if l[0] == '#': continue
                try:
                    # Parse BLAST line
                    flds = l.split('\t')
                    query = flds[0]
                    qstart, qend, qlen = map(int, flds[1:4])
                    subject = flds[4]
                    sstart, send, slen = map(int, flds[5:8])
                    alen = int(flds[8])
                    identity = float(flds[9])
                    matches, mismatches, gaps = map(int, flds[10:13])
                    scov = round(alen/slen,3)
                    annotation = flds[13]
                    if qstart > qend:
                        tmp = qstart
                        qstart = qend
                        qend = tmp
                        strand = '-'
                    else:
                        strand = '+'
                except:
                    print('Error while parsing BLAST output:')
                    print(l)
                    raise
                else: # Parse BLAST output
                    tmp_annot = []
                    for annot in annotation.split(';'):
                        # Ignore selected terms
                        if any(term in annot for term in ignore_terms):
                            continue
                        # Replace selected terms
                        annot = regex.replace_terms(annot).strip()
                        tmp_annot.append(annot)

                    if not tmp_annot:
                        continue

                    # Store annotation
                    if not query in annotations:
                        annotations[query] = []
                    threshold = alen*3*max_cov # alignment length * codon size * threshold
                    covered = False
                    for i, (start, end, notes) in enumerate(annotations[query]):
                        # Check if the new annotation region has been covered previously
                        if qend-start >= threshold and end-qstart >= threshold:
                            new_start = qstart if qstart < start else start
                            new_end = qend if qend > end else end
                            for annot in tmp_annot:
                                if not annot in notes:
                                    notes.append(annot)
                            annotations[query][i] = (new_start, new_end, notes)
                            covered = True
                            break
                    if not covered:
                        annotations[query].append((qstart, qend, tmp_annot))

        return annotations

    def reduce_redundancy(annots, sim_cutoff):
        ''' Post process annotations (remove redundancy) '''
        new_annots = []
        for annot in annots:
            found = False
            for nannot in new_annots:
                m = SequenceMatcher(a=annot, b=nannot)
                sim = sum(n for i,j,n in m.get_matching_blocks())
                sim_a = sim/len(annot)
                sim_n = sim/len(nannot)
                if sim_a >= sim_n and sim_a >= sim_cutoff:
                    # print('MERGE (%s)! %s ==> %s'%(sim_a, annot, nannot))
                    found = True
                    break
                elif sim_n > sim_a and sim_n >= sim_cutoff:
                    # print('SWITCH (%s)! %s <== %s'%(sim_n, annot, nannot))
                    new_annots[new_annots.index(nannot)] = annot
                    found = True
                    break

            if not found:
                # print('NEW! %s : %s'%(annot, new_annots))
                new_annots.append(annot)

        return new_annots

    # Get directories
    global work_dir, ref_dir, result_dir
    work_dir, ref_dir, result_dir = get_directories(get_ref_dir=True)

    # Get settings
    max_cov = settings['pcr']['annotation']['max_cov']
    sim_cutoff = settings['pcr']['annotation']['sim_cutoff']
    ignore_terms = settings['pcr']['annotation']['ignore_terms']
    replace_terms = settings['pcr']['annotation']['replace_terms']
    db = settings['pcr']['annotation']['blastx_settings']['database']
    # Validate inputs
    if dbpath == '':
        if not 'BLASTDB' in os.environ:
            sys.stderr.write('Warning: BLAST annotation was not possible. BLASTDB '
                             'environment variable not set!\n')
            return {}
        elif not os.path.exists(os.environ['BLASTDB']):
            sys.stderr.write('Warning: BLAST annotation was not possible. BLASTDB '
                             'path does not exist!\n')
            return {}
    elif not os.path.exists(dbpath):
        sys.stderr.write('Warning: BLAST annotation was not possible. Provided '
                         'BLAST DB path does not exist!\n')
        return {}
    else:
        os.environ['BLASTDB'] = dbpath

    dbs = os.listdir(os.environ['BLASTDB'])
    if not "%s.phr"%db in dbs and not "%s.00.phr"%db in dbs:
        sys.stderr.write(('Warning: BLAST annotation was not possible, %s not '
                          'found in %s!\n')%(db,dbs))
        return {}

    # assert ' ' not in reference, 'BLAST cannot handle spaces in reference path!'
    defaults = { # Name: (type, arg, default_value)
        'database':            (str,   '-db',                 'swissprot'),
        'evalue':              (float, '-evalue',                    10.0),
        'word_size':           (int,   '-word_size',                  6  ),
        'task':                (str,   '-task',             'blastx-fast'),
        'best_hit_overhang':   (float, '-best_hit_overhang',          0.4),
        'best_hit_score_edge': (float, '-best_hit_score_edge',        0.4),
        'num_threads':         (int,   '-num_threads',                4  )
    }
    if blast_settings is None: blast_settings = {}
    for s, v in blast_settings.items():
        assert s in defaults, 'unknown setting encountered (%s)!'%s
        assert isinstance(v, defaults[s][0]), ('invalid type for %s, expected %s,'
                                               ' found %s!')%(s, defaults[s][0],
                                                              type(v))

    # BLAST fasta to DB
    cmd = ['blastx', '-query', fasta, '-outfmt',
           ('7 qseqid qstart qend qlen sseqid sstart send slen length pident '
            'nident mismatch gaps stitle')]
    for s, (t, a, v) in defaults.items():
        if s in blast_settings: cmd.extend([a, str(blast_settings[s])])
        elif v is not None: cmd.extend([a, str(v)])

    # print(' '.join(cmd))
    name = os.path.basename(fasta).rsplit('.',1)[0]
    so_path = f'{work_dir}{name}.blastx.tsv'
    se_path = f'{work_dir}{name}.blastx.err'
    with open_(so_path, 'w') as so, open_(se_path, 'w') as se:
        ec = Popen(cmd, stdout=so, stderr=se).wait()
        if ec != 0:
            sys.stderr.write(str(cmd))
            raise RuntimeError('BLASTx failed during execution')
    # Extract annotations
    annotations = extract_annotations(so_path, ignore_terms, replace_terms,
                                      max_cov)

    # Post process annotations (remove redundancy)
    for query in annotations:
        for i, (start, end, annots) in enumerate(annotations[query]):
            new_annots = reduce_redundancy(annots, sim_cutoff)
            # Repeat to avoid similar parallel evolved cases
            new_annots = reduce_redundancy(new_annots, sim_cutoff)
            annotations[query][i] = (start, end, new_annots)

    return annotations
