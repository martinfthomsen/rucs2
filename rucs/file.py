# -*- coding: utf-8 -*-
"""
    File library
    ~~~~~~~~~~~~~~~~~~~~~

    Contains file handling functionalities.

"""
# (c) 2023 Martin Thomsen

import sys, os, gzip, json, shutil, glob
from subprocess import Popen, PIPE
from contextlib import closing

# CLASSES
class file_buffer():
    ''' This class wraps around a file object and allows the lines of the file
    to be stored in the RAM seemlesly. Buffering the file in RAM is an advantage
    if you want to expedite the IOops as quickly as possible to reduce IO load.
    This also makes it possible reset the line pointer to any line number.

    USAGE
        >>> import sys
        >>> with open_('test.txt', 'w') as f:
        ...    _ = f.write('1\n2\n3\n')
        ...
        >>> with open_('test.txt') as f_obj:
        ...    f = file_buffer(f_obj)
        ...    for l in f:
        ...       _ = sys.stdout.write(l)
        ...
        ...    f.seekline(2)
        ...    for l in f:
        ...       _ = sys.stdout.write(l)
        ...
        1
        2
        3
        3
    '''
    def __init__(self, file_obj, use_ram_buffer=True, pre_load_buffer=True):
        self.f = file_obj
        self.use_ram_buffer = use_ram_buffer
        self.pre_load_buffer = pre_load_buffer
        if use_ram_buffer:
            if pre_load_buffer:
                self.lines = self.f.readlines()
                self.ll = len(self.lines)
            else:
                self.lines = []
                self.ll = 0
            self.ln = 0
    def readline(self):
        ''' Return next line from file '''
        if self.use_ram_buffer:
            self.ln += 1
            if self.ln > self.ll:
                if self.pre_load_buffer: # EOF
                    return ''
                else:
                    # Get new line
                    l = self.f.readline()
                    if l: # Not EOF
                        self.lines.append(l)
                        self.ll += 1
                    return l
            else:
                # Return line
                return self.lines[self.ln-1]
        else:
            return self.f.readline()
    def readlines(self):
        ''' Return next line from file '''
        if self.use_ram_buffer:
            if not self.pre_load_buffer:
                # Read the remainder of the file
                self.lines.extend(self.f.readlines())
                self.ll = len(self.lines)
            ln = self.ln
            self.ln = self.ll
            return self.lines[ln:]
        else:
            return self.f.readlines()
    def __iter__(self):
        l = self.readline()
        while l:
            yield l
            l = self.readline()

        yield l
    def seekline(self, ln):
        ''' Go to requested line number in file '''
        if ln <= len(self.lines):
            self.ln = ln
        else:
            raise EOFError('seeked line number %s does not exist in %s!\n'%(
                ln, self.f.name))


def get_pairs(pairs_file):
    ''' Extract the pairs from the file. '''
    pairs = []
    headers = []
    with open_(pairs_file) as f:
        for l in f:
            l = l.strip()
            if l == '': continue # Skip empty rows
            if l.startswith('#'):
                try:
                    headers = l[1:].split('\t')
                    fidx = headers.index('forward_primer')
                    ridx = headers.index('reverse_primer')
                    try: pidx = headers.index('probe')
                    except: pidx = ''
                except:
                    sys.stderr.write(('Unable to locate forward_primer and/or '
                                      'reverse_primer column(s) in %s!')%pairs_file)
                    sys.exit(1)
            elif headers:
                d = l.split('\t')
                pairs.append([d[fidx].strip(), d[ridx].strip(),
                              d[pidx].strip() if d[pidx] else None])
            else:
                pairs.append(l.strip().split())

    return pairs

def parse_tsv(tsv_file, sep='\t', comment='#'):
    ''' Extract the fields from the tsv file.

    >>> my_tsv = parse_tsv('results_best.tsv')
    >>> my_tsv[0].keys()
    dict_keys(['sequence_id', 'product_size', 'unique_flags', 'sensitivity', 'specificity', 'noise', 'penalty', 'p3_penalty',  'forward_primer', 'forward_tm', 'forward_length', 'forward_gc%', 'forward_position', 'reverse_primer', 'reverse_tm', 'reverse_length',  'reverse_gc%', 'reverse_position', 'probe', 'probe_tm', 'probe_length', 'probe_gc%', 'probe_position', 'annotation'])
    '''
    tsv = []
    headers = []
    with open_(tsv_file) as f:
        for l in f:
            if l.startswith(comment):
                if not tsv: # skip comments after value have been loadet to the matrix
                    headers = l[1:].strip().split(sep)
            elif headers:
                tsv.append(dict(zip(headers, l.strip().split(sep))))

    return tsv


def get_fasta_files(inputs, subdir=''):
    ''' Expand truncated paths, download files for provided accession IDs and
    ignore non fasta files

    Test example:
    >>> import glob
    >>> with open('test.fa', 'w') as f: f.write(">test\nATGC\n")
    >>> get_fasta_files(['CP000672.1', 'test.*', 'not_a_thing', 'ARBW00000000'])
    WARNING: The following input path: "not_a_thing" was ignored, due to not being a valid file path or known accession id!
    ['CP000672.1_GCA_000016485.1_ASM1648v1_genomic.fna.gz', 'test.fa', 'ARBW00000000_GCA_000379905.1_ASM37990v1_genomic.fna.gz']
    '''
    paths = []
    if subdir:
        ref = f" from the {subdir} argument"
        dir = f'inputs/{subdir}/'
    else:
        ref = ""
        dir = 'inputs/'
    if not os.path.exists(dir):
        os.makedirs(dir)
    # Load previous inputfile if existing to dictionary
    try:
        prev_inputs = parse_tsv('inputs/inputs.tsv')
        headers = [x.lower() for x in prev_inputs[0].keys()]
        prev_inputs = {x['id']:x for x in prev_inputs}
    except (IndexError, FileNotFoundError) as e:
        prev_inputs = {}
        headers = ['id', 'filepath', 'accession', 'genus', 'species', 'strain']
    # Open or Make input list file
    with open('inputs/inputs.tsv', 'a') as f:
        if not prev_inputs:
            # Initiate input list
            f.write('#%s\n'%'\t'.join(headers))
        # Check input files
        for input in inputs:
            path_glob = glob.glob(input)
            # Catch non paths inputs
            if path_glob == []:
                # Assume fasta accession ID
                # Check if inputfile was processed previously
                if input in prev_inputs:
                    paths.append(prev_inputs[input]['filepath'])
                else:
                    # Try download file
                    with Popen(["download_genomes.sh", input], stdout=PIPE, stderr=PIPE, cwd=dir) as p:
                        stdout = p.stdout.read().decode('utf-8')
                        stderr = p.stderr.read().decode('utf-8')

                    sys.stderr.write(stdout)
                    sys.stderr.write(stderr)
                    path = glob.glob(f"{dir}{input}*")
                    path = path[0] if path != [] else ''

                    # Check if file was downloaded
                    if os.path.exists(path):
                        # Add path to paths
                        paths.append(path)
                        # Save input details for later use
                        values = ['' for h in headers]
                        values[headers.index('id')] = input
                        values[headers.index('filepath')] = path
                        values[headers.index('accession')] = input
                        f.write('\t'.join(values))
                        f.write('\n')
                    else:
                        sys.stderr.write(f'WARNING: The following input path{ref}: "{input}" '
                                         'was ignored, due to not being a valid file path or'
                                         ' known accession id!\n')
            else:
                paths.extend(path_glob)

    nonfasta = [os.path.basename(x) for x in paths if check_file_type(x) != 'fasta']
    if nonfasta:
        sys.stderr.write(f'WARNING: The following file paths{ref} was ignored, due'
                         f' to not being valid fasta format:\n{chr(10).join(nonfasta)}\n')

    return [path for path in paths if check_file_type(path) == 'fasta']


# ITERATORS
def seqs_from_file(filename, exit_on_err=False, to_upper=False,
                   use_ram_buffer=False):
    '''Extract sequences from a file

    Name:
        seqs_from_file
    Author:
        Martin CF Thomsen
    Date:
        18 Jul 2013
    Description:
        Iterator which extract sequence data from the input file
    Args:
        filename: string which contain a path to the input file
    Supported Formats:
        fasta, fastq

    USAGE:
    >>> import os, sys
    >>> # Create fasta test file
    >>> file_content = ('>head1 desc1\nthis_is_seq_1\n>head2 desc2\n'
                        'this_is_seq_2\n>head3 desc3\nthis_is_seq_3\n')
    >>> with open_('test.fa', 'w') as f: f.write(file_content)
    ...
    >>> # Parse and print the fasta file
    >>> for seq, name, desc in seqs_from_file('test.fa', use_ram_buffer=True):
    ...    print(">%s %s\n%s"%(name, desc, seq))
   ...
    >head1 desc1
    this_is_seq_1
    >head2 desc2
    this_is_seq_2
    >head3 desc3
    this_is_seq_3
    '''
    # VALIDATE INPUT
    if not isinstance(filename, str):
        msg = 'Filename has to be a string.'
        if exit_on_err: sys.exit('error: '+msg)
        else: raise IOError(msg)
    if not os.path.exists(filename):
        msg = 'File "%s" does not exist.'%filename
        if exit_on_err: sys.exit('error: '+msg)
        else: raise IOError(msg)
    if not check_file_type([filename]) in ['fasta', 'fastq']:
        msg = 'Invalid file format for %s, only fasta and fastq are supported.'%(
            filename)
        if exit_on_err: sys.exit('error: '+msg)
        else: raise IOError(msg)

    # EXTRACT DATA
    with open_(filename, "r") as f_obj:
        f = file_buffer(f_obj, use_ram_buffer, use_ram_buffer)
        queryseqsegments = []
        seq, name, desc = '', '', ''
        for line in f:
            # Skip empty lines
            if len(line.strip()) == 0: continue

            fields=line.strip().split()
            if line[0] == ">":
                # FASTA HEADER FOUND
                if queryseqsegments != []:
                    # YIELD SEQUENCE AND RESET
                    seq = ''.join(queryseqsegments)
                    if to_upper: seq = seq.upper()
                    yield (seq, name, desc)
                    seq, name, desc = '', '', ''
                    del queryseqsegments[:]
                name = fields[0][1:]
                desc = ' '.join(fields[1:])

            elif line[0] == "@":
                # FASTQ HEADER FOUND
                name = fields[0][1:]
                desc = ' '.join(fields[1:])
                try:
                    # EXTRACT FASTQ SEQUENCE
                    line = f.readline()
                    seq  = line.strip().split()[0]
                    # SKIP SECOND HEADER LINE AND QUALITY SCORES
                    line = f.readline()
                    line = f.readline()
                except:
                    break
                else:
                    # YIELD SEQUENCE AND RESET
                    if to_upper: seq = seq.upper()
                    yield (seq, name, desc)
                    seq, name, desc = '', '', ''

            elif len(fields[0])>0:
                # FASTA SEQUENCE PART FOUND
                # Add FASTA SEQUENCE as segment to queryseqsegments list
                queryseqsegments.append(fields[0])

        # CHECK FOR LAST FASTA SEQUENCE
        if queryseqsegments != []:
            # YIELD SEQUENCE
            seq = ''.join(queryseqsegments)
            if to_upper: seq = seq.upper()
            yield (seq, name, desc)

        del f


def load_commented_json(cjson):
    ''' Load a json file containing comments '''
    return json.loads(decomment_json(cjson))


def decomment_json(cjson):
    ''' Load a json file containing comments '''
    def strip_comments(line):
        return line.split('//')[0].rstrip().replace('\\','\\\\')

    with open_(cjson) as f:
        json_str = '\n'.join(map(strip_comments, f.readlines()))

    return json_str


def check_file_type(inputFiles):
    ''' Check whether the input files are in fasta format, reads format or
    other/mix formats.
    '''
    if isinstance(inputFiles, str): inputFiles = [inputFiles]
    all_are_fasta = True
    all_are_reads = True
    try:
        for file_ in inputFiles:
            with open_(file_, 'r') as f:
                fc = f.readline()[0]
                if fc != "@": all_are_reads = False
                if fc != ">": all_are_fasta = False

    except: return 'other'
    if all_are_fasta: return 'fasta'
    elif all_are_reads: return 'fastq'
    else: return 'other'

def open_(filename, mode='r', compresslevel=9):
    """Switch for both open() and gzip.open().

    Determines if the file is normal or gzipped by looking at the file
    extension.

    The filename argument is required; mode defaults to 'rb' for gzip and 'r'
    for normal and compresslevel defaults to 9 for gzip.

    >>> import gzip
    >>> from contextlib import closing
    >>> with open_(filename) as f:
    ...     f.read()
    """
    if filename[-3:] == '.gz':
        if mode == 'r': mode = 'rt'
        if mode == 'w': mode = 'wt'
        return closing(gzip.open(filename, mode, compresslevel))
    else:
        return open(filename, mode)

def save_as_fasta(sequences, filename, add_unique_id=True):
    with open_(filename, 'w') as f:
        for i, seq_obj in enumerate(sequences):
            if isinstance(seq_obj, str):
                #>INDEX
                f.write(">%s\n%s\n"%(i, seq_obj))
            elif isinstance(seq_obj, (tuple, list)) and len(seq_obj) == 3:
                seq, name, desc = seq_obj
                #>NODE_0 contig=1,position=641
                if add_unique_id:
                    f.write(">%s_%s %s\n%s\n"%(name, i, desc, seq))
                else:
                    f.write(">%s %s\n%s\n"%(name, desc, seq))


def save_as_fastq(seqs, file_, names=None):
    '''
    Stores the list of sequences in a fastq format on disk. Where the header is
    the sequence number according to the position in the list, and the
    qualities are set to 'b's.
    USAGE:
        >>> save_as_fastq(uA.keys(), 'kmers.fq')
    '''
    if names is None:
        names = list(map(str, range(len(seqs))))
    else:
        assert len(names) == len(seqs), 'the provided names does not fit with the number of sequences'
    try:
        with open_(file_, 'w') as o:
            for i, seq in enumerate(seqs):
                # write fastq line, with qualities as 'h's
                # h is the highest score on the Illumina phred scale encoding
                o.write('@%s\n%s\n+\n%s\n'%(names[i], seq, 'h'*len(seq)))
    except Exception as e:
        sys.stderr.write('Error: could not write sequences to fastq file!\n%s\n'%(e))
        raise


def create_symbolic_files(files, directory, reuse=False):
    ''' Create symbolic links to all the files in the directory and return list
    of symbolic files
    '''
    symbolic_files = []
    for fil in files:
        if not os.path.exists(fil):
            raise OSError('OSError: File "%s" does not exist.'%fil)
        name = os.path.basename(fil)
        sym_path = "%s/%s"%(directory.rstrip('/'), name)
        try:
            _ = os.lstat(sym_path)
        except FileNotFoundError:
            os.symlink(os.path.abspath(fil), sym_path)
        except OSError:
            raise
        else:
            if not reuse:
                raise UserWarning('UserError: Multiple references with identical names '
                                  '"%s" was found!'%(name))

        symbolic_files.append(sym_path)

    return symbolic_files


def which(program):
    ''' Method for finding the path to a program '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def load_file_list(file_):
    ''' Extracts file paths from a file, and stores them in a list.
    USAGE:
        >>> load_file_list('v10_files')
    '''
    files = None
    try:
        with open_(file_, 'r') as f:
            files = [l.strip() for l in f.readlines()]
    except Exception:
        sys.stderr.write('Error: could not extract files!\n')
    else:
        return files


def clean_up(ref_dir, clean_files=None):
    ''' Clean up (to reduce space usage) '''
    # Remove reference directory
    shutil.rmtree(ref_dir)
    if clean_files:
        # Remove requested files
        for fil in find_files(*clean_files):
            os.unlink(fil)


def find_files(*file_pattern_tuple):
    ''' Find files matching one or more patterns

    USAGE
        >>> import glob
        >>> files = find_files('*.pdf', '*.cpp')
    '''
    files = []
    for file_pattern in file_pattern_tuple:
        files.extend(glob.glob(file_pattern))

    return list(set(files))
