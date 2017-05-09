#!/usr/bin/env python3
''' RUCS - Rapid Identification of PCR Primers Pairs for Unique Core Sequences


Copyright (c) 2017, Martin Christen Frolund Thomsen, Technical University of Denmark
All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

'''
import sys, os, time, re, gzip, json, types, shutil, glob, bisect
import primer3
import numpy as np
from subprocess import call, Popen, PIPE, STDOUT as subprocess_stdout
from difflib import SequenceMatcher
from tabulate import tabulate

# GLOBALS
log = None
settings = None

################################################################################
################################################################################
################################################################################
# Program Entry Points
def main(positives, negatives, ref_input=None, kmer_size=None, quiet=False,
         clean_run=True, annotate=True, settings_file=None, name=None):
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
   log = LogObj(quiet)
   if settings_file is not None:
      load_global_settings(settings_file)
   
   if kmer_size is None: kmer_size = settings['ucs']['kmer_size']
   
   # Validate Input
   if positives is None or not positives:
      raise UserWarning('No Positive genomes provided!')
   if negatives is None:
      negatives = []
   
   # Create reference directory to store reference links, and BWA index files
   ref_dir = 'references'
   if not os.path.exists(ref_dir): os.mkdir(ref_dir)
   try:
      stats_file = '%sstats.log'%("%s_"%name if name is not None else '')
      pairs_json_file = '%spairs.json'%("%s_"%name if name is not None else '')
      results_file = '%sresults.tsv'%("%s_"%name if name is not None else '')
      products_file = '%sproducts.tsv'%("%s_"%name if name is not None else '')
      results_file_best = '%sresults_best.tsv'%("%s_"%name if name is not None else '')
      
      log.progress.add('main', 'Running RUCS', None)
      
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
      reference = create_symbolic_files([reference], ref_dir)[0]
      positives = create_symbolic_files(positives, ref_dir)
      negatives = create_symbolic_files(negatives, ref_dir)
      log.progress['input'].log_time()
      
      # Find unique core sequences
      cs_f, ucs_files = find_unique_core_sequences(positives, negatives,
                                                   reference, kmer_size)
      
      # Identify Primer Pairs
      contig_file = ucs_files[3]
      
      pairs, good_pp = find_validated_primer_pairs(contig_file, positives,
                                                   negatives,
                                                   annotate=annotate)
      
      # Check how many primer pairs have been tested (only tested pp in results)
      no_tested_pp = 0
      for p in pairs:
         if 'test' in p:
            no_tested_pp += 1
      
      # Check if any pairs were found
      if len(pairs) == 0:
         raise UserWarning('No valid PCR primer pairs could be found!')
      
      # Create json file containing all pairs
      with open(pairs_json_file, 'w') as f:
         json.dump(pairs[:no_tested_pp], f)
      
      # Create test summary file
      with open(products_file, 'w') as f:
         rp = range(len(positives))
         rn = range(len(negatives))
         f.write('#%s\t%s\n'%('\t'.join(('P_%s'%(x+1) for x in rp)),
                              '\t'.join(('N_%s'%(x+1) for x in rn))))
         f.write('\n'.join(('\t'.join(
            (','.join(map(str, sorted(x, key=int))) for x in
            (p['products']['pos'] + p['products']['neg'] if 'test' in p else [])
            )) for p in pairs[:no_tested_pp])))
      
      # Create tab separated summary file of good pairs
      with open(results_file_best, 'w') as f:
         f.write(present_pairs_full(good_pp))
      # Create tab separated summary file of all pairs
      with open(results_file, 'w') as f:
         f.write(present_pairs_full(pairs[:no_tested_pp]))
   except UserWarning as msg:
      # Clean up (to reduce space usage)
      clean_up(ref_dir, clean_run)
      sys.stderr.write("%s\n"%msg)
   except:
      # Clean up (to reduce space usage)
      log.progress['main'].log_time()
      log.progress.summary()
      log.stats.summary()
      clean_up(ref_dir, clean_run)
      raise
   else:
      # Clean up (to reduce space usage)
      clean_up(ref_dir, clean_run)
   finally:
      log.progress['main'].log_time()
      # Store time and stats in debug.log
      with open(stats_file, 'w') as f:
         log.progress.summary(f)
         log.stats.summary(f)

def find_primer_pairs(contig_file, positives, negatives, contig_names=None,
                      quiet=False, clean_run=True, annotate=True,
                      settings_file=None, name=None):
   ''' This methods allows you to run Primer Identification for a one or more
       fasta entries '''
   # Set Globals
   global log, settings
   log = LogObj(quiet)
   if settings_file is not None:
      load_global_settings(settings_file)
   
   try:
      stats_file = '%sstats.log'%("%s_"%name if name is not None else '')
      pairs_json_file = '%spairs.json'%("%s_"%name if name is not None else '')
      products_file = '%sproducts.tsv'%("%s_"%name if name is not None else '')
      results_file = '%sresults.tsv'%("%s_"%name if name is not None else '')
      results_file_best = '%sresults_best.tsv'%("%s_"%name if name is not None else '')
      
      log.progress.add('main', 'Running RUCS', None)
      
      # Create reference directory to store reference links, and BWA index files
      ref_dir = 'references'
      if not os.path.exists(ref_dir): os.mkdir(ref_dir)
      # Create symlinks for all reference
      log.progress.add('input',
                       'Prepare inputs: %s positive and %s negative genomes'%(
                       len(positives), len(negatives)),
                       'main')
      positives = create_symbolic_files(positives, ref_dir)
      negatives = create_symbolic_files(negatives, ref_dir)
      log.progress['input'].log_time()
      
      if contig_names is not None and isinstance(contig_names, str):
         contig_names = [contig_names]
      pairs, good_pp = find_validated_primer_pairs(contig_file, positives,
                                                   negatives,
                                                   contig_names=contig_names,
                                                   annotate=annotate)
      
      # Check how many primer pairs have been tested (only tested pp in results)
      no_tested_pp = 0
      for p in pairs:
         if 'test' in p:
            no_tested_pp += 1
      
      # Check if any pairs were found
      if len(pairs) == 0:
         raise UserWarning('No valid PCR primer pairs could be found!')
      
      # Create json file containing all pairs
      with open(pairs_json_file, 'w') as f:
         json.dump(pairs[:no_tested_pp], f)
      
      # Create test summary file
      with open(products_file, 'w') as f:
         rp = range(len(positives))
         rn = range(len(negatives))
         f.write('#%s\t%s\n'%('\t'.join(('P_%s'%(x+1) for x in rp)),
                              '\t'.join(('N_%s'%(x+1) for x in rn))))
         f.write('\n'.join(('\t'.join(
            (','.join(map(str, x)) for x in
            (p['products']['pos'] + p['products']['neg'] if 'test' in p else [])
            )) for p in pairs[:no_tested_pp])))
      
      # Create tab separated summary file of good pairs
      with open(results_file_best, 'w') as f:
         f.write(present_pairs_full(good_pp))
      # Create tab separated summary file of results
      with open(results_file, 'w') as f:
         f.write(present_pairs_full(pairs[:no_tested_pp]))
   except UserWarning as msg:
      # Clean up (to reduce space usage)
      clean_up(ref_dir, clean_run)
      sys.stderr.write("%s\n"%msg)
   except:
      # Clean up (to reduce space usage)
      log.progress['main'].log_time()
      log.progress.summary()
      log.stats.summary()
      clean_up(ref_dir, clean_run)
      raise
   else:
      # Clean up (to reduce space usage)
      clean_up(ref_dir, clean_run)
   finally:
      log.progress['main'].log_time()
      # Store time and stats in debug.log
      with open(stats_file, 'w') as f:
         log.progress.summary(f)
         log.stats.summary(f)

def virtual_pcr(references, pairs, output='products.tsv'):
   ''' Virtual PCR
   
   USAGE
      >>> refs = ['/full/path/to/ref1.fa', ...]
      >>> pairs = [['forward_primer1_seq', 'reverse_primer1_seq'], ...]
      >>> virtual_pcr(refs, pairs, output='my_pcr_results.tsv')
   '''
   # Create symlinks for all reference
   ref_dir = 'references'
   if not os.path.exists(ref_dir): os.mkdir(ref_dir)
   
   refs = create_symbolic_files(references, ref_dir, reuse=True)
   
   predict_pcr_results(refs, pairs, output=output)

def show_primer_probe_locs(wdir, contigs=None, name=None):
   ''' Show the location of primer bindingsites (red) and probebinding sites
   (green) using shell text coloring.
   '''
   seq_file    = '%s/unique_core_sequences.disscafs.fa'%(wdir)
   results_file = '%sresults.tsv'%("%s_"%name if name is not None else '')
   products_file = '%sproducts.tsv'%("%s_"%name if name is not None else '')
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
   p3_args = settings['pcr']['priming']['primer3']
   dna_conc = p3_args['PRIMER_DNA_CONC'] if 'PRIMER_DNA_CONC' in p3_args else 50.0
   dntp_conc = p3_args['PRIMER_DNTP_CONC'] if 'PRIMER_DNTP_CONC' in p3_args else 0.6
   dv_conc = p3_args['PRIMER_SALT_DIVALENT'] if 'PRIMER_SALT_DIVALENT' in p3_args else 1.5
   mv_conc = p3_args['PRIMER_SALT_MONOVALENT'] if 'PRIMER_SALT_MONOVALENT' in p3_args else 50.0
   salt_correction_method = p3_args['PRIMER_SALT_CORRECTIONS'] if 'PRIMER_SALT_CORRECTIONS' in p3_args else 1
   tm_method = p3_args['PRIMER_TM_FORMULA'] if 'PRIMER_TM_FORMULA' in p3_args else 1
   temp_c = p3_args['PRIMER_OPT_TM'] if 'PRIMER_OPT_TM' in p3_args else 60
   thal_type = 1      # type of thermodynamic alignment {1:any, 2:end1, 3:end2, 4:hairpin} (No effect)
   max_loop = 0       # Maximum size of loops in the structure.
   max_nn_length = 60 # Maximum length for nearest-neighbor calcs
   temponly = 1       # Return melting temperature of predicted structure
   dimer = 1          # if non-zero dimer structure is calculated (No effect)
   # NOTE: max_loop must be zero, to avoid Primer3 crashing
   
   # Initiate the thermodynamic analyses
   ThermoAnalysis = primer3.thermoanalysis.ThermoAnalysis
   p3_therm = ThermoAnalysis(thal_type, mv_conc, dv_conc, dntp_conc, dna_conc,
                              temp_c, max_loop, temponly, dimer, max_nn_length,
                              tm_method, salt_correction_method)
   tm_seq1 = round_sig(p3_therm.calcTm(seq1))
   if seq2 is not None:
      tm_seq2 = round_sig(p3_therm.calcTm(seq2))
      tm_dimer = round_sig(p3_therm.calcHeterodimer(seq1, seq2).tm)
   else:
      tm_seq2 = tm_seq1
      tm_dimer = round_sig(p3_therm.calcHeterodimer(seq1,
                                                    reverse_complement(seq1)).tm)
   
   return tm_seq1, tm_seq2, tm_dimer

def show_pcr_stats(forward, reverse, probe=None, template=None,
                   settings_file=None, title='PCR Stat Analysis'):
   '''
   USAGE:
          settings_file = '/Users/mcft/Work/CGE/repositories/primerfinder/settings.default.cjson'
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
      load_global_settings(settings_file)
   
   p3_args = settings['pcr']['priming']['primer3']
   dna_conc = p3_args['PRIMER_DNA_CONC'] if 'PRIMER_DNA_CONC' in p3_args else 50.0
   dntp_conc = p3_args['PRIMER_DNTP_CONC'] if 'PRIMER_DNTP_CONC' in p3_args else 0.6
   dv_conc = p3_args['PRIMER_SALT_DIVALENT'] if 'PRIMER_SALT_DIVALENT' in p3_args else 1.5
   mv_conc = p3_args['PRIMER_SALT_MONOVALENT'] if 'PRIMER_SALT_MONOVALENT' in p3_args else 50.0
   salt_correction_method = p3_args['PRIMER_SALT_CORRECTIONS'] if 'PRIMER_SALT_CORRECTIONS' in p3_args else 1
   tm_method = p3_args['PRIMER_TM_FORMULA'] if 'PRIMER_TM_FORMULA' in p3_args else 1
   temp_c = p3_args['PRIMER_OPT_TM'] if 'PRIMER_OPT_TM' in p3_args else 60
   thal_type = 1      # type of thermodynamic alignment {1:any, 2:end1, 3:end2, 4:hairpin} (No effect)
   max_loop = 0       # Maximum size of loops in the structure.
   # NOTE: max_loop must be zero, to avoid Primer3 crashing
   max_nn_length = 60 # Maximum length for nearest-neighbor calcs
   temponly = 1       # Return melting temperature of predicted structure
   dimer = 1          # if non-zero dimer structure is calculated (No effect)
   tm_threshold = settings['pcr']['priming']['threshold_tm'] if 'threshold_tm' in p3_args else 47
   max_probe_dist = settings['pcr']['priming']['max_probe_dist'] if 'max_probe_dist' in p3_args else 30
   min_primer_tm = p3_args['PRIMER_MIN_TM'] if 'PRIMER_MIN_TM' in p3_args else 57
   max_primer_tm = p3_args['PRIMER_MAX_TM'] if 'PRIMER_MAX_TM' in p3_args else 62
   min_probe_tm = p3_args['PRIMER_INTERNAL_MIN_TM'] if 'PRIMER_INTERNAL_MIN_TM' in p3_args else 67
   max_probe_tm = p3_args['PRIMER_INTERNAL_MAX_TM'] if 'PRIMER_INTERNAL_MAX_TM' in p3_args else 72
   
   # Initiate the thermodynamic analyses
   ThermoAnalysis = primer3.thermoanalysis.ThermoAnalysis
   p3_therm = ThermoAnalysis(thal_type, mv_conc, dv_conc, dntp_conc, dna_conc,
                             temp_c, max_loop, temponly, dimer, max_nn_length,
                             tm_method, salt_correction_method)
   
   n = "\x1B[0m"
   h = "\x1B[1;4m"
   r = "\x1B[31m"
   g = "\x1B[32m"
   print("\n%s%s%s\n"%(h, title, n))
   # Check PCR product size
   if template is None:
      size = 'N/A'
   else:
      start = template.index(forward)
      end = template.index(reverse_complement(reverse)) + len(reverse)
      size = "%sbp"%(end - start)
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
   tm1 = neg2mask(round(p3_therm.calcTm(forward),1))
   tm2 = neg2mask(round(p3_therm.calcTm(reverse),1))
   tm3 = neg2mask(round(p3_therm.calcTm(probe),1)) if probe is not None else 'N/A'
   m1, m2, m3 = n, n, n
   if isinstance(tm1, float): m1 = g if tm1 > min_primer_tm and tm1 < max_primer_tm else r
   if isinstance(tm2, float): m2 = g if tm2 > min_primer_tm and tm2 < max_primer_tm else r
   if isinstance(tm3, float): m3 = g if tm3 > min_probe_tm and tm3 < max_probe_tm else r
   print("Melting Tm:  %s%7s  %s%7s  %s%5s%s"%(m1, tm1, m2, tm2, m3, tm3, n))
   
   # Homo dimer temperature
   tm1 = neg2mask(round(p3_therm.calcHomodimer(forward).tm,1))
   tm2 = neg2mask(round(p3_therm.calcHomodimer(reverse).tm,1))
   tm3 = neg2mask(round(p3_therm.calcHomodimer(probe).tm,1)) if probe is not None else 'N/A'
   m1, m2, m3 = n, n, n
   if isinstance(tm1, float): m1 = r if tm1 > tm_threshold else g
   if isinstance(tm2, float): m2 = r if tm2 > tm_threshold else g
   if isinstance(tm3, float): m3 = r if tm3 > tm_threshold else g
   print("Homo dimer:  %s%7s  %s%7s  %s%5s%s"%(m1, tm1, m2, tm2, m3, tm3, n))
   
   # Hairpin temperature
   tm1 = neg2mask(round(p3_therm.calcHairpin(forward).tm,1))
   tm2 = neg2mask(round(p3_therm.calcHairpin(reverse).tm,1))
   tm3 = neg2mask(round(p3_therm.calcHairpin(probe).tm,1)) if probe is not None else 'N/A'
   m1, m2, m3 = n, n, n
   if isinstance(tm1, float): m1 = r if tm1 > tm_threshold else g
   if isinstance(tm2, float): m2 = r if tm2 > tm_threshold else g
   if isinstance(tm3, float): m3 = r if tm3 > tm_threshold else g
   print("Hairpin:     %s%7s  %s%7s  %s%5s%s\n"%(m1, tm1, m2, tm2, m3, tm3, n))
   
   print("              Fw-Rv  Fw-Pr  Rv-Pr")
   # Hetero dimer temperatur
   tm1 = neg2mask(round(p3_therm.calcHeterodimer(forward, reverse).tm,1))
   tm2 = neg2mask(round(p3_therm.calcHeterodimer(forward, probe).tm,1)) if probe is not None else 'N/A'
   tm3 = neg2mask(round(p3_therm.calcHeterodimer(reverse, probe).tm,1)) if probe is not None else 'N/A'
   m1, m2, m3 = n, n, n
   if isinstance(tm1, float): m1 = r if tm1 > tm_threshold else g
   if isinstance(tm2, float): m2 = r if tm2 > tm_threshold else g
   if isinstance(tm3, float): m3 = r if tm3 > tm_threshold else g
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
   >>> with open('test.fa', 'w') as f: f.write(file_content)
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
      line = ''
      readline = f.readline
      addsegment = queryseqsegments.append
      for line in f:
         if len(line.strip()) == 0: continue
         #sys.stderr.write("%s\n"%line)
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
               line = readline()
               seq  = line.strip().split()[0]
               # SKIP SECOND HEADER LINE AND QUALITY SCORES
               line = readline()
               line = readline()
            except:
               break
            else:
               # YIELD SEQUENCE AND RESET
               if to_upper: seq = seq.upper()
               yield (seq, name, desc)
               seq, name, desc = '', '', ''
         
         elif len(fields[0])>0:
            # EXTRACT FASTA SEQUENCE
            addsegment(fields[0])
      
      # CHECK FOR LAST FASTA SEQUENCE
      if queryseqsegments != []:
         # YIELD SEQUENCE
         seq = ''.join(queryseqsegments)
         if to_upper: seq = seq.upper()
         yield (seq, name, desc)

# CLASSES
class DependencyError(EnvironmentError):
   '''raise this when there's an unsolved program dependency error'''
   pass

class file_buffer():
   ''' This class wraps around a file object and allows the lines of the file
   to be stored in the RAM seemlesly. Buffering the file in RAM is an advantage
   if you want to expedite the IOops as quickly as possible to reduce IO load.
   This also makes it possible reset the line pointer to any line number.
   
   USAGE
      >>> import sys
      >>> with open('test.txt', 'w') as f:
      ...    _ = f.write('1\n2\n3\n')
      ... 
      >>> with open('test.txt') as f_obj:
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

class LogObj(object):
   '''
   USAGE
      >>> import sys, time
      >>> log = LogObj()
      >>> log.progress.add('main', 'Main', None)
      >>> log.progress.add('p1', 'Progress 1', 'main')
      >>> log.progress.add('p2', 'Progress 2', 'main')
      >>> log.stats.add_table('t1', 'Table 1', ['head 1', 'head 2'])
      >>> log.stats.add_row('t1', [1, 2])
      >>> log.stats.add_row('t1', [3, 4])
      >>> time.sleep(1)
      >>> log.progress['p1'].log_time()
      >>> log.progress['main'].log_time()
      >>> log.progress.summary()
      Time Analysis
      Seconds	Process
      1	Main
      1	   Progress 1
      -	   Progress 2
      
      >>> log.stats.summary()
      Table 1
      head 1	head 2
      1	2
      3	4
      
   '''
   class TableObj(object):
      '''  '''
      def __init__(self, title, headers):
         self.title = title
         self.headers = headers
         self.rows = []
         self.cols = len(headers)
      def add_row(self, row):
         '''  '''
         lrow = len(row)
         if lrow != self.cols:
            raise ValueError(('Row length does not fit the required length '
                              'of %s!')%lrow)
         self.rows.append(row)
   
   class StatObj(object):
      '''  '''
      def __init__(self):
         self.tables = {}
      def add_table(self, name, title, headers):
         '''  '''
         if name in self.tables:
            raise ValueError('Table name "%s" is in use!'%name)
         self.tables[name] = LogObj.TableObj(title, headers)
      def add_row(self, name, row):
         '''  '''
         if not name in self.tables:
            raise ValueError('Table name "%s" is unknown!'%name)
         else:
            self.tables[name].add_row(row)
      def summary(self, file_obj=sys.stdout):
         '''  '''
         for t in sorted(self.tables.values(), key=lambda x: x.title):
            file_obj.write(text_table(t.title, t.headers, t.rows))
   
   class ProgressObj(object):
      '''  '''
      def __init__(self, id_, name, msg, parent, lvl):
         self.id = id_
         self.timer = None
         self.time = time.time()
         self.name = name
         self.msg = msg
         self.parent = parent
         self.lvl = lvl
      def log_time(self):
         self.timer = time.time() - self.time
   
   class ProgressesObj(object):
      '''  '''
      def __init__(self, quiet=False):
         self.progresses = {}
         self.quiet = quiet
      def __getitem__(self, key):
         obj = self.progresses.get(key, None)
         if obj is None:
            raise KeyError('progress %s has not been initiated!'%key)
         return obj
      def add(self, name, msg, parent=None):
         '''  '''
         # Validate input
         if not isinstance(msg, str):
            raise ValueError('Invalid message, only strings allowed!')
         if not isinstance(name, str):
            raise ValueError('Invalid name, only strings allowed!')
         else:
            name = name.strip()
            if len(name) == 0:
               raise ValueError('Name "%s" is too short!'%name)
         if name in self.progresses:
            raise ValueError('Progress name "%s" is in use!'%name)
         if parent is not None and parent not in self.progresses:
            raise ValueError('Parent progress "%s" is unknown'%parent)
         # Find progress lvl
         lvl = 0
         p = parent
         while p is not None:
            lvl += 1
            p = self.progresses.get(p).parent
         
         # Print progress message to stdout
         if not self.quiet: sys.stdout.write('# %s%s...\n'%('   '*lvl, msg))
         
         # Add progress to progresses
         id_ = len(self.progresses)
         self.progresses[name] = LogObj.ProgressObj(id_, name, msg, parent, lvl)
      def summary(self, file_obj=sys.stdout):
         '''  '''
         progresses =sorted(self.progresses.values(), key=lambda p: float(p.id))
         rows =  [(round(p.timer, 1) if p.timer is not None else None,
                   "* %s%s"%('   '*p.lvl, p.msg)) for p in progresses]
         file_obj.write(text_table('Time Analysis', ['Seconds', 'Process'],
                                   rows))
   
   def __init__(self, quiet=False):
      self.progress = self.ProgressesObj(quiet)
      self.stats = self.StatObj()


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

class RegexTermReplacementObj():
   ''' For a given string, replace all keys from a given dictionary to their
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

# FUNCTIONS
def load_commented_json(cjson):
   ''' Load a json file containing comments '''
   return json.loads(decomment_json(cjson))

def decomment_json(cjson):
   ''' Load a json file containing comments '''
   def strip_comments(line):
      return line.split('//')[0].rstrip().replace('\\','\\\\')
   
   with open(cjson) as f:
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

def open_(filename, mode=None, compresslevel=9):
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
      if mode is None: mode = 'rb'
      return closing(gzip.open(filename, mode, compresslevel))
   else:
      if mode is None: mode = 'r'
      return open(filename, mode)

def save_as_fasta(sequences, filename, add_unique_id=True):
   with open(filename, 'w') as f:
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

def create_symbolic_files(files, directory, reuse=False):
   ''' Create symbolic links to all the files in the directory and return list
   of symbolic files
   '''
   symbolic_files = []
   for fil in files:
      if not os.path.exists(fil):
         raise OSError('OSError: File "%s" does not exist.'%fil)
      name = os.path.basename(fil)
      sym_path = "%s/%s"%(directory, name)
      try:
         _ = os.stat(sym_path)
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

def find_unique_core_sequences(positives, negatives, reference, kmer_size=20):
   ''' Find unique core sequences
   
   A unique core sequence is a sequence only found in the positive genomes
   '''
   min_seq_len = settings['pcr']['priming']['primer3']['PRIMER_PRODUCT_SIZE_RANGE'][0]
   if log is not None:
      # Initialize k-mer and sequences statistics logging
      log.progress.add('ucs', 'Finding unique core sequences', 'main')
      log.stats.add_table('kmer', 'k-mer analyses', ['k-mers', 'Process'])
      log.stats.add_table('seqs', 'Sequence Analyses',
                          ['Fasta file', 'Sequences', 'Size in bases',
                           'Seqs >%s'%min_seq_len, 'Size >%s'%min_seq_len])
      # Analyse Reference Sequence
      counts = analyse_genome(reference, min_seq_len)
      log.stats.add_row('seqs', [os.path.basename(reference)] + counts)
   
   # Compute intersecting k-mers from positive references
   core_kmers = find_intersecting_kmers(positives, kmer_size,
                                        min_seq_len=settings['ucs']['min_seq_len_pos'])
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
      counts = analyse_genome(cs_files[1], min_seq_len)
      log.stats.add_row('seqs', [os.path.basename(cs_files[1])] + counts)
   
   # Compute complementing k-mers to negative references
   unique_core_kmers = find_complementing_kmers(core_kmers, negatives,
                                                kmer_size,
                                                min_seq_len=settings['ucs']['min_seq_len_neg'])
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
      counts = analyse_genome(ucs_files[1], min_seq_len)
      log.stats.add_row('seqs', [os.path.basename(ucs_files[1])] + counts)
      
      log.progress['ucs'].log_time()
   
   return cs_files, ucs_files

def analyse_genome(genome, min_seq_len=300):
   ''' Analyse the genome, and annotate sequence and base data. '''
   buffer = settings['input']['use_ram_buffer']
   # counts (seqs, bases, seqs >threshold, bases >threshold)
   counts = [0, 0, 0, 0]
   for seq, n, d in seqs_from_file(genome, use_ram_buffer=buffer):
      seqlen = len(seq)
      counts[0] += 1
      counts[1] += seqlen
      if seqlen >= min_seq_len:
         counts[2] += 1
         counts[3] += seqlen
   
   return counts

def find_intersecting_kmers(files, kmer_size=20, min_seq_len=500):
   '''
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
      if log is not None:
         log.progress.add('core%s'%i, 'Processing %s'%os.path.basename(genome),
                          'core')
      # EXTRACT K-MERS FROM SEQEUNCE DATA FROM INPUT FILES
      to_upper = settings['input']['to_upper']
      buffer = settings['input']['use_ram_buffer']
      kmers_i = extract_kmers_from_file(genome, i, kmer_size, '',
                                        settings['ucs']['min_kmer_count'],
                                        settings['ucs']['rev_comp'],
                                        min_seq_len, to_upper, buffer)
      if log is not None:
         log.progress['core%s'%i].log_time()
         log.stats.add_row('kmer',
                           [len(kmers_i),
                            'Extraction of k-mers from %s'%(
                              os.path.basename(genome))
                            ])
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
                            revcom=True, min_seq_len=500, to_upper=False,
                            buffer=False):
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
   # Extract K-mers from the sequences
   kmers = {}
   seqcount = 0
   file_type = check_file_type([filename])
   for i, (seq, name, desc) in enumerate(seqs_from_file(filename,
                                                        to_upper=to_upper,
                                                        use_ram_buffer=buffer)):
      if len(seq) < min_seq_len: continue # Skip small sequences
      # Extract kmers from sequence (GenomePrefix_SequencePrefix_KmerPosition) to 'kmers'
      extract_kmers(kmers, seq, "%s_%d"%(genome_prefix, i), kmer_size,
                   kmer_prefix)
      # Extract kmers from the reverse complement sequence
      if revcom:
         extract_kmers(kmers, reverse_complement(seq), "%s_%d"%(genome_prefix, i),
                      kmer_size, kmer_prefix)
      seqcount += 1
   # File type dependend filtering
   if file_type == 'fastq':
      # Filter K-mers with low occurences
      sys.stderr.write("# Filtering K-mers with low coverage...\n")
      sum_ = 0
      count_kmers = len(kmers)
      for kmer, count in kmers.items():
         if count < fastq_kmer_count_threshold:
            del kmers[kmer]
            sum_ += count
            #sys.stderr.write("K-mer (%s) removed due to low occurrence (%s)!\n"%(kmer, count))
      count_kmers -= len(kmers)
      sys.stderr.write(("%s K-mers filtered with an average coverage of %.4f!"
                        "\n")%(count_kmers, float(sum_)/count_kmers))
   return kmers

def extract_kmers(kmers, seq, position_prefix, size=20, kmer_prefix=''):
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
      position_prefix: A prefix which is added to the position note.
      size: The size of the generated kmers. (Sequences shorter than the K-mer
         size will not generate K-mers.)
      kmer_prefix: A prefix used to filter which K-mers are stored. (Only
         K-mers starting with this prefix will be added.)
   USAGE
      >>> kmers = {}
      >>> seq = "abcdefghijklmnopqrstuvwxyz"
      >>> extract_kmers(kmers, seq, 'alphabet', 20, '')
      >>> kmers
      {'abcdefghijklmnopqrst': ['alphabet_0'],
       'defghijklmnopqrstuvw': ['alphabet_3'],
       'bcdefghijklmnopqrstu': ['alphabet_1'],
       'efghijklmnopqrstuvwx': ['alphabet_4'],
       'cdefghijklmnopqrstuv': ['alphabet_2'],
       'fghijklmnopqrstuvwxy': ['alphabet_5'],
       'ghijklmnopqrstuvwxyz': ['alphabet_6']}
   '''
   # ITERATE THE SEQUENCE TO FIND K-MERS
   kmer_prefix_len = len(kmer_prefix)
   seqiter = range(len(seq)-size+1)
   for j in seqiter:
      kmer = seq[j:j+size]
      # CHECK IF K-MER STARTS WITH SPECIFIED PREFIX
      if kmer_prefix == kmer[:kmer_prefix_len]:
         # ADD K-MER
         if kmer in kmers: kmers[kmer] += 1
         else: kmers[kmer] = 1

def reverse_complement(seq):
   ''' Compute the reverse complementary DNA string.
   >>> print(reverse_complement('ATGCTTGCGATCGTGTCAGTGCGTATGCTAGCTTCTAGGTTCGA'))
   TCGAACCTAGAAGCTAGCATACGCACTGACACGATCGCAAGCAT
   '''
   return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

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
   if op in ['union', 'uni', 'u', 'U', 'or', '|']:
      op = 'union'
      collection = set(a) | set(b)
   elif op in ['intersection', 'int', 'i', 'I', 'and', '&']:
      op = 'intersection'
      collection = set(a) & set(b)
   elif op in ['symmetric difference', 'sym', 's', 'xor', '^']:
      op = 'symmetric difference'
      collection = set(a) ^ set(b)
   elif op in ['complement', 'difference', 'diff', 'd', 'not', '-']:
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
   
   # Create BLAST DB
   name = os.path.basename(reference).rsplit('.',1)[0]
   if not os.path.exists("%s.nin"%(reference)):
      cmd = ['makeblastdb', '-title', name, '-in', reference, '-parse_seqids']
      if 'dbtype' in blast_settings:
         cmd.extend(['-dbtype', blast_settings['dbtype']])
      else:
         cmd.extend(['-dbtype', defaults['dbtype'][-1]])
      # print(' '.join(cmd))
      so_path = 'blast_db.out.txt'
      se_path = 'blast_db.err.txt'
      with open(so_path, 'w') as so, open(se_path, 'w') as se:
         ec = Popen(cmd, stdout=so, stderr=se).wait()
         if ec != 0: raise RuntimeError('BLAST formatdb failed during execution')
   del defaults['dbtype']
   
   # BLAST fasta to DB
   cmd = ['blastn', '-db', reference, '-query', fasta,
          '-outfmt', '6 qseqid sseqid sstrand qstart qend sstart send qseq sseq'] # evalue
   for s, (t, a, v) in defaults.items():
      if s in blast_settings: cmd.extend([a, str(blast_settings[s])])
      elif v is not None: cmd.extend([a, str(v)])
   
   # print(' '.join(cmd))
   so_path = 'blastn.%s.tsv'%(name)
   se_path = 'blastn.err.txt'
   with open(so_path, 'w') as so, open(se_path, 'w') as se:
      ec = Popen(cmd, stdout=so, stderr=se).wait()
      if ec != 0:
         sys.stderr.write('Failing command:\n%s\n\n'%(' '.join(cmd)))
         raise RuntimeError('BLASTn failed during execution')
   
   # Extract alignments
   primers = dict((i, seq) for i,(seq, n, d) in enumerate(seqs_from_file(fasta, use_ram_buffer=buffer)))
   alignments = {}
   with open(so_path) as f:
      for l in f:
         l = l.strip()
         if l == '': continue
         if l[0] == '#': continue
         try:
            # Parse BLAST line
            flds = l.split('\t')
            query = int(flds[0])
            contig_name = flds[1]
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
   
   return alignments

def align_to_ref(reference, fastq, bwa_settings=None, output_prefix='aln',
                 ignore_flags=4):
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
      log.progress.add('aln',
                       "Aligning %s to %s"%(os.path.basename(fastq),
                                            os.path.basename(reference)),
                       'ucs')
   # Align sequences to reference
   bam = bwa(reference, fastq, None, bwa_settings, output_prefix)
   # Extract alignment details
   alignments = extract_alignment_details(bam, ignore_flags)
   
   if log is not None:
      log.progress['aln'].log_time()
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
      'max_mismatch':       (int, '-n',             6), # edit distance
      'max_mismatch_score': (int, '-M',            10), # edit distance
      'max_gaps':           (int, '-o',             1), # max number of gaps
      'max_gap_ext':        (int, '-e',             0), # max gap extensions
      'gap_penalty_open':   (int, '-O',            10), # gap penalty open
      'gap_penalty_ext':    (int, '-E',            10), # gap penalty extension
      'seed_length':        (int, '-l',            35), # Length of the seed
      'max_mismatch_seed':  (int, '-k',             0), # edit distance in seed
      'threads':            (int, '-t',             4)  # paralelisation using multi threads
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
      with open(so_path, 'w') as so, open(se_path, 'w') as se:
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
   with open(sai_file, 'w') as so, open(se_path, 'w') as se:
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
      with open(path.sam, 'w') as sam, open(path.bwalog, 'w') as se:
         # Run BWA
         p1 = Popen(cmd1, stderr=se, stdout=sam)
         ec1 = p1.wait()
         if ec1 != 0: raise RuntimeError('bwa failed during execution (return code: %s)\n%s'%(ec1, ' '.join(cmd1)))
      with open(path.bam, 'w') as bam, open(path.bwalog, 'a') as se:
         # Convert SAM to BAM output for disk space efficiency
         p2 = Popen(cmd2, stdout=bam, stderr=se)
         # Wait for the programs to finish
         ec2 = p2.wait()
         if ec2 != 0: raise RuntimeError('samtools failed during execution (return code: %s)\n%s'%(ec2, ' '.join(cmd2)))
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
      with open(sai_file_2, 'w') as so, open(se_path_2, 'w') as se:
         ec = Popen(cmd, stdout=so, stderr=se).wait()
         if ec != 0: raise RuntimeError('bwa aln failed during execution')
      
      # Align using BWA (paired end) $ bwa sampe <REFERENCE> <SAI_1> <SAI_2> <FASTQ_1> <FASTQ_2>      -P?
      cmd = [path.bwa, 'sampe', '-n', '99', '-o', '2000', ref, sai_file, sai_file_2, fastq, paired_fastq]
      # sys.stderr.write("#CMD=%s\n"%(' '.join(cmd)))
      with open(path.bam, 'w') as so, open(path.bwalog, 'w') as se:
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

def compute_consensus_sequences(kmers, reference, kmer_size=20,
                                charspace='nATGC', name_prefix='consensus',
                                buffer=False):
   ''' Method for computing scaffolds and contigs from a '''
   scaffolds_file = "%s.scaffolds.fa"%name_prefix
   contigs_file = "%s.contigs.fa"%name_prefix
   auxiliary_file="%s.aux.tsv"%name_prefix
   dissected_scafs_file = "%s.disscafs.fa"%name_prefix
   clen = len(charspace)
   
   # Compute scaffold length
   scaffold_lengths = dict((name, len(seq))
                           for seq, name, desc in seqs_from_file(reference, use_ram_buffer=buffer))
   
   # Divide k-mers in scaffolds
   kmer_dict = {}
   for kmer in kmers:
      for contig, cigar, strand, position in kmers[kmer]:
         if not contig in kmer_dict: kmer_dict[contig] = {}
         if not kmer in kmer_dict[contig]: kmer_dict[contig][kmer] = []
         kmer_dict[contig][kmer].append(position)
   
   # Compute scaffold consensuses
   with open(auxiliary_file, 'w') as f:
      _ = f.write('# Auxiliary details on the k-mer-based scaffolds\n')
      _ = f.write('# Base\tContig\tPosition\tDepth\tConfidence\n')
      scaffolds = []
      DNAbins = np.asarray(list(charspace))[np.newaxis].T
      for scaffold in kmer_dict:
         # Build scaffold matrix
         mat = np.zeros([scaffold_lengths[scaffold],clen], dtype=np.int)
         for kmer in kmer_dict[scaffold]:
            for p in kmer_dict[scaffold][kmer]:
               dlen = scaffold_lengths[scaffold] - p
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
   dissected_scaffolds = sorted(filter(
      lambda x: x != '',
      (y.strip(charspace[0]) for x in
       map(lambda x: x[0].split(charspace[0]*1000), scaffolds) for y in x)
      ), key=lambda x: len(x)-x.count('n'), reverse=True)
   
   # Store contigs as fasta (position annotation is 0-indexed)
   save_as_fasta(dissected_scaffolds, dissected_scafs_file)
   
   return (scaffolds_file, contigs_file, auxiliary_file, dissected_scafs_file)

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

def find_complementing_kmers(kmers, files, kmer_size=20, min_seq_len=500):
   ''' '''
   if log is not None:
      log.progress.add('pan',
                       'Computing k-mer complement to the negative genomes',
                       'ucs')
   
   # Loop through set B files
   for i, genome in enumerate(files):
      if log is not None:
         log.progress.add('pan%s'%i, 'Processing %s'%os.path.basename(genome),
                          'pan')
      # Extract K-mers from the B file
      to_upper = settings['input']['to_upper']
      buffer = settings['input']['use_ram_buffer']
      negative_kmers = extract_kmers_from_file(genome, i, kmer_size, '',
                                               settings['ucs']['min_kmer_count'],
                                               settings['ucs']['rev_comp'],
                                               min_seq_len, to_upper, buffer)
      # TODO: Filter low counts?
      
      # Remove K-mers from A which are found in B
      kmers = set_op('complement', kmers, negative_kmers, method='list')
      if log is not None:
         log.progress['pan%s'%i].log_time()
         log.stats.add_row('kmer',
                           [len(negative_kmers),
                            'Extraction k-mers from %s'%(
                              os.path.basename(genome))
                            ])
   
   if log is not None:
      log.progress['pan'].log_time()
      log.stats.add_row('kmer', [len(kmers), 'Computation of k-mer complement'])
   
   return kmers

############ FUNCTIONS (FOR PCR PRIMER IDENTIFICATION)
def load_global_settings(settings_file='settings.default.cjson'):
   ''' This method is run at compile time, and will initialise global
   dependencies '''
   global settings
   if os.path.exists(settings_file):
      sf = settings_file
   else:
      # Look for the settings file relative to the script directory
      sf = "%s/%s"%(os.path.dirname(os.path.realpath(__file__)),settings_file)
      if not os.path.exists(sf):
         raise UserWarning('Settings file not found! (%s)'%(settings_file))
   settings = load_commented_json(sf)

def find_validated_primer_pairs(contig_file, p_refs, n_refs,
                                contig_names=None, annotate=True,
                                top_x_only=None):
   ''' if you do not have scaffolds, the contigs can be passed here
   
   USAGE
   >>> contig_file = '/Users/mcft/Desktop/test/old/contigs_tem.fa'
   >>> p_refs = ['/Users/mcft/Desktop/test/old/genome_pos_1.fna']
   >>> n_refs = ['/Users/mcft/Desktop/test/old/genome_neg_1.fna']
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
   p3_args = settings['pcr']['priming']['primer3']
   qpcr = p3_args['PRIMER_PICK_INTERNAL_OLIGO'] == 1
   to_upper = settings['input']['to_upper']
   filter_grade = settings['pcr']['priming']['threshold_grade']
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
   primer_pairs = []
   good_pp = []
   too_short = 0
   skipped = 0
   ignored = 0
   no_pairs = 0
   for i, (seq, name, desc) in enumerate(seqs_from_file(contig_file,
                                                        to_upper=to_upper,
                                                        use_ram_buffer=buffer)):
      if contig_names is not None and name not in contig_names:
         # Ignore unspecified contigs
         ignored += 1
         continue
      if len(seq) < p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0]:
         # Ignore sequences which are too short for Primer3 settings
         too_short += 1
         continue
      if len(seq) - seq.count('n') < 50:
         # Ignore scarce sequences content
         too_short += 1
         continue
      if len(good_pp) >= settings['pcr']['no_good_pp']:
         # Ignore contig, since the desired number of good pairs is reached
         skipped += 1
         continue
      if top_x_only is not None and i >= top_x_only:
         # Ignore contig, since the desired number of contigs have been analysed
         skipped += 1
         continue
   # COMPUTE primer pairs for sequence
      if log is not None:
         log.stats.add_table('prims%s'%i,
                             'Primer3 primer scanning of contig %s'%name,
                             ['Right', 'Left', 'Probe', 'Category'])
         log.stats.add_table('pair%s'%i,
                             'Primer3 pair scanning of contig %s'%name,
                             ['Pairs', 'Category'])
         log.progress.add('scan%s'%i,
                          'Scanning contig %s for primer pairs'%name, 'pcr')
      else:
         print(' Scanning contig %s for primer pairs'%name)
      pairs, notes = compute_primer_pairs(seq)
      if log is not None:
         log.progress['scan%s'%i].log_time()
         # Process primer3 notes
         right = dict([x.rsplit(' ', 1) for x in notes['right'].split(', ')])
         left = dict([x.rsplit(' ', 1) for x in notes['left'].split(', ')])
         probe = dict([x.rsplit(' ', 1) for x in notes.get('internal', '').split(', ') if x != ''])
         pair = dict([x.rsplit(' ', 1) for x in notes['pair'].split(', ')])
         headers = ['considered', 'too many Ns', 'GC content failed',
                    'long poly-x seq', 'low tm', 'high tm',
                    'high hairpin stability', 'ok']
         headers_pair = ['considered', 'unacceptable product size',
                         'high any compl', 'high end compl',
                         'tm diff too large', 'no internal oligo', 'ok']
         # Warn about headers unaccounted for
         miss_head = ', '.join(h for h in left if not h in headers)
         if miss_head:
            print(('Warning: One or more forward primer headers were not '
                   'reported: %s')%miss_head)
         miss_head = ', '.join(h for h in right if not h in headers)
         if miss_head:
            print(('Warning: One or more reverse primer headers were not '
                   'reported: %s')%miss_head)
         miss_head = ', '.join(h for h in probe if not h in headers)
         if miss_head:
            print(('Warning: One or more probe headers were not '
                   'reported: %s')%miss_head)
         miss_head = ', '.join(h for h in pair if not h in headers_pair)
         if miss_head:
            print(('Warning: One or more pair headers were not '
                   'reported: %s')%miss_head)
         # Log Primer3 notes as statistics
         for h in headers:
            log.stats.add_row('prims%s'%i, [right.get(h, '-'), left.get(h, '-'), probe.get(h, '-'), h])
         
         for h in headers_pair:
            log.stats.add_row('pair%s'%i, [pair.get(h, '-'), h])
         
      if len(pairs) == 0:
         no_pairs += 1
         continue
      # EXTRACT primer pairs
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
      
      # Store k-mers as fastq
      primers_fa = 'primers.fa'
      save_as_fasta(primers, primers_fa)
      if log is not None:
         log.stats.add_table('val%s'%i,
                             'Primer Prediction Overview (%s) for Contig %s'%(
                                len(primers), name),
                             ['Aligned','Grade 1', 'Grade 2', 'Grade 3',
                              'Grade 4', 'Reference'])
      
      # Create primer-reference mapping matrix
      if log is not None:
         log.progress.add('pos%s'%i, 'Aligning primers to positive references',
                          'pcr')
      else:
         print(' Aligning primers to positive references')
      for j, ref in enumerate(p_refs): # Positive references
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
                     0, os.path.basename(ref)]
         # Filter primer matches which fail the thermodynamic test
         alignments = filter_primer_alignments(ref, alignments, probes)
         # Update primer-reference matrix
         for p in primers:
            primers[p]['pos'].append(alignments[p] if p in alignments else [])
         
         if log is not None:
            log.progress["%s_pos_%s"%(i, j)].log_time()
            # Update prim_log with grade data
            for v in alignments.values():
               max_g = 0
               for a in v:
                  #(contig_name, grade, strand, position)
                  if a[1] > max_g:
                     max_g = a[1]
               
               if max_g > 0:
                  prim_log[max_g] += 1
            
            log.stats.add_row('val%s'%i, prim_log)
      
      if log is not None:
         log.progress['pos%s'%i].log_time()
         log.progress.add('neg%s'%i, 'Aligning primers to negative references',
                          'pcr')
      else:
         print(' Aligning primers to negative references')
      for j, ref in enumerate(n_refs): # Negative references
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
                     0, os.path.basename(ref)]
         # Filter primer matches which fail the thermodynamic test
         alignments = filter_primer_alignments(ref, alignments, probes)
         # Update primer-reference matrix
         for p in primers:
            primers[p]['neg'].append(alignments[p] if p in alignments else [])
         
         if log is not None:
            log.progress["%s_neg_%s"%(i, j)].log_time()
            # Update prim_log with grade data
            for v in alignments.values():
               max_g = 0
               for a in v:
                  #(contig_name, grade, strand, position)
                  if a[1] > max_g:
                     max_g = a[1]
               
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
         penalty_pos = sum([estimate_primer_rank(aln, scheme='positive')
                         for aln in primers[p]['pos']])
         penalty_neg = sum([estimate_primer_rank(aln, scheme='negative')
                         for aln in primers[p]['neg']])
         primers[p]['penalty'] = (penalty_pos, penalty_neg)
      
      # Compute primer pair rank
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
         penalty_pos = primers[p_fw]['penalty'][0] + primers[p_rv]['penalty'][0]
         penalty_neg = primers[p_fw]['penalty'][1] + primers[p_rv]['penalty'][1]
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
      
      # Sort primer pairs according to their rank and penalty
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
      
      # Add evaluated pairs to the final list of primer pairs
      for p in pairs:
         if ('test' in p and
             'penalty' in p['test']):
            primer_pairs.append(p)
      
      if log is not None:
         log.progress.add('status_%s'%(name),
                          'Current good pp: %s'%(len(good_pp)), 'pcr')
      else:
         print(' Current good pp: %s'%(len(good_pp)))
   
   # Annotate PCR products with flanks (regions)
   if annotate:
      if log is not None:
         log.progress.add('annotate', 'Annotating PCR product environment', 'pcr')
      else:
         print(' Annotating PCR product environment')
      flank_size = settings['pcr']['annotation']['flank_size']
      regions = {}
      for p in primer_pairs:
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
      env_fasta = 'pcr_products_with_skirts.fa'
      save_as_fasta(sequences, env_fasta, add_unique_id=False)
      
      # Analyse PCR product environments
      annotations = get_blast_annotations(env_fasta,
         settings['pcr']['annotation']['blastx_settings'],
         settings['pcr']['annotation']['blast_db_path'])
      for p in primer_pairs:
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
            too_short, p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0]))
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
      return (t['sensitivity'], t['specificity'], u, 1/(1+t['noise']))
   else:
      return (0,0,0,0)

def compute_primer_pairs(query_sequence):
   ''' Use the Primer 3 software to identify primer pair candidates for the
   query sequence
   
   USAGE
      >>> seq = 'TTTCCCGGCTACTTTTGACCACAGCCTTGACGTTCTCAGCTTTGC'
      >>> results = compute_primer_pairs(seq)
      >>> print('%s\t%s\t%s'%(results[0]['pair']['product_size'],
      ...                     results[0]['left']['sequence'],
      ...                     results[0]['right']['sequence']))
      45	TTTCCCGGCTACTTTTGACCA	GCAAAGCTGAGAACGTCAAGG
   
   find more information in the manual:
   http://primer3.sourceforge.net/primer3_manual.htm
   '''
   # Get Parameters
   p3_args = settings['pcr']['priming']['primer3']
   p3_seq_args = settings['pcr']['priming']['primer3_seq']
   # Filter None items
   p3_args = {k:v for k,v in p3_args.items() if v is not None}
   p3_seq_args = {k:v for k,v in p3_seq_args.items() if v is not None}
   
   # Check if sequence size is valid
   if (isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0], int) and
       len(query_sequence) < p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0]):
      return []
   
   # Set Sequence
   p3_seq_args['SEQUENCE_TEMPLATE'] = query_sequence
   
   # Find Primer Pairs
   tem_primer3_results = primer3.bindings.designPrimers(p3_seq_args, p3_args)
   
   # Extract Primer Pairs
   tem_primer3_pairs, notes = primer3_parser(tem_primer3_results)
   
   return tem_primer3_pairs, notes

def primer3_parser(primer3_results):
   ''' Parse Primer3 designPrimers output, and sort it into a hierachical
   dictionary structure of primer pairs.
   
   This method return 2 outputs, the list of primer pairs and a dictionary with
   notes (the explanatory output from Primer3).
   
   Author: Martin CF Thomsen
   '''
   primer_pairs = {}
   notes = {}
   for k in primer3_results:
      if 'PRIMER_RIGHT' == k[:12]:
         key = 'right'
         tmp = k[13:].split('_', 1)
         if tmp[0].isdigit():
            id = int(tmp[0])
            if not id in primer_pairs:
               primer_pairs[id] = {'pair': {}, 'right': {}, 'left': {},
                                   'internal': {}}
            if len(tmp) > 1:
               key2 = tmp[1].lower()
               value = primer3_results[k]
               if isinstance(value, float):
                  value = round_sig(value, sig_fig=3)
               primer_pairs[id][key][key2] = value
            else:
               primer_pairs[id][key]['position'] = primer3_results[k][0]
               primer_pairs[id][key]['length'] = primer3_results[k][1]
         elif tmp[0] == 'EXPLAIN':
            notes[key] = primer3_results[k]
         elif tmp == ['NUM','RETURNED']: pass
         else:
            sys.stderr.write("%s\n"%k)
      elif 'PRIMER_LEFT' == k[:11]:
         key = 'left'
         tmp = k[12:].split('_', 1)
         if tmp[0].isdigit():
            id = int(tmp[0])
            if not id in primer_pairs:
               primer_pairs[id] = {'pair': {}, 'right': {}, 'left': {},
                                   'internal': {}}
            if len(tmp) > 1:
               key2 = tmp[1].lower()
               value = primer3_results[k]
               if isinstance(value, float):
                  value = round_sig(value, sig_fig=3)
               primer_pairs[id][key][key2] = value
            else:
               primer_pairs[id][key]['position'] = primer3_results[k][0]
               primer_pairs[id][key]['length'] = primer3_results[k][1]
         elif tmp[0] == 'EXPLAIN':
            notes[key] = primer3_results[k]
         elif tmp == ['NUM','RETURNED']: pass
         else:
            sys.stderr.write("%s\n"%k)
      elif 'PRIMER_PAIR' == k[:11]:
         key = 'pair'
         tmp = k[12:].split('_', 1)
         if tmp[0].isdigit():
            id = int(tmp[0])
            if not id in primer_pairs:
               primer_pairs[id] = {'pair': {}, 'right': {}, 'left': {},
                                   'internal': {}}
            if len(tmp) > 1:
               key2 = tmp[1].lower()
               value = primer3_results[k]
               if isinstance(value, float):
                  value = round_sig(value, sig_fig=3)
               primer_pairs[id][key][key2] = value
            else:
               sys.stderr.write("%s, %s\n"%(k, primer3_results[k]))
         elif tmp[0] == 'EXPLAIN':
            notes[key] = primer3_results[k]
         elif tmp == ['NUM','RETURNED']: pass
         else:
            sys.stderr.write("%s\n"%k)
      elif 'PRIMER_INTERNAL' == k[:15]:
         key = 'internal'
         tmp = k[16:].split('_', 1)
         if tmp[0].isdigit():
            id = int(tmp[0])
            if not id in primer_pairs:
               primer_pairs[id] = {'pair': {}, 'right': {}, 'left': {},
                                   'internal': {}}
            if len(tmp) > 1:
               key2 = tmp[1].lower()
               value = primer3_results[k]
               if isinstance(value, float):
                  value = round_sig(value, sig_fig=3)
               primer_pairs[id][key][key2] = value
            else:
               primer_pairs[id][key]['position'] = primer3_results[k][0]
               primer_pairs[id][key]['length'] = primer3_results[k][1]
         elif tmp[0] == 'EXPLAIN':
            notes[key] = primer3_results[k]
         elif tmp == ['NUM','RETURNED']: pass
         else:
            sys.stderr.write("%s, %s\n"%(k, tmp[0]))
      else:
         sys.stderr.write("%s\n"%k)
      
   return list(map(primer_pairs.get, sorted(primer_pairs.keys()))), notes

def round_sig(number, sig_fig=3):
   ''' Round the number to the specified number of significant figures
   
   USAGE
      >>> import numpy as np
      >>> for i in range(6): print(i+1, round_sig(123.456, sig_fig=i+1))
      ...
      1 100.0
      2 120.0
      3 123.0
      4 123.5
      5 123.46
      6 123.456
   '''
   if number == 0: return number
   else: return round(number, -int(np.log10(abs(number))) -1 + sig_fig)

def filter_primer_alignments(ref, alignments, probes=[]):
   ''' Remove primer alignments where the heterodimer thermodynamics do not
   meet the threshold criteria
   
   ref is a path to the reference file
   alignments is a dict of with primer seq as keys and a list of alignments as
      values EG. {'TGATAAGGCGATGGCAGATGA': ['contig-name_+_4351']}

   # Grades
   # 0 = no priming         alignment_similarity < 0.6
   # 1 = unlikely priming   Tm < threshold_tm
   # 2 = potential priming  Tm > threshold_tm
   # 3 = probable priming   Tm = temp_c +- viable_tm_diff
   # 4 = definite priming   Tm = temp_c +- good_tm_diff
   
   USAGE:
      >>> ref = 'old/reference_tem.fa'
      >>> alignments = {'CAACATTTTCGTGTCGCCCTT': ['reference_+_1043',
      ...                                         'reference_+_300']}
      >>> filter_primer_alignments(ref, alignments)
      {'CAACATTTTCGTGTCGCCCTT': ['reference_+_1043']}
   '''
   buffer = settings['input']['use_ram_buffer']
   p3_args = settings['pcr']['priming']['primer3']
   dna_conc = p3_args['PRIMER_DNA_CONC'] if 'PRIMER_DNA_CONC' in p3_args else 50.0
   dna_conc_probe = p3_args['PRIMER_INTERNAL_DNA_CONC'] if 'PRIMER_INTERNAL_DNA_CONC' in p3_args else 50.0
   dntp_conc = p3_args['PRIMER_DNTP_CONC'] if 'PRIMER_DNTP_CONC' in p3_args else 0.6
   dv_conc = p3_args['PRIMER_SALT_DIVALENT'] if 'PRIMER_SALT_DIVALENT' in p3_args else 1.5
   mv_conc = p3_args['PRIMER_SALT_MONOVALENT'] if 'PRIMER_SALT_MONOVALENT' in p3_args else 50.0
   salt_correction_method = p3_args['PRIMER_SALT_CORRECTIONS'] if 'PRIMER_SALT_CORRECTIONS' in p3_args else 1
   tm_method = p3_args['PRIMER_TM_FORMULA'] if 'PRIMER_TM_FORMULA' in p3_args else 1
   temp_c = p3_args['PRIMER_OPT_TM'] if 'PRIMER_OPT_TM' in p3_args else 60
   thal_type = 1      # type of thermodynamic alignment {1:any, 2:end1, 3:end2, 4:hairpin} (No effect)
   max_loop = 0       # Maximum size of loops in the structure.
   max_nn_length = 60 # Maximum length for nearest-neighbor calcs
   temponly = 1       # Return melting temperature of predicted structure
   dimer = 1          # if non-zero dimer structure is calculated (No effect)
   # NOTE: max_loop must be zero, to avoid Primer3 crashing
   
   # Initiate the thermodynamic analyses
   ThermoAnalysis = primer3.thermoanalysis.ThermoAnalysis
   p3_primer = ThermoAnalysis(thal_type, mv_conc, dv_conc, dntp_conc, dna_conc,
                              temp_c, max_loop, temponly, dimer, max_nn_length,
                              tm_method, salt_correction_method)
   p3_probe = ThermoAnalysis(thal_type, mv_conc, dv_conc, dntp_conc,
                             dna_conc_probe, temp_c, max_loop, temponly, dimer,
                             max_nn_length, tm_method, salt_correction_method)
   # EXTRACT contigs from reference
   contigs = dict((name, seq) for seq, name, desc in seqs_from_file(ref, use_ram_buffer=buffer))
   
   # VALIDATE primer alignments
   validated_alignments = {}
   for primer in alignments:
      validated_alignments[primer] = []
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
               tm = p3_probe.calcHeterodimer(seq, seq2).tm
            else:
               tm = p3_primer.calcHeterodimer(seq, seq2).tm
         
         # Evaluate thermodynamic results
         grade = 1 if tm > 0 else 0
         if primer in probes:
            target_tm = p3_args['PRIMER_INTERNAL_OPT_TM']
         else:
            target_tm = temp_c
         # primer grade scheme
         prim_settings = settings['pcr']['priming']
         if tm >= prim_settings['threshold_tm']: grade += 1
         if (    tm > target_tm - prim_settings['viable_tm_diff']
             and tm < target_tm + prim_settings['viable_tm_diff']):
            grade += 1
         if (    tm > target_tm - prim_settings['good_tm_diff']
             and tm < target_tm + prim_settings['good_tm_diff']):
            grade += 1
         if grade >= prim_settings['threshold_grade']:
            validated_alignments[primer].append((contig_name, grade, strand,
                                                 position))
         # else:
         #    print('BAD:', tm, primer, alignment, ref)
   
   return validated_alignments

def align_seqs(a, b=None):
   ''' Align 2 sequences, '-' are used for gaps.
   
   USAGE
      >>> from difflib import SequenceMatcher
      >>> seq1 = 'TGGCCGCAAATGCATTCCCCA'
      >>> seq2 = 'GGCAATGCATTATGTGGGCAT'
      >>> print("%s\n%s\nSimilarity score: %s"%align_seqs(seq1, seq2))
      TGGCCGCAAATGCATTCCC----CA-
      -GGC----AATGCATTATGTGGGCAT
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

def estimate_primer_rank(alignments, scheme='positive'):
   '''  '''
   penalty = settings["pcr"]["priming"]["penalties"][scheme]
   # Get grade summary
   grade_sum = {0:0,1:0,2:0,3:0,4:0}
   for aln in alignments:
      grade_sum[aln[1]] += 1
   
   # Compute primer penalty
   penalty_score = 0
   if scheme == 'positive':
      if grade_sum[4] == 0:
         penalty_score += penalty['no_grade_4']
         if grade_sum[3] == 0:
            penalty_score += penalty['no_grade_3']
         elif grade_sum[3] > 1:
            penalty_score += penalty['multi_grade_3']
         else:
            if grade_sum[2] > 0:
               penalty_score += penalty['multi_grade_2']
            elif grade_sum[1] > 0:
               penalty_score += penalty['multi_grade_1']
      elif grade_sum[4] > 1:
         penalty_score += penalty['multi_grade_4']
      else:
         if grade_sum[3] > 0:
            penalty_score += penalty['multi_grade_3']
         elif grade_sum[2] > 0:
            penalty_score += penalty['multi_grade_2']
         elif grade_sum[1] > 0:
            penalty_score += penalty['multi_grade_1']
   elif scheme == 'negative':
      if grade_sum[4] > 0:
         penalty_score += penalty['grade_4']
      elif grade_sum[3] > 0:
         penalty_score += penalty['grade_3']
      elif grade_sum[2] > 0:
         penalty_score += penalty['grade_2']
      elif grade_sum[1] > 0:
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
               (contig_name, (strand, (grade, pos)))
               for contig_name, grade, strand, pos in primers[p_fw]['pos'][i]))
            rv_locs = AdvancedDictionary((
               (contig_name, (strand, (grade, pos)))
               for contig_name, grade, strand, pos in primers[rc_p_rv]['pos'][i]))
            if qpcr:
               probe_locs = AdvancedDictionary((
                  (contig_name, (strand, (grade, pos)))
                  for contig_name, grade, strand, pos in primers[probe]['pos'][i]))
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
            if product[-1] > lb and product[-1] < hb:
               target_match = [product[-2][0] >= mg, product[-2][1] >= mg,
                               qpcr and product[-2][2] >= mg]
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
         
         # UPDATE pair with PCR product lengths
         product_lengths = [x[-1] for x in products]
         p['products']['pos'].append(product_lengths)
         
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
               (contig_name, (strand, (grade, pos)))
               for contig_name, grade, strand, pos in primers[p_fw]['neg'][i]))
            rv_locs = AdvancedDictionary((
               (contig_name, (strand, (grade, pos)))
               for contig_name, grade, strand, pos in primers[rc_p_rv]['neg'][i]))
            if qpcr:
               probe_locs = AdvancedDictionary((
                  (contig_name, (strand, (grade, pos)))
                  for contig_name, grade, strand, pos in primers[probe]['neg'][i]))
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
            # Identify binding primers/probe
            good_binders = np.zeros(3, dtype=np.bool)
            for j, g in enumerate(product[-2]):
               if g >= mg: good_binders[j] = True
            
            if good_binders.any():
               # Update primer/probe uniqueness flags
               remove_flags = np.logical_not(good_binders)
               current_flags = np.array(num2binarray(p['test']['unique_flags'], 3))
               p['test']['unique_flags'] = binarray2num(current_flags * remove_flags)
            
            # Filter products with bad grades
            if any(g < mg for g in product[-2]):
               therm_count += 1
            else:
               filtered_products.append(product)
         
         # Update products
         products = filtered_products
         match_count += len(products)
         
         # UPDATE pair with PCR product lengths
         product_lengths = [x[-1] for x in products]
         p['products']['neg'].append(product_lengths)
         
         # Compute specificity
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
         found = False
         for gp in good_pp:
            dist = (abs(p['left']['position']  - gp['left']['position']) +
                    abs(p['right']['position'] - gp['right']['position']))
            if dist < settings['pcr']['min_pair_distance']:
               found = True
               break
         if not found:
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

def binarray2num(arr):
   ''' Convert binary array to a number. '''
   return int(np.array(arr)[::-1].dot(1 << np.arange(len(arr) - 1, -1, -1)))

def num2binarray(num, min_len=0):
   ''' Convert number to binary array, set min_len if you want a minimum array
   length. '''
   return list(x == '1' for x in reversed("{0:b}".format(num).zfill(min_len)))

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
   '''
   p3_args = settings['pcr']['priming']['primer3']
   mg = settings['pcr']['min_pcr_grade']
   qpcr = p3_args['PRIMER_PICK_INTERNAL_OLIGO'] == 1
   products = []
   total_count = 0
   small_count = 0
   therm_count = 0
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
         for g_fw, p_fw in p_fw_list['+']:
            for g_rv, p_rv in p_rv_list['+']:
               # Get PCR grade and product size
               product_size = p_rv + rv_len - p_fw
               pcr_grades = [g_fw, g_rv]
               if qpcr:
                  if probe_plus is not None:
                     probes = probe_plus[np.logical_and(probe_plus[:,1]>(p_fw+fw_len), probe_plus[:,1]<(p_rv-probe_len))]
                  else:
                     probes = np.array([])
                  pcr_grades.append(0 if probes.shape[0] == 0 else probes[:,0].max())
               # Filter bad products
               if filter_on_grade and any(g < mg for g in pcr_grades):
                  therm_count += 1
                  continue
               # Filter small products
               if product_size < p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0]:
                  small_count += 1
                  continue
               # Filter big products
               if product_size > p3_args['PRIMER_PRODUCT_SIZE_RANGE'][1]:
                  break
               # Add match to alignment_list
               products.append([contig_name, '+', p_fw,
                                contig_name, '+', p_rv,
                                pcr_grades, product_size])
      
      # Find pairs on the - strand
      if '-' in p_fw_list and '-' in p_rv_list:
         p_fw_list.sort(key=lambda x: x[1], reverse=True)
         p_rv_list.sort(key=lambda x: x[1], reverse=True)
         for g_fw, p_fw in p_fw_list['-']:
            for g_rv, p_rv in p_rv_list['-']:
               product_size = p_fw + fw_len - p_rv
               pcr_grades = [g_fw, g_rv]
               if qpcr:
                  if probe_minus is not None:
                     probes = probe_minus[np.logical_and(probe_minus[:,1]>(p_rv+rv_len), probe_minus[:,1]<(p_fw-probe_len))]
                  else:
                     probes = np.array([])
                  pcr_grades.append(0 if probes.shape[0] == 0 else probes[:,0].max())
               # Filter bad products
               if filter_on_grade and any(g < mg for g in pcr_grades):
                  therm_count += 1
                  continue
               # Filter small products
               if product_size < p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0]:
                  small_count += 1
                  continue
               # Filter big products
               if product_size > p3_args['PRIMER_PRODUCT_SIZE_RANGE'][1]:
                  break
               # Add match to alignment_list
               products.append([contig_name, '-', p_fw,
                                contig_name, '-', p_rv,
                                pcr_grades, product_size])
   
   large_count = total_count - small_count - therm_count - len(products)
   return products, (total_count, small_count, large_count, therm_count)

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
      get_blast_annotations('unique_core_sequences.contigs.fa', dbpath='/Applications/cmdline/blastdb/')
   '''
   def extract_annotations(so_path, ignore_terms, replace_terms, max_cov):
      ''' Extract BLAST annotations from tsv file '''
      annotations = {}
      
      # Prepare a regular expression from the replace terms
      regex = RegexTermReplacementObj(replace_terms)
      
      with open(so_path) as f:
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
               # Ignore selected terms
               if any(term in annotation for term in ignore_terms):
                  continue
               # Replace selected terms
               annotation = regex.replace_terms(annotation).strip()
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
                     if not annotation in notes:
                        notes.append(annotation)
                     annotations[query][i] = (new_start, new_end, notes)
                     covered = True
                     break
               if not covered:
                  annotations[query].append((qstart, qend, [annotation]))
      
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
   
   # Get settings
   max_cov = settings['pcr']['annotation']['max_cov']
   sim_cutoff = settings['pcr']['annotation']['sim_cutoff']
   ignore_terms = settings['pcr']['annotation']['ignore_terms']
   replace_terms = settings['pcr']['annotation']['replace_terms']
   # Validate inputs
   if dbpath == '':
      if not 'BLASTDB' in os.environ:
         sys.stderr.write('Annotation not possible! No BLAST DB found.\n')
         return {}
   if not 'refseq_protein.00.psq' in os.listdir(dbpath if dbpath != '' else os.environ['BLASTDB']):
      sys.stderr.write('Annotation not possible! refseq_protein.00.psq not '
                       'found.\n')
      return {}
   
   # assert ' ' not in reference, 'BLAST cannot handle spaces in reference path!'
   defaults = { # Name: (type, arg, default_value)
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
   cmd = ['blastx', '-db', dbpath+'refseq_protein', '-query', fasta,
          '-outfmt', ('7 qseqid qstart qend qlen sseqid sstart send slen length'
                      ' pident nident mismatch gaps stitle')]
   for s, (t, a, v) in defaults.items():
      if s in blast_settings: cmd.extend([a, str(blast_settings[s])])
      elif v is not None: cmd.extend([a, str(v)])
   
   # print(' '.join(cmd))
   name = os.path.basename(fasta).rsplit('.',1)[0]
   so_path = '%s.blastx.tsv'%(name)
   se_path = '%s.blastx.err'%(name)
   with open(so_path, 'w') as so, open(se_path, 'w') as se:
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

def present_pairs_full(pairs):
   ''' Presents essential aspects of the primer pairs
   USAGE
      >>> print(present_pairs_full(pairs[:3]))
      # product_size	tm	forward primer	reverse primer	test score
      515	60.0	CAACATTTTCGTGTCGCCCTT	TCGTCGTTTGGTATGGCTTCA	0
      605	60.0	CAACATTTTCGTGTCGCCCTT	AATTGTTGCCGGGAAGCTAGA	0
      636	60.0	CAACATTTTCGTGTCGCCCTT	TGCAACTTTATCCGCCTCCAT	0
   '''
   def _repr(p, h):
      return '\t'.join(["%s"]*len(h))%(
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
   h = ['sequence_id', 'product_size', 'unique_flags', 'sensitivity',
        'specificity', 'noise', 'penalty', 'p3_penalty', 'forward_primer',
        'forward_tm', 'forward_length', 'forward_gc%', 'forward_position',
        'reverse_primer', 'reverse_tm', 'reverse_length', 'reverse_gc%',
        'reverse_position', 'probe', 'probe_tm', 'probe_length',
        'probe_gc%', 'probe_position', 'annotation']
   return "#%s\n%s"%('\t'.join(h), '\n'.join((_repr(p, h) for p in pairs)))

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

def clean_up(ref_dir, clean_run=True):
   ''' Clean up (to reduce space usage) '''
   # Remove reference directory
   shutil.rmtree(ref_dir)
   if clean_run:
      # Remove requested files
      for fil in find_files(*settings["clean_files"]):
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
   
   return files

def text_table(title, headers, rows, table_format='psql'):
   ''' Create text table
   
   USAGE:
      >>> from tabulate import tabulate
      >>> title = 'My Title'
      >>> headers = ['A','B']
      >>> rows = [[1,2],[3,4]]
      >>> print(text_table(title, headers, rows))
      +-----------+
      | My Title  |
      +-----+-----+
      |   A |   B |
      +=====+=====+
      |   1 |   2 |
      |   3 |   4 |
      +-----+-----+
   '''
   # Create table
   table = tabulate(rows, headers, tablefmt=table_format)
   # Prepare title injection
   width = len(table.split('\n')[0])
   tlen = len(title)
   if tlen + 4 > width:
      # Truncate oversized titles
      tlen = width - 4
      title = title[:tlen]
   spaces = width - 2 - tlen
   left_spacer = ' '*int(spaces / 2)
   right_spacer = ' '*(spaces - len(left_spacer))
   # Update table with title
   table = '\n'.join(['+%s+'%('-'*(width-2)),
                      '|%s%s%s|'%(left_spacer, title, right_spacer),
                      table, '\n'])
   return table

def predict_pcr_results(refs, pairs, output=None, fail_on_non_match=False):
   ''' Predict the results of a PCR
   
   USAGE
      >>> predict_pcr_results(refs, pairs)
      Forward primer has %s binding locations hereof %s are good binders.
      Reverse primer has %s binding locations hereof %s are good binders.
      Potential products sizes: %s
   '''
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
   
   # Store k-mers as fastq
   primers_fa = 'primers.fa'
   save_as_fasta(dict(primers), primers_fa)
   
   # Align primers to reference
   pcr_results = [['NA' for r in refs] for p in pairs]
   for i, ref in enumerate(refs):
      alignments = blast_to_ref(ref, primers_fa,
                                settings['pcr']['priming']['blastn_settings'],
                                buffer)
      # Filter primer matches which fail the thermodynamic test
      alignments = filter_primer_alignments(ref, alignments, probes)
      for j, pair in enumerate(pairs):
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
               (contig_name, (strand, (grade, pos)))
               for contig_name, grade, strand, pos in alignments[fw]))
            rv_locs = AdvancedDictionary((
               (contig_name, (strand, (grade, pos)))
               for contig_name, grade, strand, pos in alignments[rv_rc]))
            if probe:
               probe_locs = AdvancedDictionary((
                  (contig_name, (strand, (grade, pos)))
                  for contig_name, grade, strand, pos in alignments[probe]))
            else:
               probe_locs = AdvancedDictionary({})
         except:
            if fail_on_non_match:
               print("error: Primer matching failed for PCR prediction (%s)!"%(ref))
               raise
            else:
               # Skip pairs, where one of the primers do not match anything
               pcr_results[j][i] = ''
               continue
         
         # FIND potential pcr products
         products, counts = find_pcr_products(fw_locs, rv_locs, probe_locs,
                                              len(fw), len(rv), len(probe))
         total_count = counts[0]
         small_count = counts[1]
         large_count = counts[2]
         therm_count = counts[3]
         match_count = len(products)
         
         # UPDATE pair with PCR product lengths
         products = ','.join([str(x[-1]) for x in products])
         
         pcr_results[j][i] = products
   
   # Print PCR results
   if output is not None:
      if '/' not in output or os.path.exists(os.path.dirname(output)):
         with open(output, 'w') as f:
            f.write('\t'.join(map(os.path.basename, refs)))
            f.write('\n')
            f.write('\n'.join(['\t'.join(p) for p in pcr_results]))
            f.write('\n')
      else:
         print('\t'.join(map(os.path.basename, refs)))
         print('\n'.join(['\t'.join(p) for p in pcr_results]))
         raise OSError('Error: Output path directory do not exist!')
   else:
      print('\t'.join(map(os.path.basename, refs)))
      print('\n'.join(['\t'.join(p) for p in pcr_results]))

def generate_sequence_atlas_data_file(reference_file, core_sequence_file, unique_core_sequence_file, aux_file, file_name='results.json'):
   ''' Generate Sequence Atlas data file (json)
   
   >>> os.chdir('/Users/mcft/Desktop/test/mrc-1')
   >>> generate_sequence_atlas_data_file(
   ...    'reference.fa',
   ...    'core_sequences.contigs.fa',
   ...    'unique_core_sequences.contigs.fa',
   ...    'unique_core_sequences.aux.tsv')
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
   with open(aux_file) as f:
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
   with open(file_name, 'w') as f: json.dump(data, f)

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

def get_primer_and_probe_bindsites(results_file, unique_only=False):
   ''' Extract the sequences and the primer and probe location found in the
   sequences.
   '''
   count = 0
   pribind_locs = {}
   probind_locs = {}
   with open(results_file) as f:
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
         #    print('Ignored Forward: %s %s %s'%(tmp[10], tmp[12], unique_flags))
         
         if not unique_only or unique_flags[1]:
            rv_len = int(tmp[15].strip())
            rv_pos = int(tmp[17].strip())
            pribind_locs[sequence_id].add(rv_pos-rv_len+1,rv_pos)
         # else:
         #    print('Ignored Reverse: %s %s %s'%(tmp[10], tmp[12], unique_flags))
         
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

def find_ucs(positives, negatives, ref_input=None, kmer_size=None, quiet=False,
             clean_run=True, settings_file=None, name=None):
   ''' This script computes the core sequences of the positive genomes and
   removes sequences covered by any of the negative genomes, thus creating
   unique core sequences.
   '''
   # Set Globals
   global log
   log = LogObj(quiet)
   if settings_file is not None:
      load_global_settings(settings_file)
   
   # Validate Input
   if positives is None or not positives:
      raise UserWarning('No Positive genomes provided!')
   if negatives is None:
      negatives = []
   
   stats_file = '%sstats.log'%("%s_"%name if name is not None else '')
   
   if kmer_size is None:
      kmer_size = settings['ucs']['kmer_size']
   
   min_seq_len = settings['pcr']['priming']['primer3']['PRIMER_PRODUCT_SIZE_RANGE'][0]
   
   # Create reference directory to store reference links, and BWA index files
   ref_dir = 'references'
   if not os.path.exists(ref_dir): os.mkdir(ref_dir)
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
      reference = create_symbolic_files([reference], ref_dir)[0]
      positives = create_symbolic_files(positives, ref_dir)
      negatives = create_symbolic_files(negatives, ref_dir)
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
      clean_up(ref_dir, clean_run)
      sys.stderr.write("%s\n"%msg)
   except:
      # Clean up (to reduce space usage)
      log.progress['main'].log_time()
      log.progress.summary()
      log.stats.summary()
      clean_up(ref_dir, clean_run)
      raise
   else:
      # Clean up (to reduce space usage)
      clean_up(ref_dir, clean_run)
   finally:
      log.progress['main'].log_time()
      # Store time and stats in stats.log
      with open(stats_file, 'w') as f:
         log.progress.summary(f)
         log.stats.summary(f)

# Set entry Points Methods
def full(args):
   ''' Run full diagnostic: fucs + fppp '''
   main(args.positives, args.negatives, args.reference, quiet=True,
        clean_run=True, annotate=True)

def fucs(args):
   ''' Find Unique Core Sequences '''
   find_ucs(args.positives, args.negatives, args.reference, quiet=True,
            clean_run=True)

def fppp(args):
   ''' Find PCR Primer Pairs '''
   settings['pcr']['seq_selection'] = None
   find_primer_pairs(args.template, args.positives, args.negatives,
                     quiet=True, clean_run=True, annotate=True)

def anno(args):
   ''' Annotate sequences using the refseq BLAST DB '''
   annotations = get_blast_annotations(args.template,
      settings['pcr']['annotation']['blastx_settings'],
      settings['pcr']['annotation']['blast_db_path'])
   print(annotations)

def vpcr(args):
   ''' Virtual PCR - Simulate PCR and predict PCR product for the provided
   primer pairs against the provided references. '''
   pairs = get_pairs(args.pairs)
   virtual_pcr(args.references, pairs)

def pcrs(args):
   ''' Show PCR Statistics - Analyse the temperature and more for the provided
   primer pairs. '''
   # Get pairs
   pairs = get_pairs(args.pairs)
   to_upper = settings['input']['to_upper']
   buffer = settings['input']['use_ram_buffer']
   
   # Get first template entry
   if args.template is not None:
      try:
         for seq, n, d in seqs_from_file(args.template, to_upper=to_upper,
                                         use_ram_buffer=buffer):
            template = seq
            break
      except:
         sys.stderr.write('Template: %s\n'%args.template)
         sys.stderr.write('Warning: Template not valid!\n')
         template = None
   else:
      template = None
   
   # Show PCR stats for each pair
   for i, seqs in enumerate(pairs):
      forward = seqs[0]
      reverse = seqs[1]
      probe   = seqs[2] if len(seqs) == 3 else None
      show_pcr_stats(forward, reverse, probe, template,
                     title='PCR Stat Analysis for pair %s'%(i+1))

def test(args):
   ''' Virtual PCR - Simulate PCR and predict PCR product for the provided
   primer pairs against the provided references. '''
   ''
   test_dir = os.path.dirname(os.path.realpath(__file__))
   pos = ["%s/test/bla.fa"%(test_dir)]
   neg = ["%s/test/sul.fa"%(test_dir)]
   main(pos, neg, None, quiet=True, clean_run=False, annotate=True)

def get_pairs(pairs_file):
   ''' Extract the pairs from the file. '''
   pairs = []
   headers = []
   with open(pairs_file) as f:
      for l in f:
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

#
#################################### MAIN #####################################
#
# Initialise global dependencies required
load_global_settings()

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
   args = parser.parse_args()
   
   entry_points = ['full', 'fucs', 'fppp', 'vpcr', 'anno', 'pcrs', 'test']
   args.entry_point = args.entry_point[0]
   
   if args.entry_point not in entry_points:
      sys.stderr.write('Unknown entry point provided, for help use the -h '
                       'option!\n')
      sys.exit(1)
   
   # Load Default Settings
   load_global_settings(args.settings_file)
   
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
   if args.anneal_tm is not None:
      tm = float(args.anneal_tm)
      settings['pcr']['priming']['primer3']['PRIMER_MIN_TM'] = tm - 3
      settings['pcr']['priming']['primer3']['PRIMER_OPT_TM'] = tm
      settings['pcr']['priming']['primer3']['PRIMER_MAX_TM'] = tm + 3
      settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_MIN_TM'] = tm + 8
      settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_OPT_TM'] = tm + 10
      settings['pcr']['priming']['primer3']['PRIMER_INTERNAL_MAX_TM'] = tm + 12
   if args.max_3end_gc is not None:
      settings['pcr']['priming']['primer3']['PRIMER_MAX_END_GC'] = int(args.max_3end_gc)
   if args.product_size_min is not None and args.product_size_max is not None:
      pr = [int(args.product_size_min), int(args.product_size_max)]
      settings['pcr']['priming']['primer3']['PRIMER_PRODUCT_SIZE_RANGE'] = pr
   settings['pcr']['priming']['primer3']['PRIMER_PICK_INTERNAL_OLIGO'] = 1 if args.pick_probe else 0
   if args.annotation_evalue is not None:
      settings['pcr']['annotation']['blastx_settings']['evalue'] = float(args.annotation_evalue)
   
   # Handle wildcards in positives, negatives and references
   if args.positives is not None:
      fi = [os.path.basename(x) for path in args.positives for x in glob.glob(path) if check_file_type(x) != 'fasta']
      if fi:
         sys.stderr.write(('WARNING: The following positive files were ignored,'
                           ' due to not being fasta!\n%s\n\n')%('\n'.join(fi)))
      args.positives = [x for path in args.positives for x in glob.glob(path) if check_file_type(x) is 'fasta']
   if args.negatives is not None:
      fi = [os.path.basename(x) for path in args.negatives for x in glob.glob(path) if check_file_type(x) != 'fasta']
      if fi:
         sys.stderr.write(('WARNING: The following negative files were ignored,'
                           ' due to not being fasta!\n%s\n\n')%('\n'.join(fi)))
      args.negatives = [x for path in args.negatives for x in glob.glob(path) if check_file_type(x) is 'fasta']
   if args.references is not None:
      fi = [os.path.basename(x) for path in args.references for x in glob.glob(path) if check_file_type(x) != 'fasta']
      if fi:
         sys.stderr.write(('WARNING: The following references were ignored,'
                           ' due to not being fasta!\n%s\n\n')%('\n'.join(fi)))
      args.references = [x for path in args.references for x in glob.glob(path) if check_file_type(x) is 'fasta']
   
   # Run Service
   print('Running %s'%(args.entry_point))
   locals().get(args.entry_point)(args)
