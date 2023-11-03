# -*- coding: utf-8 -*-
"""
    Primer3 Library
    ~~~~~~~~~~~~~~~~~~~~~

    Contains Primer3 functionalities.

"""
# (c) 2023 Martin Thomsen

import sys, pickle, primer3  # refactor pickle away

from rucs import settings, round_sig, get_directories
from rucs.file import open_


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
    # Get directories
    global work_dir, ref_dir, result_dir
    work_dir, ref_dir, result_dir = get_directories(get_ref_dir=True)

    # Get Parameters
    p3_args = settings['pcr']['priming']['primer3']
    p3_seq_args = settings['pcr']['priming']['primer3_seq']
    # Filter None items
    p3_args = {k:v for k,v in p3_args.items() if v is not None}
    p3_seq_args = {k:v for k,v in p3_seq_args.items() if v is not None}

    # Check if sequence size is valid
    if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'], list):
        if isinstance(p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0], list):
            product_length_limits = [lim for lims in p3_args['PRIMER_PRODUCT_SIZE_RANGE'] for lim in lims]
            min_seq_len = min(product_length_limits)
        else:
            min_seq_len = p3_args['PRIMER_PRODUCT_SIZE_RANGE'][0]
    else:
        exit('Settings error: invalid input format for PRIMER_PRODUCT_SIZE_RANGE, expected "list", got: ' + p3_args['PRIMER_PRODUCT_SIZE_RANGE'])

    if len(query_sequence) < min_seq_len:
        return []

    # Set Sequence
    p3_seq_args['SEQUENCE_TEMPLATE'] = query_sequence

    # Find Primer Pairs
    tem_primer3_results = primer3.bindings.design_primers(p3_seq_args, p3_args)

    # Extract Primer Pairs
    tem_primer3_pairs, notes = primer3_parser(tem_primer3_results)

    return tem_primer3_pairs, notes


def configure_p3_thermoanalysis():
    ''' Configure the P3 thermodynamic analysis and provide it as global variables '''
    if not "p3_primer" in globals():
        p3_args = settings['pcr']['priming']['primer3']
        global p3_primer, p3_probe
        # type of thermodynamic alignment {1:any, 2:end1, 3:end2, 4:hairpin} (No effect) # REMOVED P3_PY-2.0
        # thal_type = 1
        # Concentration of monovalent cations in mM (type: float, p3_default: 50.0)
        mv_conc = p3_args['PRIMER_SALT_MONOVALENT'] if 'PRIMER_SALT_MONOVALENT' in p3_args else 50.0
        # Concentration of divalent cations in mM (type: float, p3_default: 1.5)
        dv_conc = p3_args['PRIMER_SALT_DIVALENT'] if 'PRIMER_SALT_DIVALENT' in p3_args else 1.5
        # Concentration of dNTP in mM (type: float, p3_default: 0.6)
        dntp_conc = p3_args['PRIMER_DNTP_CONC'] if 'PRIMER_DNTP_CONC' in p3_args else 0.6
        # Concentration of DNA in mM (for primer and probe respectively) (type: float, p3_default: 50.0)
        dna_conc = p3_args['PRIMER_DNA_CONC'] if 'PRIMER_DNA_CONC' in p3_args else 50.0
        dna_conc_probe = p3_args['PRIMER_INTERNAL_DNA_CONC'] if 'PRIMER_INTERNAL_DNA_CONC' in p3_args else 50.0
        # TODO - Concentration of DMSO in percentage (type: float, p3_default: 0.0)
        dmso_conc = p3_args['PRIMER_DMSO_CONC'] if 'PRIMER_DMSO_CONC' in p3_args else 0.0
        # TODO - DMSO correction factor, default 0.6 (type: float, p3_default: 0.6)
        dmso_fact = p3_args['PRIMER_DMSO_FACTOR'] if 'PRIMER_DMSO_FACTOR' in p3_args else 0.6
        # TODO - Concentration of formamide in mol/l (type: float, p3_default: 0.8)
        formamide_conc = p3_args['PRIMER_FORMAMIDE_CONC'] if 'PRIMER_FORMAMIDE_CONC' in p3_args else 0.8
        # TODO - Actual annealing temperature of the PCR reaction in Celsius (type: float, p3_default: -10.0)
        annealing_temp_c = p3_args['PRIMER_ANNEAL'] if 'PRIMER_ANNEAL' in p3_args else -10.0
        # Simulation temperature from which hairpin structures (dG) will be calculated in Celsius (type: float, p3_default: 37.0)
        temp_c = p3_args['PRIMER_OPT_TM'] if 'PRIMER_OPT_TM' in p3_args else 60
        # Maximum hairpin loop size measured in number of base pairs. (type: int, p3_default: 30)
        max_loop = 0 # NOTE: max_loop must be zero, to avoid Primer3 crashing
        # Return melting temperature of predicted structure (type: int, p3_default: 0)
        temp_only = 1
        # if non-zero dimer structure is calculated (No effect) (type: float, p3_default: 0)
        debug = 0
        # Maximum length for nearest-neighbor calculations (type: int, p3_default: 60)
        max_nn_length = 60
        # Method used to calculate melting temperatures (Tm) (breslauer [0] or santalucia [1]) (type: int, p3_default: 1)
        tm_method = p3_args['PRIMER_TM_FORMULA'] if 'PRIMER_TM_FORMULA' in p3_args else 1
        # Method used for salt corrections applied to melting temperature calculations (schildkraut [0], santalucia [1], owczarzy [2]) (type: int, p3_default: 1)
        salt_correction_method = p3_args['PRIMER_SALT_CORRECTIONS'] if 'PRIMER_SALT_CORRECTIONS' in p3_args else 1

        # Other p3_args? What are these?
        # output_structure: bool = False
        # calc_type_wrapper = 'ANY'

        # Initiate the thermodynamic analyses
        ThermoAnalysis = primer3.thermoanalysis.ThermoAnalysis
        p3_primer = ThermoAnalysis(mv_conc=mv_conc, dv_conc=dv_conc,
            dntp_conc=dntp_conc, dna_conc=dna_conc, dmso_conc=dmso_conc,
            dmso_fact=dmso_fact, formamide_conc=formamide_conc,
            annealing_temp_c=annealing_temp_c, temp_c=temp_c, max_loop=max_loop,
            temp_only=temp_only, debug=debug, max_nn_length=max_nn_length,
            tm_method=tm_method, salt_correction_method=salt_correction_method)
        p3_probe  = ThermoAnalysis(mv_conc=mv_conc, dv_conc=dv_conc,
            dntp_conc=dntp_conc, dna_conc=dna_conc_probe, dmso_conc=dmso_conc,
            dmso_fact=dmso_fact, formamide_conc=formamide_conc,
            annealing_temp_c=annealing_temp_c, temp_c=temp_c, max_loop=max_loop,
            temp_only=temp_only, debug=debug, max_nn_length=max_nn_length,
            tm_method=tm_method, salt_correction_method=salt_correction_method)

    return p3_primer, p3_probe

def primer3_parser(primer3_results):
    ''' Parse Primer3 designPrimers output, and sort it into a hierachical
    dictionary structure of primer pairs.

    This method return 2 outputs, the list of primer pairs and a dictionary with
    notes (the explanatory output from Primer3).

    Author: Martin CF Thomsen
    '''
    primer_pairs = {}
    notes = {}
    unknown_key_error = False
    for k in primer3_results:
        if 'PRIMER_RIGHT' == k[:12]:
            key = 'right'
            tmp = k[13:].split('_', 1)
            if tmp[0] == '': # New primer3 2.0 key containing list format of everything
                pass
            elif tmp[0].isdigit():
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
            elif tmp == ['NUM','RETURNED']:
                pass
            else:
                unknown_key_error = True
                sys.stderr.write("Unknown key: '%s' (primer3_parser:PRIMER_RIGHT)\n"%k)
        elif 'PRIMER_LEFT' == k[:11]:
            key = 'left'
            tmp = k[12:].split('_', 1)
            if tmp[0] == '': # New primer3 2.0 key containing list format of everything
                pass
            elif tmp[0].isdigit():
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
            elif tmp == ['NUM','RETURNED']:
                pass
            else:
                unknown_key_error = True
                sys.stderr.write("Unknown key: '%s' (primer3_parser:PRIMER_LEFT)\n"%k)
        elif 'PRIMER_PAIR' == k[:11]:
            key = 'pair'
            tmp = k[12:].split('_', 1)
            if tmp[0] == '': # New primer3 2.0 key containing list format of everything
                pass
            elif tmp[0].isdigit():
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
                    unknown_key_error = True
                    sys.stderr.write(f"Unknown key: '{k}' '{primer3_results[k]}' (primer3_parser:PRIMER_PAIR)\n")
            elif tmp[0] == 'EXPLAIN':
                notes[key] = primer3_results[k]
            elif tmp == ['NUM','RETURNED']:
                pass
            else:
                unknown_key_error = True
                sys.stderr.write("Unknown key: '%s' (primer3_parser:PRIMER_PAIR)\n"%k)
        elif 'PRIMER_INTERNAL' == k[:15]:
            key = 'internal'
            tmp = k[16:].split('_', 1)
            if tmp[0] == '': # New primer3 2.0 key containing list format of everything
                pass
            elif tmp[0].isdigit():
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
            elif tmp == ['NUM','RETURNED']:
                pass
            else:
                unknown_key_error = True
                sys.stderr.write("Unknown key: '%s' (primer3_parser:PRIMER_INTERNAL)\n"%k)
        else:
            unknown_key_error = True
            sys.stderr.write("Unknown key: '%s' (primer3_parser)\n"%k)

    # On unknown key errors: dump input to file
    if unknown_key_error:
        with open_(f'{work_dir}primer3_results.pkl', 'wb') as f:
            pickle.dump(primer3_results, f)
    return list(map(primer_pairs.get, sorted(primer_pairs.keys()))), notes
