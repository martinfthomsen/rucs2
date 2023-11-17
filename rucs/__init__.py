# -*- coding: utf-8 -*-
"""
    Init file
    ~~~~~~~~~~~~~~~~~~~~~

    Contains global functionalities for the library.

    * load_global_settings - load customisable settings from the settings files
    * setup_directories - setup standard directories (work, references and results) used through the RUCS libraries
    * get_directories - Fetch the global defined standard directories for (work, references and results)
    * text_table - convert title, headers and a data matrix into a nice looking printable table
    * binarray2num - Convert a binary array (list of 1s and 0s) into an integer representing the array
    * num2binarray - Convert a representing a binary array into a list of 1s and 0s.
    * round_sig - round a long decimal number to a certain number of decimals to look nice for printing

    Initialise the global "settings" used through out the RUCS libraries

"""
# (c) 2023 Martin Thomsen

__all__ = ['os', 'np']

import os
import numpy as np
from tabulate import tabulate

from rucs.file import load_commented_json

def load_global_settings(settings_file=None):
    ''' This method is run at compile time, and will initialise global
    dependencies '''
    global settings
    sf = None
    if settings_file is None and 'settings' not in globals():
        # Set settings with default path
        settings_file = 'settings.default.cjson'
        if not os.path.exists(settings_file):
            # Look for the settings file relative to the script directory
            sf = "%s/%s"%(os.path.dirname(os.path.realpath(__file__)),settings_file)
        else:
            sf = settings_file
    elif os.path.exists(settings_file):
        sf = settings_file
    if sf is not None:
        if not os.path.exists(sf):
            raise UserWarning('Settings file not found! (%s)'%(settings_file))

        # (Re-)Load settings
        settings = load_commented_json(sf)

    return settings


def setup_directories(get_ref_dir=False):
    ''' Set the work, reference and result directory. Create any missing paths and return the paths '''
    global work_dir, ref_dir, result_dir

    work_dir = settings['input']['work_dir']
    if isinstance(work_dir, str) and work_dir != '':
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
        work_dir = work_dir+'/' if work_dir[-1] != '/' else work_dir
    else:
        work_dir = ''

    result_dir = settings['input']['result_dir']
    if isinstance(result_dir, str) and result_dir != '':
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)
        result_dir = result_dir+'/' if result_dir[-1] != '/' else result_dir
    else:
        result_dir = ''

    ref_dir = f"{work_dir}references"
    if get_ref_dir:
        if not os.path.exists(ref_dir):
            os.mkdir(ref_dir)
        ref_dir = ref_dir+'/' if ref_dir[-1] != '/' else ref_dir

    return (work_dir, ref_dir, result_dir) if get_ref_dir else (work_dir, result_dir)


def get_directories(get_ref_dir=False):
    ''' Set the work, reference and result directory. Create any missing paths and return the paths '''
    return (work_dir, ref_dir, result_dir) if get_ref_dir else (work_dir, result_dir)


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


def binarray2num(arr):
    ''' Convert binary array to a number. '''
    return int(np.array(arr)[::-1].dot(1 << np.arange(len(arr) - 1, -1, -1)))


def num2binarray(num, min_len=0):
    ''' Convert number to binary array, set min_len if you want a minimum array
    length. '''
    return list(x == '1' for x in reversed("{0:b}".format(num).zfill(min_len)))


def round_sig(number, sig_fig=3, lmin=-1.0e+300, lmax=1.0e+300):
    ''' Round the number to the specified number of significant figures

    USAGE
       >>> import numpy as np
       >>> round_sig(-1.7976931348623157e+308, lmin=-5000)
       -5000
       >>> round_sig(1.7976931348623157e+308, lmax=100)
       100
       >>> round_sig(0)
       0
       >>> for i in range(6): print(i+1, round_sig(123.456, sig_fig=i+1))
       ...
       1 100.0
       2 120.0
       3 123.0
       4 123.5
       5 123.46
       6 123.456
    '''
    if number == 0 or number == 1:
         return number
    elif number < lmin:
        return lmin
    elif number > lmax:
        return lmax
    else:
       return round(number, -int(np.log10(abs(number))) -1 + sig_fig)


# Initialise global dependencies required
settings = load_global_settings()
