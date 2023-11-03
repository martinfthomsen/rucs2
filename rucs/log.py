# -*- coding: utf-8 -*-
"""
    Logging library
    ~~~~~~~~~~~~~~~~~~~~~

    Contains log handling functionalities.

"""
# (c) 2023 Martin Thomsen

import sys, time
from rucs import text_table

# CLASSES
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

    def set_quiet(self, quiet=False):
        ''' change the quiet setting for the progress logging '''
        self.progress.quiet = quiet


# Initialise global dependencies required
log = LogObj()
