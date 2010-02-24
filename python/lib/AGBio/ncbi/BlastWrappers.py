#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import os
from AGBio.UtilityWrappers import *

BLAST_BIN_DIR = '/usr/local/ncbi/blast/bin/'

class BlastBaseWrapper(BaseUtilityWrapper, HasOutFile):
    """Base class for the NCBI blast tools.
    """
    def __init__(self, utilitypath='', outfile=None):
        '''
        '''
        BaseUtilityWrapper.__init__(self, utilitypath)
        HasOutFile.__init__(self, '-out ', outfile)


class BaseCmdWrapper(BlastBaseWrapper):
    '''
    '''
    def __init__(self, utilitypath=None, entry=None, info=False, db='nr', target_only=True,
                 outfile=None):
        '''
        '''
        BlastBaseWrapper.__init__(self, utilitypath, outfile)
        ## if entry:
        ##     self['entry'] = ['', ','.join(entry)]
        ## if info:
        ##     self['info'] = ['', '']
        ## if db:
        ##     self['db'] = ['', db]
        ## if target_only:
        ##     self['target_only'] = ['', '']
            
    @property
    def entry(self):
        return self['entry']

    def setentry(self, opt, value):
        if value:
            self['entry'] = [opt, ','.join(value)]
        elif 'entry' in self.keys():
            del self['entry']
        self.updateCline()
        
    @property
    def info(self):
        return self['info']
   
    def setinfo(self, opt, value):
        if value:
            self['info'] = [opt, '']
        elif 'info' in self.keys():
            del self['info']
        self.updateCline()
        
    @property
    def db(self):
        return self['db']
    
    def setdb(self, opt, value):
        if value:
            self['db'] = [opt, value]
        elif 'db' in self.keys():
            del self['db']
        self.updateCline()

    @property
    def target_only(self):
        return self['target_only']
   
    def settarget_only(self, opt, value):
        if value:
            self['target_only'] = [opt, '']
        elif 'target_only' in self.keys():
            del self['target_only']
        self.updateCline()

        
class FastaCmdWrapper(BaseCmdWrapper, HasOutFile):
    '''fastacmd wrapper.
    '''
    def __init__(self, entry=None, info=False, db='nr', target_only=True,
                 outfile=None):
        '''
        '''
        BaseCmdWrapper.__init__(self, 'fastacmd', entry, info, db, target_only,
                                outfile)
        HasOutFile.__init__(self, '-o ', outfile)
        if entry:
            print entry
            self.entry = entry
            print self.entry
        if info:
            self.info = ' '
        if db:
            self.db = db
        if target_only:
            self.target_only = ' '
        self.updateCline()

    @property
    def entry(self):
        return self['entry']
    @entry.setter   
    def entry(self, value):
        BaseCmdWrapper.setentry(self, '-s ', value)
        self.updateCline()

    @property
    def info(self):
        return self['info']
    @info.setter   
    def info(self, value):
        BaseCmdWrapper.setinfo(self, '-I ', value)
        self.updateCline()

    @property
    def db(self):
        return self['db']
    @db.setter
    def db(self, value):
        BaseCmdWrapper.setdb(self, '-d ', value)
        self.updateCline()

    @property
    def target_only(self):
        return self['target_only']
    @target_only.setter   
    def target_only(self, value):
        BaseCmdWrapper.settarget_only(self, '-t ', value)
        self.updateCline()


class BlastDbCmdWrapper(BlastBaseWrapper):
    '''blastdbcmd wrapper.
    '''
    def __init__(self, entry=None, info=False, db='nr', target_only=True,
                 outfile=None, version='plus'):
        '''
        '''
        BlastBaseWrapper.__init__(self, BLAST_BIN_DIR + 'blastdbcmd', outfile)
        if entry:
            self['entry'] = ['-entry ', ','.join(entry)]
        if info:
            self['info'] = ['-info', '']
        if db:
            self['db'] = ['-db ', db]
        if target_only:
            self['target_only'] = ['-target_only', '']
        if version == 'legacy':
            self._changeCmdOptLegacy()
        self.updateCline()

    @property
    def entry(self):
        return self['entry']
    @entry.setter   
    def entry(self, value):
        if value:
            self['entry'] = ['-entry ', ','.join(value)]
        elif 'entry' in self.keys():
            del self['entry']
        self.updateCline()
        
    @property
    def info(self):
        return self['info']
    @info.setter   
    def info(self, value):
        if value:
            self['info'] = ['-info', '']
        elif 'info' in self.keys():
            del self['info']
        self.updateCline()
        
    @property
    def db(self):
        return self['db']
    @db.setter
    def db(self, value):
        if value:
            self['db'] = ['-db ', value]
        elif 'db' in self.keys():
            del self['db']
        self.updateCline()

    @property
    def target_only(self):
        return self['target_only']
    @target_only.setter   
    def target_only(self, value):
        if value:
            self['target_only'] = ['-target_only', '']
        elif 'target_only' in self.keys():
            del self['target_only']
        self.updateCline()
