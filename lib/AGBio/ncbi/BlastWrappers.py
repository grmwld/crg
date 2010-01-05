#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import os
from AGBio.UtilityWrappers import *

BLAST_BIN_DIR = '/usr/local/ncbi/blast/bin/'

class BlastBaseWrapper(BaseUtilityWrapper, HasOutFile):
    """Base class for the NCBI blast tools.
    """
    def __init__(self, utilitypath, outfile=None):
        '''
        '''
        BaseUtilityWrapper.__init__(self, utilitypath)
        HasOutFile.__init__(self, '-out ', outfile)
        

class BlastDbCmdWrapper(BlastBaseWrapper):
    '''blastdbcmd wrapper.
    '''
    def __init__(self, entry=None, info=False, db='nr', target_only=True,
                 outfile=None):
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
        self['db'] = value

    @property
    def target_only(self):
        return self['target_only']
    @target_only.setter   
    def target_only(self, value):
        if value:
            self['target_only'] = ['target_only', '']
        elif 'target_only' in self.keys():
            del self['target_only']
        self.updateCline()


