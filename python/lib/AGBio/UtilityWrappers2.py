#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import os

TEMP = '/home/agrimaldi/temp/'
SBIN = '/soft/bin/'
HBIN = '/users/rg/agrimaldi/usr/bin/'
MBIN = '/users/rg/mmariotti/bin/'
MSCRIPT = '/users/rg/mmariotti/scripts/'
BLAST = SBIN + 'blastall'
MAFFT = SBIN + 'mafft'
TRIMAL = HBIN + 'trimal_dev'
TCOFFEE = HBIN + 't_coffee'
ADDFULLHEADERS = 'add_detail_to_titles2.py'
FILTER = MBIN + 'fetch_seq.g'
PRESELENOPROFILES = HBIN + 'prepare_alignment_selenoprofiles.py'
SELENOPROFILES = HBIN + 'selenoprofiles.py'
BLASTPARSER = 'blaster_parser.g'
DOWNLOADER = HBIN + 'get_sequence_from_gi.pl'
BSECISEARCH = 'bsecisearch_custom.pl'

NOOPT = ['', '']


def prop(fcn):
    return property(**fcn())


class BaseUtilityWrapper(dict):
    '''Base class for wrappers.

    Wrappers ease the use of common bioinformatic external tools.
    '''
    def __init__(self, utilitypath):
        '''Constructor.

        @param utilitypath : The absolute path of the tool.
        '''
        self.utilitypath = utilitypath
        self.cline = ''

    def getClineOpt(self, option):
        '''Getter of commandline options.

        @param option : The option to get
        @return : The relevent commandline option.
        '''
        return ''.join(self[option])

    def updateCline(self):
        '''Updates the commandline based on the options.
        '''
        self.cline = self.utilitypath
        for opt in self.keys():
            self.cline += (' ' + self.getClineOpt(opt))

    def run(self):
        '''Runs the commandline built with the various options.
        '''
        os.system(self.cline)

        
class HasInFile(dict):
    '''Interface providing infile input.
    '''
    def __init__(self, carg='', infile=None):
        '''Interface providing infile specification
        '''
        self._icarg = carg
        self.infile = infile

    @prop
    def infile():
        def fget(self):
            return self.getClineOpt('infile')
        def fset(self, value):
            if value:
                self['infile'] = [self._icarg, value]
            elif 'infile' in self.keys():
                del self['infile']
            self.updateCline()
        return locals()


class HasOutFile(dict):
    '''Interface providing outfile specification.
    '''
    def __init__(self, carg='', outfile=None):
        '''
        '''
        self._ocarg = carg
        self.outfile = outfile
        
    @prop
    def outfile():
        def fget(self):
            return self.getClineOpt('outfile')    
        def fset(self, value):
            if value:
                self['outfile'] = [self._ocarg, value]
            elif 'outfile' in self.keys():
                del self['outfile']
            self.updateCline()
        return locals()

