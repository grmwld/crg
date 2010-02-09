#!/usr/bin/env python2.6
# -*- coding:utf-8 -*-

from __future__ import with_statement
import os
import sys
import subprocess
sys.path.append('/users/rg/agrimaldi/usr/lib/python2.5/site-packages')
sys.path.append('/users/rg/mmariotti/libraries')
from MMlib import bbash
from AGBio.io.common import *


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


class HasInFile25(dict):
    '''Interface providing infile input.
    '''
    def __init__(self, carg='', infile=None):
        '''Interface providing infile specification
        '''
        self._icarg = carg
        self.infile = infile

    def get_infile(self):
        return self.getClineOpt('infile')
    def set_infile(self, value):
        if value:
            self['infile'] = [self._icarg, value]
        elif 'infile' in self.keys():
            del self['infile']
        self.updateCline()

    infile = property(get_infile, set_infile)


class BSeciSearchWrapper(BaseUtilityWrapper, HasInFile25):
    '''Wrapper for the bSECISearch utility.
    '''
    def __init__(self, infile, outdir='.', strands='0', frames='0'):
        BaseUtilityWrapper.__init__(self, 'bsecisearch_custom.pl')
        HasInFile25.__init__(self, '-i ', infile)
        self.outdir = outdir
        self.strands = strands
        self.frames = frames
        self.updateCline()

    def get_outdir(self):
        return self.getClineOpt('outdir')
    def set_outdir(self, value):
        self['outdir'] = ['-o ', str(value)]
        self.updateCline()

    def get_strands(self):
        return self.getClineOpt('strands')
    def set_strands(self, value):
        self['strands'] = ['-s ', str(value)]
        self.updateCline()

    def get_frames(self):
        return self.getClineOpt('frames')
    def set_frames(self, value):
        self['frames'] = ['-s ', str(value)]
        self.updateCline()

    outdir = property(get_outdir, set_outdir)
    strands = property(get_strands, set_strands)
    frames = property(get_frames, set_frames)


if __name__ == '__main__':

    bsec = BSeciSearchWrapper('/users/rg/agrimaldi/Data/testbsec/testbsecis',
                              outdir='/users/rg/agrimaldi/Data/testbsec/asds/',
                              strands=1, frames=1)

    print bsec.cline

    bsec.frames = 0

    print bsec.cline

    bsec.run()
