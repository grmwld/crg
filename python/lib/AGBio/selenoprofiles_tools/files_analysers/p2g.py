#!/usr/bin/env python2.6
# -*- coding:utf-8 -*-

from __future__ import with_statement
import os
import sys
import subprocess
sys.path.append('/users/rg/mmariotti/libraries')
from MMlib import bbash
from AGBio.io.common import *


class P2G_Parser(object):
    '''Parse a p2g output file, using appropriate parsing method
    according to where it comes from (genewise or exonerate).
    '''
    def __init__(self, filename):
        self._filename = filename
        self.parse = self._source_prog()
        self.result = P2G_ParserResult()

    def _parse_genewise(self):
        #proc = subprocess.Popen(['parse_genewise.py', '-i', self._filename],
        #                        stdout=subprocess.PIPE,
        #                        shell=True)
        op = bbash(' '.join(['parse_genewise.py', '-i', self._filename]), True)
        #raw_input()
        #os.system(' '.join(['parse_genewise.py', '-i', self._filename]))
        output = [c.strip().split() for c in op.split(';')]
        self.result.query.parse(output[0], output[5][1])
        self.result.target.parse(output[1], output[5][2])
        self.result.frameshifts = output[2][1]
        self.result.orig_pos = output[3]
        self.result.score = output[4][1]
        if len(output[6]) == 1:
            self.result.stop_codons = 0
        else: self.result.stop_codons = output[6][1]

    def _parse_exonarate(self):
        pass

    def _source_prog(self):
        with open(self._filename, 'r') as iff:
            for line in iff:
                if line.strip().startswith('genewise output'):
                     return self._parse_genewise
        return self._parse_exonarate


class P2G_ParserResult(object):
    '''Class to hold results from parsing.
    '''
    def __init__(self):
        self.score = 0
        self.introns = 0
        self.frameshifts = 0
        self.query = self._result_seq()
        self.target = self._result_seq()
        self.orig_pos = []
        self.consensus = ''

    def __str__(self):
        return '\n'.join((str(self.query), str(self.target), self.score,
                          self.frameshifts, str(self.orig_pos)))

    class _result_seq(object):
        def __init__(self):
            self.name = ''
            self.start = 0
            self.end = 0
            self.sequence = ''
        def __str__(self):
            return '\n'.join((self.name,
                              self.start + ' --- ' + self.end,
                              formatText(self.sequence, 50, 4)))
        def parse(self, info, seq):
            self.name = info[0]
            self.start = info[1]
            self.end = info[2]
            self.sequence = seq
            
