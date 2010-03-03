#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import re
from AGBio.Utilities import *
from AGBio.io.Fasta import *


class AliData(object):
    def __init__(self, filename=None):
        self.filename = filename
        self.ali = None
        self.u_redundant = None
        self.u_non_redundant = None

    def load_alignment(self, filename=None):
        if not filename:
            if not self.filename:
                raise ValueError, 'Unspecified filename'
            else:
                fname = self.filename
        else:
            fname = filename
        with open(fname, 'r') as iff:
            self.ali = Alignment(loadSequences(iff))

    def find_u_positions(self):
        if not self.ali:
            raise ValueError, 'No suitable alignment found'
        self.u_redundant = self.ali.findPositions(('U', 'C', '-'), True)
        self.u_non_redundant = self.ali.findPositions(('U', 'C', '-'), False)

    def print_stat(self, redundant):
        if redundant:
            self._print_u_redundant()
        else:
            self._print_u_non_redundant()

    def scatter_score(self, stats):
        pass

    def _print_u_redundant(self):
        '''Prints stats of redundant positions
        '''
        for i in self.u_redundant['U']:
            sys.stdout.write(str(i) + ' ')
            for j in self.u_redundant:
                try:
                    sys.stdout.write(str(j)+': ')
                    sys.stdout.write(str(len(self.u_redundant[j][i])) + ' ; ')
                except KeyError:
                    sys.stdout.write('0 ; ')
            sys.stdout.write('\n')

    def _print_u_non_redundant(self):
        '''Prints stats of non redundant positions
        '''
        for i in self.u_non_redundant['U']:
            sys.stdout.write(str(i) + ' ')
            for j in self.u_non_redundant:
                try:
                    sys.stdout.write(str(j)+': ')
                    sys.stdout.write(str(len(self.u_non_redundant[j][i])) + ' ; ')
                except KeyError:
                    sys.stdout.write('0 ; ')
            sys.stdout.write('\n')
