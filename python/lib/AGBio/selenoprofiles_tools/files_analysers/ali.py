#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import re
from AGBio.Utilities import *
from AGBio.io.Fasta import *


class AliStat(dict):
    def __init__(self, other=None):
        dict.__init__(self, other)

    def prints(self, aminoacid, out=sys.stdout):
        for pos in self[aminoacid]:
            out.write(str(pos) + ' ')
            for aa in self:
                try:
                    out.write(str(aa)+':')
                    out.write(str(len(self[aa][pos])) + ';')
                except KeyError:
                    out.write('0;')
            out.write('\n')


class AliData(object):
    def __init__(self, filename=None):
        self.filename = filename
        self.ali = None
        self.u_redundant = None
        self.u_non_redundant = None
        self.load_alignment(self.filename)

    def load_alignment(self, filename=None):
        if not filename:
            if not self.filename:
                return
            else:
                fname = self.filename
        else:
            fname = filename
        with open(fname, 'r') as iff:
            self.ali = Alignment(loadSequences(iff))

    def find_x_positions(self, aas='UC-*'):
        if not self.ali:
            raise ValueError, 'No suitable alignment found'
        self.u_redundant = AliStat(self.ali.findPositions(aas,
                                                          True))
        self.u_non_redundant = AliStat(self.ali.findPositions(aas,
                                                              False))
