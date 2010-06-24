#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
from AGBio.io.Fasta import *

def eval_coverage(ref, ali):
    with open(ref) as iff:
        rref = loadSequences(iff)
    lref = 0
    print rref
    for s in rref:
        lref += len(s.sequence)
    lal = 0
    with open(ali) as iff:
        inal = False
        for line in iff:
            if line.startswith('a score'):
                inal = True
            elif line.startswith('s ') and inal == True:
                inal = False
                lal += len(line.split()[6])
    return [lref, lal]
