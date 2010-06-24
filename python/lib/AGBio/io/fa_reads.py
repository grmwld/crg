#!/usr/bin/env python2.6

from __future__ import with_statement
import os
import sys
from AGBio.Utilities import count_N

def filter_N(chunks, max_N=0):
    for chunk in chunks:
        if count_N(chunk[1]) <= max_N and count_N(chunk[3]) <= max_N:
            yield chunk

def read_file(filename):
    for line in filename:
        if line.startswith('>'):
            chunk = [line.strip()]
            for i in range(3):
            chunk.append(filename.next().strip())
            yield chunk

def prints(chunks, output=sys.stdout):
    for chunk in chunks:
        output.write('\n'.join(chunk) + '\n')

## def load_paired_end(filename):
##     dataset = DataSet()
##     with open(filename) as iff:
##         for line in iff:
##             tmp1 = SingleRead()
##             tmp2 = SingleRead()
##             line = line.strip()
##             if line.startswith('>'):
##                 tmp1.header = line
##                 tmp1.sequence = iff.next().strip()
##                 tmp2.header = iff.next().strip()
##                 tmp2.sequence = iff.next().strip()
##                 dataset.add(PairedEndRead(tmp1, tmp2))
##     return dataset

## class SingleRead(object):
##     def __init__(self, header=None, sequence=None):
##         self.header = header
##         self.sequence = sequence

## class PairedEndRead(object):
##     def __init__(self, readA, readB):
##         if readA.header != readB.header:
##             raise Exception, 'The two reads do not match'
##         self.A = readA
##         self.B = readB
##         self.header = readA.header
##         self.sequence = [readA.sequence, readB.sequence]

## class DataSet(dict):
##     def __init__(self, dataset=None):
##         if dataset:
##             self.update(dataset)

##     def add(self, elem):
##         self[elem.header] = elem.sequence
