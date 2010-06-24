#!/usr/bin/env python2.6

from __future__ import with_statement
import os
import sys
from AGBio.Utilities import count_N

def filter_N(chunks, max_N=0):
    for chunk in chunks:
        if count_N(chunk[1]) <= max_N and count_N(chunk[5]) <= max_N:
            yield chunk

def read_file(filename):
    for line in filename:
        if line.startswith('@'):
            chunk = [line.strip()]
            for i in range(7):
                chunk.append(filename.next().strip())
            yield chunk

def prints(chunks, output=sys.stdout):
    for chunk in chunks:
        output.write('\n'.join(chunk) + '\n')

def to_fasta(chunks):
    for chunk in chunks:
        tmp = ['>'+chunk[0][1:]]
        tmp.append(chunk[1])
        tmp.append('>'+chunk[4][1:])
        tmp.append(chunk[5])
        yield tmp

## def read_file(filename):
##     for line in filename:
##         if line.startswith('@'):
##             chunk = ['>'+line[1:].strip()]
##             chunk.append(filename.next().strip())
##             for i in range(2):
##                 filename.next()
##             chunk.append('>'+filename.next()[1:].strip())
##             chunk.append(filename.next().strip())
##             yield chunk
