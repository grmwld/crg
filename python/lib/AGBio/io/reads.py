#!/usr/bin/env python2.6

from __future__ import with_statement
import os
import sys

class BaseReadDataSet(object):
    '''Base class to handle reads datasets
    '''
    def __init__(self, filename):
        self.filename = filename

    def prints(self, chunks, output=sys.stdout):
        '''Prints the dataset
        '''
        for chunk in chunks:
            output.write('\n'.join(chunk) + '\n')

    def load(self, c_size, head_char):
        '''Core chunks loader.
        A chunk is the lines containing the information of
        two paired-end reads.
        8 lines if fastq (2 header, 2 sequences, 2 headers, 2 qualities)
        4 lines if fasta (2 headers, 2 sequences)
        The set of chunks are actually a generator.
        '''
        for line in self.filename:
            if line.startswith(head_char):
                chunk = [line.strip()]
                for i in range(c_size - 1):
                    chunk.append(self.filename.next().strip())
                yield chunk

    @staticmethod
    def filter_N(chunks, p1, p2, max_N=0):
        '''Filters chunks based on the composition of their sequences
        TODO : Make it more modular so that only the result of a function
        would be required to filter out the chunks
        '''
        for chunk in chunks:
            if count_N(chunk[p1]) <= max_N and count_N(chunk[p2]) <= max_N:
                yield chunk

    @staticmethod
    def trim(chunks, p1, p2, length):
        for chunk in chunks:
            tmp = chunk[:]
            for i in (p1, p2):
                tmp[i] = tmp[i][:-length]
            yield tmp


class FastA(BaseReadDataSet):
    '''Class to handle fasta datasets
    '''
    def __init__(self, filename=None):
        BaseReadDataSet.__init__(self, filename)

    def load(self):
        return BaseReadDataSet.load(self, 4, '>')

    @staticmethod
    def filter_N(chunks, max_N=0):
        return BaseReadDataSet.filter_N(chunks, 1, 3, max_N)

    @staticmethod
    def trim(chunks, length):
        return BaseReadDataSet.trim(chunks, 1, 3, length)


class FastQ(BaseReadDataSet):
    '''Class to handle fastq datasets
    '''
    def __init__(self, filename=None):
        BaseReadDataSet.__init__(self, filename)

    def load(self):
        return BaseReadDataSet.load(self, 8, '@')

    @staticmethod
    def filter_N(self, chunks, max_N=0):
        return BaseReadDataSet.filter_N(chunks, 1, 5, max_N)

    @staticmethod
    def trim(self, chunks, length):
        return BaseReadDataSet.trim(chunks, 1, 5, length)
        
    @staticmethod
    def to_fasta(chunks):
        '''Converts the current fastq dataset to fasta format.
        '''
        for chunk in chunks:
            tmp = ['>'+chunk[0][1:]]
            tmp.append(chunk[1])
            tmp.append('>'+chunk[4][1:])
            tmp.append(chunk[5])
            yield tmp


def count_N(sequence, let='N'):
    c = 0
    for i in sequence:
        if i == let:
            c += 1
    return c
