#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import re

class SequenceList(list):
    '''Class representing a set of sequences
    '''
    def __init__(self, sequences=''):
        self.extend(sequences)

    def findPattern( sequences, pattern, mode='count' ):
        '''Finds all sequences in the list with sequence matching a pattern.
        '''
        PATTERN = re.compile( pattern )
        count = 0
        matching = SequenceList()
        for sequence in sequences:
            if re.search( PATTERN, sequence.sequence ):
                matching.append( sequence )
                count += 1
        if mode == 'full':
            return matching
        elif mode == 'header':
            return [sequence.header for sequence in matching]
        else:
            return count

    @property
    def headers(self):
        '''Returns a list containing all the headers of the sequences
        '''
        return [s.header for s in self]

    @property
    def sequences(self):
        '''Returns a list containing all the sequences of the sequences
        '''
        return [s.sequence for s in self]

    def uniques(self, method='header'):
        '''Returns a new SequenceList object with unique sequences.
        '''
        usequences = SequenceList()
        if 'header':
            for i, seq in enumerate(self):
                if seq.header not in usequences.headers:
                    usequences.append(seq)
        return usequences

    def findPositions(self, aaa=None, strict=False):
        '''Finds all positions of U in the alignment.
        '''
        aa = list(aaa)
#        aa.append('!')
        print aa
        output = {}
        tmppos = {}
        for i in aa:
            output[i] = {}
            tmppos[i] = []
        for seq in self:
            for i in aa:
                tmppos[i] = []
            for i, pos in enumerate(seq.sequence):
                if pos in aa:
                    tmppos[pos].append(i)
            for t in tmppos:
                if tuple(tmppos[t]) in output[t]:
#                    if tmppos[t]:
                        output[t][tuple(tmppos[t])].append(seq)
#                    else:
#                        output['!'][tuple(tmppos[t])].append(seq)
                else:
#                    if tmppos[t]:
                        output[t][tuple(tmppos[t])] = [seq]
#                    else:
#                        output['!'][tuple(tmppos[t])] = [seq]
        return output

    def symetric_difference(self, other, method='formated'):
        '''Find elements in either the list or other but not both.
        '''
        seqlist = SequenceList()
        tmplist = self + other
        if method == 'formated':
            for i, seq in enumerate(tmplist):
                if seq not in tmplist[:i] + tmplist[i+1:]:
                    seqlist.append(seq)
        elif method == 'raw':
            rlist = [s.rawSequence for s in tmplist]
            for i, seq in enumerate(tmplist):
                if seq.rawSequence not in rlist[:i] + rlist[i+1:]:
                    seqlist.append(seq)
        else:
            sys.exit("Comparision method sould be one of ('formated', 'raw')")
        return seqlist

    def save(self, outfile):
        '''Writes the sequences in fasta format in the specified output file.
        '''
        for sequence in self:
            outfile.write( sequence.header )
            outfile.write('\n')
            outfile.write( sequence.sequence )
            outfile.write('\n')


class Sequence(object):
    '''Class representing a sequence
    '''
    def __init__(self, header, sequence):
        self.__header = header
        self.__sequence = sequence

    def __str__(self):
        return str([self.header, self.sequence])

    @property
    def header(self):
        return self.__header

    @property
    def sequence(self):
        return self.__sequence

    def __eq__(self, other):
        ss = ''.join(self.sequence.split('\n'))
        os = ''.join(other.sequence.split('\n'))
        if ss == os:
            return True
        return False

    @property
    def rawSequence(self):
        '''Returns the sequence free of any gaps and blanj characters.
        '''
        return ''.join(self.replacePattern('-', '').sequence.split('\n'))

    def replacePattern(self, pattern, replacement):
        tmp = ''
        for c in self.sequence:
            if c == pattern:
                tmp += replacement
            else:
                tmp += c
        return Sequence(self.header, tmp)


def loadSequences( infile ):
    '''
    '''
    sequences = SequenceList()
    inseq = False
    tmp = []
    for line in infile:
        if line.startswith( '>' ):
            if inseq:
                sequences.append( Sequence( tmp[0].strip(), ''.join( tmp[1:] ).strip() ) )
                tmp = []
                inseq = False
            inseq = True
        if inseq:
            tmp.append( line )
    if tmp:
        sequences.append( Sequence( tmp[0].strip(), ''.join( tmp[1:] ).strip() ) )
    return sequences

def saveSequences( sequences, outfile ):
    '''Writes the sequences in fasta format in the specified output file.
    '''
    for sequence in sequences:
        outfile.write( sequence.header )
        outfile.write('\n')
        outfile.write( sequence.sequence )
        outfile.write('\n')

def findPattern( sequences, pattern, mode='count' ):
    '''
    '''
    PATTERN = re.compile( pattern )
    count = 0
    matching = SequenceList()
    for sequence in sequences:
        if re.search( PATTERN, sequence.sequence ):
            matching.append( sequence )
            count += 1
    if mode == 'full':
        return matching
    elif mode == 'header':
        return [sequence.header for sequence in matching]
    else:
        return count

def replacePattern( sequence, pattern, replacement ):
    '''
    '''
    tmp = ''
    for c in sequence.sequence:
        if c == pattern:
            tmp += replacement
        else:
            tmp += c
    return Sequence(sequence.header, tmp)


if __name__ == '__main__':

    with open('/users/rg/agrimaldi/Data/gos/selenoprofiles_profiles/dio_like.det.fasta', 'r') as ff:
        ss = loadSequences(ff)
    uu = ss.findPositions(('U'))

    #print uu
    count = 0
    for i in uu:
        for j in uu[i]:
            if len(uu[i][j]) > 20:
                print i, j, len(uu[i][j])
            count += len(uu[i][j])

    print count, len(ss)
    for i in uu['U'][(1687, 4679)]:
        print i.header
    for i in uu['U'][(1687, 4926)]:
        print i.header
