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
        ## for backward compatibility:
        self.save = self.prints

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

##     def save(self, outfile):
##         '''Writes the sequences in fasta format in the specified output file.
##         '''
##         for sequence in self:
##             outfile.write( sequence.header )
##             outfile.write('\n')
##             outfile.write( sequence.sequence )
##             outfile.write('\n')

    def prints(self, outfile=sys.stdout, length=80):
        for seq in self:
            seq.prints(outfile, length)


class Alignment(SequenceList):
    '''Alignment class.
    '''
    def __init__(self, sequences=''):
        '''
        '''
        SequenceList.__init__(self, sequences)
        for i, seq in enumerate(self):
            for sseq in self[i:]:
                if len(seq) != len(sseq):
                    raise Exception()
        self.length = len(self[0].sequence)
        self.profiles=[]

    def getProfiles(self):
        for pos in range():
            pass

    def findPositions(self, aaa=None, redundant=True):
        '''Finds all positions of U in the alignment.
        '''
        aa = list(aaa)
        output = {}
        if not redundant:
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
                            output[t][tuple(tmppos[t])].append(seq)
                    else:
                            output[t][tuple(tmppos[t])] = [seq]
        else:
            for i in aa:
                output[i] = {}
            for seq in self:
                for p in aa:
                    foundany = False
                    for i, pos in enumerate(seq.sequence):
                        if pos == p:
                            foundany = True
                            try:
                                if tuple([i]) in output[p]:
                                    output[p][tuple([i])].append(seq)
                                else:
                                    output[p][tuple([i])] = [seq]
                            except TypeError:
                                print pos, p, i, output[p]
                    if not foundany:
                        if () in output[p].keys():
                            output[p][()].append(seq)
                        else:
                            output[p][()] = [seq]
        return output


class Sequence(object):
    '''Class representing a sequence
    '''
    def __init__(self, header, sequence):
        self.__header = header
        self.__sequence = ''.join(sequence.split())

    def __str__(self):
        return str([self.header, self.sequence])

    def __len__(self):
        return len(self.sequence)

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

    def prints(self, outfile=sys.stdout, length=80):
        if length > 0:
            llength = length
        else:
            llength = len(self)
        outfile.write(self.header)
        outfile.write('\n')
        for buf in range(0, len(self), llength):
            outfile.write(self.sequence[buf:buf+llength])
            outfile.write('\n')

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
            tmp.append( line.strip() )
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
        aa = loadSequences(ff)
    ss = Alignment(aa)
    uu = ss.findPositions(('U','C','-'), True)

    count = 0
    for i in uu:
        for j in uu[i]:
            if len(uu[i][j]) > 200:
                print i, j, len(uu[i][j])
            count += len(uu[i][j])
    
    print count, len(ss)

    print len(uu['U'][(1687,)]), len(uu['C'][(1687,)]), len(uu['-'][(1687,)])
