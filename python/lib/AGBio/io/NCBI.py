#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import re
from Bio.Blast import NCBIStandalone, NCBIXML
from AGBio.io.common import *
import AGBio.io.Fasta as Fasta
from AGBio.Utilities import *


class PsiBlastParser(NCBIStandalone.PSIBlastParser):
    '''Class used to parse psiblast output files.
    '''
    def __init__(self):
        NCBIStandalone.PSIBlastParser.__init__(self)

    def extractData(self, psirounds, evalue=10, spattern=None, hpattern=None, expr='e | s | h' ):
        '''Parses the blast and extract essential data
        '''
        if spattern:
            sp = re.compile(spattern)
        ##if hpattern:
        ##    hp = re.compile(hpattern)
        orexpr = expr.split('|')
        andexpr = expr.split('&')
        sequences = Fasta.SequenceList()
        for psiround in psirounds.rounds:
            for alignment in psiround.alignments:
                for hsp in alignment.hsps:
                    #if spattern in str(hsp.sbjct):
                    #if spattern in hsp.sbjct or hsp.expect <= evalue:
                    if hsp.expect <= evalue:
                    #        or (spattern and (re.search(sp,
                    #                                   hsp.sbjct))):
                            ##or (hpattern and (re.search(hp,
                            ##                             alignment.title.split(']')[0]+']')))
                        sys.stderr.write(str(hsp.expect)+str(hsp.sbjct)+'\n')
                        header = '|'.join(alignment.title.split('|')[:5])
                        if header.endswith(' gi'):
                            header = header[:-3]
                        #header += ' ' + str( hsp.expect )
                        sequences.append(Fasta.Sequence(header,
                                                        hsp.sbjct))
                        #sys.stderr.write(str(repr(sequences)))
        return sequences


class PsiBlastXMLParser(NCBIXML.BlastParser):
    """Parser for psiblast outputs in xml format.
    """
    
    def __init__(self, blastoutput):
        """Constructor.
        """
        NCBIXML.BlastParser.__init__(self)
        self._blastoutput = blastoutput

    @prop
    def blastoutput():
        def fget(self):
            return self._blastoutput
        def fset(self, value):
            self._blastoutput = value
        return locals()

    @prop
    def results():
        def fget(self):
            return self._results
        def fset(self, value):
            self._results = value
        return locals()
    
    def parse(self):
        self.results = NCBIXML.parse(self.blastoutput)

    def extractData(self, evalue=10, fmt=None, ALL=False, outfile=sys.stdout,
                    includepatternsiff=None, includepatterns=None,
                    excludepatterns=None):
        """Extract data of interest from the parsed blast output.
        
        Arguments:
        - `evalue`: Minimum evalue to keep
        - `fmt`: Format the output
        - `ALL`: Get everything
        """
        ffmt = []

        if outfile == sys.stdout or type(outfile).__name__ == 'StringO':
            outf = outfile
        else:
            outf = open(outfile, 'w')
        
        if fmt:
            ffmt = self._parseFormat(fmt)
            ffmtDict = {}

        if ALL:
            ffmt = ('query', 'id', 'accession', 'evalue', 'query_seq',
                    'sbjct_seq', 'full_sbjct_seq', 'gi')

        for result in self.results:
            for al in result.alignments:
                setattr(al, 'gi', al.hit_id.split('|')[1])
                setattr(al, 'header', '>' + al.title.split('>')[0].strip())
                if 'header' in ffmt:
                    ffmtDict['header'] = al.header
                if 'id' in ffmt:
                    ffmtDict['id'] = al.hit_id
                if 'accession' in ffmt:
                    ffmtDict['accession'] = al.accession
                if 'gi' in ffmt:
                    ffmtDict['gi'] = al.gi
                for hsp in al.hsps:
                    okflag = False
                    if (hsp.expect <= evalue \
                        and (not (self._contains(al, excludepatterns) \
                             or self._contained(al, excludepatterns))) \
                        and (True if not includepatternsiff \
                             else (self._contains(al, includepatternsiff) \
                                   or self._contained(al, includepatternsiff))) \
                        or (includepatterns and self._contains(al, includepatterns)) \
                        ):
                        okflag = True
                        if 'evalue' in ffmt:
                            ffmtDict['evalue'] = hsp.expect
                        if 'sbjct_seq' in ffmt:
                            ffmtDict['sbjct_seq'] = hsp.sbjct
                if okflag:
                    line = ' '.join([str(ffmtDict[opt]) for opt in ffmt])
                    writeline(outf, line)
        outf.close()

    def _parseFormat(self, fmt):
        allowed = ('query', 'id', 'gi', 'accession', 'header', 'evalue',
                   'query_seq', 'sbjct_seq', 'full_sbjct_seq')
        ffmt = fmt.split(',')
        for i in ffmt:
            assert i in allowed
        return ffmt
    
    def _contains(self, al, patterns):
        if patterns != None:
            allowed = ('title', 'sbjct_seq', 'header')
            excluded = False
            for keyword, patterns in patterns.items():
                if keyword not in allowed: break
                if keyword not in al.__dict__:
                    raise KeyError, 'invalid keyword : '+keyword
                for pattern in patterns:
                    if pattern in getattr(al, keyword):
                        return True
        return False

    def _contained(self, al, patterns):
        if patterns != None:
            allowed = ('gi')
            excluded = False
            for keyword, patterns in patterns.items():
                if keyword not in allowed: break
                if keyword not in al.__dict__:
                    raise KeyError, 'invalid keyword : '+keyword
                if getattr(al, keyword) in patterns:
                    return True
        return False
