#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import re
from Bio.Blast import NCBIStandalone, NCBIXML
from AGBio.io.common import *
import AGBio.io.Fasta as Fasta
import AGBio.UtilityWrappers as UtilityWrappers


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

    @property
    def blastoutput(self):
        return self._blastoutput
    @blastoutput.setter   
    def blastoutput(self, value):
        self._blastoutput = value

    @property
    def results(self):
        return self._results
    @results.setter   
    def results(self, value):
        self._results = value
        
    def parse(self):
        self.results = NCBIXML.parse(self.blastoutput)

    def extractData(self, evalue=10, fmt=None, ALL=False, outfile=sys.stdout):
        """Extract data of interest from the parsed blast output.
        
        Arguments:
        - `evalue`: Minimum evalue to keep
        - `fmt`: Format the output
        - `ALL`: Get everything
        """
        ffmt = []
        
        if outfile != sys.stdout:
            outf = open(outfile, 'w')
        else:
            outf = outfile
        
        if fmt:
            ffmt = self._parseFormat(fmt)
            ffmtDict = {}

        if ALL:
            ffmt = ('query', 'id', 'accession', 'evalue', 'query_seq',
                   'sbjct_seq', 'full_sbjct_seq')

        for result in self.results:
            for al in result.alignments:
                if 'header' in ffmt:
                    ffmtDict['header'] = '>' + al.title.split('>')[0].strip()
                if 'id' in ffmt:
                    ffmtDict['id'] = '>' + al.hit_id
                if 'accession' in ffmt:
                    ffmtDict['accession'] = al.accession
                for hsp in al.hsps:
                    if hsp.expect <= evalue:
                        if 'evalue' in ffmt:
                            ffmtDict['evalue'] = hsp.expect
                line = ' '.join([str(ffmtDict[opt]) for opt in ffmt])
                writeline(outf, line)
        outf.close()

    def _getGI(self, alignment):
        return alignment.id
        

    def _parseFormat(self, fmt):
        allowed = ('query', 'id', 'accession', 'evalue', 'query_seq',
                   'sbjct_seq', 'full_sbjct_seq')
        
        ffmt = fmt.split(',')
        for i in ffmt:
            assert i in allowed

        return ffmt
