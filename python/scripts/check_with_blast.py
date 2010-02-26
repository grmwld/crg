#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import optparse
import shutil
from cStringIO import StringIO
sys.path.append('/users/rg/agrimaldi/Code/python/python/lib/')
from AGBio.ncbi.BlastWrappers import *
from AGBio.io.NCBI import PsiBlastXMLParser
from AGBio.Utilities import *


class HeadEvalueDict(BiDict):
    def __init__(self, dictionary={}):
        BiDict.__init__(self, dictionary)

    @prop
    def reverse():
        def fget(self):
            output = {}
            for key, value in self.items():
                try:
                    for v in value:
                        if v not in output:
                            output[v] = key
                        else:
                            output[v] = min(float(key), float(output[v]))
                except TypeError:
                    print v, key
                    sys.exit('beuh')
            return output
        return locals()


def main():

    parser = optparse.OptionParser()

    parser.add_option('-s', '--entry',
                      dest='gi_entry',
                      help='GI to check against the database',
                      metavar='GI')

    parser.add_option('-b', '--blast_flavour',
                      dest='blast_flavour',
                      help='what kind of blast should be performed ?',
                      metavar='BLAST')

    parser.add_option('-d', '--database_check',
                      dest='dbc',
                      help='location of the database that should be used for checking.',
                      metavar='DB')
    
    parser.add_option('-D', '--database_fetch',
                      dest='dbf',
                      help='location of the database that should be used for fetching the sequence from the gi provided.',
                      metavar='DB')

    parser.add_option('-a', '--ncore',
                      dest='ncore',
                      type='int',
                      help='number of cores to use for the blast.',
                      metavar='INT')

    parser.add_option('-n', '--num_top_hits',
                      dest='num_top_hits',
                      type='int',
                      help='number of top hits to consider.',
                      metavar='INT')

    parser.add_option('-f', '--config_file',
                      dest='config_file',
                      help='location of the file containing filters.',
                      metavar='FILE')

    parser.add_option('-T', '--temp',
                      dest='temp',
                      help='temporary folder.',
                      metavar='DIR')

    parser.set_defaults(temp = '/tmp/',
                        ncore = 1,
                        num_top_hits = 1)

    (options, args) = parser.parse_args()

    outputentryfa = os.path.join(options.temp, options.gi_entry + '.fa')
    outputblast = os.path.join(options.temp, options.gi_entry + '.xml')
    outputpf = os.path.join(options.temp, options.gi_entry + '.index')


    fetcher = FastaCmdWrapper([options.gi_entry], db=options.dbf,
                              outfile=outputentryfa)
    
    blaster = BlastAllWrapper(outputentryfa, outputblast,
                              flavour=options.blast_flavour,
                              db=options.dbc, gis=True, ncore=options.ncore)

    xmlparser = PsiBlastXMLParser(outputblast)

    print fetcher.cline
    fetcher.run()
    print blaster.cline
    blaster.run()

    with open(outputblast, 'r') as iff:
        xmlparser = PsiBlastXMLParser(iff)
        xmlparser.parse()
        xmlparser.extractData(fmt='evalue,header', outfile=outputpf)

    results = HeadEvalueDict()
    nseqs = 0

    while nseqs < options.num_top_hits:
        with open(outputpf, 'r') as iff:
            for line in iff:
                print line
                evalue = line.split()[0]
                header = ' '.join(line.split()[1:])
                if evalue not in results:
                    results[evalue] = [header]
                else:
                    results[evalue].append(header)
        nseqs = len(results.reverse)

    for k, v in results.items():     
        if float(k) < 1e-80:
            print k, v
        
    for k, v in results.reverse.items():
        print k, v
            

if __name__ == '__main__':
    main()
