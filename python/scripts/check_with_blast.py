#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import optparse
import shutil
import logging
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
                for v in value:
                    if v not in output:
                        output[v] = key
                    else:
                        output[v] = min(float(key), float(output[v]))
            return output
        return locals()


def parsefilters(filename):
    filters = []
    with open(filename, 'r') as iff:
        for line in iff:
            filters.extend([f.strip() for f in line.split(',')])
    return filters

def main():

    parser = optparse.OptionParser()

    fetchgroup = optparse.OptionGroup(parser, 'Options to work with a GI')
    blastgroup = optparse.OptionGroup(parser, 'Blast related options')

    fetchgroup.add_option('-s', '--entry',
                          dest='gi_entry',
                          help='GI to check against the database',
                          metavar='GI')

    fetchgroup.add_option('-D', '--database_fetch',
                          dest='dbf',
                          help='location of the database that should be used for fetching the sequence from the gi provided.',
                          metavar='DB')

    blastgroup.add_option('-b', '--blast_flavour',
                      dest='blast_flavour',
                      help='what kind of blast should be performed ?',
                      metavar='BLAST')

    blastgroup.add_option('-d', '--database_check',
                          dest='dbc',
                          help='location of the database that should be used for checking.',
                          metavar='DB')

    blastgroup.add_option('-a', '--ncore',
                          dest='ncore',
                          type='int',
                          help='number of cores to use for the blast.',
                          metavar='INT')

    parser.add_option('-q', '--query',
                      dest='fasta_query',
                      help='query in fasta format',
                      metavar='FILE')

    parser.add_option('-o', '--output',
                      dest='outputfile',
                      help='name of the output file. default is stdout',
                      metavar='FILE')

    parser.add_option('-n', '--num_top_hits',
                      dest='num_top_hits',
                      type='int',
                      help='number of top hits to consider.',
                      metavar='INT')

    parser.add_option('-f', '--filters_file',
                      dest='filters_file',
                      help='location of the file containing filters.',
                      metavar='FILE')

    parser.add_option('-v', '--verbosity',
                      dest='verbosity', action='count',
                      help='set verbosity level')

    parser.add_option('-T', '--temp',
                      dest='temp',
                      help='temporary folder.',
                      metavar='DIR')

    parser.add_option_group(fetchgroup)
    parser.add_option_group(blastgroup)

    parser.set_defaults(temp = '/tmp/',
                        ncore = 1,
                        num_top_hits = 1,
                        outputfile = sys.stdout)

    (options, args) = parser.parse_args()

    if len(sys.argv) == 1:
        parser.error('No options specified. check_with_blast.py --help for details.')

    log_level = logging.WARNING
    if options.verbosity == 1:
        log_level = logging.INFO
    elif options.verbosity >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level,
                        format='%(levelname)-6s:%(filename)s  %(message)s')

    if options.filters_file:
        includefilters = parsefilters(options.filters_file)
    else:
        includefilters = None
    logging.info('Filters : '+str(includefilters))

    if options.gi_entry:
        outputentryfa = os.path.join(options.temp, options.gi_entry + '.fa')
        outputblast = os.path.join(options.temp, options.gi_entry + '.xml')
        outputpf = os.path.join(options.temp, options.gi_entry + '.index')
        fetcher = FastaCmdWrapper([options.gi_entry], db=options.dbf,
                                  outfile=outputentryfa)
        blastqueryfile = outputentryfa
    elif options.fasta_query:
        outputblast = os.path.join(options.temp,
                                   os.path.basename(options.fasta_query) \
                                   + '.xml')
        outputpf = os.path.join(options.temp,
                                os.path.basename(options.fasta_query) \
                                + '.index')
        blastqueryfile = options.fasta_query
    
    blaster = BlastAllWrapper(blastqueryfile, outputblast,
                              flavour=options.blast_flavour,
                              db=options.dbc, gis=True, ncore=options.ncore)

    xmlparser = PsiBlastXMLParser(outputblast)

    if options.gi_entry:
        logging.info(' '+fetcher.cline)
        fetcher.run()
    logging.info('Running blast : '+blaster.cline)
    blaster.run()

    with open(outputblast, 'r') as iff:
        logging.info('Parsing the xml output -> '+outputpf)
        xmlparser = PsiBlastXMLParser(iff)
        xmlparser.parse()
        xmlparser.extractData(fmt='evalue,header', outfile=outputpf)

    results = HeadEvalueDict()

    with open(outputpf, 'r') as iff:
        for line in iff:
            evalue = line.split()[0]
            header = ' '.join(line.split()[1:])
            if evalue not in results:
                results[evalue] = [header]
            else:
                results[evalue].append(header)

    topindexes = results.keys()
    topindexes.sort(lambda e1, e2: cmp(float(e1), float(e2)))
    
    topseqs = [(e, results[e]) for e in topindexes[:options.num_top_hits]]

    finaloutput = []

    for eseq in topseqs:
        for header in eseq[1]:
            if options.filters_file:
                for ikw in includefilters:
                    if ikw in header:
                        finaloutput.append((eseq[0], header))
            else:
                finaloutput.append((eseq[0], header))
    logging.debug(str(finaloutput))
    try:
        if options.outputfile != sys.stdout:
            off = open(options.outputfile, 'w')
        else:
            off = sys.stdout
        off.write('# '+str(includefilters)+'\n')
        for oo in finaloutput:
            off.write(oo[0]+' : '+oo[1]+'\n')
    except Exception, (e):
        print e
    finally:
        if off == sys.stdout:
            off.close()


if __name__ == '__main__':
    main()
