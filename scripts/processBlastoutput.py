#!/usr/bin/env python2.6

import os
import sys
import time
import random
import optparse
import tempfile
import AGBio.IO.Fasta as Fasta
from AGBio.UtilityWrappers import *

TEMP = '/home/agrimaldi/temp/'

def findNotInMainIndex(infile, mainindexfile):
    '''Find gi not in the main index
    '''
    mifile = open(mainindexfile, 'r')
    tifile = open(infile, 'r')
    #tff = '|'.join([gi.split('|')[1] for gi in tifile.readlines())
    found = []
    for i in tifile:
        try:
            mifile.seek(0)
            nfound = True
            for j in mifile:
                if i.split('|')[1] == j.split('|')[1]:
                    nfound = False
                    break
            if nfound:
                found.append(i)
        except IndexError:
            sys.stderr.write(i + ' ' + j)
    mifile.close()
    tifile.close()
    return found

def findNotInMainIndex2(infile, mainindexfile):
    '''Find gi not in the main index.
    Use sets for improved efficiency.
    '''
    mset = set()
    tset = set()
    with open(mainindexfile, 'r') as mfile:
        for line in mfile:
            mset.add(line.split('|')[1])
    with open(infile, 'r') as tfile:
        for line in tfile:
            tset.add(line.split('|')[1])
    return list(tset.difference(mset))

def countHeaders(infile):
    '''Counts the number of fasta headers in a file.
    '''
    count = 0
    with open(infile, 'r') as handle:
        for line in handle:
            if line.startswith('>'):
                count += 1
    return count

def genTempfilename(dirr='./', prefix='tmp'):
    '''Generate temporary filename
    '''
    if not dirr.endswith('/'):
        dirr += '/'
    #random.seed()
    hashn = str(random.getrandbits(128))
    return ''.join((dirr, prefix, hashn))

def uniq(infile):
    '''Use remove duplicate headers, keeping the one with the lowest evalue.
    '''
    headict={}
    with open(infile) as ifile:
        for line in ifile:
            sline = line.split()
            if sline[0] in headict:
                if float(sline[1]) < float(headict[sline[0]]):
                    headict[sline[0]] = sline[1]
            else:
                headict[sline[0]] = sline[1]
    with open(infile, 'w') as ofile:
        for header, evalue in headict.items():
            ofile.write(header + ' ' + evalue + '\n')


def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='fasta file in which selenoproteins should be looked for.',
                       metavar='FILE' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='base output filename',
                       metavar='FILE' )

    parser.add_option( '-m', '--main',
                        dest='mainfilename',
                       help='base name for the main index and main fasta',
                       metavar='FILE' )
    
    parser.add_option( '-e', '--evalue',
                       dest='evalue',
                       help='e-value threshold.',
                       metavar='FLOAT' )

    parser.add_option( '-p', '--pattern',
                       dest='pattern',
                       help='pattern to look for in sequences.',
                       metavar='REGEX')

    parser.add_option( '-v', '--verbose',
                       dest='verbosity',
                       help='verbosity level : 0=none ; 1=standard ; 2=detailed ; 3=full',
                       metavar='INTEGER' )

    parser.set_defaults( verbosity = '1',
                         evalue = '10',
                         pattern = None)

    (options, args) = parser.parse_args()

    verbosity = int(options.verbosity)
    evalue = options.evalue
    pattern = options.pattern

    mainindex = ''.join(( options.mainfilename.split('.')[0], '.index' ))
    mainrepo = ''.join(( options.mainfilename.split('.')[0], '.fasta' ))
    blastindexfile = ''.join(( options.outputfilename.split('.')[0], '.index.0' ))
    blastfastafile = ''.join(( options.outputfilename.split('.')[0], '.fasta.0' ))
    tmpdlrepofilename = genTempfilename(dirr=TEMP, prefix='getseqfromgi')

    os.system(' '.join(( 'touch', blastindexfile )))
    os.system(' '.join(( 'touch', blastfastafile )))
    os.system(' '.join(( 'touch', tmpdlrepofilename )))


    blastparser = BlasterParserWrapper( options.inputfilename,
                                        blastindexfile,
                                        evalue=evalue,
                                        pattern=pattern )

    findInRepo = AddFullHeadersWrapper( mainrepo,
                                        blastfastafile,
                                        blastindexfile)

    downloader = DownloaderWrapper( None,
                                    tmpdlrepofilename )

    todownload = [] ## findNotInMainIndex( blastindexfile, mainindex )


    ## Parse the blast output file.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Parsing blast output.\n' )
    if verbosity >= 2:
            sys.stderr.write( blastparser.cline + '\n' )
    blastparser.run()

    ## Remove the space between '>' and 'gi' in the index file.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Sanitizing the index file.\n' )
    os.system(''.join(("sed -i 's/> gi|/>gi|/g' ", blastindexfile )))

    ## Only keep one copy of a header, the one with the best evalue.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Keeping only best evalues.\n' )
    uniq(blastindexfile)

    ## While the number of header differs between the index and the fasta file, keep trying to download.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Building fasta.0 file.\n' )

    try:
        while countHeaders( blastindexfile ) != countHeaders( blastfastafile ):
            tmpdlrepofilename = genTempfilename(dirr=TEMP, prefix='getseqfromgi')
            os.system(' '.join(( 'touch', tmpdlrepofilename )))
            if verbosity >= 2:
                sys.stderr.write( '\n' )
                sys.stderr.write( '    index.0 : ' + str(countHeaders( blastindexfile )) + '\n' + \
                                  '    fasta.0 : ' + str(countHeaders( blastfastafile )) + '\n')
            todownload = findNotInMainIndex( blastindexfile, mainindex )
            if verbosity >= 2:
                sys.stderr.write( '\n' )
                sys.stderr.write('>>> Downloading ' + str(len(todownload)) + ' sequences ...\n')
            for header in todownload:
                time.sleep(0.001)
                downloader.infile = header.split('|')[1]
                downloader.outfile = tmpdlrepofilename
                if verbosity >= 3:
                    sys.stderr.write( '\n' )
                    sys.stderr.write(downloader.cline)
                downloader.run()
            if verbosity >= 2:
                sys.stderr.write( '\n' )
                sys.stderr.write('>>> Adding newly downloaded sequences to the repository.\n')
            os.system(' '.join(('cat', tmpdlrepofilename, '>>', mainrepo)))
            if verbosity >= 2:
                sys.stderr.write( '\n' )
                sys.stderr.write( '>>> Sanitizing the main repository file.\n' )
            os.system(' '.join(("sed -i '/^$/d'", mainrepo )))
            os.system(' '.join(("sed -i '/>Error:/d'", mainrepo )))
            if verbosity >= 2:
                sys.stderr.write( '\n' )
                sys.stderr.write('>>> Updating headers of the main repository.\n')
            os.system(' '.join(('grep ">"', mainrepo, '>', mainindex)))
            if verbosity >= 2:
                sys.stderr.write( '\n' )
                sys.stderr.write('>>> Fetching sequences from the main repository.\n')
            if verbosity >= 3:
                sys.stderr.write( findInRepo.cline + '\n' )
            findInRepo.run()
    except KeyboardInterrupt:
        sys.exit('Manual exit. The downloaded sequences won\'t be added to the main repository.')
    finally:
        os.system(' '.join(('rm', tmpdlrepofilename)))
        pass

if __name__ == '__main__':
    main()
