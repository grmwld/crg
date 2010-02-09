#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import optparse

BINDIR = '/users/rg/agrimaldi/usr/bin/'
BLAST = 'blastpgp -F F -b 5000 -I'
TEMP = '/home/agrimaldi/temp'
PROFILE_PREP = BINDIR + 'prepare_alignment_selenoprofiles.py -temp ' + TEMP


def prepareprofile( ifile, ofile ):

    command = ' '.join( ( PROFILE_PREP, '-i', ifile, '-output', ofile, '-all' ) )
    
    os.system( command )


def runblast( ifile, ofilename, database, pssm, ncore ):

    command = ' '.join( ( BLAST, '-i', ifile, '-d', database, '-R', pssm, '-a', ncore, '-o', ofilename ) )

    os.system( command )


def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='file containing the alignments that will be used to build the PSSM using prepare_alignment_selenoprofiles.py.',
                       metavar='FILE' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='base name used for outputs',
                       metavar='NAME' )

    parser.add_option( '-d', '--database',
                       dest='database',
                       help='database used by blast.',
                       metavar='FILE' )

    parser.add_option( '-c', '--n_core',
                       dest='ncore',
                       help='number of cores to use during blast',
                       metavar='INTEGER' )

    parser.set_defaults( ncore = '1' )

    (options, args) = parser.parse_args()

    ifile = options.inputfilename

    if options.outputfilename:

        ofile = options.outputfilename

    else :

        ofile = ifile.split('.')[0]


    query = ''.join( ( ofile, '.query.fa' ) )
    pssm = ''.join( ( ofile, '.profile.chk' ) )
    blastOutput = ''.join( ('./blast/', os.path.split( ofile )[1], '.blast' ) )

    print
    print '>>> Preparing pssm of ' + ifile
    prepareprofile( ifile, ofile )

    print
    print '>>> Blasting ' + ifile
    runblast( query, blastOutput , options.database, pssm, options.ncore )

if __name__ == '__main__':

    main()
