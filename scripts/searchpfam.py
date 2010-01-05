#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import optparse

BASE_DIR = '/users/rg/agrimaldi/Data/pfam/24_0'


def parsePfamFile( id, pfamfile, optimistic=False ):

    pfamlist = []

    tmp = []

    inseq = False

    if optimistic:
        
        infam = False

    for line in pfamfile:
        
        if line.startswith( '>' ):

            if inseq:

                pfamlist.append( ''.join( tmp ).strip() )

                tmp = []
                
                inseq = False

            if line.split()[1].startswith( id ):

                inseq = True

                if optimistic:

                    infam = True

            elif optimistic and infam:
                
                break

        if inseq:
            
            tmp.append( line )
            
    pfamlist.append( ''.join( tmp ).strip() )

    return pfamlist



def resurect( id, deadfile ):

    corpse = False

    for line in deadfile:

        splitline = line.split()

        if len( splitline ) == 3 and splitline[ 1 ] == 'AC' and splitline[ 2 ] == id:

            corpse = True
            
        elif corpse and splitline[ 1 ] == 'FW':

            return splitline[ 2 ]

    return 'buried'



def main():

    parser = optparse.OptionParser()

    parser.add_option( '-f', '--file',
                       dest='filename',
                       help='pfam filename in which the family should be looked for',
                       metavar='FILE' )

    parser.add_option( '-o', '--output',
                       dest='outputfilename',
                       help='output file name',
                       metavar='FILE' )

    parser.add_option( '-i', '--id',
                       dest='accessid',
                       help='Accession number of the family that should be looked for',
                       metavar='ID' )

    parser.add_option( '-p', '--optimistic',
                       action='store_true', dest='optimistic', default=False,
                       help='Stops searching as soon as no new id matches' )

    parser.add_option( '-r', '--resurect',
                       action='store_true', dest='resurect', default=False,
                       help='If a family is dead, try to find the id in which it has been merged' )
    

    parser.set_defaults( filename=os.path.join( BASE_DIR, 'Pfam-A.fasta' ) )

    (options, args) = parser.parse_args()


    if options.resurect:

        file = open( os.path.join( BASE_DIR, 'Pfam-A.dead' ), 'r' )
        
        print resurect( options.accessid, file )

        file.close()

    else:

        ifile = open( options.filename, 'r' )
        ofile = open( options.outputfilename, 'w' )

        results = parsePfamFile( options.accessid, ifile, optimistic=options.optimistic )

        for result in results:
            
            ofile.write( result )
            ofile.write( '\n' )

        ifile.close()
        ofile.close()


if __name__ == '__main__':

    main()
