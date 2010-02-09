#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import optparse

BASE_DIR = '/users/rg/agrimaldi/Data/cog/cog/COG'
NEW_ORG = re.compile( r'\w{3}:' )


def parseWhogFile( cogID, whogfile ):

    accessIdList = []
        
    found = False

    header = True

    for line in whogfile:
        
        if line.startswith( '[' ):

            if line.split()[1] == cogID:

                found = True

            elif found:
                
                break

        elif found:

            for id in parseWhogLine( line ):

                accessIdList.append( id )

    return accessIdList


def parseWhogLine( line ):

    if re.match( NEW_ORG, line.strip() ):

        return line.split()[ 1: ]

    else:

        return line.split()


def parseMyvaFile( ids, myvafile ):

    coglist = []

    tmp = []

    inseq = False

    done = False

    for line in myvafile:
        
        if line.startswith( '>' ):

            if inseq:

                coglist.append( ''.join( tmp ).strip() )

                tmp = []
                
                inseq = False

            for id in ids:

                if line[1:].startswith( id ):

                    inseq = True

        if inseq:
            
            tmp.append( line )
            
    coglist.append( ''.join( tmp ).strip() )

    return coglist


def main():

    parser = optparse.OptionParser()

    parser.add_option( '-w', '--whogfile',
                       dest='whogfilename',
                       help='whog filename in which the family ids should be looked for',
                       metavar='WHOG_FILE' )

    parser.add_option( '-m', '--myvafile',
                       dest='myvafilename',
                       help='myva filename in which the family sequences should be looked for',
                       metavar='MYVA_FILE' )

    parser.add_option( '-o', '--output',
                       dest='outputfilename',
                       help='name of the file in which the family should be stored',
                       metavar='FILE' )

    parser.add_option( '-i', '--id',
                       dest='cogid',
                       help='COG Accession number of the family that should be looked for',
                       metavar='COG_ID' )

    parser.set_defaults( whogfilename=os.path.join( BASE_DIR, 'whog' ),
                         myvafilename=os.path.join( BASE_DIR, 'myva' ) )

    (options, args) = parser.parse_args()

    whogfile = open( options.whogfilename, 'r' )
    myvafile = open( options.myvafilename, 'r' )
    outputfile = open( options.outputfilename, 'w' )

    results = parseMyvaFile( parseWhogFile( options.cogid, whogfile ), myvafile )

    whogfile.close()
    myvafile.close()

    for sequence in results:

        outputfile.write( sequence )
        outputfile.write( '\n' )

    outputfile.close()

if __name__ == '__main__':

    main()
