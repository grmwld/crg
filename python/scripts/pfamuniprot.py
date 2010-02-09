#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import optparse

INFILE = '/users/rg/agrimaldi/Data/gos_raw/pfam/uniprotIDs'
OUTFILE = '/users/rg/agrimaldi/Data/gos_raw/pfam/uniprotIDs.unique'

def removedoblelines( infile, outfile ):

    setID = set()

    for line in infile:

        setID.add( line.strip() )

    for i in setID:

        outfile.write( i )
        outfile.write('\n')


if __name__ == '__main__':

    infile = open(INFILE, 'r')
    outfile = open(OUTFILE, 'w')

    removedoblelines(infile, outfile)

    infile.close()
    outfile.close()
