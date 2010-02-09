#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

CP = '/bin/cp '
SRC_DIR = '/users/rg/agrimaldi/Data/prk/'
DST_DIR = '/users/rg/agrimaldi/Data/gos_raw/'

prktuple = ( ('PRK11873', 'ars_s_amt'),
             ('PRK09564', 'nadh_ox') )

for t in prktuple:

    print( CP + SRC_DIR + t[0] + '.aln ' + DST_DIR + t[1] + '.faln' )
    os.system( CP + SRC_DIR + t[0] + '.aln ' + DST_DIR + t[1] + '.faln' )
