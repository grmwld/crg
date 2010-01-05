#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

CP = '/bin/cp '
SRC_DIR = '/users/rg/agrimaldi/Data/cdd/'
DST_DIR = '/users/rg/agrimaldi/Data/gos_raw/'

cdtuple = ( ('cd02195', 'seld'),
            ('cd02950', 'trx_like_3'),
            ('cd03032', 'arsc_2'),
            ('cd02953', 'trx_like_2'),
            ('cd01449', 'sulft_1') )

for t in cdtuple:

    print( CP + SRC_DIR + t[0] + '.FASTA ' + DST_DIR + t[1] + '.faln' )
    os.system( CP + SRC_DIR + t[0] + '.FASTA ' + DST_DIR + t[1] + '.faln' )
