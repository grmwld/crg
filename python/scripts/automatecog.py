#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

COMMAND = '/users/rg/agrimaldi/Code/crg/python/scripts/searchcog.py'

cogtuple = ( ('COG0526', 'trx_like_1'),
             ('COG2128', 'ahpd_like_1'),
             ('COG0243', 'fdha'),
             ('COG0737', 'usha_like'),
             ('COG0386', 'gpx'),
             ('COG0625', 'gst'),
             ('COG0607', 'sulft_2'),
             ('COG1651', 'dsbg_like'),
             ('COG2897', 'sulft_3'),
             ('COG4802', 'frx'),
             ('COG5640', 'trypsin_like'),
             ('COG2331', 'fmdb') )

for t in cogtuple:

    print ( COMMAND + ' -i ' + t[0] + ' -o ' + t[1] )
    os.system( COMMAND + ' -i ' + t[0] + ' -o ' + t[1] )
