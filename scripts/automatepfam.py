#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

COMMAND = '/users/rg/agrimaldi/Code/crg/python/scripts/searchpfam.py'

pfamtuple = ( ('PF10262', 'selw_like' ),
              ('PF00578', 'prx_like_1' ),
              ('PF07355', 'prdb' ),
              ('PF04592', 'prx_like_2' ),
              ('PF03960', 'arsc_1' ),
              ('PF02627', 'ahpd_like_2' ),
              ('PF08534', 'prx_like_3' ),
              ('PF01323', 'dsba' ),
              ('PF00462', 'grx_like' ),
              ('PF00837', 'dio_like' ),
              ('PF05237', 'moeb' ),
              ('PF02635', 'dsre_like' ),
              ('PF01625', 'msra' ),
              ('PF07355', 'grdb' ),
              ('PF02566', 'osmc_like' ),
              ('PF02508', 'nadh_uo_e' ),
              ('PF02662', 'frhd' ),
              ('PF07610', 'os_hp2' ),
              ('PF04723', 'grda' ),
              ('PF00122', 'atpase_e1_e2' ),
              ('PF01035', 'mett' ),
              ('PF00070', 'hdra' ),
              ('PF04945', 'yhs') )

for t in pfamtuple:

    print ( COMMAND + ' -i ' + t[0] + ' -o ' + t[1] )
    os.system( COMMAND + ' -i ' + t[0] + ' -o ' + t[1] )
