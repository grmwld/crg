#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import optparse
import shutil
import traceback
from cStringIO import StringIO
sys.path.append('/users/rg/agrimaldi/Code/python/python/lib/')
from AGBio.io.common import *
from AGBio.Utilities import *
from AGBio.selenoprofiles_tools.results_analyser import *
from AGBio.selenoprofiles_tools.files_analysers.p2g import *


def main():

    parser = optparse.OptionParser()

    parser.add_option('-f', '--folder',
                      dest='sp_folder',
                      help='selenoprofiles output folder',
                      metavar='DIR')

    parser.add_option('-c', '--check_with_blast',
                      action='store_true', default=False,
                      dest='check_with_blast',
                      help='chack candidates with blast')

    parser.add_option('-o', '--output',
                      dest='outputfile',
                      help='name of the output file. default is stdout',
                      metavar='FILE')

    parser.add_option('-T', '--temp',
                      dest='temp',
                      help='temporary folder',
                      metavar='DIR')

    parser.set_defaults(outputfile = sys.stdout,
                        temp = '/tmp/')

    (options, args) = parser.parse_args()

    organisms_folders = [os.path.join(options.sp_folder, f) \
                         for f in os.listdir(options.sp_folder)]

    try:
        tmp_run_folder = FolderWrapper(genTempfilename(options.temp,
                                                       'check_'))
        tmp_run_folder.create()
        for f in organisms_folders:
            tmp_org_folder = tmp_run_folder.sub_folder(os.path.basename(f))
            tmp_org_folder.create()
            sp_parser = GenomeFolderParser(f)
            sp_parser.parse(sec=True, cys=True, thr=True, arg=True, uga=True)
            sp_parser.parseResultFiles(p2g=True, bsecisearch=False)
            for p2g in sp_parser.p2g:
                tmp_query_file = tmp_org_folder.sub_file()
                tmp_query_file.create()
                p2g.result.target.fasta().prints()
        print '=========='
    except Exception:
        print traceback.print_exc()
    finally:
        tmp_run_folder.delete(recursive=True)
    
if __name__ == '__main__':
    main()
