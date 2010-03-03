#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import optparse
import shutil
import traceback
import logging
import subprocess
from cStringIO import StringIO
sys.path.append('/users/rg/agrimaldi/Code/python/python/lib/')
from AGBio.io.common import *
from AGBio.Utilities import *
from AGBio.selenoprofiles_tools.results_analyser import *
from AGBio.selenoprofiles_tools.files_analysers.p2g import *

def check_with_blast(organisms_folders=None, temp=None):
    try:
        tmp_run_folder = FolderWrapper(genTempfilename(temp, 'check_'))
        logging.debug('Creating run folder : '+tmp_run_folder.abspath)
        tmp_run_folder.create()
        for f in organisms_folders:
            tmp_org_folder = tmp_run_folder.sub_folder(os.path.basename(f))
            logging.debug('Creating organism folder : '+tmp_org_folder.abspath)
            tmp_org_folder.create()
            sp_parser = GenomeFolderParser(f)
            sp_parser.parse(sec=True, cys=True, thr=True, arg=True, uga=True)
            sp_parser.parseResultFiles(p2g=True, bsecisearch=False)
            for p2g in sp_parser.p2g:
                tmp_query_file = tmp_org_folder.sub_file()
                logging.debug('Creating query file : '+tmp_query_file.abspath)
                tmp_query_file.create()
                logging.info(p2g.filename)
#                logging.verbose1(p2g.result.target.fasta())
                with open(tmp_query_file.abspath, 'w') as off:
                    p2g.result.target.fasta().prints(off)
                cmd = ' '.join(['check_with_blast.py', '-v',
                                '-q', tmp_query_file.abspath,
                                '-b', 'blastp',
                                '-d', '/seq/databases/nr_uncompressed/nr',
                                '-a', '4',
                                '-n', '10',
                                '-T', temp])
                subprocess.call(cmd, shell=True)
    except Exception:
        print traceback.print_exc()
    finally:
        logging.info('Cleaning temporary files'.rjust(10))
        tmp_run_folder.delete(recursive=True)

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

    parser.add_option('-v', '--verbosity',
                      dest='verbosity', action='count',
                      help='set verbosity level')

    parser.set_defaults(outputfile = sys.stdout,
                        temp = '/tmp/')

    (options, args) = parser.parse_args()

##     VERBOSE_1 = 21
##     VERBOSE_2 = 22

##     logging.addLevelName(VERBOSE_1, 'verbose1')
##     logging.addLevelName(VERBOSE_2, 'verbose2')
    log_level = logging.WARNING
    if options.verbosity == 1:
        log_level = logging.INFO
##     elif options.verbosity == 2:
##         log_level = logging.VERBOSE_1
##     elif options.verbosity == 3:
##         log_level = logging.VERBOSE_2
    elif options.verbosity >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level,
                        format='%(levelname)-7s%(name)s  %(message)s')

    organisms_folders = [os.path.join(options.sp_folder, f) \
                         for f in os.listdir(options.sp_folder)]

    if options.check_with_blast:
        check_with_blast(organisms_folders, options.temp)
    
if __name__ == '__main__':
    main()
