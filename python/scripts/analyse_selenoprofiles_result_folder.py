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
import cStringIO
from cStringIO import StringIO
sys.path.append('/users/rg/agrimaldi/Code/python/python/lib/')
from AGBio.io.common import *
from AGBio.Utilities import *
from AGBio.selenoprofiles_tools.results_analyser import *
from AGBio.selenoprofiles_tools.files_analysers.p2g import *
from AGBio.selenoprofiles_tools.files_analysers.ali import *


def check_with_blast(sp_parser=None, org_folder=None, temp=None):
    sp_parser = sp_parser
    sp_parser.parse_p2gs()
    for p2g in sp_parser.p2g:
        query_file = org_folder.sub_file(genTempfilename(prefix=''))
        logging.info('Creating query file : '+query_file.abspath)
        query_file.create()
        output_file = org_folder.sub_file(p2g.filename)
        print output_file, output_file.abspath
        raw_input()
        with open(query_file.abspath, 'w') as off:
            p2g.result.target.fasta().prints(off)
        cmd = ' '.join(['check_with_blast.py', '-v',
                        '-o', output_file.abspath,
                        '-q', query_file.abspath,
                        '-b', 'blastp',
                        '-d', '/seq/databases/nr_uncompressed/nr',
                        '-a', '4',
                        '-n', '10',
                        '-T', org_folder.abspath])
        subprocess.call(cmd, shell=True)

def compute_u_stats(sp_parser=None, org_folder=None):
    sp_parser = sp_parser
    sp_parser.ali_stats()
    logging.debug(str([s.filename for s in sp_parser.ali]))
    for ali in sp_parser.ali:
        try:
            prot_name = '.'.join(os.path.basename(ali.filename).split('.')[:-1])
            tmp_prot_folder = org_folder.sub_folder('.'.join([prot_name, 'rad']))
            logging.info('Creating temporary protein folder : ' \
                         +tmp_prot_folder.abspath)
            tmp_prot_folder.create()
        except OSError, (e):
            if e.errno == 17:
                pass
        tmp_out_file = tmp_prot_folder.sub_file('ali.stats')
        logging.info('Creating temporary output file : '+tmp_out_file.abspath)
        logging.info('Statistics computed for file '+ali.filename)
        with open(tmp_out_file.abspath, 'w+') as off:
            ali.u_redundant.prints('U', off)
            off.seek(0)
            for line in off:
                logging.info(line.strip())

def write_to_result_folder():
    pass

def main():

    parser = optparse.OptionParser()

    parser.add_option('-f', '--folder',
                      dest='sp_folder',
                      help='selenoprofiles output folder',
                      metavar='DIR')

    parser.add_option('-c', '--check_with_blast',
                      action='store_true', default=False,
                      dest='check_with_blast',
                      help='check candidates with blast')

    parser.add_option('-u', '--u_stats',
                      action='store_true',
                      dest='u_stats',
                      help='Generate statistics about U dispertion')

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

    log_level = logging.WARNING
    if options.verbosity == 1:
        log_level = logging.INFO
    elif options.verbosity >= 2:
        log_level = logging.DEBUG
    log_format = '%(levelname)-7s%(asctime)s  %(filename)s  %(message)s'
    logging.basicConfig(level=log_level, format=log_format)

    organisms_folders = [os.path.join(options.sp_folder, f) \
                         for f in os.listdir(options.sp_folder)]

    tmp_run_folder = FolderWrapper(genTempfilename(options.temp, 'spa'))
    logging.info('Creating run folder : '+tmp_run_folder.abspath)
    tmp_run_folder.create()
    
    delete_tmp = False
    try:
        for f in organisms_folders:
            tmp_org_folder = tmp_run_folder.sub_folder(os.path.basename(f))
            logging.info('Creating temporary organism folder : ' \
                          +tmp_org_folder.abspath)
            tmp_org_folder.create()
            sp_parser = GenomeFolderParser(f)
            sp_parser.parse(sec=True, cys=True, thr=True,
                            arg=True, uga=True)
            if options.check_with_blast:
                check_with_blast(sp_parser, tmp_org_folder)
            if options.u_stats:
                logging.info('Computing statistics on '+sp_parser.rootdir)
                compute_u_stats(sp_parser, tmp_org_folder)
    except KeyboardInterrupt:
        sys.stderr.write('\nManual quit during processing.\n')
        delete_tmp = raw_input('Delete temporary files ? [y/n] ').lower()\
                     in ['y', 'yes']
    finally:
        if delete_tmp:
            tmp_run_folder.delete(recursive=True)
            
if __name__ == '__main__':
    main()
