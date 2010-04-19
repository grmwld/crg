#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import multiprocessing
import optparse
import logging
import time

def worker(cmd):
    oc = os.system(cmd)
    log_emitter((cmd, oc))
    return oc

def init_jobs(filename):
    jobs = []
    with open(filename, 'r') as iff:
        for line in iff:
            jobs.append(line.strip())
    return jobs

def emitter_builder(logger):
    def emitter(result):
        if result[1] == 0:
            logger.info(result[0])
        else:
            logger.error(result[0])
    return emitter


def main():
    
    parser = optparse.OptionParser()

    parser.add_option('-f', '--file',
                      dest='filecmd',
                      help='bash file containing the commands to execute.',
                      metavar='FILE')

    parser.add_option('-n', '--num_proc',
                      dest='num_proc', type='int',
                      help='Number of simultaneous processes.',
                      metavar='INT')

    parser.add_option('-L', '--logfile',
                      dest='logfile',
                      help='File to use for logging',
                      metavar='FILE')
    
    parser.set_defaults(num_proc = multiprocessing.cpu_count(),
                        logfile = None)

    (opt, args) = parser.parse_args()

    filecmd = os.path.abspath(os.path.expanduser(opt.filecmd))
    logfile = opt.logfile
    if not logfile:
        logfile = '.'.join([filecmd, 'LOG'])

    logger = logging.getLogger(filecmd)
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(logfile)
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)

    global log_emitter
    log_emitter = emitter_builder(logger)

    pool = multiprocessing.Pool(processes=opt.num_proc)
    logger.info('======  START  ======')
    logger.info(' '.join(sys.argv))
    try:
        result = pool.map_async(worker, init_jobs(filecmd)).get(999999)
    except KeyboardInterrupt:
        pool.terminate()
    finally:
        logger.info('======  END  ======')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
