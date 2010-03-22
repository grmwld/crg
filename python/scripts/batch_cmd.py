#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import multiprocessing
import optparse
import time

def worker(cmd):
    os.system(cmd)
    return

def init_jobs(filename):
    jobs = []
    with open(filename, 'r') as iff:
        for line in iff:
            jobs.append(line.strip())
    jobs.reverse()
    return jobs

def main():
    
    parser = optparse.OptionParser()

    parser.add_option( '-f', '--file',
                       dest='filecmd',
                       help='bash file containing the commands to execute.',
                       metavar='FILE' )

    parser.add_option( '-n', '--num_proc',
                       dest='num_proc', type='int',
                       help='Number of simultaneous processes.',
                       metavar='INT' )
    
    parser.set_defaults(num_batch = 2)

    (opt, args) = parser.parse_args()

    pool = multiprocessing.Pool(processes=opt.num_proc)
    result = pool.map(worker,
                      init_jobs(opt.filecmd))

if __name__ == '__main__':
    main()
