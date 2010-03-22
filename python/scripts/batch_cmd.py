#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

    parser.add_option( '-n', '--num_batch',
                       dest='num_batch', type='int',
                       help='Number of simultaneous steps allowed to be performed',
                       metavar='INT' )
    
    parser.set_defaults(num_batch = 2)

    (opt, args) = parser.parse_args()

    pool = multiprocessing.Pool(processes=2)
    result = pool.map(worker,
                      ['blastp -query ~/tt -out ~/tt.out -db nr',
                       'blastp -query ~/ttttt -out ~/ttttt.out -db nr'])

if __name__ == '__main__':
    main()
