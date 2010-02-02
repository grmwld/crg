#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import optparse
from ete2 import Tree

def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--input_file',
                      dest='infilename',
                      help='input filename, newick tree.',
                      metavar='FILE')
    
    parser.add_option('-g', '--gui',
                      action='store_true', dest='gui', default=False,
                      help='Use GUI for printing.')
    

    (options, args) = parser.parse_args()

    
    
    t = Tree(options.infilename)

    if options.gui:
        t.show()
    else:
        print t


if __name__ == '__main__':
    try:
        main()
    except Exception, (e):
        print e
        
