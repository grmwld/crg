#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import optparse
from ete2 import Tree, faces

def long2short(longname):
    ln = longname.split()
    return ln[0][0] + ln[1][:3]

def layout(node):
    if node.is_leaf():
        node.img_style['size'] = 10
        shortNameFace = faces.TextFace(long2short(node.name))
        longNameFace = faces.TextFace(node.name)
#        faces.add_face_to_node(longNameFace, node, column=0)
        faces.add_face_to_node(shortNameFace, node, column=0)
        
        if node.name.split()[1].startswith('m'):
            node.img_style["fgcolor"] = "#FFFFFF"
    else:
        node.img_style['size'] = 0
        

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
        t.show(layout)
    else:
        print t


if __name__ == '__main__':
    try:
        main()
    except Exception, (e):
        print e
        
