#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import optparse

def load_nodes(nodes):
    data = {}
    c = None
    p = None
    with open(nodes, 'r') as iff:
        for line in iff:
            c, p = [i.strip() for i in line.split('|')[:2]]
            if p in data:
                data[p].append(c)
            else:
                data[p] = [c]
    return data

def get_sub_taxids(tree, taxid, only_leaves=False):
    '''Get children taxids in a recursive manner.
    '''
    def _do_leaves(tree, taxid, subs=set()):
        '''Get children leaves taxids in a recursive manner
        '''
        lsubs = subs
        try:
            for c in tree[taxid]:
                lsubs.update(_do_leaves(tree, c, lsubs))
        except KeyError, (e):
            return [taxid]
        return lsubs
    def _do_all(tree, taxid, subs=set()):
        '''Get all children taxids in a recursive manner
        '''
        lsubs = subs
        try:
            for c in tree[taxid]:
                lsubs.update(_do_all(tree, c, lsubs))
            lsubs.add(taxid)
        except KeyError, (e):
            return [taxid]
        return lsubs
    _do = _do_leaves if only_leaves else _do_all
    tmp = _do(tree, taxid)
    return tmp


def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--taxid',
                      dest='taxid',
                      help='taxid from wich all children nodes or parents nodes should be returned',
                      metavar='ID',
                      type='int')

    parser.add_option('-l', '--only_leaves',
                      dest='only_leaves', action='store_true', default=False,
                      help='only return those children nodes that are leaves')

    parser.add_option('-n', '--nodes_file',
                      dest='nodes_file',
                      help='file containing relations between nodes from NCBI',
                      metavar='FILE')

    parser.set_defaults(nodes_file = '/users/rg/mmariotti/libraries/nodes.dmp')

    (opt, args) = parser.parse_args()

    if not opt.taxid:
        parser.error('You must at least provide a taxid as an input')

    tree = load_nodes(opt.nodes_file)

    sub_nodes = get_sub_taxids(tree, str(opt.taxid),
                               only_leaves=opt.only_leaves)
    
    for i in sub_nodes:
        sys.stdout.write(i+'\n')


if __name__ == '__main__':
    main()
