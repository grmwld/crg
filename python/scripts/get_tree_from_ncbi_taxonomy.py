#!/usr/bin/env python 

import os
import sys 

sys.path.append('/users/rg/mmariotti/libraries/')

from string import strip
from ete2 import TreeNode, Tree

help_msg="""
This program was built modifying the script NCBI_tree_of_life.py by Jaime Huerta Cepas.

It reads the NCBI taxonomy database (text file format) and reconstructs the whole taxonomy tree of life. 
Then, it prunes all the nodes that do not correspond to the species in input.
In the end, you have a tree of all species in input according to ncbi taxonomy.

Note that reconstructing the tree from the raw NCBI format may take some
minutes.

usage:

get_tree_from_ncbi_taxonomy.py species_file

species_file must contain the wanted species, one per line.
Options:

-s	after computing the final tree, open the interactive ete environment
-w [filename]	after computing the final tree, it outputs it in newick format
"""

def main():
    from MMlib import bash, command_line
    def_opt= {  'i':'species_file', 'names':'/users/rg/mmariotti/libraries/names.dmp', \
                'nodes':'/users/rg/mmariotti/libraries/nodes.dmp', 'alt':0, 'w':0, 's':0} 
    opt=command_line(def_opt, help_msg, 'i' )
    species_file=opt['i']
    names_file  =opt['names']
    nodes_file  =opt['nodes']


    t, id2node, all_wanted_species  = load_NCBI(species_file, names_file, nodes_file )
    all_wanted_nodes = get_wanted_nodes(t, id2node, all_wanted_species)
    t= prune_tree(t,  all_wanted_nodes)


    if opt['s']:
            t.show()
    if opt['w']:
            if opt['w']==1:
                    opt['w']=raw_input('Please specify a name for the newick output file:').strip()
            open(opt['w'], "w").write( t.write())


    

def load_NCBI(species_file, names_file, nodes_file ):
    if not os.path.isfile(species_file):
            print "ERROR "+species_file+' can\'t be read. Exiting... '
            sys.exit(8)

    all_wanted_species={} # species_name:   taxid (string)

    print "Reading wanted species from file: "+species_file
    ifile=open(species_file, 'r')
    for iline in ifile:
            species_name=iline.strip()
            all_wanted_species[species_name]=-1
    ifile.close()


    # This sets Unbuffered stdout/auto-flush
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    id2node= {}
    node2parentid = {}
    all_ids = set([])
    all_nodes = []
    id2name= {}

    # Loads info from NCBI taxonomy files
    if os.path.exists(nodes_file):
        NODESFILE = open(nodes_file)
    elif os.path.exists(nodes_file+"bz2"):
        import bz2
        NODESFILE = bz2.BZ2File(nodes_file+'.bz2')
    else:
        print nodes_file+' file is missing. '
        sys.exit(8)

    if os.path.exists(names_file):
        NAMESFILE = open(names_file)
    elif os.path.exists(names_file+"bz2"):
        import bz2
        NAMESFILE = bz2.BZ2File(names_file+'.bz2')
    else:
        print names_file +' file is missing. '
        sys.exit(8)


    # Reads taxid/names transaltion
    print 'Loading species names from "names.dmp" file...',
    for line in NAMESFILE:
        # lines are redundant. synonyms are on different lines defined by the same id. So, we store only lines with "scientific name".
        line = line.strip()
        fields = map(strip, line.split("|"))
        nodeid, name = fields[0], fields[1]

        if all_wanted_species.has_key(name):
            all_wanted_species[name]=nodeid

        if fields[3]=='scientific name':
            #storing name that will appear afterwards in the ete2 node
            id2name[nodeid] = name

    print len(id2name)

    any_species_is_missing=0
    for species_name in all_wanted_species:
            if all_wanted_species[species_name]==-1:
                    print "ERROR the species name \""+species_name+"\" was not found!"
                    any_species_is_missing=1
    if any_species_is_missing:                    
        sys.exit(9)


    # Reads node connections in nodes.dmp
    print 'Loading node connections from "nodes.dmp" file...', 
    for line in NODESFILE:
        line = line.strip()
        fields = map(strip, line.split("|"))
        nodeid, parentid = fields[0], fields[1]
        if nodeid =="" or parentid == "":
            raw_input("Wrong nodeid!")

        # Stores node connections
        all_ids.update([nodeid, parentid])

        # Creates a new TreeNode instance for each new node in file
        n = TreeNode()
        # Sets some TreeNode attributes
        n.add_feature("name", id2name[nodeid])
        n.add_feature("taxid", nodeid)

        # updates node list and connections
        node2parentid[n]=parentid
        id2node[nodeid] = n

    print len(id2node)


    # Reconstruct tree topology from previously stored tree connections
    print 'Reconstructing tree topology...'
    for node in id2node.itervalues():
        parentid = node2parentid[node]
        parent = id2node[parentid]
        # node with taxid=1 is the root of the tree
        if node.taxid == "1":
            t = node
        else:
            parent.add_child(node)
    return t, id2node, all_wanted_species

def get_wanted_nodes(t, id2node, all_wanted_species):
    #computing list to give to the pruning method
    all_wanted_nodes= set()
    for species_name in all_wanted_species:
        all_wanted_nodes.add(id2node[all_wanted_species[species_name]])
    return all_wanted_nodes
    
def prune_tree(t, nodes_to_keep):
    """ Fixed (and faster) prunning algorithm. Use this until I fix
    the problem within the main ETE branch.

    'nodes_to_keep' must be the list of node instances that you want
    to keep in the final tree. All nodes must be leaves, if not, they
    are automatically converted into leaves by removing their
    children.

    So far, this function is quite verbose. Printing slows down a bit
    the process, but you can follow the progress...
    """ 
    print "Getting tree path..."
    # Converts to set to speed up searches
    if type(nodes_to_keep) == set:
        to_keep=nodes_to_keep
    else:
        to_keep=set(nodes_to_keep)

    print "Checking that all nodes are leaves..."
    not_leaves = [n for n in nodes_to_keep if not n.is_leaf()]
    if len(not_leaves)>0:
        print "\nFixing", len(not_leaves), "non-leaf nodes..."
        # Converts all internal species nodes into leaves by removing all
        # their sub-species or strains
        for nl in not_leaves: 
            for c in nl.get_children():
                c.detach()
    to_detach = []
    
    print "Selecting unused nodes"
    counter = 0
    for node in t.traverse("postorder"):
        print "\r", counter,
        counter +=1
        for c in node.children:
            if c in to_keep:
                to_keep.add(node)
                break
        if node not in to_keep:
            to_detach.append(node)
            for c in node.children:
                to_detach.remove(c)
    print "\nDetaching", len(to_detach), "nodes"
    counter = 0
    for node in to_detach:
        print "\r", counter,
        counter +=1
        node.detach()
    print "\nFixing", len(to_keep), "orphan nodes"
    counter = 0
    for node in to_keep:
        print "\r", counter,
        counter +=1
        if len(node.children) == 1:
            node.delete()
    return t


if __name__=="__main__":
    main()
