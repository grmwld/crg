#! /usr/bin/python
import sys
import os
import optparse
import shutil
import logging
from cStringIO import StringIO
sys.path.append('/users/rg/agrimaldi/Code/python/python/lib/')
from AGBio.ncbi.BlastWrappers import *
from AGBio.Utilities import *
sys.path.append('/users/rg/mmariotti/libraries')
sys.path.append('/users/rg/mmariotti/scripts')


def main():

    parser = optparse.OptionParser()

    parser.add_option('-d', '--database_check',
                      dest='dbc',
                      help='location of the database that should be used for checking.',
                      metavar='DB')

    parser.add_option('-i', '--input',
                      dest='fasta_query',
                      help='Full path to a .seq file inside a selenoprofiles output folder.',
                      metavar='FILE')

    parser.add_option('-a', '--ncore',
                      dest='ncore',
                      type='int',
                      help='number of cores to use for the blast.',
                      metavar='INT')

    parser.add_option('-b', '--blast_flavour',
                      dest='blast_flavour',
                      help='what kind of blast should be performed ?',
                      metavar='BLAST')

    parser.add_option('-v', '--verbosity',
                      dest='verbosity', action='count',
                      help='set verbosity level')

    parser.add_option('-T', '--temp',
                      dest='temp',
                      help='temporary folder.',
                      metavar='DIR')

    parser.set_defaults(temp = '/home/agrimaldi/temp',
                        ncore = 1,
                        blast_flavour = 'blastp')

    (opts, args) = parser.parse_args()

    sp_res = opts.fasta_query

    blast_output = os.path.join(opts.temp, 'tmpblast.xml')

    blaster = BlastAllWrapper(sp_res, blast_output,
                              flavour=opts.blast_flavour,
                              db=opts.dbc, gis=True, ncore=opts.ncore)

    print blaster.cline

    blaster.run()

    gis_file = os.path.join(opts.temp, 'tmpgis')
    os.system("read_blast_xml -n 50 -i " + blast_output + " | grep S_ID | gawk -F'|' '{print $2}' > " + gis_file)

    blast_res = os.path.join(opts.temp, 'blast_res.fasta')
    os.system("fastacmd -i "+ gis_file +" -o "+ blast_res+ " -d " + opts.db)

    blast_res_aln = os.path.join(opts.temp, 'blast_res.faln')
    os.system("custom_align.py -i "+blast_res + " -o " + blast_res_aln + " -m tcoffee -r gi -a 4 -u U:C -t " + opts.temp)


if __name__ == '__main__':
    main()
