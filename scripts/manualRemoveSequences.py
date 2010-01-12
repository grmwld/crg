#!/usr/bin/env python2.6

import os
import sys
import optparse
import AGBio.io.Fasta as Fasta

def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='fasta file in which selenoproteins should be looked for.',
                       metavar='FILE' )

    parser.add_option( '-a', '--alignmentfile',
                       dest='alfilename',
                       help='alignment file used when details are requested.',
                       metavar='FILE' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='base output filename',
                       metavar='FILE' )

    parser.add_option( '-f', '--keep_prefilter',
                       dest='keepprefilter',
                       help='prefilters all sequences that have the given pattern in their name and keep them.',
                       metavar='PATTERN' )
    
    parser.add_option( '-F', '--throw_prefilter',
                       dest='throwprefilter',
                       help='prefilters all sequences that have the given pattern in their name and throw them.',
                       metavar='PATTERN' )

    parser.set_defaults( keepprefilter = False,
                         throwprefilter = False,
                         alfilename = False )

    (options, args) = parser.parse_args()

    with open(options.inputfilename, 'r') as inf:
        sequences = Fasta.loadSequences(inf)

    if options.alfilename:
        with open(options.alfilename, 'r') as alf:
            alignment = Fasta.Alignment(Fasta.loadSequences(alf))
        nrdetail = alignment.findPositions(('U','C','-'), False)
        rdetail = alignment.findPositions(('U','C','-'), True)
        
    if options.keepprefilter:
        kpatterns = options.keepprefilter.split(',')
    if options.throwprefilter:
        tpatterns = options.throwprefilter.split(',')

    kept_seq = []
    thrown_seq = []
    man_check_list = []

    for seq in sequences:
        kept = False
        thrown = False
        if options.keepprefilter:
            for pattern in kpatterns:
                if pattern in seq.header:
                    kept_seq.append(seq)
                    kept = True
        if options.throwprefilter:
            for pattern in tpatterns:
                if pattern in seq.header:
                    thrown_seq.append(seq)
                    thrown = True
        if not kept and not thrown:
            man_check_list.append(seq)

    for idx, seq in enumerate(man_check_list):
        gi = seq.header.split('|')[1]
        choice = 'r'
        decided = False
        print seq.header
        while not decided:
            choice = raw_input('# '+str(idx)+' / '+str(len(man_check_list))+' -- Keep ? ')
            if choice == 'y':
                kept_seq.append(seq)
                decided = True
            elif choice == 'n':
                thrown_seq.append(seq)
                decided = True
            elif choice == 's':
                os.system('fetch_seq.g -v TITLE="'+gi+'" -v ALL=1 '+options.inputfilename )
            elif choice == 'd' and options.alfilename:
                print
                print 'General Detail :'
                for pos in rdetail['U']:
                    sys.stdout.write('    '+str(pos) + ' ')
                    for xpos in rdetail:
                        try:
                            sys.stdout.write(str(xpos)+': ')
                            sys.stdout.write(str(len(rdetail[xpos][pos])) + ' ; ')
                        except KeyError:
                            sys.stdout.write('0 ; ')
                    sys.stdout.write('\n')
                print
                for seqal in alignment:
                    if seqal.header == seq.header:
                        tmpseq = seqal
                tmppos = [i for i, x in enumerate(tmpseq.sequence) if x == 'U']
                print 'In the alignment provided :'
                print '    U :', tmppos
                print '    U in those positions :', [len(rdetail['U'][(l,)]) for l in tmppos]
                print '    C in those positions:', [len(rdetail['C'][(l,)]) for l in tmppos]
                print '    - in those positions:', [len(rdetail['-'][tuple((l,))]) for l in tmppos]
            elif choice == 'q':
                while cc not in ('y', 'n'):
                    cc = raw_input('Manual quit. Would you like to save your changes ? ')
                    if cc in 'y':
                        pass
                    if cc in 'n':
                        pass
            else:
                print 'Wrong command'

    with open(options.outputfilename, 'w') as of:
        kept_seq.prints(of, 80)


if __name__ == '__main__':

    main()
