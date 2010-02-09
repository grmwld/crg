import sys
import random

def genmer( size ):

    mer = ''

    for i in range( size ):
    
        nt = random.randint( 0, 3 )

        if nt == 0:
            mer += 'a'
        elif nt == 1:
            mer += 'c'
        elif nt == 2:
            mer += 'g'
        else:
            mer += 't'

    return mer


def main():

#    if sys.argc != 3:
#        sys.exit(1)

    for i in range( int( sys.argv[1] ) ):

        print genmer( int( sys.argv[2] ) )


if __name__ == '__main__':

    main()
        

    
