import sys, re, string, argparse

parser = argparse.ArgumentParser(description='Transform output of bedtools getfasta\n')
parser.add_argument('-revRev', action='store_true',help='reverse reverseness, i.e. write reverse complement for the plus strand (+), and the sense strand for the negative strand (-)')
parser.add_argument('-noRev', action='store_true',help='do not reverse complement sequences')
args = parser.parse_args()

complement = string.maketrans("ACGTU","TGCAA")

while True:
    line = sys.stdin.readline()
    if not line:
        break
    line=line.rstrip()
    if not re.search('^>',line):
        raise RuntimeError("line \"%s\" doesnt start with >, but expected descriptor" % (line))
    line2=sys.stdin.readline()
    if not line2:
        raise RuntimeError('descriptor line without sequence')
    line2=line2.rstrip()
    line=re.sub('::.*$','',line) # apparently some versions of 'bedtools getfasta' add '::' and then the coordinates, and some don't.  we just remove them
    doOutput=False
    doRevComp=False
    if re.search('[-]$',line):
        doOutput=True
        doRevComp=True # that'd be the default
    elif re.search('[+]$',line): 
        doOutput=True
        doRevComp=False
    elif re.search('nostrand$',line):
        # neither strand
        # don't output anything
        doOutput=False
    else:
        raise RuntimeError("could not classify line \"%s\" in terms of strandness" % (line))
    if doOutput:
        if args.revRev:
            doRevComp=not doRevComp
        if args.noRev and doRevComp: # user doesn't want this
            doRevComp=False
        if doRevComp:
            # reverse-complement sequence
            print(line)
            print(line2.translate(complement)[::-1])
        else:
            # sense sequence
            print(line)
            print(line2)
