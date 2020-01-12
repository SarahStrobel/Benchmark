import sys, re, string

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
    if re.search('[-]$',line):
        # reverse-complement sequence
        print line
        print line2.translate(complement)[::-1]        
    elif re.search('[+]$',line): 
        # sense sequence
        print line
        print line2
    elif re.search('nostrand::',line):
        # neither strand
        # don't output anything
        pass
    else:
        raise RuntimeError("could not classify line \"%s\" in terms of strandness" % (line))
