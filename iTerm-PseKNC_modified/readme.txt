The files in this directory are taken from the iTerm-PseKNC software package available at http://lin-group.cn/server/iTerm-PseKNC/download.php

The iTerm-PseKNC.py file was modified by Sarah Strobel so that it was easier to integrate into the benchmark's scripts.

The iTerm-PseKNC software package was described in this paper: https://www.ncbi.nlm.nih.gov/pubmed/30247625

For iTerm-PseKNC to run the user must have LibSVM installed.
Additionally it uses Python3 and the libraries sys, subprocess, os, itertools, pandas and numpy.

The original content of the readme.txt file from the iTerm-PseKNC software package follows:

1. Install the svm environment by the following intruction in current path:
make
2. The predicted result of the sequence can be obtained by running the following intruction:
    python iTerm-PseKNC.py inputfile
the inputfile is in FASTA format and stores the sequence to be predicted, and then the predict
result is writen in "result.txt"
