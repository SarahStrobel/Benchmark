The files in this directory are taken from the iTerm-PseKNC software package available at http://lin-group.cn/server/iTerm-PseKNC/download.php

The iTerm-PseKNC.py file was modified by Sarah Strobel and Zasha Weinberg so that it was easier to integrate into the benchmark's scripts.

The iTerm-PseKNC software package was described in this paper: https://www.ncbi.nlm.nih.gov/pubmed/30247625

For iTerm-PseKNC to run the user must have LibSVM installed.
Additionally it uses Python3 and the libraries sys, subprocess, os, itertools, pandas, argparse and numpy.

The command-line format is explained by the script with
python iTerm-PseKNC_modified.py -h

and is invoked by our bash script Termi_Benchmark.sh
