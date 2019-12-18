# Benchmark for terminator prediction using RNA-Seq and Term-Seq data


Download [Termi script](https://github.com/SarahStrobel/Benchmark/blob/master/Termi.sh) and run with `sh Termi.sh` <br/>
Download [Benchmark script](https://github.com/SarahStrobel/Benchmark/blob/master/Termi_Benchmark.sh) and run with `sh Termi_Benchmark.sh` <br/>

Attention:
For the script to run you need the RNA-Seq FASTQ files by Warrier et al. that were kindly shared with us and put them in the folder "yourWorkingDirectory/Termi/RNASeq" that is made by the Termi.sh script.

### Used Programs:<br/>

* [git](https://git-scm.com/)<br/>
<t/>version-control system for tracking changes in source code<br/>
* [SRA-Toolkit-2.9.6-ubuntu64](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)<br/>
<t/>enables reading ("dumping") of sequencing files from the SRA database<br/>
* [fastp-0.20.0](https://github.com/OpenGene/fastp)<br/>
<t/>all-in-one preprocessing for FastQ files<br/>
* [Novoalign-V3.02.07.Linux3.0](http://www.novocraft.com/products/novoalign/)<br/>
<t/>mapping of short reads onto a reference genome<br/>
* [Samtools-1.9](http://www.htslib.org/download/)<br/>
<t/>set of utilities for interacting with and post-processing short sequence read alignments<br/>
* [Bedtools-v2.27.1](https://bedtools.readthedocs.io/en/latest/index.html)<br/>
<t/>tools for a wide-range of genomics analysis tasks<br/>
* [BLAST-2.9.0+-x64-linux](https://blast.ncbi.nlm.nih.gov/Blast.cgi)<br/>
<t/>Basic Local Alignment Search Tool finds regions of local similarity between sequences<br/>
* [Infernal-0.81](http://eddylab.org/infernal/)<br/>
<t/>Infers RNA Alignment using covariance models (CMs)<br/>
* [Easel](http://eddylab.org/infernal/)<br/>
<t/>Integrated in Infernal-1.1.2<br/>
* [RNIE](https://github.com/ppgardne/RNIE)<br/>
<t/>Terminator prediction software<br/>
* [RNAmotif-3.1.1](http://casegroup.rutgers.edu/casegr-sh-2.5.html)<br/>
<t/>Terminator prediction software<br/>
* [iTerm-PseKNC](http://lin-group.cn/server/iTerm-PseKNC/download.php)<br/>
<t/>Terminator prediction software<br/>
* [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/)<br/>
<t/>Integrated software for support vector classification<br/>

### Used Python2 packages:<br/>

`pip install packagename`<br/>
will be installed to /usr/bin/ on Linux or in your Python installation on Windows<br/>
<br/>
`pip install packagename --user`<br/>
will be installed to ~/.local/bin/ on Linux or to %AppData\Python\Scripts\ on Windows<br/>

* [argparse](https://docs.python.org/2/library/argparse.html)<br/>
* [bisect](https://docs.python.org/2/library/bisect.html)<br/>
* [collections](https://docs.python.org/2/library/collections.html)<br/>
* [glob](https://docs.python.org/2/library/glob.html)<br/>
* [itertools](https://docs.python.org/2/library/itertools.html)<br/>
* [math](https://docs.python.org/2/library/math.html)<br/>
* [matplotlib](https://matplotlib.org/)<br/>
* [numpy](https://numpy.org/)<br/>
* [operator](https://docs.python.org/2/library/operator.html)<br/>
* [os](https://docs.python.org/2/library/os.html)<br/>
* [re](https://docs.python.org/2/library/re.html)<br/>
* [sys](https://docs.python.org/2/library/sys.html)<br/>
* [tabulate](https://pypi.org/project/tabulate/)<br/>

### Python3 packages used by iTerm-PseKNC_modified.py:<br/>

`pip3 install packagename`<br/>
will be installed to /usr/bin/ on Linux or in your Python installation on Windows<br/>
<br/>
`pip3 install packagename --user`<br/>
will be installed to ~/.local/bin/ on Linux or to %AppData\Python\Scripts\ on Windows<br/>

* [itertools](https://docs.python.org/3/library/itertools.html)<br/>
* [numpy](https://numpy.org/)<br/>
* [os](https://docs.python.org/3/library/os.html)<br/>
* [pandas](https://pandas.pydata.org/)<br/>
* [subprocess](https://docs.python.org/3/library/subprocess.html)<br/>
* [sys](https://docs.python.org/3/library/sys.html)<br/>

### B.subtilis, E.faecalis and L.monocytogenes RNA-Seq / Term-Seq Data (Dar et al., 2016):<br/>
* [Paper](https://www.ncbi.nlm.nih.gov/pubmed/27120414)<br/>
* [fastq files](https://www.ncbi.nlm.nih.gov/sra?term=ERP014057)<br/>


### S.pneumoniae Term-Seq Data (Warrier et al., 2018):<br/>
* [Paper](https://www.ncbi.nlm.nih.gov/pubmed/30517198)<br/>
* [fastq files](https://www.ncbi.nlm.nih.gov/sra/?term=SRP136114)<br/>



### Additional Programs:<br/>

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<br/>
* [BamQC](https://github.com/s-andrews/BamQC)<br/>
* [IGV](https://software.broadinstitute.org/software/igv/)<br/>
* [RNAfold](http://rna.tbi.univie.ac.at/)<br/>
* [mfold](http://unafold.rna.albany.edu/?q=mfold)<br/>
* [Segemehl](https://www.bioinf.uni-leipzig.de/Software/segemehl/)<br/>
* [RNAmotif](http://casegroup.rutgers.edu/casegr-sh-2.5.html)<br/>


### Additional Python3 packages:<br/>

* [scikit-learn](https://scikit-learn.org/stable/)<br/>
* [mlxtend](http://rasbt.github.io/mlxtend/api_subpackages/mlxtend.plotting/)<br/>




