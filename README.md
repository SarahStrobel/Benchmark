# Benchmark for terminator prediction using RNA-Seq and Term-Seq data



Download [Benchmark script](https://github.com/SarahStrobel/Benchmark/blob/master/Benchmark.sh) and run with `sh Benchmark.sh` <br/>


### Used Programs:<br/>

* [git](https://git-scm.com/)<br/>
<t/>version-control system for tracking changes in source code<br/>
* [SRA-Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)<br/>
<t/>enables reading ("dumping") of sequencing files from the SRA database<br/>
* [fastp](https://github.com/OpenGene/fastp)<br/>
<t/>all-in-one preprocessing for FastQ files<br/>
* [Novoalign](http://www.novocraft.com/products/novoalign/)<br/>
<t/>mapping of short reads onto a reference genome<br/>
* [Samtools](http://www.htslib.org/download/)<br/>
<t/>set of utilities for interacting with and post-processing short sequence read alignments<br/>
* [Bedtools](https://bedtools.readthedocs.io/en/latest/index.html)<br/>
<t/>tools for a wide-range of genomics analysis tasks<br/>
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)<br/>
<t/>Basic Local Alignment Search Tool finds regions of local similarity between sequences<br/>
* [Infernal](http://eddylab.org/infernal/)<br/>
<t/>Infers RNA Alignment using covariance models (CMs)<br/>
* [RNIE](https://github.com/ppgardne/RNIE)<br/>
<t/>Terminator prediction software<br/>


### Used Python3 packages:<br/>

`pip install packagename`<br/>
will be installed to /usr/bin/ on Linux of or in your Python installation on Windows<br/>
<br/>
`pip install packagename --user`<br/>
will be installed to ~/.local/bin/ on Linux or to %AppData\Python\Scripts\ on Windows<br/>

* [glob](https://docs.python.org/3/library/glob.html)<br/>
* [numpy](https://numpy.org/)<br/>
* [matplotlib](https://matplotlib.org/)<br/>
* [argparse](https://docs.python.org/3/library/argparse.html)<br/>
* [re](https://docs.python.org/3/library/re.html)<br/>
* [collections](https://docs.python.org/3/library/collections.html)<br/>
* [sys](https://docs.python.org/3/library/sys.html)<br/>
* [os](https://docs.python.org/3/library/os.html)<br/>
* [operator](https://docs.python.org/3/library/operator.html)<br/>
* [tabulate](https://pypi.org/project/tabulate/)<br/>
* [bisect](https://docs.python.org/3.0/library/bisect.html)<br/>



### B.subtilis, E.faecalis and L.monocytogenes RNA-Seq / Term-Seq Data (Dar et al., 2016):<br/>
* [Paper](https://www.ncbi.nlm.nih.gov/pubmed/27120414)<br/>
* [fastq files](https://www.ncbi.nlm.nih.gov/sra?term=ERP014057)<br/>



### S.pneumoniae Term-Seq Data (Warrier et al., 2018):<br/>
* [Paper](https://www.ncbi.nlm.nih.gov/pubmed/30517198)<br/>
* [fastq files](https://www.ncbi.nlm.nih.gov/sra/?term=SRP136114)<br/>



### Additional Programs:<br/>

* [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<br/>
* [IGV](https://software.broadinstitute.org/software/igv/)<br/>
* [RNAfold](http://rna.tbi.univie.ac.at/)<br/>
* [mfold](http://unafold.rna.albany.edu/?q=mfold)<br/>
* [Segemehl](https://www.bioinf.uni-leipzig.de/Software/segemehl/)<br/>
* [RNAmotif](http://casegroup.rutgers.edu/casegr-sh-2.5.html)<br/>



### Additional Python3 packages:<br/>

* [scikit-learn](https://scikit-learn.org/stable/)<br/>
* [mlxtend](http://rasbt.github.io/mlxtend/api_subpackages/mlxtend.plotting/)<br/>




