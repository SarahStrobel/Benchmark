# Benchmark for terminator prediction using RNA-Seq and Term-Seq data

## Replicating the benchmark from source

The scripts necessary to replicate the benchmark data are included here.  To run these, you will need:
* To use a UNIX system.  We have only tested this with CentOS (Linux)
* To install required software (see below)
* run the script <tt>Termi.sh</tt> (in this directory) to generate the benchmark data
* run the script <tt>Termi_Benchmark.sh</tt> (in this directory) to run the tests for RNIE, iTerm-PseKNC and the RNAmotif pattern by Lesnik et al

Attention: At the moment you'll also need the RNA-Seq FASTQ files by Warrier et al. that were kindly shared with us and put them in the folder "yourWorkingDirectory/Termi/RNASeq" that is made by the Termi.sh script.  We will remove this requirement once the files are available in the NCBI'a SRA.

## Required software

In order to run the benchmark scripts, you'll need the following software.  We have given the version of the software we used.  We expect that other versions of the software will likely give the same results, but, of course, we cannot guarantee this.

The case where the version number is definitely important is with the Infernal software (see below). We use the <tt>cmsearch</tt> and <tt>esl-shuffle</tt> commands from this software.  We used version 1.0.2, because the RNIE script is dependant on it.  If you use a newer version of Infernal, you will not be able to replicate the benchmark.  The appropriate versions of the <tt>cmsearch</tt> and <tt>esl-shuffle</tt> commands should be available in the <tt>$PATH</tt>.

The first think that the <tt>Termi.sh</tt> script does is to check that the necessary software is available.  The script will also check that the cmsearch and esl-shuffle commands are from the appropriate version of Infernal.  In other cases, the script does not check version information.

Here is the software and version numbers:
* bash shell
* Python 2 and 3 (Python 2 is required for exactly replicated the Benchmark results.  Similar results can be obtained using only Python 3.  If you desire this, edit the <tt>Termi.sh</tt> and <tt>Termi_Benchmark.sh</tt> scripts as follows.  Set <tt>PYTHON=python3</tt> and set <tt>SCRIPTS=Scripts-python2</tt>.  Note: the reason for the requirement for Python 2 is that the order of keys in the dictionary has changed in Python 3.  As a result, when multiple maximum values are identical, which one is selected will depend on the dictionary keys' order.  Although these differences are not significant, they do lead to non-identical results.)
* Perl. (This is required for the RNIE script)
* [git](https://git-scm.com/)<br/>
<t/>version-control system for tracking changes in source code<br/>
* [SRA-Toolkit-2.9.6-ubuntu64](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)<br/>
<t/>enables reading ("dumping") of sequencing files from the SRA database.  Note: the currently available package lacks the <tt>fasterq-dump</tt> program, which we need.  If the command is missing, the <tt>Termi.sh</tt> script will report that fact.  If you install a package and the <tt>fasterq-dump</tt> command is not available, we recommend that you install NCBI's pre-compiled verseion (available at the above URL) <br/>
* [fastp-0.20.0](https://github.com/OpenGene/fastp)<br/>
<t/>all-in-one preprocessing for FastQ files<br/>
* [Novoalign-V3.02.07.Linux3.0](http://www.novocraft.com/products/novoalign/)<br/>
<t/>mapping of short reads onto a reference genome.  We used the academic version, which is free to academic users.<br/>
* [Samtools-1.9](http://www.htslib.org/download/)<br/>
<t/>set of utilities for interacting with and post-processing short sequence read alignments<br/>
* [Bedtools-v2.27.1](https://bedtools.readthedocs.io/en/latest/index.html)<br/>
<t/>tools for a wide-range of genomics analysis tasks<br/>
* [BLAST-2.9.0+-x64-linux](https://blast.ncbi.nlm.nih.gov/Blast.cgi)<br/>
<t/>Basic Local Alignment Search Tool finds regions of local similarity between sequences<br/>
* [Infernal-1.0.2](http://eddylab.org/infernal/)<br/>
<t/>Infers RNA Alignment using covariance models (CMs).  <br/>
* [Easel](http://eddylab.org/infernal/)<br/>
<t/>Integrated in Infernal-1.0.2.  We require only the <tt>esl-shuffle</tt> command.  There is no pre-packaged software package available that includes the esl-shuffle command.  To install, download Infernal version 1.0.2 (at the above URL), run the <tt>configure</tt> script, run <tt>make</tt>, then cd into the <tt>easel</tt> subdirectory, and run <tt>make install</tt>, or <tt>cp</tt> the <tt>esl-shuffle</tt> command from <tt>easel/miniapps</tt> into a directory within your <tt>$PATH</tt>.<br/>

The following software is required for the <tt>Termi_Benchmark.sh</tt> script, but is not required for <tt>Termi.sh</tt>
* [RNAmotif-3.1.1](http://casegroup.rutgers.edu/casegr-sh-2.5.html)<br/>
<t/>Terminator prediction software<br/>
* [LIBSVM-3.23-3](https://www.csie.ntu.edu.tw/~cjlin/libsvm/)<br/>
<t/>Integrated software for support vector classification<br/>

The repository comes with the following software that is obtained from other sources.  You do not need to install it.
* [RNIE](https://github.com/ppgardne/RNIE)<br/>
<t/>Terminator prediction software<br/>
* [iTerm-PseKNC](http://lin-group.cn/server/iTerm-PseKNC/download.php)<br/>
<t/>Terminator prediction software.  We have modified the original Python script slightly so that it integrates better with our scripts.<br/>


### Required Python 2 packages:<br/>

The following Python 2 packages are required to exactly replicate the benchmark results.  If you wish to use Python 3, see the comment above.

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

### Required Python 3 packages<br/>

These Python 3 packages are required for the <tt>iTerm-PseKNC_modified.py</tt> script.

* [itertools](https://docs.python.org/3/library/itertools.html)<br/>
* [numpy](https://numpy.org/)<br/>
* [os](https://docs.python.org/3/library/os.html)<br/>
* [pandas](https://pandas.pydata.org/)<br/>
* [subprocess](https://docs.python.org/3/library/subprocess.html)<br/>
* [sys](https://docs.python.org/3/library/sys.html)<br/>

<!--
### B.subtilis, E.faecalis and L.monocytogenes RNA-Seq / Term-Seq Data (Dar et al., 2016):<br/>
* [Paper](https://www.ncbi.nlm.nih.gov/pubmed/27120414)<br/>
* [fastq files](https://www.ncbi.nlm.nih.gov/sra?term=ERP014057)<br/>


### S.pneumoniae Term-Seq Data (Warrier et al., 2018):<br/>
* [Paper](https://www.ncbi.nlm.nih.gov/pubmed/30517198)<br/>
* [fastq files](https://www.ncbi.nlm.nih.gov/sra/?term=SRP136114)<br/>



### Additional Programs:<br/>

We used the following software in our project, but these software are not necessary in order to run the benchmark scripts.
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<br/>
* [BamQC](https://github.com/s-andrews/BamQC)<br/>
* [IGV](https://software.broadinstitute.org/software/igv/)<br/>
* [RNAfold](http://rna.tbi.univie.ac.at/)<br/>
* [mfold](http://unafold.rna.albany.edu/?q=mfold)<br/>
* [Segemehl](https://www.bioinf.uni-leipzig.de/Software/segemehl/)<br/>
* [RNAmotif](http://casegroup.rutgers.edu/casegr-sh-2.5.html)<br/>


### Additional Python3 packages:<br/>

We used the following Python3 packages in our project, but these packages are not necessary to run the benchmark scripts.
* [scikit-learn](https://scikit-learn.org/stable/)<br/>
* [mlxtend](http://rasbt.github.io/mlxtend/api_subpackages/mlxtend.plotting/)<br/>

-->