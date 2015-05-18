# Supplementary Data for Boseon Kim et al. (2015)

by Hyeshik Chang &lt;hyeshik@snu.ac.kr&gt;

See Article: [B. Kim, M. Ha and L. Loeff et al. (2015) TUT7 controls the fate of precursor microRNAs by using three different uridylation mechanisms. *EMBO Journal.* *in press*](http://emboj.embopress.org/content/early/2015/05/15/embj.201590931.long) 

<figure>
<img src="https://raw.githubusercontent.com/hyeshik/bskim-2015-pre-miRNA/master/images/uridylating-tuts.png" width="142px">
<figurecaption>(c) Boseon Kim</figurecaption>
</figure>

This repository includes scripts and interactive notebooks for the data analysis of the high-throughput
sequencing of pre-miRNAs in HeLa cells. See the supplementary experimental procedures included in
[the manuscript](http://emboj.embopress.org/content/early/2015/05/15/embj.201590931.long)
for the detailed descriptions of the analytic processes.


## Saved notebook sessions

1. [Modification frequency throughout the 3' end positions of templated nucleotides](https://github.com/hyeshik/bskim-2015-pre-miRNA/blob/master/notes/plot-modifications-endpos-frequency-circles.ipynb)
1. [Bird's eye view for 3' end and uridylation patterns in control and TUT247 knockdown samples](https://github.com/hyeshik/bskim-2015-pre-miRNA/blob/master/notes/plot-uridylation-rate-change-by-position.ipynb)
1. [Global change of uridylation by knockdown of TUT2/4/7](https://github.com/hyeshik/bskim-2015-pre-miRNA/blob/master/notes/plot-global-uridylation-changes.ipynb)
1. [Relationships between 3' end trimming and uridylation?](https://github.com/hyeshik/bskim-2015-pre-miRNA/blob/master/notes/plot-trimming-and-uridylation.ipynb)

## Availability of the sequencing data

You can download the FASTQ files enclosing the reads from pre-miRNA sequencing
from the NCBI Gene Expression Omnibus (GEO) with an
[accession number GEO64482](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64482).

## Prerequisites

* [Python (≥ 3.4)](https://www.python.org)
* [FASTX_Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [cutadapt](https://code.google.com/p/cutadapt/)
* [snakemake](https://bitbucket.org/johanneskoester/snakemake)
* [numpy](http://www.numpy.org)
* [pandas](http://pandas.pydata.org)
* [IPython Notebook or Jupyter+IPython kernel](http://ipython.org)
* [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html)
* [BEDTools](https://github.com/arq5x/bedtools2)
* [samtools](http://www.htslib.org)
* [Jim Kent's Genome Browser utilities (faSize only)](https://genome.ucsc.edu/util.html)

## Instructions

1. Check out the repository

   ```
   $ git clone https://github.com/hyeshik/bskim-2015-pre-miRNA.git
   $ cd bskim-2015-pre-miRNA
   ```

1. Run the pipeline.

   ```
   $ snakemake
   ```

## The MIT License

Copyright (c) 2014 Hyeshik Chang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
