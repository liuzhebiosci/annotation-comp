# annotation-comp
This README file discusses the Python code and datasets associated with the paper:
Article title: A semi-automated genome annotation comparison and integration scheme
Authors: Zhe Liu, Hongwu Ma and Igor Goryanin

This directory contains supporting python codes and datasets to perform annotation comparison and determination tasks. 

Datasets
The ‘datasets’ folder contains the data used in our research, in particular, the contents of data used are listed below:
1.	This file ‘pfam_mapped_dic.txt’ contains the Pfam dictionary constructed from Pfam database [1];
2.	This file ‘gene_info_dic.txt’ contains the gene symbol dictionary constructed from Entrez gene database [2].
3.	The ‘database data’ folder contains the automated annotation results and the OMA data for each specific genome.

Python Code
The ‘code’ folder contains supporting Python code to pre-process, compare and determine genome annotations from multiple annotation sources. It should be noted that BLAST software should be installed prior to the usage of this code. BLAST can be found on NCBI BLAST website (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). ‘pre-treat.py’ and ‘pre.py’ are used to format and retrieve the OGA annotations from OMA database. ‘annoComp.py’, ‘annoDeter.py’ and ‘pre.py’ are used to compare the annotations and derive a consensus annotation result based on the comparison result. The detailed descriptions of how to use the code can be found in the comments.

Reference
1.	Pfam data ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam26.0/database_files/.
2.	NCBI data ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz.
3.	Schneider A, Dessimoz C, Gonnet GH: OMA Browser—Exploring orthologous relations across 352 complete genomes. Bioinformatics 2007, 23(16):2180-2182.


