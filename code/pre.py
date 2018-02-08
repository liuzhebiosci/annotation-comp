from __future__ import division

import numpy as np
import re
import operator
import numpy
import math
import nltk
import pdb
import os
import subprocess
import pickle
import shutil

from datetime import datetime
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import BlastallCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from numpy import zeros,dot
from numpy.linalg import norm
from nltk.tokenize import *
from nltk.stem.wordnet import WordNetLemmatizer
from nltk.corpus import wordnet

def parseQueryFile(fasta_f, consensus_anno, query_f):
    '''
    This function is used to extract the fasta file for a specific organism
    from the input file
    1. fasta_f: query genome fasta file
    2. consensus_anno: consensus annotation file, comprising three/four genome
    annotation results.
    3. query_f: output file 
    '''

    fasta=open(fasta_f, 'r')
    fasta_lines=fasta.readlines()
    fasta_len=len(fasta_lines)

    for i in xrange(fasta_len):
        fasta_lines[i]=re.sub("\n", "", fasta_lines[i])

    fasta_format=open(fasta_f+"2", 'w')
    fasta_format.write("".join(fasta_lines[1:]))
    fasta_format.close()

    fasta_format=open(fasta_f+"2", 'r')

    annofile=open(consensus_anno, 'r')
    annolines=annofile.readlines()
    annofile.close()

    seqfile=open(query_f, "w")
    gene_id=0

    for annoline in annolines:
        annoline=re.sub("\n", "", annoline)
        annolist=annoline.split("\t")
        annolist=[str(x) for x in annolist]
    ##    if annolist[5]=='peg' or annolist[9]=='CDS':
        gene_id=gene_id+1
        seqfile.write(">")

        seqfile.write(str(gene_id)+'|')
        seqfile.write('|'.join(annolist[:5]))
        seqfile.write('\n')
        seq_len=abs(int(annolist[2])-int(annolist[1]))+1

        if annolist[0]=="+":
            if (int(annolist[1])-1)>0:
                fasta_format.seek(int(annolist[1])-1)  
                dna=Seq(fasta_format.read(seq_len), generic_dna)
                seqfile.write(str(dna))

        if annolist[0]=="-":
            if (int(annolist[1])-1)>0:
                fasta_format.seek(int(annolist[1])-1)     
                dna=Seq(fasta_format.read(seq_len), generic_dna).reverse_complement()
                seqfile.write(str(dna))
                
        seqfile.write("\n")

    fasta_format.close()
    seqfile.close()

def parseDbFile(fasta_f, anno_f, seq_f):
    '''
    This function is used to extract the fasta file for a specific organism
    from the input file
    anno_f: organism specific annotation file, it can be found in
    the 'annotation files' folder
    fasta_f: database genome sequence .fna file.
    seq_f: database genome output file 
    '''
    
    fasta=open(fasta_f, 'r')
    fasta_lines=fasta.readlines()
    fasta_len=len(fasta_lines)

    for i in xrange(fasta_len):
        fasta_lines[i]=re.sub("\n", "", fasta_lines[i])

    fasta_format=open(fasta_f+"2", 'w')
    fasta_format.write("".join(fasta_lines[1:]))
    fasta_format.close()

    fasta_format=open(fasta_f+"2", 'r')
       
    anno_f=open(anno_f, 'r')
    fasta_f=open(fasta_f, "r")
    org_fasta_lines=fasta_f.readlines()[1:]

    org_fasta=""

    for org_fasta_line in org_fasta_lines:
        org_fasta=org_fasta+re.sub("\n", "", org_fasta_line)
        
    seq_f=open(seq_f, "w")

    annolines=anno_f.readlines()
    anno_f.close()

    for annoline in annolines:
        annoline=re.sub("\n", "", annoline)
        annolist=annoline.split("\t")
        annolist=[str(x) for x in annolist]

        coordinate=re.sub("complement\(", "", annolist[3])
        coordinate=re.sub("\)", "", coordinate)
        
        if not coordinate[:4]=='join':
            l_coordinate=coordinate.split("..")[0]
            r_coordinate=coordinate.split("..")[1]

            seq_f.write(">")
            seq_f.write(annolist[0][5:]+'|')
            seq_f.write('|'.join(annolist[1:]))
            seq_f.write('\n')
            seq_len=abs(int(l_coordinate)-int(r_coordinate))+1

            if l_coordinate<r_coordinate:
                fasta_format.seek(int(l_coordinate)-1)     
                dna=Seq(fasta_format.read(seq_len), generic_dna)
                seq_f.write(str(dna))
                
            if l_coordinate>r_coordinate:
                fasta_format.seek(int(l_coordinate)-1)     
                dna=Seq(fasta_format.read(seq_len), generic_dna).reverse_complement()
                seq_f.write(str(dna))
            
            seq_f.write("\n")
    seq_f.close()

    anno_f.close()
    fasta_f.close()

def launchBLAST(blast_dir, blast_db, blast_query, blastRes, db_fname):
   '''
   This function is used to extract the BLAST search 
   annotations for a specific organism
   In order to use this function, you need to provide:
   1. local BLAST software directory 
   2. BLAST database file
   3. BLAST query file
   4. BLAST search results output file
   '''
   
   #launch BLAST search
   os.chdir(blast_dir)
   result_handle=subprocess.Popen(args="makeblastdb -in "+blast_db+" -dbtype nucl -parse_seqids -input_type fasta", shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   result_handle=subprocess.Popen(args="blastn -query "+blast_query+" -db "+blast_db+" -evalue 0.000001"+" -outfmt 5", shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)

   #parse BLAST search result and write them out
   blast_records = NCBIXML.parse(result_handle.stdout)

   f=open(blastRes, 'w')
   header1=["gene_id", "direction", "left_end", "right_end", "length", "annotation"]
   header2=["match_line", "match_gene_id", "contig_start", "contig_end", "node_start", "node_end", "evalue", "alignment_length", "score", "bits", "identities", "positives", "gaps\n"]

   f.write("\t".join(header1+header2))
   
   for blast_record in blast_records:
      ##write out the query data
      f.write("\t".join(blast_record.query.split("|"))+"\t")
      if blast_record.alignments:
         alignment=blast_record.alignments[0]
         ##the line number of DB file
         f.write(str(int(alignment.title.split("|")[0]))+"\t")
         ##write the gene ID
         f.write(db_fname+"%05d" % int(alignment.title.split("|")[0])+"\t")
         hsp=alignment.hsps[0]

         f.write(str(hsp.query_start)+"\t")
         f.write(str(hsp.query_end)+"\t")
         f.write(str(hsp.sbjct_start)+"\t")
         f.write(str(hsp.sbjct_end)+"\t")
         f.write(str(hsp.expect)+"\t")
         l=len(hsp.query)

         f.write(str(l)+"\t")
         f.write(str(hsp.score)+"\t")
         f.write(str(hsp.bits)+"\t")
         f.write(str(hsp.identities)+"\t")
         f.write(str(hsp.positives)+"\t")
         f.write(str(hsp.gaps)+"\n")
      else:
         f.write("\n")
   f.close()
   
def extractAnno(parseDir, clusterID, annoDir, blastRes, extractFile, blastDir):
    '''
    This function is used to extract the one to one orthologous
    annotations for a specific organism
    In order to use this function, you need to provide:
    1. cur_folder is the folder containing the pair wise orthology results
       and the annotation files. 
    2. clusterID is the matched organism ID file
    3. annoDir is the annotation file directory
    4. blastRes is the the matched organism blast results 
    5. extractFile is the extracted cluster annotatioin file 
    '''   

    os.chdir(parseDir)
    
    #open organism ID file
    org_cl_file=open(clusterID, "r")
    cl_lines=org_cl_file.readlines()
    for cl_line in cl_lines:
       cl_list=cl_line.split("\t")
    
    #construct the gene id dictionary from the gene id cluster file
    gene_cl_dict=dict()
    
    l=0
    for line in cl_lines:
       gene_cl_dict[line.split("\t")[0]]=l
       l=l+1
       
    #open extraction result file
    cl_anno_file=open(blastDir+extractFile, "w")
    
    #open blast result file
    
    f=open(blastDir+blastRes, "r")
    blast_lines=f.readlines()
    
    ##extract the cluster annotation data from the annotation file and
    ##id cluster file
    
    k=0
    
    for blast_line in blast_lines[1:]:
       k=k+1
       print 'extract anntoation cluster '+str(k)
       blast_list=blast_line.split("\t")
       if len(blast_list)>8 and not blast_list[7]=="":
          if gene_cl_dict.has_key(blast_list[7]):
             cl_list=cl_lines[gene_cl_dict[blast_list[7]]].split("\t")[:-1][:30]
             cl_list_len=len(cl_list)
             
             cl_anno_file.write(cl_list[0])
             cl_anno_file.write("\t")
    
             for i in xrange(0, cl_list_len-1):
                 if cl_list[i]=="":
                     pass
                 else:
                     file_name=cl_list[i][0:5]
                     gene_id=cl_list[i][5:]
                     
                     anno_file=open(annoDir+file_name+".txt", "r")
                     anno_lines=anno_file.readlines()
                   
                     if int(anno_lines[int(gene_id)-1][5:10])==int(gene_id):
                         if not anno_lines[int(gene_id)-1].split("\t")[4]=="":
    ##                        print anno_lines[int(gene_id)-1].split("\t")[4]
                            cl_anno_file.write(anno_lines[int(gene_id)-1].split("\t")[4])
                            cl_anno_file.write("\t")
                     anno_file.close()
    
       cl_anno_file.write("\n")
          
    cl_anno_file.close()
    org_cl_file.close()

##def format_blast_db(blast_db):
##   db_file=open(blast_db, 'r')
##   db_lines=db_file.readlines()
##   db_format_fname=blast_db+"format"
##   db_format_f=open(db_format_fname, 'w')
##
##   for db_line in db_lines:
##      if db_line[0]=='>':
##         db_format_f.write('\n'+db_line)
##      else:
##         db_format_f.write(re.sub('\n', '', db_line))
##
##   db_format_f.close()
##   db_file.close()
##   return(db_format_fname)

