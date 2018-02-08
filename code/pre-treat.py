import sys
import os
import re
import random
import gc
import operator
import numpy
import math
import nltk
import pickle
import m
import pre
import array

from datetime import datetime
from numpy import zeros,dot
from datetime import datetime

if __name__ == '__main__':
    '''
    This module aims to use BLAST to establish the orthology and then use these
    othology to extract the orthologous gene anntoations from OMA database. The
    query sequence refers to the genome that you would like to compare the
    annotation, the database sequence refers to the genome that you would like
    to use to build up the OGA from OMA database.

    In order to use this module, you should do the following things:
    1. Install local BLAST from NCBI database.
    2. Point the variable 'blastDir' to the BLAST folder.
    3. Download the query and database genome .fna file
    4. Select database genome annotation file from the OMA database.
    5. The gene recognition integration results of the multiple annotations.
    '''

    tstart = datetime.now()
    blastDir="C:/Program Files/NCBI/blast-2.2.27+/"

    org_abb_list=['ecol', 'ctra', 'hpyl', 'mgen', 'mytu', 'rpro']
    org_cap_list=['ECOLI', 'CHLTR', 'HELPY', 'MYCGE', 'MYCTU', 'RICPR']

    for i in xrange(len(org_abb_list)):
        org_abb=org_abb_list[i]
        org_cap=org_cap_list[i]
        
        ##parse the fna file to BLAST db file
        os.chdir(blastDir)
        fasta_f='.\\genome data\\'+org_cap+'.fna'
        anno_f='.\\genome data\\'+org_cap+'.txt'
        db_f='.\\db\\'+org_cap+'_db.fa'
        pre.parseDbFile(fasta_f, anno_f, db_f)
        
        print 'parse BLAST db file done'
        print datetime.now()-tstart

        ##parse the fna file to BLAST query file format
        consensus_anno=".\\consensus annotation\\"+org_abb+"_output.txt"
        query_f=".\\query\\"+org_cap+"_query.txt"
        pre.parseQueryFile(fasta_f, consensus_anno, query_f)

        print 'parse BLAST query file done'
        print datetime.now()-tstart

        ##launch BLAST search
        blastRes=".\\results\\blast_res_"+org_abb+".txt"
        db_fname=org_cap
        ##   blast_db_f=format_blast_db('db/'+blast_db)
        
        pre.launchBLAST(blastDir,
          db_f, 
          query_f,
          blastRes,
          db_fname)

        print 'BLAST search done'
        print datetime.now()-tstart

        ##extract gene cluster annotations 
        parseDir="C:/Users/Liu Zhe/Desktop/New folder/arc1/New folder/one to one ortholog/"
        clusterID='.\\pairs files\\prok only cluster\\'+org_cap+'.txt'
        annoDir=".\\annotation files\\"
        extractFile='.\\results\\'+org_cap+' to '+org_cap+' cluster annotation.txt'
        
        pre.extractAnno(parseDir, clusterID, annoDir, blastRes, extractFile, blastDir)
        
        print 'extract gene cluster annotations done'
        print datetime.now()-tstart

