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

def extract_anno(anno_dir, anno_names):
    '''
    This function aims to extract the following files from the annotation files:
    1. word_dic_list: entire annotation word dictionary 
    2. anno_list: annotation list 
    3. anno_proc_list: processed annotation list
    '''
    word_dic_list=list()
    anno_list=list()
    anno_proc_list=list()

    for i in xrange(len(anno_names)):
        anno_file = open(anno_dir+anno_names[i], "r")
        anno_lines=anno_file.readlines()
        anno_file.close()
        
        anno_list.append([])
        
        anno_proc_list.append([])
        word_dic_list.append(dict())
        for anno_line in anno_lines:
            element=anno_line.split('\t')
            element_proc=proc(preproc(element[0]))
            
            anno_match=[re.split(r'\;|\s|\||\,|\-|\/|\(|\)', preproc(x)) for x in element]

            anno_list_l=re.sub('\n', '', anno_line).split('\t')
            parse_tag=True
            sub_info=sub_EC(anno_list_l[0])
            if len(anno_list_l)>2:
                if anno_list_l[1]=='':
                    if not sub_info[1]=='':
                        anno_list_l[0]=sub_info[0]
                        anno_list_l[1]=sub_info[1]
                else:
                    if not anno_list_l[1]==sub_info[1] and not sub_info[1]=='':
                        anno_list_l[0]=sub_info[0]
                        anno_list_l[1]=anno_list_l[1]+','+sub_info[1]
                    else:
                        anno_list_l[1]=anno_list_l[1]

                anno_list[i].append('\t'.join(anno_list_l)) 
                anno_proc_list[i].append(element_proc)
                
                for x in anno_match:
                    for y in x:
                        if not word_dic_list[i].has_key(y) and not y=="":
                            word_dic_list[i][y]=1
    
    return(word_dic_list, anno_list, anno_proc_list)

def comp_baseline(anno_list_1, anno_list_2):
    '''
    This function performs the baseline comparison between two annotation results.
    '''
    res_list=list()
    Ng=len(anno_list_1)

    for i in xrange(0, Ng): 
        res_list.append("0")

    for l in xrange(0, Ng):
        if res_list[l]=='0' and not anno_list_1[l]=="" and not anno_list_2[l]=="":
            anno_list_l_1=re.sub("\n", "", anno_list_1[l]).split("\t")
            anno_list_l_2=re.sub("\n", "", anno_list_2[l]).split("\t")
            
            anno_l_1=anno_list_l_1[0]
            anno_l_2=anno_list_l_2[0]
            
            ec_list_l_1=anno_list_l_1[1]
            ec_list_l_2=anno_list_l_2[1]

            gs_list_l_1=anno_list_l_1[2]
            gs_list_l_2=anno_list_l_2[2]

            ##One of the annotations is missing  
            if anno_l_1=='' and anno_l_2=='':
                res_list[l]='base_both miss'
            
            if anno_l_1=='' and not anno_l_2=='':
                res_list[l]='base_one miss'
            
            if not anno_l_1=='' and anno_l_2=='':
                res_list[l]='base_one miss'
            
            #construct the EC list
            ec_set1=set(ec_list_l_1.split(';'))
            ec_set2=set(ec_list_l_2.split(';'))
            
            ####EC are the same
            if re.search('-', ec_list_l_1)==None and re.search('-', ec_list_l_2)==None:
                if not ec_set1==set(['']) and not ec_set2==set(['']):
                    if len(ec_set1.symmetric_difference(ec_set2))==0:
                        res_list[l]='base_ec'

            ####EC are different
            if re.search('-', ec_list_l_1)==None and re.search('-', ec_list_l_2)==None:
                if not ec_set1==set(['']) and not ec_set2==set(['']):
                    if len(set(ec_set1)&set(ec_set2))==0:
                        res_list[l]='base_ec-d'

            #compare the annotation via gene symbol
            if anno_list_l_1[3]==anno_list_l_2[3] and not anno_list_l_1[3]=="":
                res_list[l]="base_tigr"
                
            #compare the annotation via gene symbol
            if gs_list_l_1.lower()==gs_list_l_2.lower() and not gs_list_l_1=="":
                res_list[l]="base_gs"

            #compare annotations with HP and non-HP
            if anno_l_1=="hypothetical protein" or anno_l_1=="conserved hypothetical protein":
                if not anno_l_2=="hypothetical protein" and not anno_l_2=="conserved hypothetical protein":
                    res_list[l]="base_one HP"
                        
            if anno_l_2=="hypothetical protein" or anno_l_2=="conserved hypothetical protein":
                if not anno_l_1=="hypothetical protein" and not anno_l_1=="conserved hypothetical protein":
                    res_list[l]="base_one HP"
                        
            #both of the annotations are the same
            if anno_l_1==anno_l_2 and not anno_l_2=="":
                res_list[l]="base_lit-s"

    return res_list

def comp_id_anno(res_list, anno_proc_1, anno_proc_2):
    '''
    This function compares database ID and gene annotation
    '''
    Ng=len(res_list)

    for l in xrange(0, Ng):
        if res_list[l]=='0':
            ##compare hypothetical protein
            if not anno_proc_1[l]==[""] and not anno_proc_2[l]==[""]:
                cmp_res=cmpAnno(anno_proc_1[l], anno_proc_2[l])
                hp_res=cmpHPAnno(anno_proc_1[l], anno_proc_2[l])
                
            ##compare hypothetical protein
            if not hp_res=='0':
                res_list[l]=hp_res
                
            ##compare the same annotations
            if not cmp_res=='0':
                res_list[l]='lit-s'

    return res_list

def cmpAnno(anno_1, anno_2):
    #compare the annotation in the bracket separately
    ignore_list=['like', 'chain', 'possible', '+', 'subfamily', 'system', 'component', 'membrane', 'uncharacterised', 'subunit', 'homolog']
    
    res='0'

    #token difference between two annotations
    set_dif=set(anno_2).symmetric_difference(set(anno_1))
    sym_dif_len=len(set_dif)

    #automated annotation contains more information
    if set_dif==[]:
        res='1'
    else:            
        if len(set_dif.difference(ignore_list))==0:
            res='1'

    return res

def cmpHPAnno(anno_1, anno_2):
    res='0'
    
    #deal with different annotations
    if anno_1=='hypothetical':
        if not anno_2=='hypothetical':
            res='one HP'

    if not anno_1=='hypothetical':
        if anno_2=='hypothetical':
            res='one HP'

    return res

def preproc(anno_str):
    '''
    This function processed gene annotation by
    1. removing some database ID
    2. uninformative annotations
    '''
    anno_ini=anno_str
    ##eliminate the useless information in annotations
    anno_str=re.sub("^\"", "", anno_str)
    anno_str=re.sub("\"$", "", anno_str)
    anno_str=re.sub("TIGR[0-9]{5}", "", anno_str)
    anno_str=re.sub("DUF[0-9]{0,4}", "", anno_str)
    anno_str=re.sub("FOG\:", "", anno_str)
    anno_str=re.sub("FIG[0-9]{1,10}\:", "", anno_str)
    anno_str=re.sub("FIGfam[0-9]{1,10}", "", anno_str)
    anno_str=re.sub("UPF[0-9]{4}", "", anno_str)
    anno_str=re.sub("COG[0-9]{4}\:", "", anno_str)
    anno_str=re.sub("\(TC\s.*\\)", "", anno_str)
    anno_str=re.sub("\(EC\s[0-9.-]{0,15}\)", "", anno_str)

    ##replace the annotation as hypothetical protein
    anno_str=re.sub("clustering\swith.*|clustered\swith.*", "hypothetical protein", anno_str)
    anno_str=re.sub("coding\sregion.*", "hypothetical protein", anno_str)
    anno_str=re.sub("conserved\sdomain\sprotein|[Uu]ncharacterized\sconserved\sprotein|[Uu]ncharacterized\sprotein\sconserved\sin\sbacteria", "hypothetical protein", anno_str)
    anno_str=re.sub("conserved\shypothetical\sprotein", "hypothetical protein", anno_str)
    anno_str=re.sub("synthetase", "synthase", anno_str)

    if not anno_ini=="" and anno_str=="":
        anno_str="hypothetical protein"

    return anno_str

def proc(anno_str):
    '''
    This function further processed gene annotation by
    1. lowcase the annotation text
    2. lemmatize the annotation
    3. eliminate the stop word
    '''
    
    anno_ini=anno_str
    
    token_list=list()

    #seperate tokens
    token_list=nltk.word_tokenize(anno_str)

    #change to lower case
    token_lowcase=list()
    token_lowcase=[x.lower() for x in token_list]

    #Lemmatize the annotation
    penn_tag=[x[1] for x in nltk.pos_tag(token_lowcase)]
    morphy_tag = {'NN':wordnet.NOUN,'JJ':wordnet.ADJ,'VB':wordnet.VERB,'RB':wordnet.ADV}

    ##, 'IN':wordnet.NOUN, 'CC':wordnet.NOUN, 'CD':wordnet.NOUN, 'DT':wordnet.NOUN, 'EX':wordnet.NOUN, 'FW':wordnet.NOUN, 'IN':wordnet.NOUN, 'LS':wordnet.NOUN, 'MD':wordnet.NOUN, 'PD':wordnet.NOUN, 'PO':wordnet.NOUN, 'PR':wordnet.NOUN, 'RP':wordnet.NOUN, 'SY':wordnet.NOUN, 'TO':wordnet.NOUN, 'UH':wordnet.NOUN, 'WD':wordnet.NOUN, 'WP':wordnet.NOUN, 'WP':wordnet.NOUN, 'WR':wordnet.NOUN, '-N':wordnet.NOUN, '``':wordnet.NOUN
    token_lemma=list()
    #####automatic verb recognition
    wnl=nltk.WordNetLemmatizer()
    morphy_tag_res=list()
    
    for i in xrange(0, len(token_lowcase)):
        if morphy_tag.has_key(penn_tag[i][:2]):
            morphy_tag_res.append(morphy_tag[penn_tag[i][:2]])
        else:
            morphy_tag_res.append(wordnet.NOUN)
            
    token_lemma=[wnl.lemmatize(token_lowcase[i], morphy_tag_res[i]) for i in xrange(0, len(token_lowcase))]

    #eliminate stop word
    stops=['n-terminal', 'uncharacterised', 'subfamily', 'relate', '-', 'cluster', 'uncharacterized', 'biosynthetic', 'biosynthesis', 'molecular', 'pathway', 'system', 'c-terminal', 'activity', 'binding', 'containing', 'possible', 'have', 'probable', 'for', 'region', 'unknown', 'to', 'bind', 'bacterial', 'ec', 'in', 'unkonwn', 'superfamily', 'subunit', 'component', 'and', 'enzyme', 'associate', 'type', 'involve', 'function', 'consist', 'domain-containing', 'putative', 'uncharacterise', 'probable', 'predict', 'family', 'translate', 'translation', 'transcription', 'signal', 'protein', '[', ']', '(', ')', '--', ';', ':', ',', '/', 'or', 'with', 'the', 'contain', 'of', 'product', 'related', 'domain', 'an', '|', '.','conserve', '#', 'response', 'chain']
    ##regulation, uncharacterised, regulator, transport, regulatory, cell, catalytic, group
    s = set(stops)
    token_stop = [x for x in token_lemma if x not in s]
    
    if not anno_ini=='' and len(token_stop)==0:
        token_stop=['hypothetical']

    return token_stop

def comp_pfam_dic(f_pfam_dir, res_list, anno_list_1, anno_list_2, anno_proc_1, anno_proc_2, word_dic_1, word_dic_2):
    '''
    ##This function aims to compare two annotations with pfam dictionary
    ##When
    ##1. Both of the annotations are the same as the pfam annotaitons
    ##2. Both of the anntoations have all of the pfam keywords in common
    ##They are assigned as the same annotations.
    '''
    
    Ng=len(res_list)
    word_dic_list=[word_dic_1, word_dic_2]

    pfam_file=open(f_pfam_dir, 'r')
    pfam_dic=dict()

    ##find the pfam keywords in the annotation and construct a pfam dictionary
    while True:
        pfam_line=pfam_file.readline()
        
        if pfam_line=='':
            break
        
        pfam_list=re.sub('\n', '', pfam_line).split('\t')

        for ind in xrange(2):
            if word_dic_list[ind].has_key(pfam_list[0]):
                pfam_dic[pfam_list[0]]=pfam_list[1]

    pfam_file.close()

    for l in xrange(0, Ng):
        if res_list[l]=='0' and not anno_list_1[l]=="" and not anno_list_2[l]=="":
            pfam_s=set([])
            pfam_lab=list()
            
            anno_match_l_1=re.split(r'\;|\s|\||\,|\-|\/|\(|\)', preproc(anno_list_1[l].split("\t")[0]))
            anno_match_l_2=re.split(r'\;|\s|\||\,|\-|\/|\(|\)', preproc(anno_list_2[l].split("\t")[0]))

            anno_1ist=[anno_list_1, anno_list_2]
            anno_match_list=[anno_match_l_1, anno_match_l_2]

            anno_proc=[anno_proc_1[l], anno_proc_2[l]]
            
            ##extract all of the pfam data
            pfam_list=pfam_dic.keys()

            for anno in anno_match_list:
                for x in anno:
                    if len(x)>=3 and x[0].isupper():
                        if pfam_dic.has_key(x):
                            pfam_s=pfam_s|set([x])
            
            if len(pfam_s)==1:
                #check if pfam annotation is the same as another anntoation
                for ind in xrange(2):
                    for y in pfam_dic[list(pfam_s)[0]].split('|'):
                        pfam_res=cmpAnno(anno_proc[ind], proc(preproc(y)))
                        if not pfam_res=='0' and not ind in pfam_lab:
                            pfam_lab.append(ind)
                        
            #check if all of the pfam IDs appear in both of the annotations
            if len(pfam_s)>=1:
                if set(pfam_s).issubset(set(anno_match_list[0])):
                    if set(pfam_s).issubset(set(anno_match_list[1])):
                        res_list[l]='pfam-sh'

            if len(set(pfam_lab))>1:
                res_list[l]='pfam-sh'

    return res_list

def comp_gs_dic(f_gs_dir, res_list, anno_list_1, anno_list_2, anno_proc_1, anno_proc_2, word_dic, anno_list):
    '''
    ##This function aims to compare two annotations with gene symbol dictionary
    ##When
    ##1. Both of the annotations are the same as the gene symbol annotaitons
    ##2. Both of the anntoations have the gene symbol in common
    ##They are assigned as the same annotations.
    '''

    gslist_dic=dict()
    gs_dic=dict()

    gs_anno=dict()
    gs_list=list()

    anno_match_list=list()
    Ng=len(res_list)
    
    ##extract all of the gene symbols
    for i in xrange(len(anno_list)):
        for l in xrange(0, Ng):
            gs_list=anno_list[i][l].split("\t")[2]
            
            if not gslist_dic.has_key(gs_list):
                gslist_dic[gs_list]=1
                
        #if a word is four character length and ending by a capital letter,
        #we consider it as gene symbol as well.
        for k in word_dic[i].keys():
            if len(k)==4 and k[-1].isupper():        
                if not gslist_dic.has_key(k):
                    gslist_dic[k]=1
                if not gslist_dic.has_key(k):
                    gslist_dic[k]=1

    ##construct gene symbols dictionary
    gs_file=open(f_gs_dir, 'r')
    
    while True:
        gs_line=gs_file.readline()
        
        if gs_line=='':
            break
        
        gs=gs_line.split('\t')[0]
        gs_func=re.sub(r'\n', '', gs_line.split('\t')[3])

        if gslist_dic.has_key(gs):
            gs_dic[gs]=gs_func

    ##compare annotation via gs
    for l in xrange(0, Ng):
        if res_list[l]=='0' and not anno_list_1[l]=="" and not anno_list_2[l]=="":
            gs_s=set([])
            gs_lab=list()

            gs_list_l_1=anno_list_1[l].split("\t")[2]
            gs_list_l_2=anno_list_2[l].split("\t")[2]
            gs_list=[gs_list_l_1, gs_list_l_2]

            anno_match_1=re.split(r'\;|\s|\||\,|\-|\/|\(|\)', preproc(anno_list_1[l].split("\t")[0]))
            anno_match_2=re.split(r'\;|\s|\||\,|\-|\/|\(|\)', preproc(anno_list_2[l].split("\t")[0]))
            anno_match=[anno_match_1, anno_match_2]
            anno_proc=[anno_proc_1[l], anno_proc_2[l]]
            
            ##match any gene symbol in the dictionary
            for ind in xrange(2):
                if not gs_list[ind]=='' and gs_dic.has_key(gs_list[ind]):
                    gs_s=gs_s|set([gs_list[ind]])
                    gs_lab.append(ind)
                ##match any word from the annotation in the dictionary            
                for x in anno_match[ind]:
                    if not x =='' and gs_dic.has_key(x):
                        gs_tmp1=x
                        gs_tmp2=gs_dic[x]
                        gs_s=gs_s|set([x])
                        gs_lab.append(ind)
                                
            ##check if gene annotation is the same as another annotation
            if len(gs_s)>0:
                for gs_ind in gs_s:
                    gs_dic_anno=gs_dic[gs_ind].split('|')
                    for gs_dic_anno_ind in gs_dic_anno:
                        gs_dic_anno_proc=nltk.word_tokenize(gs_dic_anno_ind.lower())
                        gs_res1=cmpAnno(gs_dic_anno_proc, anno_proc[0])
                        gs_res2=cmpAnno(gs_dic_anno_proc, anno_proc[1])
                        
                        if not gs_res1=='0':
                            gs_lab.append(0)
                        if not gs_res2=='0':
                            gs_lab.append(1)
                            
                    ##check if gene symbol appears in another annotation
                    for ind in xrange(2):
                        if gs_ind in anno_match[ind] and not ind in gs_lab:
                            gs_lab.append(ind)
                            
                if len(set(gs_lab))>1:
                    res_list[l]='gs-sh'
    return(res_list)

def comp_orth(f_orth_dir, res_list, anno_list_1, anno_list_2, anno_proc_1, anno_proc_2):
    '''
    compare using orthologous gene dataset
    '''
    f_orth=open(f_orth_dir, 'r')
    Ng=len(res_list)

    for l in xrange(0, Ng):
        if not res_list[l]=="0":
            orth_line=f_orth.readline().lower()
        else:
            orth_lab=list()
            orth_line=f_orth.readline().lower()
            orth_list=orth_line.split('\t')[1:]
            
            anno_preproc_1=preproc(anno_list_1[l].lower().split("\t")[0])
            anno_preproc_2=preproc(anno_list_2[l].lower().split("\t")[0])
            
            anno_token_1=nltk.word_tokenize(anno_preproc_1)
            anno_token_2=nltk.word_tokenize(anno_preproc_2)
            
            anno_preproc=[anno_preproc_1, anno_preproc_2]
            anno_token=[anno_token_1, anno_token_2]
            
            anno_proc=[anno_proc_1[l], anno_proc_2[l]]

            for ind in xrange(2):
                ##compare annotation by exact match and ensure that
                ##the annotation is not hypothetical protein
                if anno_preproc[ind] in orth_list and not anno_proc[ind]==['hypothetical']:
                    orth_lab.append(ind)
                    
                #compare annotations using preprocessing and without doing fully preprocessing
                for x in orth_list:
                    x_proc=nltk.word_tokenize(preproc(x))
                    orth_res=cmpAnno(x_proc, anno_token[ind])
                    if not orth_res=='0' and not anno_proc[ind]==['hypothetical']:
                        orth_lab.append(ind)

            if len(set(orth_lab))>1:
                res_list[l]='orth'
    f_orth.close()
    return res_list

def constQueDict(anno_proc_list, doc_freq):
    '''
    construction of query word dictionary
    '''
    que_idx=dict()
    que_len=list()

    que_freq=doc_freq

    for anno_proc in anno_proc_list:
        for anno in anno_proc:
            que_len.append(len(anno))
            for word in anno:
                que_freq.setdefault(word,1)
                que_freq[word]+=1

    keys=que_freq.keys()
    keys.sort()

    for i in xrange(len(que_freq)):
        que_idx[keys[i]]=(i, que_freq[keys[i]])
        
    num_que=0
    
    for x in que_len:
        if not x==0:
            num_que+=1
    
    return(que_idx, num_que, que_len)

def constDocDict(anno_proc_list, blast_dict):
    '''
    construction of document word dictionary
    '''
    word_freq=dict()

    for anno_proc in anno_proc_list:
        for anno in anno_proc:
            for word in anno:
                word_freq.setdefault(word,0)
        
    for blast in blast_dict:
        for word in blast:
            word_freq.setdefault(word,0)

    doc_len=[len(blast) for blast in blast_dict]
    
    num_doc=0

    for x in doc_len:
        if not x==0:
            num_doc+=1

    doc_freq=dict()
    doc_freq=word_freq

    for blast in blast_dict:
        for word in blast:
            doc_freq.setdefault(word,1)
            doc_freq[word]+=1

    keys=doc_freq.keys()
    keys.sort()
    doc_idx=dict()

    for i in xrange(len(doc_freq)):
        doc_idx[keys[i]]=(i, doc_freq[keys[i]])

    return(doc_idx, num_doc, doc_len, word_freq)

def add_word(word,d):
    """
    Add words into a dictionary
    """
    d.setdefault(word,0)
    d[word]+=1
    
def doc_vec(doc,key_idx):
    """
    Construct a vector containing words in ranked order
    """
    vec=np.zeros(len(key_idx))
    
    for word in doc:
        keydata=key_idx.get(word, None)
        if keydata: vec[keydata[0]]+=1
        
    return vec

def comp_que_w(query, que_idx, num_que, que_len, que_len_list):
    """
    function used to compute query weight
    """
    ele_num=len(que_idx)
    que_con_t=zeros(ele_num)
    query_weight=zeros(ele_num)
    fir_part=zeros(ele_num)
    sec_part=zeros(ele_num)
    
    for word in query:
        keydata=que_idx.get(word, None)
        if keydata: que_con_t[keydata[0]]=1

    que_con_t=que_con_t+0.00001
    que_tf=doc_vec(query, que_idx)

    fir_part=que_tf/que_len
    sec_part=np.log10(num_que/que_con_t)
    query_weight=fir_part*sec_part
    
    sum_weight=np.sqrt(np.sum(np.square(query_weight)))
    query_weight=query_weight/sum_weight

    return query_weight

def comp_doc_w(doc, doc_idx, num_doc, doc_len, doc_len_list):
    """
    function used to compute document weight
    """
    ele_num=len(doc_idx)

    doc_weight=zeros(ele_num)
    doc_con_t=zeros(ele_num)
    fir_part=zeros(ele_num)
    sec_part=zeros(ele_num)

    avg_doc_len=np.sum(doc_len)/num_doc

    for w in doc:
        keydata=doc_idx.get(w, None)
        if keydata: doc_con_t[keydata[0]]=1
        
    doc_con_t=doc_con_t+0.00001
    doc_tf=doc_vec(doc,doc_idx)
    
    fir_part=doc_tf/(doc_tf+1.5*(doc_len/avg_doc_len))
    sec_part=np.log10(num_doc/doc_con_t)
    
    doc_weight=fir_part*sec_part
    sum_weight=np.sqrt(np.sum(np.square(doc_weight)))
    doc_weight=doc_weight/sum_weight

    return doc_weight

def tfidf(res_list, DocRes, QueRes, anno_proc_list, bla_list):
    """
    function used to compute the tf-idf score between query and document
    """
    tfidf_list=list()
    doc_len=DocRes[2]
    Ng=len(res_list[0])
    c_len=len(res_list)

    doc_idx=DocRes[0]
    num_doc=DocRes[1]
    doc_len=DocRes[2]
    word_freq=DocRes[3]

    que_idx=QueRes[0]
    num_que=QueRes[1]
    
    for l in xrange(0, Ng):
        res_lab=0
        for i in xrange(c_len):
            if res_list[i][l]=="0":
                res_lab=1
                break
            
        if res_lab==1:
            bla_list_l=bla_list[l]
            anno_proc=anno_proc_list[l]
            
            que_len=QueRes[2][l]
            anno_proc_len=len(anno_proc)
            
            if not doc_len[l]==0:
                doc_weight=comp_doc_w(bla_list_l, doc_idx, num_doc, len(bla_list_l), doc_len)
                if not anno_proc==['hypothet'] and not anno_proc==[]:
                    query_weight=comp_que_w(anno_proc, que_idx, num_que, anno_proc_len, que_len)
                    sim_score=np.sum(doc_weight*query_weight)
                    tfidf_list.append(sim_score)
                else:
                    tfidf_list.append(0)
            else:
                tfidf_list.append(0)
        else:
            tfidf_list.append(0)

    return(tfidf_list)

def cmp_tfidf(res_list, tfidf_list1, tfidf_list2, threshold):
    '''
    This function compare the two annotations based on tf-idf score
'''
    Ng=len(res_list)

    for l in xrange(0,Ng):
        if res_list[l]=='0':
            if tfidf_list1[l]>threshold and tfidf_list2[l]>threshold:
                res_list[l]='tf_idf'
                    
    return res_list

def map_res(res_list, anno_proc_list, dif_tag):
    '''
    This function map the same annotations 
    '''
    Ng=len(res_list[0])
    res_len=len(anno_proc_list)
    
    ind_list=list()
    
    for i in xrange(0, res_len):
        for j in xrange(i+1, res_len):
            ind_list.append([i, j])

    ind_len=len(ind_list)

    for l in xrange(0, Ng):
        tmp_list=list()
        s_res=list()

        for ind in xrange(0,ind_len):
            ind1=ind_list[ind][0]
            ind2=ind_list[ind][1]
            
            if tmp_list==[]:
                if not res_list[ind][l] in dif_tag:
                    tmp_list.append(ind1)
                    tmp_list.append(ind2)
            else:
                if not res_list[ind][l] in dif_tag:
                    if ind1 in tmp_list or ind2 in tmp_list:
                        tmp_list.append(ind1)
                        tmp_list.append(ind2)

        tmp_set=list(set(tmp_list))
        tmp_len=len(tmp_set)
        
        for i in xrange(tmp_len):
            for j in xrange(i+1, tmp_len):
                if not anno_proc_list[i][l]==['hypothetical'] and not anno_proc_list[j][l]==['hypothetical']:
                    map_ind=ind_list.index([tmp_set[i], tmp_set[j]])
                    if res_list[map_ind][l]=="0":
                        res_list[map_ind][l]="map"

    return (res_list)

##def writeResult(res_list, resfile):
##    resfile=open("results/baseline.txt", 'w')
##    resfile.write("----------------------------------------------------------------")
##
##    for l in xrange(0, Ng):
##        for res_ele in res_list[l]:
##            resfile.write(str(res_ele)+"\t")
##        resfile.write("\n")
##
##    resfile.write("----------------------------------------------------------------")
##    resfile.close()

def sub_EC(anno_tmp):
    '''
    This function extract the EC ID
    '''
    annolist=list()
    annolist=['', '']

    #extract and eliminate EC
    
    ec_findall=re.findall(r'\s\(EC\s([\.0-9\s\-]{7,31})\)', anno_tmp)
    
    if not ec_findall==None:
        annolist[1]=', '.join(ec_findall)
        annolist[0]=re.sub('\s\(EC\s([\.0-9\s\-]{7,31})\)', '', anno_tmp)
    else:
        annolist[0]=anno_tmp
        annolist[1]=''

    return annolist

def printResult(res_list, res_sta_file):
    '''
    This function print the result
    '''
    times=dict()
    
    time_seq=['base_both miss', 'base_one miss', 'base_ec', 'base_ec-d', 'base_tigr', 'base_gs', 'base_one HP', 'base_lit-s', 'lit-s', 'pfam-sh', 'gs-sh', 'orth', 'tf_idf', 'map', '0']
    
    res_sta_file=open(res_sta_file, 'a')
    
    Ng=len(res_list)
    
    for l in xrange(0, Ng):
        if not times.has_key(res_list[l]):
            times[res_list[l]]=1
            continue
        if times.has_key(res_list[l]):
            times[res_list[l]]=times[res_list[l]]+1

    print "----------------------------------"
    for x in time_seq:
        if times.has_key(x):
            print str(x)+"\t"+str(times[x])+"\t"+"{0:.2f}%".format(times[x]/Ng*100)
            res_sta_file.write(str(x)+"\t"+str(times[x])+"\t"+"{0:.2f}%".format(times[x]/Ng*100)+"\n")
        else:
            print x+"\t"+"0"+"\t"+'0'
            res_sta_file.write(str(x)+"\t"+'0'+"\t"+'0'+"\n")
    
    print 'Entire gene number\t'+str(Ng)
    res_sta_file.write('Entire gene number\t'+str(Ng)+"\n")
    res_sta_file.write("----------------------------------\n")
    print "----------------------------------"

    res_sta_file.close()
    
def printResult2(res_list, plot_file):
    '''
    This function print the result step by step
    '''
    times2=dict()
    
    time_seq=['base_both miss', 'base_one miss', 'base_ec', 'base_ec-d', 'base_tigr', 'base_gs', 'base_one HP', 'base_lit-s', 'lit-s', 'pfam-sh', 'gs-sh', 'orth', 'tf_idf', 'map', '0']
    baseline=['base_both miss', 'base_one miss', 'base_ec', 'base_ec-d', 'base_tigr', 'base_gs', 'base_one HP', 'base_lit-s']
    improvement=['lit-s', 'pfam-sh', 'gs-sh', 'orth', 'tf_idf', 'map']
    
    Ng=len(res_list)
    
    for l in xrange(0, Ng):
        if res_list[l] in baseline:
            if not times2.has_key('baseline'):
                times2['baseline']=1
                continue

            if times2.has_key('baseline'):
                times2['baseline']=times2['baseline']+1
                
        if res_list[l] in improvement:
            if not times2.has_key('improvement'):
                times2['improvement']=1
                continue

            if times2.has_key('improvement'):
                times2['improvement']=times2['improvement']+1

    print 'baseline\t'+str(times2['baseline'])+'\t'+"{0:.2f}%".format(times2['baseline']/Ng*100)    
    print 'improvement\t'+str(times2['improvement'])+'\t'+"{0:.2f}%".format(times2['improvement']/Ng*100)

    plot_file=open(plot_file, 'a')

    for x in ['baseline', 'improvement']:
        if times2.has_key(x):
            plot_file.write(str(times2[x])+"\t"+"{0:.2f}%".format(times2[x]/Ng)+"\t")
        else:
            plot_file.write('0'+"\t"+'0')

    plot_file.write("\n")
    plot_file.close()

def annoDeter(anno_list, anno_proc_list, res_list, annoDeter_f_dir, annoDeter_sta_f):
    '''
    The function is used to determine the annotations based on
    annotation comparison results. The inputs are:
    1. anno_list: annotation list
    1. anno_proc_list: processed annotation list
    1. res_list: annotation comparison result list
    1. annoDeter_f: the output file for annotation determination results
    '''
    ##extract the annotation information in the anno_proc_list and res_list
    infoRes=extractInfo(anno_proc_list, res_list)
        
    nonhp_list=infoRes[0]
    lab_res=infoRes[1]
    
    res_len=len(res_list)
    cmp_len=int(res_len*(res_len-1)/2)

    annoDet_list=list()
    Ng=len(res_list[0])

    for l in xrange(0,Ng):
        annoDet_list.append(["","",""])
        
    for l in xrange(0, Ng):
        lab_res_l=lab_res[l]
        nonhp_len=nonhp_list[l]
            
        for i in xrange(0,res_len):
            anno_proc_list[i][l]=list(set(anno_proc_list[i][l]))
            
        if lab_res_l==[]:
            ##all of the annotations are hypothetical protein
            if nonhp_list[l]==0:
                if annoDet_list[l][0]=="":
                    annoDet_list[l][0]="all-hp"
                    annoDet_list[l][1]="hypothetical protein"
                    annoDet_list[l][2]="HP\t"
                    
            ##only one of the annotations has informative annotation
            if nonhp_list[l]==1:
                if annoDet_list[l][0]=="":
                    annoDet_list[l][0]="majority"
                    for i in xrange(0,res_len):
                        if not anno_proc_list[1][l]==[] and not anno_proc_list[1][l]==["hypothetical"]:
                            annoDet_list[l][1]='\t'.join(anno_list[1][l].split('\t')[:2])
                        else:
                            annoDet_list[l][1]='\t'.join(anno_list[i][l].split('\t')[:2])
        else:
            ##all of the annotations are the same
            if nonhp_list[l]<=len(lab_res_l):
                if annoDet_list[l][0]=="":
                    annoDet_list[l][0]="all-same"
                    for i in xrange(0,res_len):
                        if not anno_proc_list[1][l]==[] and not anno_proc_list[1][l]==["hypothetical"]:
                            annoDet_list[l][1]='\t'.join(anno_list[1][l].split('\t')[:2])
                        else:
                            annoDet_list[l][1]='\t'.join(anno_list[i][l].split('\t')[:2])

            ##majority of the annotations are the same
            if 0.5*nonhp_len<len(lab_res_l):
                if annoDet_list[l][0]=="":
                    annoDet_list[l][0]="majority"
                    for i in xrange(0,res_len):
                        if not anno_proc_list[1][l]==[] and not anno_proc_list[1][l]==["hypothetical"]:
                            annoDet_list[l][1]='\t'.join(anno_list[1][l].split('\t')[:2])
                        else:
                            annoDet_list[l][1]='\t'.join(anno_list[i][l].split('\t')[:2])
                                        
        if annoDet_list[l][0]=="":
            annoDet_list[l][0]="man-rest"
            annoDet_list[l][1]="MANUAL-REST\t"

    print "already determined"

    ##print the results
    annos=dict()
    times=dict()

    anno_seq=["all-hp", "all-same", "majority", "man-rest"]
    anno_source=["HP", "blast null", "blast", "need check"]

    annoDeter_f=open(annoDeter_f_dir, 'w')
    annoDeter_sta_f=open(annoDeter_sta_f, 'w')
    
    for l in xrange(0, Ng):
        tmp=annoDet_list[l][0]
        if not tmp in times:
            times[tmp]=1
            continue
        if tmp in times:
            times[tmp]=times[tmp]+1

    for x in anno_seq:
        if times.has_key(x):
            annoDeter_sta_f.write(str(x)+"\t"+str(times[x])+"\t"+str(times[x]/Ng)+'\n')
            print str(x)+"\t"+str(times[x])+"\t"+str(times[x]/Ng)
        else:
            annoDeter_sta_f.write(x+"\t"+"0"+"\t"+"0"+'\n')
            print x+"\t"+"0"+"\t"+"0"

    print "----------------------------------"

    ##write the results to the file

    for l in xrange(0, Ng):
        for i in xrange(cmp_len):
            annoDeter_f.write(str(res_list[i][l])+"\t")
            
        for y in annoDet_list[l]:
            annoDeter_f.write(str(y)+"\t")
        annoDeter_f.write(str(nonhp_list[l])+"\t")
        annoDeter_f.write(", ".join([str(x+1) for x in lab_res[l]]))
        annoDeter_f.write("\n")

    annoDeter_f.close()
    annoDeter_sta_f.close()
    
##    ##write the results to the log file
##    logfile=open('results/'+org_abb+' anno deter hist.txt', 'a')
##    logfile.write("------------------------------------------------")
##    for l in xrange(0, Ng):
##        for i in xrange(cmp_len):
##            logfile.write(str(res_list[i][l])+"\t")
##            
##        for y in annoDet_list[l]:
##            logfile.write(str(y)+"\t")
##        logfile.write(str(nonhp_list[l])+"\t")
##        logfile.write(", ".join([str(x+1) for x in lab_res[0]]))
##        logfile.write("\n")
##    logfile.write("------------------------------------------------")
##    logfile.close()

def extractInfo(anno_proc_list, res_list):
    #####################loading non-hp data
    dif_tag_entire=set(['ec-d', 'both miss', '0', 'base_ec-d', 'base_both miss', 'base_one miss', 'base_one HP'])

    nonhp_list=list()
    lab_res=list()
    Ng=len(anno_proc_list[0])
    res_len=len(res_list)

    ind_list=list()
    cmp_len=int(res_len*(res_len-1)/2)
    
    for i in xrange(0, cmp_len):
        for j in xrange(i+1, cmp_len):
            ind_list.append([i, j])
            
    for l in xrange(0, Ng):
        k=0
        t=0
        lab_res.append([])
        lab_tmp=list()

        for i in xrange(0, cmp_len):        
            if not res_list[i][l] in dif_tag_entire:
                lab_tmp.append(ind_list[i][0])
                lab_tmp.append(ind_list[i][1])
                
        for i in xrange(0, res_len):
            if not anno_proc_list[i][l]==["hypothetical"] and not anno_proc_list[i][l]==[]:
                k=k+1

        lab_res[l]=list(set(lab_tmp))
        nonhp_list.append(k)
        
    return(nonhp_list, lab_res)
