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
import array

from datetime import datetime
from numpy import zeros,dot
from datetime import datetime

##if _DEBUG == True:
##    import pdb 
##    pdb.set_trace()

if __name__ == '__main__':
    '''
    This module aims to compare the annotations from multiple sources, in order
    to use this module, you need to set up the following things:
    1. 'pfam_mapped_dic.csv', 'gene_info_dic.csv' files that can be found in
    datasets folder.
    '''
    tstart = datetime.now()
    gc.disable()
    print os.getcwd()
    os.chdir("..")
    
    myfolder="tmp/"
    
    org_abb_list=['ecol', 'ctra', 'hpyl', 'mgen', 'mytu', 'rpro']
    org_cap_list=['ECOLI', 'CHLTR', 'HELPY', 'MYCGE', 'MYCTU', 'RICPR']
    
    for i in xrange(len(org_abb_list)):
        
        org_abb=org_abb_list[i]
        org_cap=org_cap_list[i]
        print 'running '+org_abb
        
        #########################extract annotation results############################
        input_file=[org_abb+"_igs.txt", org_abb+"_img.txt", org_abb+"_rast.txt"]

        extract_res=m.extract_anno("database data/", input_file)
        anno_dir="database data/"
        anno_names=input_file
        
        word_dic=extract_res[0]
        anno_list=extract_res[1]
        anno_proc_list=extract_res[2]
        print 'extracting annotation done'
        
        #########################baseline comparison results############################
        pkl_anno_proc=open(myfolder+org_abb+'_proc_list'+'.pkl','wb')
        pickle.dump(anno_proc_list, pkl_anno_proc)
        pkl_anno_proc.close()
        
        pkl_word_dic=open(myfolder+org_abb+'_word_dic'+'.pkl','wb')
        pickle.dump(word_dic, pkl_word_dic)
        pkl_word_dic.close()
        
        pkl_anno_list=open(myfolder+org_abb+'_anno_list'+'.pkl','wb')
        pickle.dump(anno_list, pkl_anno_list)
        pkl_anno_list.close()

        ##read in the blast result data
        blast_file=open('database data/'+org_cap+' to '+org_cap+' cluster annotation.txt', 'r')
        i=0
        bla_col=list()

        while True:
            blast_line=blast_file.readline()
            if blast_line=='':
                break
            bla_col.append(m.proc(m.preproc(' '.join(blast_line.split("\t")[1:]))))

            i=i+1
            print i
        
        blast_file.close()
        print 'loading blast annotation done'
        pkl_bla_col=open(myfolder+org_abb+'_bla_col'+'.pkl','wb')
        pickle.dump(bla_col, pkl_bla_col)
        pkl_bla_col.close()
        
        pkl_bla_col=open('tmp/'+org_abb+'_bla_col'+'.pkl','rb')
        bla_col=pickle.load(pkl_bla_col)
        pkl_bla_col.close()
        
        pkl_word_dic=open('tmp/'+org_abb+'_word_dic'+'.pkl','rb')
        word_dic=pickle.load(pkl_word_dic)
        pkl_word_dic.close()

        pkl_anno_proc=open('tmp/'+org_abb+'_proc_list'+'.pkl','rb')
        anno_proc_list=pickle.load(pkl_anno_proc)
        pkl_anno_proc.close()
                        
        pkl_anno_list=open('tmp/'+org_abb+'_anno_list'+'.pkl','rb')
        anno_list=pickle.load(pkl_anno_list)
        pkl_anno_list.close()
        
        print "loading data"
        print datetime.now()-tstart

        #########################baseline comparison results############################
        ind_list=list()
        file_len=len(input_file)
        
        for i in xrange(0, file_len):
            for j in xrange(i+1, file_len):
                ind_list.append([i, j])

        res_list=list()
        tfidf_list=list()

        r_len=len(anno_list)
        c_len=len(anno_list[0])
        
        for r in xrange(file_len*(file_len-1)/2):
            res_list.append([])
            tfidf_list.append([])
            for c in xrange(c_len):
                res_list[r].append("0")
                tfidf_list[r].append(0)

        DocRes=m.constDocDict(anno_proc_list, bla_col)
        QueRes=m.constQueDict(anno_proc_list, DocRes[3])

        pkl_DocRes=open(myfolder+org_abb+'_DocRes'+'.pkl','wb')
        pickle.dump(DocRes, pkl_DocRes)
        pkl_DocRes.close()
        
        pkl_QueRes=open(myfolder+org_abb+'_QueRes'+'.pkl','wb')
        pickle.dump(QueRes, pkl_QueRes)
        pkl_QueRes.close()

        print "construct dictionary"
        print datetime.now()-tstart

        for ind in xrange(0,file_len*(file_len-1)/2):        
            ind1=ind_list[ind][0]
            ind2=ind_list[ind][1]

            ##baseline annotation comparison 
            res_list[ind]=m.comp_baseline(anno_list[ind1], anno_list[ind2])
            print "baseline done"
            print datetime.now()-tstart
        
        for ind in xrange(0,file_len*(file_len-1)/2):
            ind1=ind_list[ind][0]
            ind2=ind_list[ind][1]
            
            ##basic annotation ID comparison
            res_list[ind]=m.comp_id_anno(res_list[ind], anno_proc_list[ind1], anno_proc_list[ind2])
            print "basic comparison done"
            print datetime.now()-tstart

            ##use pfam dictionary to compare annotations
            f_pfam_dir="datasets/pfam_mapped_dic.csv"
            res_list[ind]=m.comp_pfam_dic(f_pfam_dir, res_list[ind], anno_list[ind1], anno_list[ind2], anno_proc_list[ind1], anno_proc_list[ind2], word_dic[ind1], word_dic[ind2])
            print "pfam comparison done"
            print datetime.now()-tstart
            
            ##compare using gene symbol dataset
            f_gs_dir="datasets/gene_info_dic.csv"
            res_list[ind]=m.comp_gs_dic(f_gs_dir, res_list[ind], anno_list[ind1], anno_list[ind2], anno_proc_list[ind1], anno_proc_list[ind2], word_dic, anno_list)
            print "gene symbol dictionary based comparison done"
            print datetime.now()-tstart

            ##use orthologous gene annotation to compare annotations
            f_orth_dir='database data/'+org_cap+' to '+org_cap+' cluster annotation.txt'
            res_list[ind]=m.comp_orth(f_orth_dir, res_list[ind], anno_list[ind1], anno_list[ind2], anno_proc_list[ind1], anno_proc_list[ind2])
            print "orthologs comparison done"
            print datetime.now()-tstart

        for ind in xrange(0,file_len):
            tfidf_list[ind]=m.tfidf(res_list, DocRes, QueRes, anno_proc_list[ind], bla_col)

        ##output the tfidf results
        tfidf_res_f='results/'+org_abb+' tfidf res.txt'
        tfidf_res_file=open(tfidf_res_f, "w")
        tfidf_rlen=len(tfidf_list[0])
        tfidf_clen=len(tfidf_list)
        
        for l in xrange(tfidf_rlen):
            for i in xrange(tfidf_clen):
                tfidf_res_file.write("%1.3f" % tfidf_list[i][l])
                tfidf_res_file.write("\t")
            tfidf_res_file.write("\n")
            
        tfidf_res_file.close()
        
        threshold=0.7
        
        for ind in xrange(0, file_len*(file_len-1)/2):
            ind1=ind_list[ind][0]
            ind2=ind_list[ind][1]

            ##use vector space model to compare annotations
            res_list[ind]=m.cmp_tfidf(res_list[ind], tfidf_list[ind1], tfidf_list[ind2], threshold)
            print "tf-idf comparison done"
            print datetime.now()-tstart

        dif_tag_entire=set(['ec-d', 'both miss', '0', 'base_ec-d', 'base_both miss', 'base_one miss', 'base_one HP'])

        ##map the same annotations from baseline result
        res_list=m.map_res(res_list, anno_proc_list, dif_tag_entire)
        print "mapping between the entire results done"
        print datetime.now()-tstart

        #print out the comparison results step by step
        for ind in xrange(0,file_len*(file_len-1)/2):
            m.printResult(res_list[ind], 'results/'+org_abb+'_sta_res.txt')

        #print out the baseline and improvement results
        for ind in xrange(0,file_len*(file_len-1)/2):
            m.printResult2(res_list[ind], 'results/'+org_abb+'_plot_res.txt')
            
        ##output the comparison results
        res_f='results/'+org_abb+' pairwise comp.txt'
        res_file=open(res_f, "w")
        l_len=len(res_list[0])
        res_len=len(res_list)
        
        for l in xrange(l_len):
            for i in xrange(file_len):
                res_file.write('\t'.join(anno_list[i][l].split('\t')[:3]))
                res_file.write("\t")
            for i in xrange(res_len):
                res_file.write(res_list[i][l])
                res_file.write("\t")
            res_file.write("\n")
            
        res_file.close()

        pkl_tfidf_list=open('results/'+org_abb+'_tfidf_list'+'.pkl','wb')
        pickle.dump(tfidf_list, pkl_tfidf_list)
        pkl_tfidf_list.close()

        pkl_res_list=open('results/'+org_abb+'_res_list'+'.pkl','wb')
        pickle.dump(res_list, pkl_res_list)
        pkl_res_list.close()
