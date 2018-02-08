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

from datetime import datetime


##_DEBUG=True#False#
##
##if _DEBUG == True:
##    import pdb 
##    pdb.set_trace()


####################################Loading data################################
if __name__ == '__main__':
    '''
    The module aims to use the annotation comparison results to derive a consensus
    annotation result. 
    '''
    tstart=datetime.now()
    gc.disable()
    print os.getcwd()
    os.chdir("..")

    org_abb_list=['ecol', 'ctra', 'hpyl', 'mgen', 'mytu', 'rpro']
    org_cap_list=['ECOLI', 'CHLTR', 'HELPY', 'MYCGE', 'MYCTU', 'RICPR']

    for i in xrange(len(org_abb_list)):        
        org_abb=org_abb_list[i]
        org_cap=org_cap_list[i]
        
        ############################Loading data################################
        pkl_anno_proc_list=open('tmp/'+org_abb+'_proc_list'+'.pkl','rb')
        anno_proc_list=pickle.load(pkl_anno_proc_list)
        pkl_anno_proc_list.close()
                        
        pkl_anno_list=open('tmp/'+org_abb+'_anno_list'+'.pkl','rb')
        anno_list=pickle.load(pkl_anno_list)
        pkl_anno_list.close()

        pkl_res_list=open('results/'+org_abb+'_res_list'+'.pkl','rb')
        res_list=pickle.load(pkl_res_list)
        pkl_res_list.close()

        print "loading data done"
        print datetime.now()-tstart
        
        ############################Annotation determination################################
        annoDeter_res_f='results/'+org_abb+'_anno_deter.txt'
        annoDeter_sta_f='results/'+org_abb+'_sta_anno_deter.txt'
        
        m.annoDeter(anno_list, anno_proc_list, res_list, annoDeter_res_f, annoDeter_sta_f)
        print "Annotation determination done"
        print datetime.now()-tstart
