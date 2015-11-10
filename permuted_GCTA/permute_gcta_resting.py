#!/usr/bin/python

import sys, os
import numpy as np
import copy


N = 10000

in_tab = np.loadtxt("phenotype_key.txt", dtype='str')
inds = range(1,49)

GCTA_OUT_BASE='LUPI_maf_48_resting'

#def permute(tab, col=2, skiprows=1): #permute the second column and skip the header.
def permute(tab, col=3, skiprows=1): #permute the second column and skip the header.
    tmp = tab[skiprows:,col].copy()
    np.random.shuffle(tmp)
    tab[skiprows:,col] = tmp
    return(tab)

def get_her():
    fp = open(GCTA_OUT_BASE+".hsq")
    lines = fp.readlines()
    fp.close()
    return(lines[4].split()[1])


re = np.empty(N)

tab = copy.deepcopy(in_tab)
for i in range(N):
    tab = permute(tab)
    print tab
    np.savetxt("permutes_resting/tmp_{}.pheno".format(i), tab[1:], fmt='%s')
    cmd = "./gcta64  --reml  --grm  LUPI_maf  --pheno " + "permutes_resting/tmp_{}.pheno".format(i) + " --mpheno 2 --reml-alg 2  --prevalence 0.43  --out LUPI_maf_48_resting --reml-maxit 10000 --thread-num 20"
    os.system(cmd)
    her = get_her()
    re[i] = her


np.savetxt(GCTA_OUT_BASE+"_heritability_10k.txt", re, fmt='%s')

