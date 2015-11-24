#!/usr/bin/python

import matplotlib as MPL
MPL.use('agg') # no X (so show won't work)

from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
#from matplotlib import rc   #for adding italics. Via latex style
#rc('text', usetex=True)
import pylab as P
import math
import numpy
import commands
import sys
from scipy import stats

DATA_DIR='/mnt/lanzarobas/home/bradmain/arabiensis/VCFs/'
FST_LIM = [-0.05, 0.25]
DSTAT_LIM = [-40, 60]
#FST_COLOR = 'b'
FST_SIG_COLOR = 'b'
DSTAT_COLOR = 'r'
INV_HEIGHT=0.05
#TITLE="Sequence Differentiation Between Homozygous 2Rb Inversion States (PCA3 Split)"
TITLE="Genome-wide FST (sliding windows)\nbetween PCA Clusters"

LEG_LINES = []
LEG_LABELS = []


#input windowed FST from vcftools

fig, axes = P.subplots(ncols=2,nrows=3)
fig.subplots_adjust(left=None, bottom=None, right=None, top=None,
                    wspace=.01, hspace=None)
((chrX,N), (chr2R, chr2L), (chr3R, chr3L)) = axes
#N.axis('off')

"""
#Second Y axis
chrXd = chrX.twinx()
chr2Rd = chr2R.twinx()
chr2Ld = chr2L.twinx()
chr3Rd = chr3R.twinx()
chr3Ld = chr3L.twinx()
"""

def smoothListGaussian(list,strippedXs=False,degree=5):  
     window=degree*2-1  
     weight=numpy.array([1.0]*window)  
     weightGauss=[]  
     for i in range(window):  
	i=i-degree+1  
	frac=i/float(window)  
	gauss=1/(numpy.exp((4*(frac))**2))  
	weightGauss.append(gauss)  
	weight=numpy.array(weightGauss)*weight  
	smoothed=[0.0]*(len(list)-window)  
	for i in range(len(smoothed)):  
		smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
	return smoothed 

inversions=["2Rc","2Rb","2La"]
## plot inversions
inv={}
for line in open("/mnt/lanzarobas/home/bradmain/gambiae/gene_flow/pest_M/An_gambiae_karyotype.gtf"):
    i=line.strip().split()
    chr=i[0]
    l=int(i[3])
    r=int(i[4])
    name=i[9].strip(";").strip('"')
    if name not in inversions:
        continue
    num=int(i[-1].strip(";").strip('"'))
    if chr not in inv:
        inv[chr]={}
    if name not in inv[chr]:
        inv[chr][name]={}
    inv[chr][name][num]=[l/1.0e6,r/1.0e6]
outer=[inv["2R"]["2Rb"][1][0],inv["2R"]["2Rb"][2][1]]
inner=[inv["2R"]["2Rb"][1][1],inv["2R"]["2Rb"][2][0]]
Couter=[inv["2R"]["2Rc"][1][0],inv["2R"]["2Rc"][2][1]]
Cinner=[inv["2R"]["2Rc"][1][1],inv["2R"]["2Rc"][2][0]]
outer2La=[inv["2L"]["2La"][1][0],inv["2L"]["2La"][2][1]]
inner2La=[inv["2L"]["2La"][1][1],inv["2L"]["2La"][2][0]]
#for N in inv["2R"]["2Rb"]:
#    outer.append(inv["2R"]["2Rb"][N][1])
#    inner.append(inv["2R"]["2Rb"][N][0])

print 'outer',outer
print 'inner',inner
#chr2R.plot(outer,[INV_HEIGHT,INV_HEIGHT],'k-',linewidth=5,alpha=0.5)
#chr2R.plot(inner,[INV_HEIGHT,INV_HEIGHT],'y-',linewidth=5)
#chr2R.plot(Couter,[INV_HEIGHT,INV_HEIGHT],'k-',linewidth=5,alpha=0.5)
#chr2R.plot(Cinner,[INV_HEIGHT,INV_HEIGHT],'g-',linewidth=5)
#chr2L.plot(outer2La,[INV_HEIGHT,INV_HEIGHT],'k-',linewidth=5,alpha=0.5)
#chr2L.plot(inner2La,[INV_HEIGHT,INV_HEIGHT],'y-',linewidth=5)
#chr3R.plot([12.5,38],[INV_HEIGHT,INV_HEIGHT],'y-',linewidth=5,alpha=0.5)
chr2R.plot(outer,[INV_HEIGHT,INV_HEIGHT],'k-',linewidth=15,alpha=0.5)
chr2R.plot(inner,[INV_HEIGHT,INV_HEIGHT],'y-',linewidth=15,label='2Rb inversion')
chrX.plot(inner,[INV_HEIGHT+1000,INV_HEIGHT+1000],'y-',linewidth=15,label='2Rb inversion') #just plotting out of range on X for legend purposes
chr2R.text(numpy.mean(inner)-.5,INV_HEIGHT-0.01,'b',fontweight='bold',fontsize=14)
chr2R.plot(Couter,[INV_HEIGHT,INV_HEIGHT],'k-',linewidth=15,alpha=0.5)
chr2R.plot(Cinner,[INV_HEIGHT,INV_HEIGHT],'g-',linewidth=15,label='2Rc inversion')
chrX.plot(Cinner,[INV_HEIGHT+1000,INV_HEIGHT+1000],'g-',linewidth=15,label='2Rc inversion') #just plotting out of range on X for legend purposes
chr2R.text(numpy.mean(Cinner)-.5,INV_HEIGHT-0.01,'c',fontweight='bold',fontsize=14)
#chr2L.plot(outer2La,[INV_HEIGHT,INV_HEIGHT],'k-',linewidth=15,alpha=0.5)
#chr2L.plot(inner2La,[INV_HEIGHT,INV_HEIGHT],'y-',linewidth=15)
chr3R.plot([12.5,38],[INV_HEIGHT,INV_HEIGHT],'r-',linewidth=15,alpha=0.5,label='3Ra inversion')
chrX.plot([12.5+1000,38+1000],[INV_HEIGHT,INV_HEIGHT],'r-',linewidth=15,alpha=0.5,label='3Ra inversion') #just plotting out of range on X for legend purposes
chr3R.text(numpy.mean([12.5,38]),INV_HEIGHT-0.01,'a',fontsize=14,fontweight='bold')
#chr3R.legend()
or7=[22.849252,22.858650]
or40=[22.823983,22.825656]
gr53=[24.694665,24.698605]
gr13=[24.811173,24.812613]
or39=[24.850239,24.851846]
or38=[24.857474,24.859095]

def fst_plotter(fst_files,FST_COLOR,style,newLEGEND):    
    fstD={}
    fstmean={}
    leg_done = False
    for file in fst_files:
        for line in open(file):
            i=line.strip().split()
            chr=i[0]
            #skip unknown and Y chromosomes
            if chr=="CHROM" or chr=="UNKN" or chr=="Y_unplaced":
                continue
            if chr not in fstD:
                fstD[chr]={}
                fstmean[chr]={}
            pos=int(i[1])+24999   #moves x position to middle of 50kb bin
            if i[2]=="-nan":
                continue
            fst=float(i[4]) #i[4] is the weighted fst
            fstM=float(i[5]) #i[5] is the mean fst
            if pos not in fstD[chr]:
                fstD[chr][pos]=fst
                fstmean[chr][pos]=fstM
    F=[]
    Fs=[]
    for CHROM in fstD:
        x=numpy.array(sorted(fstD[CHROM]))
        xmean=sorted(fstmean[CHROM])
        y=[]
        ymean=[]
        for i in x:
            F.append(fstD[CHROM][i])
            y.append(fstD[CHROM][i])
            ymean.append(fstmean[CHROM][i])

        ax = globals()['chr'+CHROM]
        #tmp, = ax.plot(x/1.0e6, y, '-', color=FST_COLOR, linewidth=1.5)
        tmp, = ax.plot(x/1.0e6, y, style, color=FST_COLOR, linewidth=1.5,label=newLEGEND)
        #if( not leg_done ):
        #    LEG_LINES.append(tmp)
        #    LEG_LABELS.append(r"$F_{\mathrm{ST}}$ pre- vs post-2006 $A. coluzzii$")
        #    leg_done = True
    chrX.legend()
    #LEG_LINES.append(leg_fst_sig)
    #LEG_LABELS.append(r"$F_{\mathrm{ST}}$ 99.9 percentile level")





# actually plot fst (on top)
#fst_plotter([DATA_DIR+"pca1_pca2.windowed.weir.fst"],'b','--', "PCA1 vs PCA2")
#fst_plotter([DATA_DIR+"pca3_pca2.windowed.weir.fst"],'k','-', "PCA3 vs PCA2")
#fst_plotter(["pca1_pca2.windowed.weir.fst"],'b','--', "PCA1 vs PCA2")
#fst_plotter(["pca3_pca2.windowed.weir.fst"],'k','-', "PCA3 vs PCA2")
fst_plotter(["pca1_pca2.windowed.weir.fst"],'b','--', "Left PCA cluster vs middle")
fst_plotter(["pca3_pca2.windowed.weir.fst"],'k','-', "Right PCA cluster vs middle")



# chromosome names
for C in ['X', '2R', '2L', '3R', '3L']:
    ax = globals()['chr'+C]
    if( C[-1] == 'L' ):
        x = 0.975
        ha = 'right'
    else:
        x = 0.025
        ha = 'left'
        #ax.text(x, 0.95, r'\textbf{'+C+'}', size='xx-large', ha=ha, va='top', transform=ax.transAxes)
    ax.text(x, 0.95, C, size='xx-large', ha=ha, va='top', transform=ax.transAxes)


chrX.set_ylabel("$F_{\mathrm{ST}}$",color='b',fontsize=24)
chr2R.set_ylabel("$F_{\mathrm{ST}}$",color='b',fontsize=24)
chr3R.set_ylabel("$F_{\mathrm{ST}}$",color='b',fontsize=24)
chr3R.set_xlabel(r"position [Mb]",fontsize=24)
chr3L.set_xlabel(r"position [Mb]",fontsize=24)

chr2L.get_yaxis().set_visible(False)
chr3L.get_yaxis().set_visible(False)

chrX.set_ylim(FST_LIM)
chrX.set_xlim(0,22)
chr2L.set_ylim(FST_LIM)
chr2R.set_ylim(FST_LIM)
chr3L.set_ylim(FST_LIM)
chr3R.set_ylim(FST_LIM)


#P.show()
chrX.set_title(TITLE, y=1.04, fontsize=24)

##################### PCA PLOT

human=[line.strip() for line in open("../pca/allhumanfed.txt")]
cattle=[line.strip() for line in open("../pca/allcattlefed.txt")]


cattlex=[]
cattley=[]
humanx=[]
humany=[]
for line in open("../pca/LUPI_maf_pca.eigenvec"):
    i=line.strip().split()
    pc1=i[2]
    pc2=i[4]
    if i[1] in human:
        humanx.append(pc1) 
        humany.append(pc2) 
    elif i[1] in cattle:
        cattlex.append(pc1) 
        cattley.append(pc2) 
    else:
        print "not human or cattle-fed:", line.strip()
        gamx.append(pc1) 
        gamy.append(pc2) 
        ###P.text(pc1,pc2,i[1],color='g',fontsize=14)

ax = N
ax.set_xlim(-.4,.3)
ax.set_ylim(-.35,.45)
pos = ax.get_position()
pts = pos.get_points()
w = pts[1,0]-pts[0,0]
h = pts[1,1]-pts[0,1]
nw = w*0.6
nh = h*0.8
x0 = pts[0,0]+(w-nw)/2.0
y0 = pts[0,1]+0.01 #+(h-nh)
print pts, w, h
ax.set_position([x0, y0, nw, nh])

ax.plot(cattlex,cattley,'bo',label="cattlefed")
ax.plot(humanx,humany,'ro',label="humanfed")
#P.text(-.38,-.3,"P<0.01; humanfed vs cattlefed 2x3 Fisher Exact")
ax.set_xlabel("PCA1")
ax.set_ylabel("PCA2")
ax.set_xlim(-.4,.3)
ax.set_ylim(-.35,.45)
leg = ax.legend(numpoints=1, ncol=2, loc=8, bbox_to_anchor=(0.5, 1.01))
leg.get_frame().set_alpha(0.5)
#P.title(r"PCA on all \textit{An. arabiensis} SNPs",fontsize=20)
ax.set_title("PCA on Genome-wide SNPs",fontsize=24, y=1.34)


################ Final adjustments and save

fig.set_size_inches(14.4, 9.6)
#P.show()
P.savefig('pca_based_fst.svg', dpi=300)


