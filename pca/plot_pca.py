import pylab as P
#from matplotlib import rc   #for adding italics. Via latex style
#rc('text', usetex=True)

human=[line.strip() for line in open("allhumanfed.txt")]
cattle=[line.strip() for line in open("allcattlefed.txt")]


cattlex=[]
cattley=[]
humanx=[]
humany=[]
for line in open("LUPI_maf_pca.eigenvec"):
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
P.plot(cattlex,cattley,'bo',label="cattlefed")
P.plot(humanx,humany,'ro',label="humanfed")
#P.text(-.38,-.3,"P<0.01; humanfed vs cattlefed 2x3 Fisher Exact")
P.xlabel("PCA1")
P.ylabel("PCA2")
P.xlim(-.4,.3)
P.ylim(-.35,.45)
P.legend()
#P.title(r"PCA on all \textit{An. arabiensis} SNPs",fontsize=20)
P.title("PCA on Genome-wide SNPs",fontsize=20)
P.savefig('/home/bradmain/public_html/2d_PCA.pdf', dpi=300)
P.show()
