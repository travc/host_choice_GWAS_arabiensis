import pylab as P

Hvalue_from_real_data = 0.963332

x=[]
count=0
pvalue=0
for line in open("LUPI_maf_48_heritability_10k.txt"):
    count=count+1
    if float(line.strip()) > Hvalue_from_real_data:
        pvalue=pvalue+1
    x.append(float(line.strip()))
n, bins, patches = P.hist(x, 500)
print float(pvalue)/float(count) #proportion of permutations with a heritability > than the real data. ~ permuted pvalue
P.show()
