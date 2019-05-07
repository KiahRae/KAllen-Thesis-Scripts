# This gives the percent GC content and then see if the averages between mapped and unmapped differ
# Then change the script to the distribution of GC content read by read
#then run one-way ANOVA or t-test to determine if there is a significant difference


fx = open ('Subsample_Hybrid_Mapped.fasta')
cont = fx.readlines()
fx.close

cont2=cont[0:]

# Separating Nucleotides
x=0
ls=[]
for x in cont2:
    if '>' not in x:
        ls.append(x)

nwls=[]
for y in ls:
    toto = y.replace('\n', '')
    nwls.append(toto)

#print nwls [1:5]

MappedGC=[]
for read in nwls:

    GC = read.count('G')
    GC2 = GC + read.count('C')
    ratio1 = ((float(GC2) / len(read)))*100

    MappedGC.append(ratio1)

fk = open('SUBSAMPLE_UnMap.fasta')
ct = fk.readlines()
fk.close

ct2=ct[0:]
z=0
lis=[]
for z in ct2:
    if '>' not in z:
        lis.append(z)

nwlis=[]
for p in lis:
    tyty = p.replace('\n', '')
    nwlis.append(tyty)

UnmappedGC=[]
for rd in nwlis:
    GeC = rd.count('G')
    GeC2 = GeC + rd.count('C')
    ratio = ((float(GeC2) / len(rd)))*100
    UnmappedGC.append(ratio)

#print(GCls)

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

fig = plt.figure()


MappedGC.sort()
UnmappedGC.sort()
GClistmean = np.mean(MappedGC)
GCliststd = np.std(MappedGC)


plt.title('Per Sequence Guanine-Cytosine (GC) Content')
plt.xlabel('Mean GC content (%)')
plt.hist(MappedGC,bins=50,label = "mapped reads")
plt.hist(UnmappedGC,bins=50,alpha = 0.5,label = "unmapped reads")
plt.legend()



#One-way ANOVA
fvalue, pvalue = stats.f_oneway(MappedGC, UnmappedGC)
print(fvalue, pvalue)

#Equal variance
w, pvalue = stats.levene(MappedGC, UnmappedGC)
print(w, pvalue) #the pvalue < so there is unequal variance

#tot= stats.ttest_rel(MappedGC, UnmappedGC)

#print(tot)

print('t-test')

from scipy.stats import ttest_ind
t, p = ttest_ind(MappedGC, UnmappedGC)
print(p, t)

plt.text(60,900,"t-test = " + str(t.round(2)) + ", $P$ = " + str(p))

plt.show()

plt.savefig('GCcontent.png', dpi = 600)