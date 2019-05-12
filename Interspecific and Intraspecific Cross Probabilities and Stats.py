
############# probability that the observed genotypes are interspecific crosses ###########

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt


fig = plt.figure()
fu = open ('Comandrae_distribution.csv')
ct = fu.readlines()
fu.close

cte=ct[0:]

########## Probability lists ##########
a=0
Comandrae=[]
for a in cte:
    if 'no data' not in a:
        Comandrae.append(a.replace('\n',''))
#print(Comandrae)


fp = open ('Bowser_distribution.csv')
cont = fp.readlines()
fp.close

cont1=cont[0:]

q=0
Bowser=[]
for q in cont1:
    if 'no data' not in q:
        Bowser.append(q.replace('\n',''))
#print(Bowser)


fx = open ('Rocky_distribution.csv')
cont = fx.readlines()
fx.close

cont1=cont[0:]

o=0
Rocky=[]
for o in cont1:
    if 'no data' not in o:
        Rocky.append(o.replace('\n',''))
#print(Rocky)


fy = open ('hybrid_distribution.csv')
ct = fy.readlines()
fy.close

ls=ct[0:]

x=0
Hybrid=[]
for x in ct:
    if 'no data' not in x:
        Hybrid.append(x.replace('\n',''))
#print(Hybrid)

fk = open ('Mckendrick_distribution.csv')
ckt = fk.readlines()
fk.close

liste =ckt[0:]

p=0
Mckendrick=[]
for p in liste:
    if 'no data' not in p:
        Mckendrick.append(p.replace('\n',''))
#print(Mckendrick)

convertedR = [float(num) for num in Rocky]
convertedM = [float(num) for num in Mckendrick]
convertedH = [float(num) for num in Hybrid]
convertedB = [float(num) for num in Bowser]
convertedC = [float(num) for num in Comandrae]

########## Mean and standard deviation of each list #############

Rockymean = np.mean(convertedR)
Rockystd = np.std(convertedR)

Hybridmean = np.mean(convertedH)
Hybridstd = np.std(convertedH)

Mckendrickmean = np.mean(convertedM)
Mckendrickstd = np.std(convertedM)

Bowsermean = np.mean(convertedB)
Bowserstd = np.std(convertedB)

Comandraemean = np.mean(convertedC)
Comandraestd = np.std(convertedC)

print(Rockymean, Mckendrickmean, Hybridmean, Bowsermean, Comandraemean)

print(len(Rocky))
Rocky.sort()
Hybrid.sort()
Mckendrick.sort()
Bowser.sort()
Comandrae.sort()

### Probability Distribution Plots ##################

ax1 = plt.subplot2grid((5,2), (0, 0), colspan=2)
ax2 = plt.subplot2grid((5,2),(1, 0), colspan=2) #define your graph first
ax3 = plt.subplot2grid((5,2),(2, 0),colspan=2)
ax4 = plt.subplot2grid((5,2),(3, 0),colspan=2)
ax5 = plt.subplot2grid((5,2),(4, 0),colspan=2)

ax1.hist(Hybrid,bins=30,color='lightgreen')
ax2.hist(Rocky,bins=30)
ax3.hist(Mckendrick,bins=30, color='grey')
ax4.hist(Bowser,bins=30,color='skyblue')
ax5.hist(Comandrae,bins=30,color='darkred')


ax1.legend(["$\it{Cronartium}$ $\it{x}$ $\it{flexili}$"])
ax2.legend(['Rocky Mountain sites'])
ax3.legend(['Mckendrick Pass'])
ax4.legend(['Bowser'])
ax5.legend(['Prince George sites'])

ax5.set_xlabel('Probability $\it{P(x)}$')
ax1.set_xticks([0,87])
ax2.set_xticks([0,59])
ax3.set_xticks([0,62])
ax4.set_xticks([0,66])
#ax3.set_ylabel('# Individual Genotypes')
ax1.text(31, 300, 'mean = 0.3175')
ax2.text(21, 4000, 'mean = 0.0491')
ax3.text(22, 2000, 'mean = 0.0523')
ax4.text(24,1700, 'mean=0.0536')
ax5.text(31,3000, 'mean=0.0417')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.spines['top'].set_visible(False)

plt.xticks([0,42,85])
ax1.set_title('Probability Observed Genotype is from an Interspecific Cross')
plt.subplots_adjust(left=0.15, right=0.9, top=0.96, bottom=0.07)
#plt.show()
#fig.savefig('LociProbability.png', dpi = 600)

####### Pairwise T-test
# t-test for Bowser
print('T-test Bowser v Hybrid')

print(stats.ttest_ind_from_stats(Bowsermean, Bowserstd, 7660, Hybridmean, Hybridstd, 1208, equal_var=True))
print(stats.ttest_ind_from_stats(Bowsermean, Bowserstd, 7660, Mckendrickmean, Mckendrickstd, 9144, equal_var=True))
print(stats.ttest_ind_from_stats(Bowsermean, Bowserstd, 7660, Comandraemean, Comandraestd, 7990, equal_var=True))
print(stats.ttest_ind_from_stats(Bowsermean, Bowserstd, 7660, Rockymean, Rockystd, 16486, equal_var=True))

# t-test for Mckendrick
print('T-test Mckendrick v Hybrid')
print(stats.ttest_ind_from_stats(Mckendrickmean, Mckendrickstd, 9144, Hybridmean, Hybridstd, 1208, equal_var=True))
print(stats.ttest_ind_from_stats(Mckendrickmean, Mckendrickstd, 9144,Comandraemean, Comandraestd, 7990, equal_var=True))
print(stats.ttest_ind_from_stats(Mckendrickmean, Mckendrickstd, 9144, Rockymean, Rockystd, 16486, equal_var=True))

# t-test for Rocky
print('T-test Rocky v Hybrid')
print(stats.ttest_ind_from_stats(Rockymean, Rockystd, 16486, Hybridmean, Hybridstd, 1208, equal_var=True))
print(stats.ttest_ind_from_stats(Rockymean, Rockystd, 16486, Comandraemean, Comandraestd, 7990, equal_var=True))

# t-test for Comandrae
print('T-test Comandrae v Hybrid')
print(stats.ttest_ind_from_stats(Comandraemean, Comandraestd, 7990, Hybridmean, Hybridstd, 1208, equal_var=True))


# mannwhitneyu pairwise t-test for non-normally distributed data
from scipy.stats import mannwhitneyu
print('Mannwhitney Bowser')
print(mannwhitneyu(Bowser,Hybrid))
print(mannwhitneyu(Bowser,Mckendrick))
print(mannwhitneyu(Bowser,Rocky))
print(mannwhitneyu(Bowser,Comandrae))

print('Mannwhitney Mckendrick')
print(mannwhitneyu(Mckendrick,Hybrid))
print(mannwhitneyu(Mckendrick, Rocky))
print(mannwhitneyu(Mckendrick, Comandrae))

print('Mannwhitney Rocky')
print(mannwhitneyu(Rocky, Hybrid))
print(mannwhitneyu(Rocky, Comandrae))

print('Mannwhitney Comandrae')
print(mannwhitneyu(Comandrae,Hybrid))



######## one way ANOVA #################

print(stats.f_oneway(Hybrid,Mckendrick,Rocky,Bowser,Comandrae))

groups = (Hybrid,Mckendrick,Rocky,Bowser,Comandrae)
import statsmodels.api as sm
from statsmodels.formula.api import ols
mod = ols((Hybrid,Mckendrick,Rocky,Bowser,Comandrae).fit())

aov_table = sm.stats.anova_lm(mod, typ=2)
print(aov_table)


######### shapiro-wilks test for normality

print(stats.shapiro(Hybrid))
print(stats.shapiro(Mckendrick))
print(stats.shapiro(Rocky))

##### Kruskal-Wallis H-test tests the null hypothesis that the population median of all of the groups are equal. It is a non-parametric version of ANOVA

print(stats.kruskal(Hybrid,Mckendrick,Rocky, Bowser, Comandrae))

