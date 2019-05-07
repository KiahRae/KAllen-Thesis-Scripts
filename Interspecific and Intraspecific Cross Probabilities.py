'''
############ for the genotype 0/1 where 1 = com and 0 = rib ##########

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

fig = plt.figure()
fx = open ('Bowser_Prob_0Rib_1Com.csv')
cont = fx.readlines()
fx.close

Bowser=cont[0:]


fy = open ('Rocky_Prob_0rib_1com.csv')
ct = fy.readlines()
fy.close

Rocky=ct[0:]

fk = open ('Mcken_Prob_0Rib_1Com.csv')
ckt = fk.readlines()
fk.close

Mckendrick =ckt[0:]

ax1 = plt.subplot2grid((2,2), (0, 0))
ax2 = plt.subplot2grid((2,2),(0, 1)) #define your graph first
ax3 = plt.subplot2grid((2,2),(1, 0),colspan=2)

Bowser.sort()
Rocky.sort()
Mckendrick.sort()


ax1.hist(Rocky,bins=30)
ax2.hist(Bowser,bins=30)
ax3.hist(Mckendrick,bins=30)



ax1.set_title('Rocky Mountain')
ax2.set_title('Bowser')
ax3.set_title('Mckendrick')
ax3.set_xlabel('Probability observed genotype (0/1) is from an interspecific cross')
ax3.set_ylabel('# of Individuals')
plt.show()

############ for the genotype 0/1 where 0 = com and 1 = rib ##########
'''
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

fig = plt.figure()
fx = open ('Bowser_Prob_1Rib_0Com.csv')
cont = fx.readlines()
fx.close

Bowser=cont[0:]


fy = open ('Rocky_Prob_1Rib_0Com.csv')
ct = fy.readlines()
fy.close

Rocky=ct[0:]

fk = open ('Mcken_Prob_1Rib_0Com.csv')
ckt = fk.readlines()
fk.close

Mckendrick =ckt[0:]

ax1 = plt.subplot2grid((2,2), (0, 0))
ax2 = plt.subplot2grid((2,2),(0, 1)) #define your graph first
ax3 = plt.subplot2grid((2,2),(1, 0),colspan=2)

Bowser.sort()
Rocky.sort()
Mckendrick.sort()


ax1.hist(Rocky,bins=30)
ax2.hist(Bowser,bins=30)
ax3.hist(Mckendrick,bins=30)



ax1.set_title('Rocky Mountain')
ax2.set_title('Bowser')
ax3.set_title('Mckendrick')
ax3.set_xlabel('Probability observed genotype (0/1) is from an interspecific cross')
ax3.set_ylabel('# of Individuals')
plt.show()



######## one way ANOVA #################
print(stats.f_oneway(Bowser,Mckendrick,Rocky))

######### shapiro-wilks test for normality
print(stats.shapiro(Bowser))
print(stats.shapiro(Mckendrick))
print(stats.shapiro(Rocky))

##### Kruskal-Wallis H-test tests the null hypothesis that the population median of all of the groups are equal. It is a non-parametric version of ANOVA
print(stats.kruskal(Bowser,Mckendrick,Rocky))

import statsmodels.api as sa
import scikit_posthocs as sp
