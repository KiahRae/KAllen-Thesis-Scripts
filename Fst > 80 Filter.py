# This script filters out the SNPs with Fst > 80

fx = open ('Prob_McKendrick.csv')
cont = fx.readlines()
fx.close

cont2=cont[0:]
#print(cont2)

diallelic=[]
for line in cont2:
    if '#DIV/0!' not in line:
        new = line.replace('\n','').split(',')
        diallelic.append(new)
print(diallelic[99])


list=[]
for a in diallelic[2:]:
    if float(a[10]) > 0.80:
        list.append(a)

print(list)


######### INDICES FOR EACH SPECIES #########
#Ribicola 1:5
#Comandrae 6:10
#Hybrid 11:15
#print(list)

#Rib(0|0) Rib(1|0)	Rib(0|1) Rib(1|1) Rib(.|.) Com(0|0)	Com(1|0) Com(0|1) Com(1|1) Com(.|.)	Flex (0|0) Flex(1|0) Flex(0|1) Flex1|1) Flex (.|.)

newlist=[]
for h in list:
    newlist.append(h[0] + '\t' + h[1] + '\t' +  h[2] + '\t' +  h[3] + '\t' + h[4] + '\t' + h[5] + '\t' + h[6] + '\t' + h[7] + '\t' +  h[8] + '\t' + h[9] + '\t' + h[10])
print(newlist)


import numpy as np
np.savetxt("BowsMck_Fst_80.txt", newlist, delimiter=",", fmt='%s')

