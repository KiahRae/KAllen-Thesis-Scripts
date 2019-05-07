# This script uses the output of the Fast-GBS filtered SNPs as input
# It pulls out each sample individually and then groups them by sample site and species
# It then writes a CSV with the data for the SNPs in each sample site

fx = open ('ribi_SNPs_diallelic.recode.vcf')
cont = fx.readlines()
fx.close

cont2=cont[0:]

cont2=cont[47:]


diallelic=[]
for c in cont2:
    c=c.replace(' ', '')
    if ',' not in c[25:28]:
        diallelic.append(c)

line=0
ls=[]
for line in diallelic:
    new = line.replace('"','').split('\t')
    ls.append(new)
#print(ls[0])

######Lists to separate the different species by isolates ######## ##### the rows of data are separated by ','
Ribicola = []
for R in ls:
    Ribicola.append(R[1] + ':' + R[20]+':' + R[21] + ':' + R[22] +':' + R[23] + ':' + R[24] + ':' + R[30] + ':' + R[41] + ':' + R[46] + ':'+ R[47] + ':' + R[48] + ':'+ R[49] + ':' + R[50] + ':' + R[51] + ':' + R[52]+ ':' + R[53]+ ':' + R[54] + ':' + R[55] + ':' + R[63] + ':' + R[69] + ':' + R[70] + ':'+ R[71] + ':' + R[72] + ':' + R[73] + ':' + R[74] + ':' + R[75] + ':' + R[76] + ':' + R[77] + ':' + R[78] + ':' +  R[85] + ':' +R[96] + ':' + R[107]+ ':' + R[36]+ ':' + R[37] + ':'+ R[38] + ':' + R[39] + ':' + R[40]+ ':' + R[42] + ':' + R[43] + ':'+ R[44] + ':' + R[45] + ':' + R[129] +':' + R[140] + ':' +R[151] + ':' + R[162] + ':' + R[173] + ':' + R[195] + ':' +R[184] + ':' + R[206] + ':' + R[316] + ':' + R[327] + ':' + R[333] + ':' + R[334]+ ':' + R[335] + ':' + R[336]+ ':' + R[337] + ':' + R[338]+ ':' + R[339] + ':' + R[340]+ ':' + R[341] + ':' + R[342]+ ':' + R[343] + ':' + R[344]+ ':' + R[345] + ':' + R[346] + ':' + R[347] + ':' + R[348] + ':' + R[391] + ':' + R[392] + ':' + R[393] + ':' + R[68])
#print(Ribicola[0])

Comandrae = []
for C in ls:
    Comandrae.append(C[1] + ':' + C[26]+':'+C[27]+':'+C[28]+':'+ C[29] + ':' + C[31] + ':' + C[32] + ':' + C[33] + ':' + C[34] + ':' + C[35] + ':' + C[56] + ':' + C[57] + ':' + C[58] + ':' + C[59]+ ':' + C[60] + ':' + C[61] + ':' + C[62] + ':' + C[64] + ':' + C[79] + ':'+ C[349] + ':' + C[350] + ':' + C[351] + ':' + C[352])
print(Comandrae[0])

Hybrid = []
for H in ls:
    Hybrid.append(H[1] + ':' + H[66] + ':' + H[67])
#print(Hybrid)

########### Lists for each Hybrid Zone/ Positive Control Separately ##########
McKendrick=[]
for M in ls:
    McKendrick.append(M[1] + ':' + M[36] + ':' + M[37] + ':'+ M[38] + ':' + M[39] + ':' + M[40]+ ':' + M[42] + ':' + M[43] + ':'+ M[44] + ':' + M[45] + ':' + M[129] +':' + M[140] + ':' +M[151] + ':' + M[162] + ':' + M[173] + ':' + M[195] + ':' +M[184] + ':' + M[206])
#print(McKendrick[1])

Bowser=[]
for B in ls:
    Bowser.append(B[1] + ':' + B[19] + ':' + B[30] + ':' + B[41] + ':' + B[46] + ':'+ B[47] + ':' + B[48] + ':'+ B[49] + ':' + B[50] + ':' + B[51] + ':' + B[52]+ ':' + B[53]+ ':' + B[54] + ':' + B[55] + ':' + B[63] + ':' + B[74] + ':' +  B[85] + ':' + B[96] + ':' + B[107])
#print(Bowser[0])

RockyMtn=[]
for RM in ls:
    RockyMtn.append(RM[1] + ':' + RM[10]+':' + RM[11] + ':' + RM[12] +':' + RM[13] + ':' + RM[14] + ':' + RM[15]+':' + RM[16] + ':' + RM[17] +':' + RM[18] + ':' + RM[316] + ':' + RM[327] + ':' + RM[333] + ':' + RM[334]+ ':' + RM[335] + ':' + RM[336]+ ':' + RM[337] + ':' + RM[338]+ ':' + RM[339] + ':' + RM[340]+ ':' + RM[341] + ':' + RM[342]+ ':' + RM[343] + ':' + RM[344]+ ':' + RM[345] + ':' + RM[346] + ':' + RM[347] + ':' + RM[348]+ ':' + RM[217]+ ':' + RM[228] + ':' + RM[239] + ':' + RM[250] + ':' + RM[261] + ':' + RM[272] + ':' + RM[283] + ':' + RM[296])
#print(RockyMtn[0])


#so I need to keep all of the columns so I need to split it by ':' and replace '"'

geno=[]
for isolate in ls:
    genotype=str(isolate).split(':')
    geno.append(genotype)
#print(geno)


data=[]
for character in geno:
    cleanedup=str(character).replace('"','')
    data.append(cleanedup)
#print(data)


####Now I need to calculate the Fst Values for each parental species across each line
##### To do this I need the 0/0, 1/1, 0/1 and 1/0 counts for each species for each row
'''
###### RIBICOLA ######

ribigene=[]
for isolate1 in Ribicola:
    genotype1=str(isolate1).split(':')
    ribigene.append(genotype1)
#print(ribigene)

data1=[]
for character1 in ribigene:
    cleanedup1=str(character1).replace('"','').replace("'", '').replace(' ', '').replace('[','').replace(']','')
    data1.append(cleanedup1)
#print(data1[5])


data2=[]
for t in data1:
    toto = t.split(',')
    data2.append(toto)
#print(data2[5])

data3=[]
for up in data1:
    tata = up.replace(',','\t')
    data3.append(tata)

#import numpy as np
#np.savetxt("fixed_loci_ribi.txt", data3, delimiter=",", fmt='%s')

'''
##########COMANDRAE##########

comagene=[]
for isolate1 in Comandrae:
    genotype1=str(isolate1).split(':')
    comagene.append(genotype1)
#print(comagene)

data1=[]
for character1 in comagene:
    cleanedup1=str(character1).replace('"','').replace("'", '').replace(' ', '').replace('[','').replace(']','')
    data1.append(cleanedup1)
print(data1[5])


data2=[]
for t in data1:
    toto = t.split(',')
    data2.append(toto)
#print(data2[5])

data3=[]
for up in data1:
    tata = up.replace(',','\t')
    data3.append(tata)

#import numpy as np
#np.savetxt("fixed_loci_coma.txt", data3, delimiter=",", fmt='%s')

'''
###########HYBRID#############
hybridgene=[]
for isolate1 in Hybrid:
    genotype1=str(isolate1).split(':')
    hybridgene.append(genotype1)
#print(hybridgene)

data1=[]
for character1 in hybridgene:
    cleanedup1=str(character1).replace('"','').replace("'", '').replace(' ', '').replace('[','').replace(']','')
    data1.append(cleanedup1)
print(data1[5])


data2=[]
for t in data1:
    toto = t.split(',')
    data2.append(toto)
#print(data2[5])

data3=[]
for up in data1:
    tata = up.replace(',','\t')
    data3.append(tata)


########### McKENDRICK PASS #############
mckendgene=[]
for isolate1 in McKendrick:
    genotype1=str(isolate1).split(':')
    mckendgene.append(genotype1)
#print(mckendgene)

data1=[]
for character1 in mckendgene:
    cleanedup1=str(character1).replace('"','').replace("'", '').replace(' ', '').replace('[','').replace(']','')
    data1.append(cleanedup1)
print(data1[5])


data2=[]
for t in data1:
    toto = t.split(',')
    data2.append(toto)
#print(data2[5])

data3=[]
for up in data1:
    tata = up.replace(',','\t')
    data3.append(tata)

########### Bowser #############
bowsergene=[]
for isolate1 in Bowser:
    genotype1=str(isolate1).split(':')
    bowsergene.append(genotype1)
#print(bowsergene)

data1=[]
for character1 in bowsergene:
    cleanedup1=str(character1).replace('"','').replace("'", '').replace(' ', '').replace('[','').replace(']','')
    data1.append(cleanedup1)
print(data1[5])


data2=[]
for t in data1:
    toto = t.split(',')
    data2.append(toto)
#print(data2[5])

data3=[]
for up in data1:
    tata = up.replace(',','\t')
    data3.append(tata)


########### Rocky Mountain #############

RockyMtngene=[]
for isolate1 in RockyMtn:
    genotype1=str(isolate1).split(':')
    RockyMtngene.append(genotype1)
#print(RockyMtngene)

data1=[]
for character1 in RockyMtngene:
    cleanedup1=str(character1).replace('"','').replace("'", '').replace(' ', '').replace('[','').replace(']','')
    data1.append(cleanedup1)
#print(data1[3])


data2=[]
for t in data1:
    toto = t.split(',')
    data2.append(toto)
#print(data2[0])
'''
fx = open ('Bowser_Mck_POS_ID.csv')
coy = fx.readlines()
fx.close

print(len(coy))

pos=[]
for id in coy:
    id2 = id.replace('\n','')
    pos.append(id2)
#print(pos[0])

poslis=[]
i=0
for lnie in data2:
    if data2[i][0] in pos:
       poslis.append(lnie)
    i=i+1
print(len(poslis))


data4=[]
for tt in poslis:
    titi = str(tt).replace(',', '\t').replace("'",'')
    data4.append(titi)

########### This gives me the raw Excel Files -- The rest of the calculations are done in Excel
#import numpy as np
#np.savetxt("fixed_PG1.txt", data4, delimiter=",", fmt='%s')
