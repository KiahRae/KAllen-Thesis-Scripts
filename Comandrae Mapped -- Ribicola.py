fx = open ('COMmapped.fasta')
cont = fx.readlines()
fx.close

cont2=cont[:146894]
#print(cont2[:200])


Ribi=cont2[:25000]

########### this writes a fasta file of the comandrae reads that mapped to ribicola #############
'''
ofile = open("Ribi_new.fasta", "w")
for x in range(len(Ribi)):
    ofile.write(Ribi[x])
#do not forget to close it
ofile.close()
'''

############## this code is to check the GC content of these reads ########################

#print(Ribi)

list=[]
for c in Ribi:
    c=c.replace('\n', '')
    if '>J0W1A' not in c:
        list.append(c)

print(list)

######### #Joining the nucleotides into one string
DNA = ''.join(list)

#print(DNA)

# Counting the GC content and calculating the percentage

GC = DNA.count('G')
GC2 = GC + DNA.count('C')

percent = (float(GC2)/len(DNA))*100

print(percent)

####### Counting the number of duplicates in the list #########################

def count_duplicates(list):

#takes as argument a sequence and returns the number of duplicate elements

    return len(list) - len(set(list)) ######### a set is basically a list of elements that cannot include duplicates ############

res = count_duplicates(list)
print(res)

print(len(list))

duplicat=[]
for u in list:
    if 'TGCAGCCCCAATAGCGCCTCGAATTCCTCCACCACCTACCCAGCTATTCCCACTTGACATCACCTGCTCATTGGCAATCGTGACAATGAGCGTCCG' in u:
        duplicat.append(u)

print(len(duplicat))
