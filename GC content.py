# Checking the contamination of the Hybrid GBS data
# Takes the two files from mapping the species - from the sppIDer output mapped and unmapped reads
# Subsample of 10000 reads per file

fx = open ('COMunmapped.fasta')
cont = fx.readlines()
fx.close

cont2=cont[0:]

# Separating Nucleotides
x=0
ls=[]
for x in cont2:
    if '>J0W1A' not in x:
        ls.append(x)

# Cleaning up the list
nwls=[]
for y in ls:
    toto = y.replace('\n', '')
    nwls.append(toto)

#Joining the mapped nucleotides into one string

DNA = ''.join(nwls)

#print(DNA)

# Counting the GC content and calculating the percentage

GC = DNA.count('G')
GC2 = GC + DNA.count('C')

percent = (float(GC2)/len(DNA))*100

print(percent)

# This gives the percent GC content and then see if the averages between mapped and unmapped differ
# Then change the script to the distribution of GC content read by read