# This script determines:
# the size of the genome assembly
# the number of scaffolds:
# the N50 number


fx = open('GCA_000464975.1_Cron_com_1.0_genomic.fasta')
cont = fx.read().split('>')
fx.close
#print cont
cont2 = cont[1:]

length=[] #empty list where we will put the length of each scaffold (to be used to find the total length of the genome)
scaffnum=[] #empty list where we will append one for each of the scaffolds (to be used to find the total number of scaffolds)
for x in cont2:
    toto=x.split('\n',1)[1].replace('\n','') #we are only interested in the ATGC's so we split it away from the name information using the first '\n'.  Because there are a bunch of other '\n's floating around we want to get rid of them becuase the would mess up our bp count so we replace them with nothing
    #print tata
    length.append(len(toto)) #adding the total length of each scaffold to the lsit 'length'
    scaffnum.append(1) #for each scaffold we add one to the list 'scaffnum'


length.append(len(toto))  # adding the total length of each scaffold to the list 'length'
scaffnum.append(1)  # for each scaffold we add one to the list 'scaffnum'

totalleng = sum(length)  # this is the total length of all the scaffolds and therefore the genome
print('totalleng')
print(totalleng)

totalscaff = sum(scaffnum)  # this is the total number of scaffolds
print('totalscaff')
print(totalscaff)

# Now we want to use the leng list a bit more to figure out the N50 value
length.sort()
length2 = length[::-1]
'''
print 'length2'
print length2  # this is the length of all the scaffolds in order from the largest to the smallest
'''
half = float(totalleng) / 2  # N50 uses half of the total genome size in the N50 value
print('half total genome size')
print(half)

# Now we want to find out, as we add together one scaffold at a time, the length of the final scaffold that makes the total reach or pass half of the total genome length

gettingtohalf = []  # empty list that we will add each of the scaffolds one at a time to
for x in length2:
    # print x
    if sum(gettingtohalf) < half:  # if the sum of 'gettingtohalf' is less than half the genome size then..
        gettingtohalf.append(float(x))  # ..append the length of the current scaffold to the list 'gettingtohalf'

print('gettingtohalf')
print(gettingtohalf[-1])


# length2 is the sorted list of contig lengths
# another loop to find N50

n = 0
i = 0
while i < len(length2):
    n = n + length2[i]
    if n >= float(sum(length2)) / 2:
        Nx = length2[i]
        i = len(length2)
    i = i + 1  # this only happens once when n is not >= float(sum(length2))/2

print('N50=')
print(Nx)
