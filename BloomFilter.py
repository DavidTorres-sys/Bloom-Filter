#!/usr/bin/env python
# coding: utf-8

# In[8]:


get_ipython().system('pip install bitmap')
get_ipython().system('pip install primesieve')


# In[1]:


from math import log
#class to open the files
class ReadFile:
    #read the kmer file
    def read_kmers(file_path, k):
        kmers = []
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                line = line.strip().upper()
                for i in range(len(line) - k + 1):
                    kmers.append(line[i:i+k])
        return kmers

    #read a fastQ file
    def read_fastq(file_path):
        #"a read is divide by 4 lines" 
        #"the first line give us the name and general info"
        #"the second line is the sequences of the reads"
        #"the third line is really no important"
        #"the fourth line is the base qualities"
        #we wanna store the sequencies and the base qualities of the reads
        sequences = []
        qualities = []
        #open de fastq file
        with open(file_path) as f:
            while True:
                #read the first line but no store
                f.readline()
                #store the base sequencies
                seq = f.readline().rstrip()
                #read the third line but no store
                f.readline()
                #read the base qualities
                qual = f.readline().rstrip()
                #check the end of the file and break the loop
                if len(seq) == 0:
                    break

                sequences.append(seq)
                qualities.append(qual)

        return sequences, qualities
        


# In[2]:


from bitmap import BitMap
from primesieve import n_primes
import numpy as np

def my_Bloom_filter(name, capacity):
    
    class BloomFilter:
        def __init__(self, name, capacity):
            self._buckets = n_primes(1, capacity+1)[0] 
            #optimal number of hash functions
            self._k = optimum_hash_function
            #Create a bit table (BitMap) with the calculated number of buckets
            self._table = BitMap(self._buckets)
            #Create a list of length k to store the hash functions
            self._functions = [None] * self._k 
            #generate k random numbers for a, b
            for i in range(self._k):
                self._a = np.random.randint(1, self._buckets-1)                               
                self._b = np.random.randint(0, self._buckets-1)
                self._functions[i] = (self._a, self._b)
            self._on = 0
            
        #Function to get the hash functions and the number of buckets used
        def getFunctions_(self):
            '''return constant a,b and bucket number.'''
            return self._functions, self._buckets 

        def hash_(self, key):
            '''return a list of the positions in the bit table'''
            return [((a * hash(key)+ b) % self._buckets) for (a, b) in self._functions] # return bucket number 

        #Function to insert an element into the filter
        def setItem_(self,pos):
            '''return the position to get into the hash table'''
            for i in range (self._k):
                ##Set each position corresponding to the key in the bit table to 1
                self._table.set(pos[i])
            #self._on = self._table.count()
            return pos

        # Function to check if a key is present in the filter
        def findItem_(self,key):
            '''return True if all positions are set to 1, False otherwise'''
            # Apply the hash function to the key
            pos = self.hash_(key)
            answer = True
            #Check if all positions corresponding to the key are set to 1 in the bit table
            for i in range(self._k):
                if self._table[pos[i]] == 0:
                    answer = False
                    break
            return answer
        
        # Function to calculate the filter load factor
        def load_factor_(self):
            return self._table.count()/ self._buckets
        
    return BloomFilter("table-1", capacity)
    


# In[3]:


def calculate_number_buckets(average_false_positive):
    '''return the optimum number of bucekts with respect to the expected number of false positives'''
    return round(-(10_000_000 * log(average_false_positive)) / (log(2)**2))

def calculate_hash_function(optimum_number_buckets):
    '''return the optimum number of hash function to use'''
    return round((optimum_number_buckets / 10_000_000) * log(2))


# In[4]:


#variable to handle the probability of false positives
average_false_positive = float(input('Ingrese la probabilidad de falsos positivos esperada: '))

#file path of the FastQ file
fastq_file_path = input('Ingrese la ruta del archivo FastQ: ')

#length of the k-mer
kmer_length = int(input('Ingrese el tamaño de los K-Mer: '))


# In[5]:


#variable to store the optimum number of buckets
optimum_number_buckets = calculate_number_buckets(average_false_positive)

#variable to store the number of hash function to use
optimum_hash_function = calculate_hash_function(optimum_number_buckets)


# In[6]:


#read the kmer file
kmers_file = ReadFile.read_kmers('kMer-BloomFilter.txt', kmer_length)

#store the sequences and qualities of the reads
sequences, qualities = ReadFile.read_fastq(fastq_file_path)

print(f'Reads {sequences[:5]}, K-Mer {kmers_file[:5]}')
print(f'numero de funciones hash a usar: { optimum_hash_function }, cantidad optima de cubetas: { optimum_number_buckets }')


# In[7]:


try:
    val = int(optimum_number_buckets)
except ValueError:
    print("No.. input string is not an Integer. It's a string")
    sys.exit("program is aborted....")
#instance of the class BloomFilter    
member = my_Bloom_filter("table-1", val)


# In[8]:


functions, n = member.getFunctions_()
for i in range(len(functions)):
    print(f'function is ( { functions[i][0] } * key + {functions[i][1]} ) module { n }')


# In[1]:


#add all the kmers into the hash table
for kmer in kmers_file:
    member.setItem_(member.hash_(kmer))

#loop through each sequence in the fastq file and calculate the kmers it contains
#total_sequences = len(sequences)
sequences_in_bloomfilter = 0

for seq in sequences:
    kmers = ''
    for i in range(len(seq) - kmer_length + 1):
        kmers.append(seq[i:i+kmer_length])
    
    
    #check if each calculated kmer is found 
    found = False
    for kmer in kmers:
        if member.findItem_(kmer):
            found = True
            break
    #count the number of sequences in the fastq file that have at least one kmer present in the BloomFilter
    if found: 
        sequences_in_bloomfilter += 1
        
    total_sequences = len(kmers)
        
print(total_sequences)

#percentage = (sequences_in_bloomfilter / total_sequences) * 100
#print(f"Porcentaje de secuencias del archivo fastq que hay respecto a las del kMer-BloomFilter.txt: {percentage:.2f}%")


# In[ ]:


#APPROXIMATE MATCHING


# In[92]:


#EDIT DISTANCE


# In[136]:


get_ipython().run_cell_magic('time', '', "#using recursion\ndef edit_distance_recursive(a, b):\n    if len(a) == 0: \n        return len(b)\n    elif len(b) == 0: \n        return len(a)\n    else:\n        dist_horz = edit_distance_recursive(a[:-1], b) + 1\n        dist_vert = edit_distance_recursive(a, b[:-1]) + 1\n        #if the char in the last position is diferent\n        #in the 2 string then delta = 1, otherwise delta = 0\n        if a[-1] == b[-1]:\n            dist_diag = edit_distance_recursive(a[:-1], b[:-1])\n        else:\n            dist_diag = edit_distance_recursive(a[:-1], b[:-1]) + 1\n        return min(dist_horz, dist_vert, dist_diag)\n__x = 'ATTGCTATTC'\n__y = 'ATATGACTTAC'\nedit_distance_recursive(__x, __y)")


# In[135]:


get_ipython().run_cell_magic('time', '', "#using dynamic programming\ndef edit_distance_dynamic(a, b):\n    '''return the edit distance between the strings'''\n    d = []\n    #initialize the matrix with 0\n    for i in range(len(a)+1):\n        d.append([0] * (len(b)+1))\n        \n    #the epsilon column and row is fill with the edit distance between the empty string \n    #and the position i of the string a and b.\n    for i in range(len(a)+1):\n        d[i][0] = i    \n    for i in range(len(b)+1):\n        d[0][i] = i\n    \n    \n    for i in range(1, len(a)+1):\n        for j in range(1, len(b)+1):\n            dist_horz = d[i][j-1] + 1\n            dist_vert = d[i-1][j] + 1\n            #check if the characters match \n            #if they match the edit distance wont increase\n            if a[i-1] == b[j-1]:\n                dist_diag = d[i-1][j-1]\n            #otherwise the edit distance increase    \n            else:\n                dist_diag = d[i-1][j-1] + 1\n            \n            #take the minimum value \n            d[i][j] = min(dist_horz, dist_vert, dist_diag)\n        \n    return d[-1][-1]\n\n__a = 'ATTGCTATTC'\n__b = 'ATATGACTTAC'\nedit_distance_dynamic(__a, __b)")


# In[34]:


#GLOBAL ALIGNMENT


# In[164]:


alphabet = ['A', 'C', 'G', 'T']
#penalty matrix
score = [[0, 4, 2, 4, 8], 
         [4, 0, 4, 2, 8], 
         [2, 4, 0, 4, 8], 
         [4, 2, 4, 0, 8], 
         [8, 8, 8, 8, 8]]


# In[196]:


#TODO: todavia no funciona correctamente falta mirar a detalle
def global_alignment(a, b):
    '''return the edit distance between the strings'''
    d = []
    #initialize the matrix with 0
    for i in range(len(a)+1):
        d.append([0] * (len(b)+1))
        
    #the epsilon column and row is fill with the edit distance between the empty string 
    #and the position i of the string a and b.
    #for i in range(1, len(a)+1):
        #d[i][0] = d[i-1][0] = score[alphabet.index(a[i-1])][-1]    
    #for i in range(1, len(b)+1):
        #d[0][i] = d[0][i-1] = score[-1][alphabet.index(b[i-1])]    
    
    
    for i in range(1, len(a)+1):
        for j in range(1, len(b)+1):
            dist_horz = d[i][j-1] + score[-1][alphabet.index(b[j-1])]
            dist_vert = d[i-1][j] + score[alphabet.index(a[i-1])][-1]
            #check if the characters match 
            #if they match the edit distance wont increase
            if a[i-1] == b[j-1]:
                dist_diag = d[i-1][j-1]
            #otherwise the edit distance increase    
            else:
                dist_diag = d[i-1][j-1] + score[alphabet.index(a[i-1])][alphabet.index(b[j-1])]
            
            #take the minimum value 
            d[i][j] = min(dist_horz, dist_vert, dist_diag)
            
    return d[-1][-1]

__u = 'ATATC'
__v = 'ATATT'

print(global_alignment(__u, __v))


# In[112]:


#GRAPHICS


# In[130]:


import matplotlib.pyplot as plt

# Datos
labels = ['Secuencias encontradas', 'Secuencias no encontradas']
sizes = [percentage, 100 - percentage]

# Crear el gráfico circular
fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
ax1.axis('equal')

# Añadir título al gráfico
plt.title('Porcentaje de secuencias encontradas respecto a k-mer BloomFilter')

# Mostrar el gráfico
plt.show()


# In[131]:


def phred33_to_q(qualities):
    """Turn Phred+33 ASCII-encoded quality into Q"""
    return ord(qualities) - 33

def create_histogram_qualities(qualities):
    hist = [0] * 45
    for qual in qualities:      
        for phred in qual:
            q = phred33_to_q(phred)
            hist[q] += 1
    return hist

h = create_histogram_qualities(qualities)


# In[132]:


#create a histogram of the qualities values
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.bar(range(len(h)), h)
plt.show()

