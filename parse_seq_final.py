import numpy as np
import pandas as pd
from time import time
from Bio import SeqIO
from Bio.Seq import Seq

sequences = [record.seq for index, record in enumerate(SeqIO.parse("ferrochelatase.fasta", "fasta"))]

print """
This is an implementation of the Needleman-Wunsch algorithm made by Nicolas A. Gort Freitas
as a Final Project for the course CS110: Solving Problems with Algorithms taught by
Prof. Philip Sterne in the Fall Semester of 2016.

This script finds an prints the optimal alignment of two biological sequences using dynamic
programming. Due to the nature and rigor of the algorithm, it has a considerable complexity
that makes it unsuitable for the practical needs of a bioinformatician. Nonetheless, it is
still a good didactic example of a problem with an optimal sub-structure and overlapping
subproblems, and of the capabilities of dynamic programming overall.

To show the performance of this algorithm in practical time, this script will align only
the first 50 nucleotides of the sample sequences. The sample dataset corresponds to the
sequence of the enzyme Ferrochelatase in three species: Homo sapiens (HsFECH), Arabidopsis
thaliana (AtFC-I) and Escherichia coli (EcHEMH) respectively. This enzyme catalyzes the
last step in the biosynthesis of the Heme group, and these pathways are very similar even
across biological kingdoms.\n"""

for index, record in enumerate(SeqIO.parse("ferrochelatase.fasta", "fasta")):
    print("index %i, ID = %s, length %i, with %i features"  % (index+1, record.id, len(record.seq), len(record.features)))

seq1 = ''
seq2 = ''
choice1 = 0
choice2 = 0

while (choice1 not in [1,2,3]) or (choice2 not in [1,2,3]) or (choice1==choice2):
    print '\nPlease select a valid pair of sequences from the prior list to align.\n'    
    try:
        choice1 = int(raw_input('What is the first sequence you would like to align?\n'))
        choice2 = int(raw_input('What is the second sequence you would like to align?\n'))
    except:
        pass
    print '\n'

a = sequences[choice1-1][:50]
b = sequences[choice2-1][:50]

start = time()

score_matrix = np.zeros((len(a)+1,len(b)+2))
score_matrix = pd.DataFrame(score_matrix).set_index(0)
score_matrix.columns = ['gap']+[x+str(i+1) for i,x in enumerate(list(b))]
score_matrix.index = ['gap']+[x+str(i+1) for i,x in enumerate(list(a))]

traceback = np.zeros((len(a)+1,len(b)+2))
traceback = pd.DataFrame(traceback).set_index(0)
traceback.columns = ['gap']+[x+str(i+1) for i,x in enumerate(list(b))]
traceback.index = ['gap']+[x+str(i+1) for i,x in enumerate(list(a))]

def align(a,b):
    for i in range(1,score_matrix.shape[0]):
        for j in range(1,score_matrix.shape[1]):
            score(i,j)

    position = (len(a),len(b))

    alig1 = ''
    alig2 = ''

    while position != (0,0):
        if traceback.ix[position] == 'diag':
            alig2+=score_matrix.columns[position[1]][0]
            alig1+=score_matrix.index[position[0]][0]
            position = (position[0]-1,position[1]-1)
        elif traceback.ix[position] == 'left':
            alig2+=score_matrix.columns[position[1]][0]
            alig1+='-'
            position = (position[0],position[1]-1)
        elif traceback.ix[position] == 'up':
            alig2+='-'
            alig1+=score_matrix.index[position[0]][0]
            position = (position[0]-1,position[1])

    alig1 = alig1[::-1]
    alig2 = alig2[::-1]

    print alig1
    print alig2

def score(i,j):
    
    gap = -5
    miss = -1
    match = 1
    
    score_matrix.ix[0,] = [float(gap*x) for x in range(score_matrix.shape[1])]
    score_matrix.ix[:,0] = [float(gap*x) for x in range(score_matrix.shape[0])]

    left = score_matrix.ix[i,j-1]
    up = score_matrix.ix[i-1,j]
    diag = score_matrix.ix[i-1,j-1]
    val = score_matrix.ix[i,j]
    
    if score_matrix.index[i][0] == score_matrix.columns[j][0]:
        score_matrix.ix[i,j] = max(left+gap,up+gap,diag+match)
    else:
        score_matrix.ix[i,j] = max(left+gap,up+gap,diag+miss)
    
    if score_matrix.ix[i,j] == left+gap:
        traceback.ix[i,j] = 'left'
    elif score_matrix.ix[i,j] == up+gap:
        traceback.ix[i,j] = 'up'
    elif score_matrix.ix[i,j] == diag+miss:
        traceback.ix[i,j] = 'diag'
    elif score_matrix.ix[i,j] == diag+match:
        traceback.ix[i,j] = 'diag'
        


align(a,b)
        
        
print '\nSuccessfully aligned the sequences in '+str(time()-start)+' seconds.'


