import os
import pandas as pd
from Bio import SeqIO

def get_gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

directory = 'analysis/unclassified_out/'
sequences = []

for filename in os.listdir(directory):
    if filename.endswith('_u_1.fq'):
        for record in SeqIO.parse(os.path.join(directory, filename), 'fastq'):
            sequences.append(str(record.seq))

df = pd.DataFrame({'sequence': sequences})
df['count'] = df.groupby('sequence')['sequence'].transform('count')
df = df.drop_duplicates()
df['GC%'] = df['sequence'].apply(get_gc_content)
df['length'] = df['sequence'].apply(len)
df = df[['sequence', 'count', 'GC%', 'length']]

df.to_csv('data/unclassified_stats.csv', index=False, sep=' ')

# get number of classified seqs
def get_sequences_(directory):
    sequences = []
    for filename in os.listdir(directory):
        if filename.endswith('_c_1.fq'):
            for record in SeqIO.parse(os.path.join(directory, filename), 'fastq'):
                sequences.append(str(record.seq))
    return len(sequences)