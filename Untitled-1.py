import os
os.chdir('/Users/zachr/Desktop')

# Open SNP_LIST.txt and read it into a df

import pandas as pd
df = pd.read_csv('SNP_LIST.txt', sep='\t')

#print first 5
print(df.head())

#name the col SNP

df.columns = ['SNP']

#Now create a text file, where each line is "dbsnp", then a tab, and then the original value of df.SNP, then staring next line

df['SNP'] = df['SNP'].astype(str)
df['SNP'] = 'dbsnp' + '\t' + df['SNP']

df.to_csv('SNP_LIST.txt', sep='\t', index=False, header=False)