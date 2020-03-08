## Code constructed with help from Fritz Pichardo

import sys
import os
import pandas as pd

refs = sys.argv[1]
query = sys.argv[2]
output = sys.argv[3]
output2 = sys.argv[4]


if os.path.exists(output):
	os.remove(output)
else:
	pass

df = pd.read_csv(refs, sep='\t')

query_df = pd.read_csv(query, header=None)
query_list = query_df[0].to_list()

mask = df['#gene_id'].isin(query_list)

filtered_df = df[mask]
filtered_df.to_csv(output, index=False, sep='\t')

#print(filtered_df)

df_symbols = filtered_df.copy()
df_symbols["sprotsymbols"] = filtered_df['sprot_Top_BLASTX_hit'].str.split('^').apply(lambda x: x[0])
mask2 = df_symbols['sprotsymbols'] != '.'
#print(mask2)
df_symbols = df_symbols[mask2]
df_symbols = df_symbols[['#gene_id', 'sprotsymbols']].groupby('#gene_id').first().reset_index()
df_symbols.to_csv(output2, index=False)

