import pandas as pd
import glob

pickles = sorted(glob.glob('out/*.pkl'))

dataframes = []
for filename in pickles:
    temp_df = pd.read_pickle(filename)
    if not temp_df.empty:
        dataframes.append(temp_df)

df = pd.concat(dataframes, ignore_index=True)
df.to_pickle('out/final.pkl')
