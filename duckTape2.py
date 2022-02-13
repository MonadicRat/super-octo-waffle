import pandas as pd

data = pd.read_csv("./genes_and_sites.csv", sep = ',')

for i in range(1, len(data)):
    if pd.isna(data.iloc[i]['Gene']):
        pass
    else:
        data.at[i, 'Gene'] = data.at[i, 'Gene'].split(";")[0]
    print(i/470043, data.iloc[i]['Gene'])
data.to_csv("./genes_and_sites.noreps.csv", na_rep='NA')
