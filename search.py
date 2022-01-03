from pandas import DataFrame, read_csv, Series

ref = read_csv("./meth2.csv")
sites = read_csv("./coeffs.sample.txt", sep='\t')

ref = ref.sort_values('Name')
print('Values FUCKING SORTED AGAIN, YOU STUPID MORON! I FUCKED WITH THIS SHIT FOR 2 DAYS, DAMN IT!')
nref = len(ref)

output = DataFrame(
   columns = ['Name', 'Chromosome', 'Gene', 'Island', 'Feature', 'Rho', 'P']) 
print('All preparations done. Beginning')

for cg in sites.iterrows():
    cgname = int(cg[1]['Name'][2:])
    ind = ref['Name'].searchsorted(cgname)
    if ind != 0 and ind <= nref:
        print(cg[1]['Name'], ref['UCSC_RefGene_Name'].iloc[ind])
        output = output.append(
            {'Name'      : cg[1]['Name'],
            'Chromosome' : ref['Chromosome_36'].iloc[ind],
            'Gene'       : ref['UCSC_RefGene_Name'].iloc[ind],
            'Island'     : ref['UCSC_CpG_Islands_Name'].iloc[ind],
            'Feature'    : ref['Regulatory_Feature_Group'].iloc[ind],
            'Rho'        : cg[1]['Rho'],
            'P'          : cg[1]['p']}, ignore_index=True)

output.to_csv("genes_and_sites.csv", na_rep='NA')
