from pandas import DataFrame, read_csv, Series

ref = read_csv("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/meth2.csv")
sites = read_csv("/media/yanus/3ce7ab71-4bf7-4250-894c-508fb7a28109/Meth/coeffs.p2.txt", sep='\t')
print("Data loaded")

ref = ref.sort_values('Name')
print('Values sorted')
nref = len(ref)

output = DataFrame(
   columns = ['Name', 'Chromosome', 'Gene', 'Island', 'Feature', 'Rho', 'P']) 
print('All preparations done. Beginning')

for cg in sites.iterrows():
    print(cg[1]['Name'])
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
