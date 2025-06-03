import pandas as pd
from upsetplot import generate_counts
from upsetplot import UpSet
from upsetplot import plot
from matplotlib import pyplot as plt
from upsetplot import from_contents
from upsetplot import from_memberships
from matplotlib_venn import venn2
from matplotlib_venn import venn3

df_enserink_drugs = pd.read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/enserink_lab_drug_information_pubchem.csv', index_col=False)
df_enserink_drugs = df_enserink_drugs.rename(columns={'org_drug_name': 'enserink_drug_name'})
print(df_enserink_drugs)
beat_aml_inhibitor = pd.read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/beataml_drug_information_pubchem.csv', index_col=False)
beat_aml_inhibitor = beat_aml_inhibitor.rename(columns={'org_drug_name': 'beat_aml_drug_name'})
fimm_drug_library = pd.read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_drug_library_pubchem.csv', index_col=False)
fimm_drug_library = fimm_drug_library.rename(columns={'org_drug_name': 'fimm_drug_name'})
karolinska_drug_library = pd.read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_drug_library_pubchem.csv', index_col=False)
karolinska_drug_library = karolinska_drug_library.rename(columns={'org_drug_name': 'karolinska_drug_name'})

common_drug_names = pd.merge(df_enserink_drugs, beat_aml_inhibitor, on = 'pubchem_drug_name', how='inner')
common_drug_names = pd.merge(common_drug_names,fimm_drug_library, on = 'pubchem_drug_name', how='inner')
common_drug_names = pd.merge(common_drug_names,karolinska_drug_library, on = 'pubchem_drug_name', how='inner')
print(len(common_drug_names['pubchem_drug_name']))
print(common_drug_names)
common_drug_names.to_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv', index=False)

df = from_contents({'Oslo': df_enserink_drugs['pubchem_drug_name'], 'BeatAML': beat_aml_inhibitor['pubchem_drug_name'], 'Helsinki': fimm_drug_library['pubchem_drug_name'], 'Karolinska': karolinska_drug_library['pubchem_drug_name']})
print(df)


upset = UpSet(df, show_counts="%d", show_percentages=True, subset_size="count")
upset.style_subsets(present=["BeatAML", "Oslo", "Helsinki", "Karolinska"], facecolor="blue")
upset.plot()
#plt.suptitle("With counts and % of common drug names from pubchem")
plt.show()
plt.savefig('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/V3/overlapping_drugs_pubchem.png')
	

 
""" plt.figure(figsize=(10,10))
plt.title("Venn Diagram For Number Of Overlapping Drug Names")
 
venn3([set(df_enserink_drugs['pubchem_drug_name'].to_list()), set(beat_aml_inhibitor['pubchem_drug_name'].to_list()), 
       set(fimm_drug_library['pubchem_drug_name'].to_list())],
       set_labels=('Enserink','BeatAML', 'FIMM')
     )
 
plt.show()


from venny4py.venny4py import *

#dict of sets
sets = {
    'Enserink': set(df_enserink_drugs['pubchem_drug_name'].to_list()),
    'BeatAML': set(beat_aml_inhibitor['pubchem_drug_name'].to_list()),
    'FIMM': set(fimm_drug_library['pubchem_drug_name'].to_list()), 
    'Karolinska': set(karolinska_drug_library['pubchem_drug_name'].to_list())}

colors = ["blue", "red", "green", "yellow"]
print(colors[1])
venny4py(sets=sets)

#venn.set_region_color('100', color='blue')   # Only in Set A
#venn.get_patch_by_id('BeatAML').set_color('#8dd3c7')  # Only in Set 2
#venn.get_patch_by_id('Enserink').set_color('#fdb462') # Intersection of Set 1 and Set 2
#venn.get_patch_by_id('Karolinska').set_color('#80b1d3') # Intersection of Set 1 and Set 2


plt.show() """
