"""
Script: drug_names.py

Overview:
This script gets the most common drug names in pubchem found for each given 
drug name for the different studies.
It imports drug name information from each study and runs 
function get_drug_synonyms to get the most common drug names.
Finnaly an upset plot and venn diagram of the number of overlapping drug names 
between the 4 studies is visualized.

Main Steps:
- Load data
- Pubchem API to get most common drug names
- Visualize overlap in drug names between studies

Author: Katarina Willoch
Date: 2025-06-01
"""

#Importing libraries
import pubchempy as pcp
import pandas as pd
from upsetplot import generate_counts
from upsetplot import UpSet
from upsetplot import plot
from matplotlib import pyplot as plt
from upsetplot import from_contents
from matplotlib_venn import venn2
from matplotlib_venn import venn3
from venny4py.venny4py import *

#Reading files as dataframes
initial_cleaning_output_loc = '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/'
df_enserink_drugs = pd.read_csv(fr'{initial_cleaning_output_loc}enserink_lab_drug_information.csv')
beat_aml_inhibitor  = pd.read_csv(fr'{initial_cleaning_output_loc}beat_aml_inhibitor.csv')
fimm_drug_library = pd.read_csv(fr'{initial_cleaning_output_loc}fimm_drug_library.csv')
karolinska_drug_library = pd.read_csv(fr'{initial_cleaning_output_loc}karolinska_institute_dss2.csv')

#Function to get the most common drug name in pubchem using the pubchem API
def get_drug_synonyms(df, col_name, output_loc):
    """
    Get the most common drug names from pubchem based on the ones privided by the stidues. 

    Args:
        series (pd.dataframe): The dataframe to get the original drug names from.
        col_name (String): The originally provided drug name.
        output_loc (String): Location to save the results. 

    Returns:
        pd.dataframe: The most common drug names. 
    """
    output_df = []
    for drug_name in df[col_name].unique():
        synonyms = pcp.get_synonyms(drug_name, 'name')
        if not synonyms:
            top_synonym = drug_name
            print(top_synonym)
        else:
            top_synonym = synonyms[0].get('Synonym')[0]
            print(top_synonym)
        output_df.append((drug_name,top_synonym))

    final_df = pd.DataFrame(output_df, columns=['org_drug_name', 'pubchem_drug_name'])
    final_df.to_csv(output_loc, index=False)
    return final_df

#Use the function to get the most common drug name for each dataset
df_enserink_pubchem_drugs = get_drug_synonyms(df_enserink_drugs, 'Item Name', fr'{initial_cleaning_output_loc}enserink_lab_drug_information_pubchem.csv') #'/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/enserink_lab_drug_information_pubchem.csv'
df_enserink_pubchem_drugs = df_enserink_pubchem_drugs.rename(columns={'org_drug_name': 'enserink_drug_name'})

df_beat_aml_pubchem_drugs = get_drug_synonyms(beat_aml_inhibitor, 'inhibitor', fr'{initial_cleaning_output_loc}beataml_drug_information_pubchem.csv') #'/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/beataml_drug_information_pubchem.csv'
df_beat_aml_pubchem_drugs= df_beat_aml_pubchem_drugs.rename(columns={'org_drug_name': 'beat_aml_drug_name'})

df_fimm_pubchem_drugs = get_drug_synonyms(fimm_drug_library, 'Preferred name', fr'{initial_cleaning_output_loc}fimm_drug_library_pubchem.csv') #'/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/fimm_drug_library_pubchem.csv'
df_fimm_pubchem_drugs = df_fimm_pubchem_drugs.rename(columns={'org_drug_name': 'fimm_drug_name'})

df_karolinska_pubchem_drugs = get_drug_synonyms(karolinska_drug_library, 'drug name', fr'{initial_cleaning_output_loc}karolinska_drug_library_pubchem.csv') #'/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/karolinska_drug_library_pubchem.csv'
df_karolinska_pubchem_drugs = df_karolinska_pubchem_drugs.rename(columns={'org_drug_name': 'karolinska_drug_name'})

#merging result dfs
common_drug_names = pd.merge(df_enserink_pubchem_drugs, df_beat_aml_pubchem_drugs, on = 'pubchem_drug_name', how='inner')
common_drug_names = pd.merge(common_drug_names,df_fimm_pubchem_drugs, on = 'pubchem_drug_name', how='inner')
common_drug_names = pd.merge(common_drug_names,df_karolinska_pubchem_drugs, on = 'pubchem_drug_name', how='inner')
print(len(common_drug_names['pubchem_drug_name']))
print(common_drug_names)
#saving merged dfs to csv
common_drug_names.to_csv(fr'{initial_cleaning_output_loc}common_drugnames_pubchem.csv', index=False)

#Creating dict of the dfs
df = from_contents({'Oslo': df_enserink_pubchem_drugs['pubchem_drug_name'], 'Beat AML': df_beat_aml_pubchem_drugs['pubchem_drug_name'], 'Helsinki': df_fimm_pubchem_drugs['pubchem_drug_name'], 'Karolinska': df_karolinska_pubchem_drugs['pubchem_drug_name']})
print(df)
#Creating an upset plot of the number of overlapping common drug names between the different studies
upset = UpSet(df, show_counts="%d", show_percentages=True)
upset.style_subsets(present=["Beat AML", "Oslo", "Helsinki", "Karolinska"], facecolor="blue")
upset.plot()
plt.suptitle("With counts and % of common drug names from pubchem")
plt.show()
#plt.savefig('/Users/katarinawilloch/Desktop/UiO/Project 1/Figures/overlapping_drugs_pubchem.png')

#Dict of sets
sets = {
    'Oslo': set(df_enserink_pubchem_drugs['pubchem_drug_name'].to_list()),
    'Beat AML': set(df_beat_aml_pubchem_drugs['pubchem_drug_name'].to_list()),
    'Helsinki': set(df_fimm_pubchem_drugs['pubchem_drug_name'].to_list()), 
    'Karolinska': set(df_karolinska_pubchem_drugs['pubchem_drug_name'].to_list())}
    
#Venn diagram of the number of overlapping common drug names
venny4py(sets=sets)