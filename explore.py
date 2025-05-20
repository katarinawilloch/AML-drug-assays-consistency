import pandas as pd 
import matplotlib.pyplot as plt

common_drugs = pd.read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/common_drugnames_pubchem.csv')
beat_aml = pd.read_csv('~/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_calc.csv')
beat_aml = beat_aml.loc[:, ~beat_aml.columns.str.contains('^Unnamed')]
beat_aml_merged = pd.merge(common_drugs, beat_aml, left_on='beat_aml_drug_name', right_on='inhibitor', how='inner')

# Counting the number of unique values in col1 for each unique value in col2
result = beat_aml_merged.groupby('min_conc')['pubchem_drug_name'].nunique().reset_index()

print(result)

plt.bar(result['pubchem_drug_name'], result['min_conc'])
plt.xlabel('Nr drugs')
plt.ylabel('Min Concentration')
plt.title('Nr of drugs that has each min concentration')
plt.show()

# Counting the number of unique values in col1 for each unique value in col2
result_max = beat_aml_merged.groupby('max_conc')['pubchem_drug_name'].nunique().reset_index()

print(result)

plt.bar(result_max['pubchem_drug_name'], result_max['max_conc'])
plt.xlabel('Nr drugs')
plt.ylabel('Max Concentration')
plt.title('Nr of drugs that has each max concentration')
plt.show()