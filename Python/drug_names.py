import pubchempy as pcp
import pandas as pd

df_enserink_drugs = pd.read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/enserink_lab_drug_information.csv')
beat_aml_inhibitor  = pd.read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/beat_aml_inhibitor.csv')
fimm_drug_library = pd.read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/fimm_drug_library.csv')
karolinska_drug_library = pd.read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/karolinska_institute_dss2.csv')

def get_drug_synonyms(df, col_name, output_loc):
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

#get_drug_synonyms(df_enserink_drugs, 'Item Name', '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/enserink_lab_drug_information_pubchem.csv')
#get_drug_synonyms(beat_aml_inhibitor, 'inhibitor', '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/beataml_drug_information_pubchem.csv')
#get_drug_synonyms(fimm_drug_library, 'Preferred name', '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/fimm_drug_library_pubchem.csv')
get_drug_synonyms(karolinska_drug_library, 'drug name', '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/karolinska_drug_library_pubchem.csv')

