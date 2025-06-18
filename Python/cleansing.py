"""
Script: cleansing.py

Overview:
This script reformats and gathers the original downloaded data from the publications 
and stores the cleansed data in a new output location.  
Different reformatting procedures are used for the different studies. 
Data downloaded from: 
Beat AML: https://biodev.github.io/BeatAML2/
Helsinki: https://zenodo.org/records/7274740
Oslo: https://github.com/Enserink-lab/DSCoxTools/tree/main
Karolisnka: recieved from collaborator Tom Erkers

Main Steps:
- Load data
- Reformatting Beat AML data
- Reformatting Helsinki data
- Reformatting Oslo data 
- Reformatting Karolinska data 

Author: Katarina Willoch
Date: 2025-06-01
"""

#Loading libraries
import pandas as pd
from ydata_profiling import ProfileReport
import sweetviz as sv
import ast
import os
import glob
import re 
import math
pd.set_option('display.min_rows', 30)
pd.set_option('display.max_rows', 100)

#File locations
##Input locations
original_downloaded_data_location = '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/'
karolinska_files = glob.glob('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Karolinska_DSRT/breeze_files_all/' + '*.xlsx')
##output locations 
initial_cleaning_output_loc = '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Second run/' #'/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/'
breeze_loc = '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Breeze/'


#Functions
##Saving dataframe as csv
def save_df(df, loc):
    df.to_csv(fr'{loc}.csv')

##Reading file based on input format
##Return dataframe
def get_df(input):
    type = input.split('.')[1]
    print(type)
    df = pd.DataFrame()
    if type == 'txt':
        df = pd.read_table(input, delimiter='\t', header=(0))
    elif type == 'csv': 
        df = pd.read_csv(input)
    elif type == 'xlsx':
        df = pd.read_excel(input)
    return df

def make_breeze_input_file(df, output_loc, column_map):
    print(df)
    df = df.rename(columns=column_map, errors="raise")
    df['SCREEN_NAME'] = df['SCREEN_NAME'].astype('str')
    print('rename',df)
    df = df[['DRUG_NAME', 'CONCENTRATION', 'SCREEN_NAME', 'PERCENT_INHIBITION']]
    #df['CONCENTRATION'] = df['CONCENTRATION'].astype('int64')
    print(df.dtypes)
    print(df)
    df.to_csv(output_loc, index=False)

#Beat AML 
##Creating copies of the original data and storing them in the common output location
beat_aml_calc = get_df(fr'{original_downloaded_data_location}BeatAML/beataml_probit_curve_fits_v4_dbgap.txt')
save_df(beat_aml_calc, initial_cleaning_output_loc + 'beat_aml_calc')

beat_aml_waves = get_df(fr'{original_downloaded_data_location}BeatAML/beataml_waves1to4_counts_dbgap.txt')
save_df(beat_aml_waves, initial_cleaning_output_loc + 'beat_aml_waves')

beat_aml_waves_norm = get_df(fr'{original_downloaded_data_location}BeatAML/beataml_waves1to4_norm_exp_dbgap.txt')
save_df(beat_aml_waves_norm, initial_cleaning_output_loc + 'beat_aml_waves_norm')

beat_aml_mutations = get_df(fr'{original_downloaded_data_location}BeatAML/beataml_wes_wv1to4_mutations_dbgap.txt')
save_df(beat_aml_mutations, initial_cleaning_output_loc + 'beat_aml_mutations')

beat_aml_inhibitor = get_df(fr'{original_downloaded_data_location}BeatAML/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt')
save_df(beat_aml_inhibitor, initial_cleaning_output_loc + 'beat_aml_inhibitor')
make_breeze_input_file(beat_aml_inhibitor, '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Breeze/beat_aml_breeze_input.csv', column_map={'run_index':'PLATE', 'inhibitor':'DRUG_NAME', 'well_concentration': 'CONCENTRATION', 'normalized_viability':'PERCENT_INHIBITION', 'dbgap_subject_id': 'SCREEN_NAME'})

beat_aml_sample_mapping = get_df(fr'{original_downloaded_data_location}BeatAML/beataml_waves1to4_sample_mapping.xlsx')
save_df(beat_aml_sample_mapping, initial_cleaning_output_loc + 'beat_aml_sample_mapping')

beat_aml_clinical = get_df(fr'{original_downloaded_data_location}BeatAML/beataml_wv1to4_clinical.xlsx')
save_df(beat_aml_clinical, initial_cleaning_output_loc + 'beat_aml_clinical')

beat_aml_vg_cts = get_df(fr'{original_downloaded_data_location}BeatAML/beataml2_manuscript_vg_cts.xlsx')
save_df(beat_aml_vg_cts, initial_cleaning_output_loc + 'beat_aml_vg_cts')

#Helsinki
##Creating copies of the original data and storing them in the common output location
fimm_dose_response_raw = get_df(fr'{original_downloaded_data_location}Functional_Precision_Medicine_Tumor_Board_AML/cm125_ctg_fo4b_fo5a_FIMM_125_AML_patients_oslo_14052024.csv')
fimm_drug_concentrations = pd.read_excel(fr'{original_downloaded_data_location}Functional_Precision_Medicine_Tumor_Board_AML/170515 FO5A Annotations_JS-1.xlsx', sheet_name='FO5A Annotations')

##Remove drug names 'empty' and 'Cytarabine/Idarubicin'
fimm_drug_concentrations = fimm_drug_concentrations[fimm_drug_concentrations['Preferred name'] != 'Cytarabine/Idarubicin']
fimm_drug_concentrations = fimm_drug_concentrations[fimm_drug_concentrations['Preferred name'] != 'empty']
##Make concentration type float
fimm_drug_concentrations['Final Conc'] = fimm_drug_concentrations['Final Conc'].astype('float')
#Merge raw dose responses and concentrations
fimm_df = pd.merge(fimm_dose_response_raw, fimm_drug_concentrations, left_on = "drug",right_on="Preferred name",how="inner")

fimm_df['doseResponses'] = fimm_df['doseResponses'].apply(lambda x: [float(i) for i in ast.literal_eval(x)])
fimm_df['doseResponses'] = fimm_df.apply(lambda row: row['doseResponses'][int(row['Dilution'])-1], axis=1)

print("\nDataFrame after exploding 'values' column:")
print(fimm_df)
save_df(fimm_df, fr'{initial_cleaning_output_loc}fimm_raw_dose_responses.csv')

fimm_sample_annotation = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_0_Common_sample_annotation_252S.xlsx')
save_df(fimm_sample_annotation, initial_cleaning_output_loc + 'fimm_sample_annotation')

fimm_clinical_summary = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_1_1_Clinical_summary_186_Patients.xlsx')
save_df(fimm_clinical_summary , initial_cleaning_output_loc + 'fimm_clinical_summary')

fimm_drug_library = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_2_Drug_library_515D.xlsx')
save_df(fimm_drug_library , initial_cleaning_output_loc + 'fimm_drug_library')

#fimm_drug_response = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_3_Drug_response_sDSS_164S_17Healthy.xlsx')
#save_df(fimm_drug_response , initial_cleaning_output_loc + 'fimm_drug_response')

fimm_assay_dets = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_4_DSRT_assay_details_164S_DM.xlsx')
save_df(fimm_assay_dets , initial_cleaning_output_loc + 'fimm_assay_dets')

fimm_mutation_variant_allele_frequency = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_5_Mutation_Variant_Allele_Frequency_225S_340G.csv')
save_df(fimm_mutation_variant_allele_frequency , initial_cleaning_output_loc + 'fimm_mutation_variant_allele_frequency')

fimm_binary_mutations = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_6_Binary_mutation_225S_57G.xlsx')
save_df(fimm_binary_mutations , initial_cleaning_output_loc + 'fimm_binary_mutations')

fimm_rna_seq_healthy = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_7_RNA_seq_CPM_163S_4Healthy.csv')
save_df(fimm_rna_seq_healthy , initial_cleaning_output_loc + 'fimm_rna_seq_healthy')

fimm_rna_raw_reads = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv')
save_df(fimm_rna_raw_reads , initial_cleaning_output_loc + 'fimm_rna_raw_reads')

fimm_seq_sample_annotation = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_9_RNA-seq_sample_annotation_167S.csv')
save_df(fimm_seq_sample_annotation , initial_cleaning_output_loc + 'fimm_seq_sample_annotation')


#Oslo
##Creating copies of the original data and storing them in the common output location
enserink_lab_dose_response = get_df(fr'{original_downloaded_data_location}Enserink-lab/Dose response data and Hill curve fits_2023-07.csv')
save_df(enserink_lab_dose_response, initial_cleaning_output_loc + 'enserink_lab_dose_response')

enserink_lab_drug_sensitivity = get_df(fr'{original_downloaded_data_location}Enserink-lab/Drug sensitivity metrics_2023-07.csv')
save_df(enserink_lab_drug_sensitivity, initial_cleaning_output_loc + 'enserink_lab_drug_sensitivity')

enserink_lab_hill_drug_sensitivity = get_df(fr'{original_downloaded_data_location}Enserink-lab/Drug sensitivity metrics_Hill rAUCs_2023-07.csv')
save_df(enserink_lab_hill_drug_sensitivity, initial_cleaning_output_loc + 'enserink_lab_hill_drug_sensitivity')

enserink_lab_drug_information = get_df(fr'{original_downloaded_data_location}Enserink-lab/Selleck_Drug information.xlsx')
save_df(enserink_lab_drug_information, initial_cleaning_output_loc + 'enserink_lab_drug_information')

enserink_lab_drug_information = get_df(fr'{original_downloaded_data_location}Enserink-lab/Selleck_Drug information.xlsx')
save_df(enserink_lab_drug_information, initial_cleaning_output_loc + 'enserink_lab_drug_information')

enserink_lab_sample_annotation = get_df(fr'{original_downloaded_data_location}Enserink-lab/enserink_saple_annotation.xlsx')
save_df(enserink_lab_sample_annotation, initial_cleaning_output_loc + 'enserink_lab_sample_annotation')



#Karoslinska
##Gathering the files
karolinska_reponse = pd.DataFrame()
for files in karolinska_files:
    print(files)
    basename = os.path.basename(files)
    parts = basename.split('_')
    match = re.search(r'[A-Z]+_\d+', basename)

    if len(parts) >= 1:
        part1 = parts[0]
    if match:
        p_id = match.group(0)
    
    print(f"Part 1: {part1}, Patient_id: {p_id}")
    df = pd.read_excel(fr'{original_downloaded_data_location}karolinska_DSRT/breeze_files_all/{basename}', sheet_name='EC50')
    df['sample'] = parts[0]
    df['Patient.num'] = p_id
    print(df)
    karolinska_reponse = pd.concat([karolinska_reponse, df], ignore_index=True)

karolinska_reponse.to_csv(fr'{initial_cleaning_output_loc}karolinska_normalised_reponse_breeze.csv')

#Making concentration of type float
karolinska_reponse['Max.Conc.tested'] = karolinska_reponse['Max.Conc.tested'].astype('float64')
karolinska_reponse['x'] = 10
# Melting the DataFrame to long format
karolinska_reponse_melted = karolinska_reponse.melt(id_vars=['Min.Conc.tested', 'Max.Conc.tested', 'Patient.num', 'DRUG_NAME', 'x', 'sample'], 
                    value_vars=['D1', 'D2', 'D3', 'D4', 'D5'], 
                    var_name='D', 
                    value_name='Normalised_reponse')
print('DF_MELTED')
print(karolinska_reponse_melted)

karolinska_reponse_melted['CONCENTRATION'] = karolinska_reponse_melted['Max.Conc.tested']
karolinska_reponse_melted['CONCENTRATION'][karolinska_reponse_melted['D'] == 'D1'] = karolinska_reponse_melted[karolinska_reponse_melted['D'] == 'D1'].apply(lambda row: row['Max.Conc.tested']/math.pow(row['x'],4), axis=1)
karolinska_reponse_melted['CONCENTRATION'][karolinska_reponse_melted['D'] == 'D2'] = karolinska_reponse_melted[karolinska_reponse_melted['D'] == 'D2'].apply(lambda row: row['Max.Conc.tested']/math.pow(row['x'],3), axis=1)
karolinska_reponse_melted['CONCENTRATION'][karolinska_reponse_melted['D'] == 'D3'] = karolinska_reponse_melted[karolinska_reponse_melted['D'] == 'D3'].apply(lambda row: row['Max.Conc.tested']/math.pow(row['x'],2), axis=1)
karolinska_reponse_melted['CONCENTRATION'][karolinska_reponse_melted['D'] == 'D4'] = karolinska_reponse_melted[karolinska_reponse_melted['D'] == 'D4'].apply(lambda row: row['Max.Conc.tested']/row['x'], axis=1)
karolinska_reponse_melted['CONCENTRATION'][karolinska_reponse_melted['D'] == 'D5'] = karolinska_reponse_melted[karolinska_reponse_melted['D'] == 'D5'].apply(lambda row: row['Max.Conc.tested'], axis=1)


karolinska_reponse_melted.to_csv(fr'{initial_cleaning_output_loc}karolinska_normalised_reponse_raw.csv')
karolinska_reponse_melted['patient_sample'] = karolinska_reponse_melted['Patient.num'] + '_'+ karolinska_reponse_melted['sample']
print(karolinska_reponse_melted[karolinska_reponse_melted['Patient.num'] == 'AML_001'].sort_values(by=['DRUG_NAME', 'sample']))

karolinska_sample_annotation = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Karolinska_DSRT/clinical_info_Tero_AML_samples_DSRT.csv')
save_df(karolinska_sample_annotation, initial_cleaning_output_loc + 'karolinska_sample_annotation')
