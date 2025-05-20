import pandas as pd
from ydata_profiling import ProfileReport
import sweetviz as sv

initial_cleaning_output_loc = '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial Cleansing/'
breeze_loc = '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Breeze/'
report_output_loc = '/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Python/'
def save_df(df, loc):
    df.to_csv(fr'{loc}.csv')


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

def create_summary_report(df, output_loc, title):
    report = ProfileReport(
            df,
            title=title)
    report.to_file(output_loc)
    #sv_report = sv.analyze(df)
    #sv_report.show_html(output_loc)

def get_summary(df):
    with pd.option_context('display.max_rows', 10, 'display.max_columns', None):
        print(df)
    print(df.dtypes)
    print(df.describe)

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


#Karolinska institute
karolinska_institute = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Karolinska/DSS_mastersheet_20240117_KI-1.xlsx')
get_summary(karolinska_institute)
save_df(karolinska_institute, initial_cleaning_output_loc + 'karolinska_institute_dss2')

#Enserink-lab 
enserink_lab_dose_response = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/Dose response data and Hill curve fits_2023-07.csv')
get_summary(enserink_lab_dose_response)
save_df(enserink_lab_dose_response, initial_cleaning_output_loc + 'enserink_lab_dose_response')
#worked :)
#create_summary_report(enserink_lab_dose_response, '/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Python/enserink_lab_dose_response.html', 'enserink_lab_dose_response')
print('#################################')
print('ENSERINK BREEe')
make_breeze_input_file(enserink_lab_dose_response, '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Breeze/enserink_breeze_input.csv', column_map={'Patient.num': 'PLATE', 'drug':'DRUG_NAME', 'dose': 'CONCENTRATION', 'pred_err.breeze':'PERCENT_INHIBITION', 'Patient.ID': 'SCREEN_NAME'})
print('#################################')
enserink_lab_drug_sensitivity = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/Drug sensitivity metrics_2023-07.csv')
get_summary(enserink_lab_drug_sensitivity)
save_df(enserink_lab_drug_sensitivity, initial_cleaning_output_loc + 'enserink_lab_drug_sensitivity')
#Worked :)
#create_summary_report(enserink_lab_drug_sensitivity, '/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Python/enserink_lab_drug_sensitivity.html', 'enserink_lab_drug_sensitivity')


enserink_lab_hill_drug_sensitivity = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/Drug sensitivity metrics_Hill rAUCs_2023-07.csv')
get_summary(enserink_lab_hill_drug_sensitivity)
save_df(enserink_lab_hill_drug_sensitivity, initial_cleaning_output_loc + 'enserink_lab_hill_drug_sensitivity')
#Worked :)
#create_summary_report(enserink_lab_hill_drug_sensitivity, '/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Python/enserink_lab_hill_drug_sensitivity.html', 'enserink_lab_hill_drug_sensitivity')


enserink_lab_drug_information = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/Selleck_Drug information.xlsx')
get_summary(enserink_lab_drug_information)
save_df(enserink_lab_drug_information, initial_cleaning_output_loc + 'enserink_lab_drug_information')
#Worked :)
#create_summary_report(enserink_lab_drug_information, '/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Python/enserink_lab_drug_information.html', 'enserink_lab_drug_information')



#BeatAML
beat_aml_calc = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/BeatAML/beataml_probit_curve_fits_v4_dbgap.txt')
get_summary(beat_aml_calc)
save_df(beat_aml_calc, initial_cleaning_output_loc + 'beat_aml_calc')
#not
#create_summary_report(beat_aml_calc.drop('status', axis='columns'), '/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Python/beat_aml_calc.html', 'beat_aml_calc')

beat_aml_waves = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/BeatAML/beataml_waves1to4_counts_dbgap.txt')
get_summary(beat_aml_waves)
save_df(beat_aml_waves, initial_cleaning_output_loc + 'beat_aml_waves')
#not
#create_summary_report(beat_aml_waves, report_output_loc + 'beat_aml_waves.html', 'beat_aml_waves')

beat_aml_waves_norm = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/BeatAML/beataml_waves1to4_norm_exp_dbgap.txt')
get_summary(beat_aml_waves_norm)
save_df(beat_aml_waves_norm, initial_cleaning_output_loc + 'beat_aml_waves_norm')
#not
#create_summary_report(beat_aml_waves_norm, report_output_loc + 'beat_aml_waves_norm.html', 'beat_aml_waves_norm') 

beat_aml_mutations = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/BeatAML/beataml_wes_wv1to4_mutations_dbgap.txt')
get_summary(beat_aml_mutations)
save_df(beat_aml_mutations, initial_cleaning_output_loc + 'beat_aml_mutations')
#Worked :)
#create_summary_report(beat_aml_mutations, report_output_loc + 'beat_aml_mutations.html', 'beat_aml_mutations')

beat_aml_inhibitor = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/BeatAML/beataml_wv1to4_raw_inhibitor_v4_dbgap.txt')
get_summary(beat_aml_inhibitor)
save_df(beat_aml_inhibitor, initial_cleaning_output_loc + 'beat_aml_inhibitor')
print(beat_aml_inhibitor)
print('#################################')
print('ENSERINK BREEe')
make_breeze_input_file(beat_aml_inhibitor, '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Breeze/beat_aml_breeze_input.csv', column_map={'run_index':'PLATE', 'inhibitor':'DRUG_NAME', 'well_concentration': 'CONCENTRATION', 'normalized_viability':'PERCENT_INHIBITION', 'dbgap_subject_id': 'SCREEN_NAME'})
print('#################################')
#Worked :)
#create_summary_report(beat_aml_inhibitor, report_output_loc + 'beat_aml_inhibitor.html', 'beat_aml_inhibitor')

beat_aml_sample_mapping = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/BeatAML/beataml_waves1to4_sample_mapping.xlsx')
get_summary(beat_aml_sample_mapping)
save_df(beat_aml_sample_mapping, initial_cleaning_output_loc + 'beat_aml_sample_mapping')
#Worked :)
#create_summary_report(beat_aml_sample_mapping, report_output_loc + 'beat_aml_sample_mapping.html', 'beat_aml_sample_mapping')

beat_aml_clinical = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/BeatAML/beataml_wv1to4_clinical.xlsx')
get_summary(beat_aml_clinical)
save_df(beat_aml_clinical, initial_cleaning_output_loc + 'beat_aml_clinical')
#Worked :)
#create_summary_report(beat_aml_clinical, report_output_loc + 'beat_aml_clinical.html', 'beat_aml_clinical')

beat_aml_vg_cts = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/BeatAML/beataml2_manuscript_vg_cts.xlsx')
get_summary(beat_aml_vg_cts)
save_df(beat_aml_vg_cts, initial_cleaning_output_loc + 'beat_aml_vg_cts')
#Worked :)
#create_summary_report(beat_aml_vg_cts, report_output_loc + 'beat_aml_vg_cts.html', 'beat_aml_vg_cts')


#FIMM
fimm_sample_annotation = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_0_Common_sample_annotation_252S.xlsx')
get_summary(fimm_sample_annotation)
save_df(fimm_sample_annotation, initial_cleaning_output_loc + 'fimm_sample_annotation')
#Worked :)
#create_summary_report(fimm_sample_annotation, '/Users/katarinawilloch/Desktop/UiO/Project 1/Reports/Python/Python/fimm_sample_annotation', 'fimm_sample_annotation')


fimm_clinical_summary = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_1_1_Clinical_summary_186_Patients.xlsx')
get_summary(fimm_clinical_summary)
save_df(fimm_clinical_summary , initial_cleaning_output_loc + 'fimm_clinical_summary')
#Worked :)
#create_summary_report(fimm_clinical_summary, report_output_loc + 'fimm_clinical_summary', 'fimm_clinical_summary') 

fimm_drug_library = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_2_Drug_library_515D.xlsx')
get_summary(fimm_drug_library)
save_df(fimm_drug_library , initial_cleaning_output_loc + 'fimm_drug_library')
#Worked :)
#create_summary_report(fimm_drug_library, report_output_loc + 'fimm_drug_library', 'fimm_drug_library') 

fimm_drug_response = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_3_Drug_response_sDSS_164S_17Healthy.xlsx')
get_summary(fimm_drug_response)
save_df(fimm_drug_response , initial_cleaning_output_loc + 'fimm_drug_response')
#Worked :)
#create_summary_report(fimm_drug_response, report_output_loc + 'fimm_drug_response', 'fimm_drug_response') 

fimm_assay_dets = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_4_DSRT_assay_details_164S_DM.xlsx')
get_summary(fimm_assay_dets)
save_df(fimm_assay_dets , initial_cleaning_output_loc + 'fimm_assay_dets')
#Worked :)
#create_summary_report(fimm_assay_dets, report_output_loc + 'fimm_assay_dets', 'fimm_assay_dets') 

fimm_mutation_variant_allele_frequency = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_5_Mutation_Variant_Allele_Frequency_225S_340G.csv')
get_summary(fimm_mutation_variant_allele_frequency)
save_df(fimm_mutation_variant_allele_frequency , initial_cleaning_output_loc + 'fimm_mutation_variant_allele_frequency')
#Worked :)
#create_summary_report(fimm_mutation_variant_allele_frequency, report_output_loc + 'fimm_mutation_variant_allele_frequency', 'fimm_mutation_variant_allele_frequency') 

fimm_binary_mutations = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_6_Binary_mutation_225S_57G.xlsx')
get_summary(fimm_binary_mutations)
save_df(fimm_binary_mutations , initial_cleaning_output_loc + 'fimm_binary_mutations')
#Worked :)
#create_summary_report(fimm_binary_mutations, report_output_loc + 'fimm_binary_mutations', 'fimm_binary_mutationa´s') 

fimm_rna_seq_healthy = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_7_RNA_seq_CPM_163S_4Healthy.csv')
get_summary(fimm_rna_seq_healthy)
save_df(fimm_rna_seq_healthy , initial_cleaning_output_loc + 'fimm_rna_seq_healthy')
#Worked :)
#create_summary_report(fimm_rna_seq_healthy, report_output_loc + 'fimm_rna_seq_healthy', 'fimm_rna_seq_healthy') 

fimm_rna_raw_reads = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv')
get_summary(fimm_rna_raw_reads)
save_df(fimm_rna_raw_reads , initial_cleaning_output_loc + 'fimm_rna_raw_reads')
#Worked :)
#create_summary_report(fimm_rna_raw_reads, report_output_loc + 'fimm_rna_raw_reads', 'fimm_rna_raw_reads') 

fimm_seq_sample_annotation = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/File_9_RNA-seq_sample_annotation_167S.csv')
get_summary(fimm_seq_sample_annotation)
save_df(fimm_seq_sample_annotation , initial_cleaning_output_loc + 'fimm_seq_sample_annotation')
#Worked :)
#create_summary_report(fimm_seq_sample_annotation, report_output_loc + 'fimm_seq_sample_annotation', 'fimm_seq_sample_annotation') 

make_breeze_input_file(beat_aml_inhibitor, '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Breeze/beat_aml_breeze_input.csv', column_map={'run_index':'PLATE', 'inhibitor':'DRUG_NAME', 'well_concentration': 'CONCENTRATION', 'normalized_viability':'PERCENT_INHIBITION', 'dbgap_subject_id': 'SCREEN_NAME'})


