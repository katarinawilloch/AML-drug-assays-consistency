import pandas as pd
from ydata_profiling import ProfileReport
import sweetviz as sv
import ast
print("Pandas version:", pd.__version__)

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

fimm_dose_response_raw = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/cm125_ctg_fo4b_fo5a_FIMM_125_AML_patients_oslo_14052024.csv')
print(fimm_dose_response_raw)


fimm_drug_concentrations = pd.read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/170515 FO5A Annotations_JS-1.xlsx', sheet_name='FO5A Annotations')
print(fimm_drug_concentrations)

fimm_drug_concentrations = fimm_drug_concentrations[fimm_drug_concentrations['Preferred name'] != 'Cytarabine/Idarubicin']
fimm_drug_concentrations = fimm_drug_concentrations[fimm_drug_concentrations['Preferred name'] != 'empty']
fimm_drug_concentrations['Final Conc'] = fimm_drug_concentrations['Final Conc'].astype('float')
print('I AM HERE')
print(fimm_drug_concentrations[fimm_drug_concentrations['Preferred name'] == 'Cytarabine/Idarubicin'])
print(fimm_drug_concentrations.sort_values('Preferred name'))

df = pd.merge(fimm_dose_response_raw, fimm_drug_concentrations, left_on = "drug",right_on="Preferred name",how="inner")
print('PRINTING DF')
print(fimm_dose_response_raw)
print(df['sample'].nunique())



df['doseResponses'] = df['doseResponses'].apply(lambda x: [float(i) for i in ast.literal_eval(x)])
df['doseResponses'] = df.apply(lambda row: row['doseResponses'][int(row['Dilution'])-1], axis=1)
print(df['sample'].nunique())

#df['doseResponses'] = df.apply(lambda row: row['doseResponses'][5 - row['Final Conc']], axis=1)

# print("\nTypes in 'values' column:")
# df['doseResponses'] = df['doseResponses'].apply(lambda x: [float(i) for i in ast.literal_eval(x)])
# print(df['doseResponses'].apply(type))
# df_exploded = df.explode('doseResponses').reset_index(drop=True)



print("\nDataFrame after exploding 'values' column:")
#print(df)
save_df(df, fr'{initial_cleaning_output_loc}fimm_raw_dose_responses_test.csv')
stop
fimm_drug_raw_reponse = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_raw_dose_responses.csv.csv')
print(fimm_drug_raw_reponse)
#Raw reponse find Cytarabine/Idarubicin
print(fimm_drug_raw_reponse[fimm_drug_raw_reponse['sample'] == 'TG100-115_FH_7532_29042019_9999_BM'])

#

#breeze = df_exploded[['drug', 'concentration', 'sample', 'doseResponses']].rename(columns={'drug': 'DRUG_NAME', 'concentration': 'CONCENTRATION','sample':'SCREEN_NAME', 'doseResponses': 'PERCENT_INHIBITION'})
#breeze['PERCENT_INHIBITION'] = breeze['PERCENT_INHIBITION'].astype(float)
#print(breeze.dtypes)

#breeze.to_csv(fr'{breeze_loc}fimm_input_breexe.csv', index=False)