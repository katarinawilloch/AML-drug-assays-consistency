import pandas as pd
import math
import matplotlib.pyplot as plt
pd.set_option('display.min_rows',30)
pd.set_option('display.max_rows', 100)
df = pd.read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_breese_output.csv')


df['Max.Conc.tested'] = df['Max.Conc.tested'].astype('float64')
print(df)
df['x'] = 10
print(df['x'])


# Melting the DataFrame to long format
df_melted = df.melt(id_vars=['Min.Conc.tested', 'Max.Conc.tested', 'Patient.num', 'DRUG_NAME', 'x', 'sample'], 
                    value_vars=['D1', 'D2', 'D3', 'D4', 'D5'], 
                    var_name='D', 
                    value_name='Normalised_reponse')
print('DF_MELTED')

print(df_melted)

df_melted['CONCENTRATION'] = df_melted['Max.Conc.tested']
print(df_melted['CONCENTRATION'][(df_melted['D'] == 'D1') & (df_melted['DRUG_NAME'] != 'ARV-825')])
df_melted['CONCENTRATION'][df_melted['D'] == 'D1'] = df_melted[df_melted['D'] == 'D1'].apply(lambda row: row['Max.Conc.tested']/math.pow(row['x'],4), axis=1)
df_melted['CONCENTRATION'][df_melted['D'] == 'D2'] = df_melted[df_melted['D'] == 'D2'].apply(lambda row: row['Max.Conc.tested']/math.pow(row['x'],3), axis=1)
df_melted['CONCENTRATION'][df_melted['D'] == 'D3'] = df_melted[df_melted['D'] == 'D3'].apply(lambda row: row['Max.Conc.tested']/math.pow(row['x'],2), axis=1)
df_melted['CONCENTRATION'][df_melted['D'] == 'D4'] = df_melted[df_melted['D'] == 'D4'].apply(lambda row: row['Max.Conc.tested']/row['x'], axis=1)
df_melted['CONCENTRATION'][df_melted['D'] == 'D5'] = df_melted[df_melted['D'] == 'D5'].apply(lambda row: row['Max.Conc.tested'], axis=1)




#negative = df_melted[['Patient.num', 'DRUG_NAME', 'sample']][(df_melted['Normalised_reponse'] < 0) | (df_melted['Normalised_reponse'] > 150)]
#negative = negative.set_index(['Patient.num', 'DRUG_NAME', 'sample']).index
#filtered_df = df_melted[~df_melted.set_index(['Patient.num', 'DRUG_NAME', 'sample']).index.isin(negative)].reset_index(drop=True)
#print(filtered_df)

#df_melted['Normalised_reponse'][df_melted['Normalised_reponse'] < 0] = 0
#df_melted['Normalised_reponse'][df_melted['Normalised_reponse'] > 100] = 100


df_melted.to_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_normalised_reponse_raw.csv')
df_melted['patient_sample'] = df_melted['Patient.num'] + '_'+ df_melted['sample']
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print(df_melted[df_melted['Patient.num'] == 'AML_001'].sort_values(by=['DRUG_NAME', 'sample']))

df_breese = df_melted[['DRUG_NAME', 'CONCENTRATION', 'patient_sample', 'Normalised_reponse']]

df_breese = df_breese.rename(columns={'patient_sample': 'SCREEN_NAME', 'Normalised_reponse': 'PERCENT_INHIBITION'})
print(df_breese.dtypes)
df_breese.to_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Breeze/karolinska_normalised_reponse_raw.xlsx', index=False, header=True)

