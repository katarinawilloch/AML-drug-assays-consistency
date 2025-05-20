import pandas as pd
import math
import matplotlib.pyplot as plt
pd.set_option('display.min_rows',30)
pd.set_option('display.max_rows', 100)
df = pd.read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Functional_Precision_Medicine_Tumor_Board_AML/6549394/cd-21-0410_supplementary_tables_suppst1.xlsx', sheet_name= "Table 7", skiprows= 2)

print(df)

# Melting the DataFrame to long format
df_melted_c = df.melt(id_vars=['Sample_ID', 'Chemical_compound'], 
                    value_vars=['C_1', 'C_2', 'C_3', 'C_4', 'C_5'], 
                    var_name='C', 
                    value_name='CONCENTRATION')
print('DF_MELTED')

df_melted_c['V'] = df_melted_c['C'].str.replace('C', 'V')
print(df_melted_c)

df_melted_v = df.melt(id_vars=['Sample_ID', 'Chemical_compound'], 
                    value_vars=['V_1', 'V_2', 'V_3', 'V_4', 'V_5'], 
                    var_name='V', 
                    value_name='Normalised_reponse')

print(df_melted_v)


df_melted_merged = pd.merge(df_melted_v[['Sample_ID', 'Chemical_compound', 'Normalised_reponse', 'V']], df_melted_c[['Sample_ID', 'Chemical_compound', 'CONCENTRATION', 'C', 'V']], on=['Sample_ID', 'Chemical_compound', 'V'])
print(df_melted_merged)

df_melted_merged.to_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/fimm_normalised_reponse_raw_supplemental.csv')


