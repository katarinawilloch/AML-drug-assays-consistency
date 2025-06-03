"""
Script: cleansing.py

Overview:
This script gets the most common drug names in pubchem found for each given 
drug name for the different studies.
It imports drug name information from each study and runs 
function get_drug_synonyms to get the most common drug names.

experimental data from Excel, performs PCA and statistical tests
(Wilcoxon, correlation), and generates plots including density, boxplots,
and volcano plots to explore associations between variables.

Main Steps:
- Load data
- Pubchem API to get most common drug names
- Visualize overlap in drug names between studies

Author: Katarina Willoch
Date: 2025-06-01
"""

#Loading libraries
import pandas as pd
from ydata_profiling import ProfileReport
import sweetviz as sv
import ast

#File locations

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

