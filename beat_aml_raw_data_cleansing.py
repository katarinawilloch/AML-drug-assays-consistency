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

beat_aml_dose_response_raw = get_df('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/beat_aml_inhibitor.csv')
print(beat_aml_dose_response_raw)

