import pandas as pd 
import os
import glob
import re 

karolinska_reponse = pd.DataFrame()
for files in glob.glob('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Karolinska/all_files/' + '*.xlsx'):
    print(files)
    basename = os.path.basename(files)
    parts = basename.split('_')
    match = re.search(r'[A-Z]+_\d+', basename)

    if len(parts) >= 1:
        part1 = parts[0]
    if match:
        p_id = match.group(0)

    
    print(f"Part 1: {part1}, Patient_id: {p_id}")
    df = pd.read_excel(files, sheet_name='EC50')
    df['sample'] = parts[0]
    df['Patient.num'] = p_id
    print(df)
    karolinska_reponse = pd.concat([karolinska_reponse, df], ignore_index=True)

karolinska_reponse.to_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/karolinska_breese_output.csv')