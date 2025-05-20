library(readr)
library("ggplot2")
library(dplyr)
library(readxl)
library(openxlsx)
library(stringr) # for string operations

pt30r_repeat_260517_D1_sheet2 <- read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/Pt30r_repeat/pt30r_repeat_260517_D1.xls', sheet=2)
#plate_ref <- read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/selleck_dispensing.xls')

clean_plate_val_long <- function(input_df, platename, patient_id){
  input_df <- input_df[-(1:5), ]
  input_df <- input_df[-nrow(input_df), ]
  colnames(input_df) <- seq_len(ncol(input_df))
  input_df$row <- LETTERS[seq_len(nrow(input_df))]
  input_df[,-ncol(input_df)] <- lapply(input_df[,-ncol(input_df)], as.numeric)
  
  # Transform the dataframe to the desired format
  input_df_long <- input_df %>%
    pivot_longer(cols = -row,  
                 names_to = "col",       # Create new column for the index (X1, X2, etc.)
                 values_to = "Value")      # Create new column for the values
  
  # Combine the Label with the Index to get A1, A2, etc.
  input_df_long$WELLID <- paste(input_df_long$row, input_df_long$col, sep = "")
  input_df_long$plate <- platename
  input_df_long$patient_id <- patient_id
  return(input_df_long)
}

get_file_data <- function(main_folder, patient_id){
  # List all Excel files in the folder and subfolders
  excel_files <- list.files(path = main_folder, 
                            pattern = "\\.xlsx$|\\.xls$", 
                            recursive = TRUE, 
                            full.names = TRUE)
  
  raw_data <- data.frame()
  
  for (file in excel_files) {
    df_file <- read_excel(file, sheet = 2)
    
    file_name <- basename(file)
    platename <- str_split_i(file_name, '_',4)
    platename <- str_split_i(platename, '.xls',1)
    #print(file_name)
    print(platename)
    
    file_raw_df <- clean_plate_val_long(df_file, platename = platename, patient_id = patient_id)
    raw_data <- rbind(raw_data, file_raw_df)
  }
  return(raw_data)
}
# Set your main folder path
main_folder <- "/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/Pt30r_repeat"  # Replace this with the actual folder path
pt30r_repeat_df <- get_file_data(main_folder, "Pt30r_repeat")
main_folder <- "/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/Pt33_frozen"  # Replace this with the actual folder path
pt33_frozen_df <- get_file_data(main_folder, "Pt33_frozen")
main_folder <- "/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/Pt34r_frozen"  # Replace this with the actual folder path
pt34r_frozen_df <- get_file_data(main_folder, "Pt34r_frozen")

all_data_raw <- rbind(pt30r_repeat_df, pt33_frozen_df, pt34r_frozen_df)

reference_conc_drug <- read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/reference_drug_conc.xlsx')

final_data_raw <- merge(all_data_raw, reference_conc_drug, by.x = c("plate", "WELLID"), by.y = c("Destination Plate Barcode", "Destination Well"))

names(final_data_raw)

final_data_raw <- final_data_raw %>%
  dplyr::rename(WELL = WELLID, PLATE = plate, DRUG_NAME = SUPPLIER_REF, CONCENTRATION = final_conc, SCREEN_NAME = patient_id, WELL_SIGNAL = Value)
names(final_data_raw)
openxlsx::write.xlsx(final_data_raw[,c("WELL", "PLATE", "DRUG_NAME", "CONCENTRATION" , "SCREEN_NAME", "WELL_SIGNAL")], '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/breeze_input.xlsx')
