library(readr)
library("ggplot2")
library(dplyr)
library(readxl)
library(openxlsx)

drug_name_ref <- read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/Selleck_CBP_DrugNames.xlsx')
plate_ref <- read_excel('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/selleck_dispensing.xls')
names(drug_name_ref)
names(plate_ref)
reference_plate_layout <- merge(drug_name_ref, plate_ref, by.x = 'FORMATTED_BATCH_ID', by.y = 'Sample Name', all=TRUE)
replacement_values <- c(`B004` = 10,
                        `B004-I1` = 0.1,
                        `B004-I2` = 0.001)
reference_plate_layout$stock_conc <- replacement_values[reference_plate_layout$`Source Plate Barcode`]

calculate_final_conc <- function(stock_conc, transfer_vol, final_vol = 25) {
  return ((stock_conc * transfer_vol) / final_vol)
}

reference_plate_layout$final_conc <- mapply(calculate_final_conc, reference_plate_layout$stock_conc, reference_plate_layout$`Transfer Volume`)
reference_plate_layout$final_conc <- reference_plate_layout$final_conc/0.001

reference_plate_layout <- reference_plate_layout %>%
  mutate(FORMATTED_BATCH_ID = ifelse(FORMATTED_BATCH_ID == "Control", 'BzCl', FORMATTED_BATCH_ID))
reference_plate_layout <- reference_plate_layout %>%
  mutate(FORMATTED_BATCH_ID = ifelse(FORMATTED_BATCH_ID == "DMSO", 'dmso', FORMATTED_BATCH_ID))
reference_plate_layout <- reference_plate_layout %>%
  mutate(SUPPLIER_REF = ifelse(is.na(SUPPLIER_REF), FORMATTED_BATCH_ID, SUPPLIER_REF))

write.xlsx(reference_plate_layout, '/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Original data/Enserink-lab/repeated from frozen/reference_drug_conc.xlsx')

