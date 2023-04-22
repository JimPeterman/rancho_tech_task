
library(readxl)
library(janitor)
library(dplyr)
library(here)

# Import the different datasets and clean up the column names.
df_patient <- read_xlsx(
  here::here("Technical Test - Data Wrangling - 20230413 (3).xlsx"), 
  sheet = "Patient_clinical_data"
  )
df_patient <- janitor::clean_names(df_patient)

df_tissue_data <- read_xlsx(
  here::here("Technical Test - Data Wrangling - 20230413 (3).xlsx"), 
  sheet = "Tissue Sample Metadata"
  )
df_tissue_data <- janitor::clean_names(df_tissue_data)

df_serum_data <- read_xlsx(
  here::here("Technical Test - Data Wrangling - 20230413 (3).xlsx"), 
  sheet = "Serum Protein data"
)
df_serum_data <- janitor::clean_names(df_serum_data)

# Don't clean these columns so they can match sample ID's.
df_rna_seq <- read_xlsx(
  here::here("Technical Test - Data Wrangling - 20230413 (3).xlsx"), 
  sheet = "RNA-seq (RPKM)"
)


# Create an empty dataframe to write the results to.
df_summary <- data.frame(
  Study_ID = character(),
  Patient_ID = numeric(),
  Unique_Patient_ID = character(),
  Sex = character(),
  Age = integer(),
  Sample_ID = character(),
  Sample_General_Pathology = character(),
  Material_type = character(),
  Gene_Symbol = character(),
  Result = numeric(),
  Result_Units = character(),
  Status = character()
)


# Go through different datasets and write results to df_summary.
# The primary loop is through the patient ID dataset involving:
# 1. Gathering patient info
# 2. Nested loop gathering any tissue data for that patient 
#   (and another nested loop to gather any gene data for that sample)
#     - Write results to df_summary
# 3. Nested loop gathering any serum data for that patient
#     - Write results to df_summary

for (i in 1:length(df_patient$patient_number)){
  # Collect the patient information.
  study_id <- toupper(df_patient$study_id[i])
  patient_id <- df_patient$patient_number[i]
  unique_patient_id <- paste(
    study_id, as.character(patient_id), sep = "_"
  )
  sex <- ifelse(df_patient$sex[i] == "M", "MALE", 
               ifelse(df_patient$sex[i] == "F", "FEMALE", "NA"))
  age <- as.integer(df_patient$age[i])

  # Add in the data from the tissue dataset.
  if (patient_id %in% df_tissue_data$patient_number){
    # Temporary tissue dataset containing only the 1 patient.
    temp_df_tissue <- filter(df_tissue_data, patient_number == patient_id)

    for (j in 1:length(temp_df_tissue$patient_number)){
      # Collect the tissue information.
      sample_id <- toupper(temp_df_tissue$sample[j])
      sample_general_pathology <- ifelse(
        temp_df_tissue$sample_type[j] == "Normal", "NORMAL",
        ifelse(temp_df_tissue$sample_type[j] == "Metastic Lung", "METASTATIC",
               ifelse(temp_df_tissue$sample_type[j] == "Liver Tumor", "PRIMARY",
                      "NA")))
      
      # Add in the data from the RNA-seq dataset.
      if (sample_id %in% colnames(df_rna_seq)){
        # Collect the gene information if present and add to summary dataframe.
        for (gene_symbol in df_rna_seq$GeneID){
          result <- df_rna_seq[df_rna_seq$GeneID == gene_symbol,][sample_id][[1]]
          result <- ifelse(is.numeric(result), result, NaN)
          status <- ifelse(is.na(result), "ERROR", "NA")
          
          df_summary <- add_row(df_summary,
                                Study_ID = study_id,
                                Patient_ID = patient_id,
                                Unique_Patient_ID = unique_patient_id,
                                Sex = sex,
                                Age = age,
                                Sample_ID = sample_id,
                                Sample_General_Pathology = sample_general_pathology,
                                Material_type = "RNA",
                                Gene_Symbol = gene_symbol,
                                Result = result,
                                Result_Units = "RPKM",
                                Status = status
                                )
        }
      } else {
        # No gene information but still write to dataframe (indicating missing).
        df_summary <- add_row(df_summary,
                              Study_ID = study_id,
                              Patient_ID = patient_id,
                              Unique_Patient_ID = unique_patient_id,
                              Sex = sex,
                              Age = age,
                              Sample_ID = sample_id,
                              Sample_General_Pathology = sample_general_pathology,
                              Material_type = "RNA",
                              Gene_Symbol = "NA",
                              Result = NaN,
                              Result_Units = "NA",
                              Status = "NOT DONE"
                              )
      }
    }
  }
  
  # Add in the data from the serum dataset.
  if (patient_id %in% df_serum_data$patient){
    # Temporary serum dataset containing only the 1 patient.
    temp_df_serum <- filter(df_serum_data, patient == patient_id)
    
    for (l in 1:length(temp_df_serum$sample)){
      # Collect the sample information.
      sample_id <- toupper(temp_df_serum$sample[l])
      sample_general_pathology <- "NA"
      
      for (m in 1:2){
        # Collect the gene symbol information.
        if (m == 1){
          gene_symbol <- "IL6"
          result <- suppressWarnings(as.numeric(temp_df_serum$serum_il_6_g_l[l]))
          result <- ifelse(is.na(result), NaN, result)
          result_units <- "g/L"
  
        } else {
          gene_symbol <- "IL6R"
          result <- suppressWarnings(as.numeric(
            temp_df_serum$serum_il_6_receptor_mg_l[l]))
          # Convert to g/L from mg/L
          result <- ifelse(is.na(result), NaN, result/1000)
          result_units <- "g/L"
        }
        
        # Determine the status information.
        status <- ifelse(is.na(result), "ERROR", "NA")
        # Add results to summary dataframe.
        df_summary <- add_row(df_summary, 
                              Study_ID = study_id,
                              Patient_ID = patient_id,
                              Unique_Patient_ID = unique_patient_id,
                              Sex = sex,
                              Age = age,
                              Sample_ID = sample_id,
                              Sample_General_Pathology = sample_general_pathology,
                              Material_type = "SERUM",
                              Gene_Symbol = gene_symbol,
                              Result = result,
                              Result_Units = result_units,
                              Status = status
                              )
      }
    }
  }
}

# Per instructions, make all missing/NULL values "NA"
df_summary <- df_summary %>% replace(is.na(.), "NA")

# Save file.
write.csv(df_summary, here::here("curated_data.csv"), row.names = FALSE)

