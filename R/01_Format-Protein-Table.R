build_protein_table_function <- function(path_proteins, path_protein_info, path_protein_keys, path_bridge)
  
{
  
  #----------------------------------------------------------#
  #---------------Read in files -----------------------------#
  #----------------------------------------------------------#
  
  #Intensity table
  proteins_raw <- data.table::fread(path_proteins)
  
  #Info table
  protein_info <- data.table::fread(path_protein_info)
  
  #Info table
  protein_keys <- data.table::fread(path_protein_keys)
  
  #Bridging file
  bridge <- data.table::fread(path_bridge) |>
    dplyr::select(`SHARE ID Number`, `MESA Participant ID`) |>
    dplyr::rename(sidno = `SHARE ID Number`, idno = `MESA Participant ID`) |>
    dplyr::mutate(sidno = factor(sidno),
                  idno = factor(idno))
  
  raw_protein_N_all <- length(unique(proteins_raw$OlinkID))
  
  raw_protein_N_proteins <- proteins_raw |>
    dplyr::filter(OlinkID %in% protein_keys$OlinkID) |>
    dplyr::summarise(unique_count = dplyr::n_distinct(OlinkID)) 
  
  raw_N_all <- length(unique(proteins_raw$SampleID))
  
  Non_bridging_samples <- protein_info |>
    dplyr::filter(Sample_Type=="SAMPLE" & is.na(Note))
  
  raw_N_nobridge <- proteins_raw |>
    dplyr::filter(SampleID %in% Non_bridging_samples$SampleID) |>
    dplyr::summarise(unique_count = dplyr::n_distinct(SampleID)) 
    
  
 
  
  #------------------------------------------------------------#
  #---------------Format proteins -----------------------------#
  #------------------------------------------------------------#
    
  formatted_proteins <- proteins_raw |>
    #remove non-proteins (assays used for QC)
    dplyr::filter(OlinkID %in% protein_keys$OlinkID) |>
    #remove bridging samples |>
    dplyr::filter(SampleID %in% Non_bridging_samples$SampleID) |>
    #Remove QC failures
    dplyr::filter(is.na(QC_Warning) | QC_Warning !="EXCLUDED") |>
    dplyr::select(SampleID, OlinkID, NPX) |>
    tidyr::pivot_wider(id_cols=SampleID, names_from = OlinkID, values_from = NPX)|>
    #merge in exam and id info, using only info in mapping file (no protein assays)
    dplyr::left_join(dplyr::select(Non_bridging_samples, SampleID, sidno, Exam, Batch), 
                     dplyr::join_by(SampleID)) |>
    dplyr::mutate(sidno = factor(sidno),
                  Exam = factor(Exam),
                  Batch = factor(Batch)) |>
    dplyr::left_join(bridge, dplyr::join_by(sidno))
  
  sidno_by_id <- as.data.frame(table(formatted_proteins$sidno, formatted_proteins$Exam)) |>
    dplyr::rename(sidno = Var1,
                  exam = Var2) |>
    dplyr::arrange(Freq)
  
  dupe_ids <- sidno_by_id |>
    dplyr::filter(Freq>1) 
    
  #------------------------------------------------------------#
  #---------------Info for README -----------------------------#
  #------------------------------------------------------------#
  
  QC_info <- list(
    filenames = list(raw_proteins_file_name = path_proteins,
                     raw_protein_info_file_name = path_protein_info,
                     raw_bridgingfile_file_name = path_bridge,
                     raw_protein_keys_file_name = path_protein_keys
    ),
    
    raw_inputs = list(raw_protein_N_all = raw_protein_N_all,
                      raw_protein_N_proteins = raw_protein_N_proteins,
                      raw_N_all = raw_N_all,
                      raw_N_nobridge = raw_N_nobridge
    ),
    formatted_protein_file_info = list(
      unique_sampleids = length(unique(formatted_proteins$SampleID)),
      duplicated_sampleids = sum(duplicated(formatted_proteins$SampleID)),
      protein_N = dim(formatted_proteins)[2],
      sample_N = dim(formatted_proteins)[1]
      ),
    
    initial_dupes =  sidno_by_id
   
  )
  
    #----------------Outputs -----------------------------#
    
    list(
      
      QC_info_out = QC_info,
      Formatted_proteins_out = formatted_proteins

    )
    
  
}