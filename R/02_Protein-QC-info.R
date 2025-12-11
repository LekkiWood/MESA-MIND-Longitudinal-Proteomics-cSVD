QC_proteins_function <- function(path_proteins, path_protein_info, path_protein_keys, Proteins_long)
  
{
  #Intensity table
  proteins_raw <- data.table::fread(path_proteins)
  
  protein_keys <- data.table::fread(path_protein_keys)
  
  #Info table
  protein_info <- data.table::fread(path_protein_info)
  
  #----------------------------------------------------------#
  #---------------Protein QC info----------------------------#
  #----------------------------------------------------------#
  qc_cols_proteins <- protein_info %>%
    dplyr::filter(Sample_Type == "CONTROL") %>%
    dplyr::pull(SampleID)
  
  message(length(qc_cols_proteins), " pooled QC IDs identified: ") 
  
  #temp <- subset(proteins_raw, OlinkID=="OID20419")
  #temp2 <- subset(temp, SampleID %in% qc_cols_proteins)
  #temp3 <- subset(temp, !SampleID %in% qc_cols_proteins)
  
  protein_QC <- proteins_raw |>
    #remove non-proteins (assays used for QC)
    dplyr::filter(OlinkID %in% protein_keys$OlinkID) |>
    #Select vars for making wide format
    dplyr::select(SampleID, OlinkID, NPX) |>
    #Make people columns and proteins rows
    tidyr::pivot_wider(id_cols=OlinkID, names_from = SampleID, values_from = NPX) |>
    #Select protein names and pooled samples
    dplyr::select(OlinkID, dplyr::all_of(qc_cols_proteins)) |>
    #Rowise calcs across proteins
    dplyr::rowwise() |>
    #Calcs
    dplyr::mutate(
      sigma = sd(dplyr::c_across(dplyr::all_of(qc_cols_proteins)), na.rm = TRUE),
      cv_percent = sqrt(2^(sigma^2) - 1) * 100
    ) |>
    dplyr::ungroup() |>
    dplyr::select(OlinkID, cv_percent) |>
    dplyr::left_join(protein_keys, dplyr::join_by(OlinkID)) |> 
    dplyr::group_by(Assay) |>
    dplyr::mutate(Retain = dplyr::case_when(dplyr::n() == 1 ~ 0,
                                            cv_percent== min(cv_percent, na.rm = TRUE) ~ 1,
                                            TRUE ~ 2 ) )  |>
    dplyr::ungroup() |>
    dplyr::mutate(Retain = factor(Retain, labels=c("Unique Protein", "Duplicated Protein with lowest CV", "Duplicated Protein and not lowest CV"))) 
  
  #check <- subset(Protein_QC, Retain != "Unique Protein")
  

  
  list(final_proteins_mapping = protein_QC)
  
}