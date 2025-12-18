#Proteins_long = tar_read(Proteins_long)
#QC_fil = tar_read(Protein_mapping_file)

final_proteins_function <- function(Proteins_long, QC_file){
  
  Included_proteins <- QC_file |>
    dplyr::filter(Retain=="Unique Protein" | Retain=="Duplicated Protein with lowest CV") |>
    dplyr::pull(OlinkID)
  
  Final_proteins <- Proteins_long |>
    dplyr::select(sidno, idno, SampleID, Exam, Batch, dplyr::all_of(Included_proteins))
  
  sidno_by_exam <- Final_proteins |>
    dplyr::group_by(Exam) |>
    dplyr::summarise(N_Pps = dplyr::n_distinct(sidno))
  
  list(
    Final_proteins = Final_proteins,
    N_by_exam = sidno_by_exam
  )
  
  
  }