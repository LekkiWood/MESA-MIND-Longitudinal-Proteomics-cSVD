
rm(list=ls())
library(targets)
library(dplyr)

path_bridge = "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESA-SHARE_IDList_Labeled.csv"
path_E1_covs = "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe1FinalLabel02092016.dta"
#E6_covs
path_E6_covs = "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe6_FinalLabel_20220513.dta"
path_E5_covs = "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAe5_FinalLabel_20140613.dta"
#afib
path_afib = "/media/RawData/MESA/MESA-Phenotypes/MESA-SHARe-Phenos/SHARe_MesaEventsThruYear2020_AF_DS.txt"
path_apoe = "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESA_ApoE_03102014.sas7bdat"
path_cvd = "/media/RawData/MESA/MESA-Phenotypes/MESA-Website-Phenos/MESAEvThru2020AllCohort_20241120.dta"
path_mb = "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESAe6as253as301_BMRICMB_08052025.csv"
path_evps = "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESAe6as253as301_BMRIPVS_20250310.csv"
path_wmh = "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/MESAe6anyFIRST_BMRIWMHVol_20240422.csv"
path_icv = "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/SHARe_AncilMesaAF_BMRIROIVol_DS.txt"
path_wmfa = "/media/RawData/MESA/MESA-Phenotypes/MESA-MIND/mesae6anyfirst_bmriTotalFAMUSE_20250828.csv"

cleaned_proteins <- targets::tar_read(Proteins_long_clean)