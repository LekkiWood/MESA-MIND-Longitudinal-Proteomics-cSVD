build_traits_function <- function(path_E1_covs, path_E5_covs, path_E6_covs, path_afib, path_bridge, path_apoe, 
                                  path_cvd, path_mb, path_evps, path_wmh, path_icv, path_wmfa, cleaned_proteins)
  
{
  
  #Bridging file
  bridge <- data.table::fread(path_bridge) |>
    dplyr::select(`SHARE ID Number`, `MESA Participant ID`) |>
    dplyr::rename(sidno = `SHARE ID Number`, idno = `MESA Participant ID`)
  

  
  #---------------------------Build outcomes file-----------------------#
  

  ###############
  #Microbleeds
  ###############
  
  MB_info <- read.csv(path_mb) |>
    dplyr::mutate(mb_present_all = dplyr::case_when(mb_n_total > 0 ~ 1,
                                                    mb_n_total == 0 ~ 0,
                                                    TRUE ~ NA_real_)) |>
    dplyr::mutate(mb_present = dplyr::case_when(qsm_swi_image_quality == 4 ~ NA_real_,
                                                TRUE ~ mb_present_all)) |>
    dplyr::select(idno, mb_present_all, mb_present, qsm_swi_image_quality)
  
  ###############
  #EVPS
  ###############
  
  EVPS_info <- read.csv(path_evps) |>
    dplyr::mutate(epvs = dplyr::case_when(pvs_exclude == 1 ~ NA_real_,
                                          TRUE ~ epvs_wholebrain_vol)) |>
    dplyr::select(idno, pvs_exclude, epvs)
  
  ###############
  #WMH
  ###############
  
  WMH_info <- read.csv(path_wmh) |>
    dplyr::mutate(wmh = dplyr::case_when(is.na(wmh_exclude) ~ NA_real_,
                                         TRUE ~ wmh_wm/1000)) |>
    dplyr::select(idno, wmh, wmh_exclude)
  

  ###############
  #WMFA
  ###############
  
  WMFA_info <- read.csv(path_wmfa)|>
    dplyr::mutate(fa = dplyr::case_when(fa_exclude==0 ~ fa_wm,
                                        TRUE ~ NA_real_)) |>
    dplyr::select(idno, qc_code, fa_exclude, fa)
  
  
  ###############
  #Traits-file
  ###############
 
  Traits <- MB_info |>
    dplyr::full_join(EVPS_info, dplyr::join_by(idno)) |>
    dplyr::full_join(WMH_info, dplyr::join_by(idno)) |>
    dplyr::full_join(WMFA_info, dplyr::join_by(idno)) |>
    dplyr::left_join(bridge, dplyr::join_by(idno)) |>
    dplyr::filter(!is.na(mb_present) | !is.na(epvs) | !is.na(epvs) | !is.na(wmh) | !is.na(fa)) 
    
  
  
  mind_ids <- Traits |>
    dplyr::select(idno, sidno)
  
  #---------------------------N with protein data-----------------------#

  
  protein_and_MIND_ids <- cleaned_proteins |>
    #dplyr::full_join(cleaned_proteins, dplyr::join_by(idno, sidno)) |>
    dplyr::filter(idno %in% mind_ids$idno) |>
    dplyr::select(idno, sidno, Exam)
  
  
  Traits <- Traits |>
    dplyr::filter(idno %in% protein_and_MIND_ids$idno)
  

  
  #------------------------------Build covariates-------------------------------
  

  
  #########Time invariant
  
  ###############
  #ICV
  ###############
  
  ICV_info <- read.table(path_icv, header=TRUE, sep="\t")|>
    dplyr::mutate(sidno = subject_id) |>
    dplyr::select(sidno, icv)
  
  ###############
  #Afib info
  ###############
  afib_info <- read.table(path_afib, header=TRUE, sep="\t") |>
    dplyr::mutate(AFprevalent = ifelse(af2020 == 1, 1, 0)) |>
    dplyr::mutate(sidno = subject_id) |>
    dplyr::select(sidno, AFprevalent)
  
  
  ###############
  #CVD
  ###############
  
  CVD_info <- foreign::read.dta(path_cvd) |>
    dplyr::mutate(MIprevalent = dplyr::case_when(mi == "No" ~ 0,
                                                 mi == "Yes" ~ 1,
                                                 TRUE ~ NA_real_),
                  CHFprevalent = dplyr::case_when(chf == "No" ~ 0,
                                                  chf == "Yes" ~ 1,
                                                 TRUE ~ NA_real_)) |>
    dplyr::select(idno, MIprevalent,CHFprevalent) 
  
  ###############
  #Apoe
  ###############
  
  ApoE_info <- haven::read_sas(path_apoe) |>
    dplyr::mutate(E4 = dplyr::case_when(ApoE %in% c(24, 34, 44) ~ 1, 
                                        ApoE %in% c(22, 23, 33) ~ 0,
                                        TRUE ~ 2)) |>
    dplyr::mutate(E4 = as.factor(E4)) |>
    dplyr::select(idno, E4)

  
  E1_always <- foreign::read.dta(path_E1_covs) |>
    dplyr::mutate(edu = dplyr::case_when(
      is.na(educ1) ~ NA_real_ ,
      educ1 == "0: NO SCHOOLING" ~ 0,
      educ1 == "1: GRADES 1-8" ~ 0,
      educ1 == "2: GRADES 9-11" ~ 0,
      TRUE ~ 1),
      gender = dplyr::case_when(gender1=="0: FEMALE" ~ 0,
                                gender1=="1: MALE" ~ 1,
                                TRUE ~ NA_real_ ),
      race = dplyr::case_when(race1c=="1: white, CAUCASIAN" ~ 1,
                              race1c=="2: CHINESE-AMERICAN" ~ 2,
                              race1c=="3: black, AFRICAN-AMERICAN" ~ 3,
                              race1c=="4: HISPANIC" ~ 4,
                              TRUE ~ NA_real_),
      BL_age = age1c
    ) |>
    dplyr::left_join(bridge, dplyr::join_by(idno)) |>
    dplyr::select(idno, sidno, edu, gender, race, BL_age)
  
  E6_always <- foreign::read.dta(path_E6_covs) |>
    dplyr::mutate(site = dplyr::case_when(site6c=="WFU" ~ 0,
                                          site6c=="COL" ~ 1,
                                          site6c=="JHU" ~ 2,
                                          site6c=="UMN" ~ 3,
                                          site6c=="NWU" ~ 4,
                                          site6c=="UCLA" ~ 6,
                                          TRUE ~ NA_real_), 
                  ldl = ldl6,
                  sbp = sbp6c,
                  htnmeds = dplyr::case_when(htnmed6c =="NO" ~ 0,
                                             htnmed6c =="YES" ~ 1,
                                             TRUE ~ NA_real_)
                  )|>
    dplyr::left_join(bridge, dplyr::join_by(idno)) |>
    dplyr::left_join(ICV_info, dplyr::join_by(sidno)) |>
    dplyr::left_join(ApoE_info, dplyr::join_by(idno)) |>
    dplyr::left_join(afib_info, dplyr::join_by(sidno)) |>
    dplyr::left_join(CVD_info, dplyr::join_by(idno)) |>
    dplyr::select(idno, sidno, site, ldl, sbp, icv, E4, AFprevalent, MIprevalent, CHFprevalent, htnmeds)
                    
  
  
  
  #########By exam
  
  ###############
  #E1
  ###############
  
  E1_traits <- foreign::read.dta(path_E1_covs) |>
    dplyr::mutate(age = age1c,
                  BMI = bmi1c,
                  diabetes = dplyr::case_when(dm031c=="NORMAL" ~ 0,
                                              dm031c=="IFG" ~ 0,
                                              dm031c=="Untreated DIABETES" ~ 1,
                                              dm031c=="Treated DIABETES" ~ 1,
                                              TRUE ~ NA_real_),
                  egfr = egfr1c,
                  smoking = dplyr::case_when(cig1c=="0: NEVER" ~ 1,
                                             cig1c=="1: FORMER" ~ 2,
                                             cig1c=="2: CURRENT" ~ 3,
                                             TRUE ~ NA_real_),
                  time = 0,
                  )|>
    dplyr::select(idno, age, BMI, diabetes, egfr, smoking, time) |>
    dplyr::left_join(E1_always, dplyr::join_by(idno)) |>
    dplyr::left_join(E6_always, dplyr::join_by(idno, sidno)) |>
    dplyr::left_join(Traits, dplyr::join_by(idno, sidno)) |>
    dplyr::mutate(Exam = 1) |>
    dplyr::select(idno, sidno, Exam, time, BL_age, age, gender, site, edu, race, BMI, 
                  smoking, ldl, sbp, diabetes, htnmeds, AFprevalent, MIprevalent, CHFprevalent, 
                  E4, egfr, icv, fa, wmh, epvs, mb_present) 
  
  
  ###############
  #E5
  ###############
  
  E5_traits <- foreign::read.dta(path_E5_covs) |>
    dplyr::mutate(age = age5c,
                  BMI = bmi5c,
                  diabetes = dplyr::case_when(dm035c=="NORMAL" ~ 0,
                                              dm035c=="IMPAIRED FASTING GLUCOSE" ~ 0,
                                              dm035c=="UNTREATED DIABETES" ~ 1,
                                              dm035c=="TREATED DIABETES" ~ 1,
                                              TRUE ~ NA_real_),
                  egfr = egfr5c,
                  smoking = dplyr::case_when(cig5c=="Never" ~ 1,
                                             cig5c=="Former" ~ 2,
                                             cig5c=="Current" ~ 3,
                                             TRUE ~ NA_real_),
                  time = e15dyc
    )|>
    dplyr::select(idno, age, BMI, diabetes, egfr, smoking, time) |>
    dplyr::left_join(E1_always, dplyr::join_by(idno)) |>
    dplyr::left_join(E6_always, dplyr::join_by(idno, sidno)) |>
    dplyr::left_join(Traits, dplyr::join_by(idno, sidno)) |>
    dplyr::mutate(Exam = 5) |>
    dplyr::select(idno, sidno, Exam, time, BL_age, age, gender, site, edu, race, BMI, 
                  smoking, ldl, sbp, diabetes, htnmeds, AFprevalent, MIprevalent, CHFprevalent, 
                  E4, egfr, icv, fa, wmh, epvs, mb_present) 
  
  ###############
  #E6
  ###############
  
  E6_traits <- foreign::read.dta(path_E6_covs) |>
    dplyr::mutate(age = age6c,
                  BMI = bmi6c,
                  diabetes = dplyr::case_when(dm036c=="NORMAL" ~ 0,
                                              dm036c=="IMPAIRED FASTING GLUCOSE" ~ 0,
                                              dm036c=="UNTREATED DIABETES" ~ 1,
                                              dm036c=="TREATED DIABETES" ~ 1,
                                              TRUE ~ NA_real_),
                  egfr = egfr6c,
                  smoking = dplyr::case_when(cig6c=="NEVER" ~ 1,
                                             cig6c=="FORMER" ~ 2,
                                             cig6c=="CURRENT" ~ 3,
                                             TRUE ~ NA_real_),
                  time = e16dyc
    )|>
    dplyr::select(idno, age, BMI, diabetes, egfr, smoking, time) |>
    dplyr::left_join(E1_always, dplyr::join_by(idno)) |>
    dplyr::left_join(E6_always, dplyr::join_by(idno, sidno)) |>
    dplyr::left_join(Traits, dplyr::join_by(idno, sidno)) |>
    dplyr::mutate(Exam = 6) |>
    dplyr::select(idno, sidno, Exam, time, BL_age, age, gender, site, edu, race, BMI, 
                  smoking, ldl, sbp, diabetes, htnmeds, AFprevalent, MIprevalent, CHFprevalent, 
                  E4, egfr, icv, fa, wmh, epvs, mb_present) 
  
  #----------------------N with MIND, protein and covs -----------------------------------
  
  missing_cov_info <- dplyr::bind_rows(E1_traits, E5_traits, E6_traits) |>
    dplyr::filter(idno %in% protein_and_MIND_ids$idno) 
  
  full_db <- dplyr::bind_rows(E1_traits, E5_traits, E6_traits) |>
    dplyr::filter(idno %in% protein_and_MIND_ids$idno)
  
  
  #|>
  #dplyr::filter(!is.na(age) & !is.na(time) & !is.na(gender) & !is.na(site) & !is.na(edu) & !is.na(race) & !is.na(BMI) & !is.na(smoking) & 
  #!is.na(ldl) & !is.na(sbp) & !is.na(diabetes) & !is.na(htnmeds) & !is.na(AFprevalent) & !is.na(MIprevalent) & 
  #!is.na(CHFprevalent) & !is.na(E4) & !is.na(egfr) & !is.na(icv))
  

  
  
  # protein_and_MIND_and_cov_ids <- full_db |>
  # dplyr::select(idno, sidno, Exam)
  
  #------------------------------------------------------------#
  #---------------Info for README -----------------------------#
  #------------------------------------------------------------#
  
  QC_info <- list(
    filenames = list(E1_covs_filename = path_E1_covs,
                     E5_covs_filename = path_E5_covs,
                     E6_covs_filename = path_E6_covs,
                     afib_filename = path_afib,
                     apoe_filename = path_apoe,
                     cvd_filename = path_cvd,
                     mb_filename = path_mb,
                     evps_filename = path_evps,
                     wmh_filename = path_wmh,
                     icv_filename = path_icv,
                     wmfa_filename = path_wmfa)
    )
  
  #----------------Outputs -----------------------------#
  
  list(
    
    QC_info_out = QC_info,
    Traits_table = full_db,
    missing_cov_info= missing_cov_info,
    mind_ids = mind_ids,
    protein_and_MIND_ids = protein_and_MIND_ids #,
    #protein_and_MIND_and_cov_ids = protein_and_MIND_and_cov_ids
    )
  
  
}