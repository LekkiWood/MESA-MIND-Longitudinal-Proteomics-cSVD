build_traits_function <- function(path_E1_covs, path_E5_covs, path_E6_covs, path_afib, path_bridge, path_apoe, 
                                  path_cvd, path_mb, path_evps, path_wmh, path_icv, path_wmfa, cleaned_proteins)
  
{
  #Bridging file
  bridge <- data.table::fread(path_bridge) |>
    dplyr::select(`SHARE ID Number`, `MESA Participant ID`) |>
    dplyr::rename(sidno = `SHARE ID Number`, idno = `MESA Participant ID`)
  
  ##########################################################################################
  #.                  Select those with protein data                                       #
  ##########################################################################################
  
  protein_ids <- cleaned_proteins |>
    dplyr::select(idno, Exam)
  
  
  ##########################################################################################
  #.                  Select those with Outcomes                                           #
  ##########################################################################################
  
  
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
    dplyr::select(idno, epvs, pvs_exclude)
  
  ###############
  #WMH
  ###############
  
  WMH_info <- read.csv(path_wmh) |>
    dplyr::mutate(wmh = dplyr::case_when(is.na(wmh_exclude) ~ NA_real_,
                                         TRUE ~ wmh_wm/1000)) |>
    dplyr::select(idno, wmh, wmh_exclude)
  
  ###############
  #ICV
  ###############
  
  ICV_info <- read.table(path_icv, header=TRUE, sep="\t")|>
    dplyr::mutate(sidno = subject_id) |>
    dplyr::select(sidno, icv)
  
  ###############
  #WMFA
  ###############
  
  WMFA_info <- read.csv(path_wmfa)|>
    dplyr::mutate(fa = dplyr::case_when(fa_exclude==0 ~ fa_wm,
                                        TRUE ~ NA_real_)) |>
    dplyr::select(idno, qc_code, fa_exclude, fa)
  
  
  ###############
  #Merge
  ###############
  
  Traits <- bridge |>
    dplyr::full_join(MB_info, dplyr::join_by(idno)) |>
    dplyr::full_join(EVPS_info, dplyr::join_by(idno)) |>
    dplyr::full_join(WMH_info, dplyr::join_by(idno)) |>
    dplyr::full_join(ICV_info, dplyr::join_by(sidno)) |> 
    dplyr::full_join(WMFA_info, dplyr::join_by(idno)) 
  
  protein_and_MIND_ids <- Traits |>
    dplyr::full_join(cleaned_proteins, dplyr::join_by(idno)) |>
    dplyr::filter(idno %in% protein_ids$idno) |>
    dplyr::filter(!is.na(fa) | !is.na(wmh) | !is.na(epvs) | !is.na(mb_present)) |>
    dplyr::select(idno, Exam)
  
  
  ##########################################################################################
  #.                  Select those with Covariates.                                        #
  ##########################################################################################
  

  
  
  
  E1_covs <- foreign::read.dta(path_E1_covs)
  E5_covs <- foreign::read.dta(path_E5_covs)
  E6_covs <- foreign::read.dta(path_E6_covs) |>
    dplyr::select(-race1c, -gender1, -age1c, -agecat1c, -mexican1, -dominic1, -puert1, -cuban1,
  -othhisp1, -lang1)
  
  ApoE_info <- haven::read_sas(path_apoe) 
  
  afib_info <- read.table(path_afib, header=TRUE, sep="\t") |>
    dplyr::mutate(sidno = subject_id) |>
    dplyr::select(-subject_id)
  
  CVD_info <- foreign::read.dta(path_cvd) |>
    dplyr::select(idno, mi, chf)
  
  ######################################
  #.   Time invariant                  #
  ######################################
  
  #idno, sidno
  #education, gender, race/ethnicity, E4, site, ldl, AFprevalent, MIprevalent, 
  #CHFprevalent, systolic blood pressure 
  
  
  Covs <- E1_covs |>
    dplyr::left_join(bridge, dplyr::join_by(idno)) |>
    dplyr::left_join(ApoE_info, dplyr::join_by(idno)) |>
    dplyr::left_join(E6_covs, dplyr::join_by(idno)) |>
    dplyr::left_join(CVD_info, dplyr::join_by(idno)) |>
    dplyr::left_join(afib_info, dplyr::join_by(sidno)) |>
    dplyr::mutate(
      edu = dplyr::case_when(
        is.na(educ1) ~ NA_real_ ,
        educ1 == "0: NO SCHOOLING" ~ 0,
        educ1 == "1: GRADES 1-8" ~ 0,
        educ1 == "2: GRADES 9-11" ~ 0,
        TRUE ~ 1
        ),
      gender = dplyr::case_when(
        gender1=="0: FEMALE" ~ 0,
        gender1=="1: MALE" ~ 1,
        TRUE ~ NA_real_ 
        ), 
      race=dplyr::case_when(
        race1c=="1: white, CAUCASIAN" ~ 1,
        race1c=="2: CHINESE-AMERICAN" ~ 2,
        race1c=="3: black, AFRICAN-AMERICAN" ~ 3,
        race1c=="4: HISPANIC" ~ 4,
        TRUE ~ NA_real_
        ),
      E4 = dplyr::case_when(
        ApoE %in% c(24, 34, 44) ~ 1, 
        ApoE %in% c(22, 23, 33) ~ 0,
        TRUE ~ 2
        ),
      site = dplyr::case_when(
        site6c=="WFU" ~1,
        site6c=="COL" ~2,
        site6c=="JHU" ~3,
        site6c=="UMN" ~4,
        site6c=="NWU" ~5,
        site6c=="UCLA" ~6,
        TRUE ~ NA_real_
      ),
      ldl = ldl6,
      AFprevalent = dplyr::case_when(
        af2020 == 0 ~ 0,
        af2020 == 1 ~ 1,
        TRUE ~ NA_real_
        ),
      MIprevalent = dplyr::case_when(
        mi == "No" ~ 0,
        mi == "Yes" ~ 1,
        TRUE ~ NA_real_
      ),
      CHFprevalent = dplyr::case_when(
        chf == "No" ~ 0,
        chf == "Yes" ~ 1,
        TRUE ~ NA_real_
      ),
      htnmed = dplyr::case_when(
        htnmed6c == "NO" ~ 0,
        htnmed6c == "YES" ~ 1,
        TRUE ~ NA_real_
      ),
      sbp = sbp6c,
      ) |>
    dplyr::mutate(
      E4 = as.factor(E4),
      race = as.factor(race),
      gender = as.factor(gender),
      site = as.factor(site)
      ) |>
    dplyr::select(idno, sidno, edu, gender, race, E4, site, ldl, AFprevalent, MIprevalent, CHFprevalent, sbp, htnmed) |>
    dplyr::left_join(Traits, dplyr::join_by(idno, sidno))
                   
    
        
  ######################################
  #.   Time variant                    #
  # Age, eGFR, BMI, smoking, diabetes  #
  ######################################  
  
  ##E1
  
Traits_E1 <- E1_covs |>
    dplyr::mutate(
      age = age1c,
      egfr = cepgfr1c,
      BMI = bmi1c,
      smoking = dplyr::case_when(
        cig1c=="0: NEVER" ~ 0,
        cig1c=="1: FORMER" ~ 1,
        cig1c=="2: CURRENT" ~ 2,
        TRUE ~ NA_real_
      ),
      diabetes = dplyr::case_when(
        dm031c =="NORMAL" ~0, 
        dm031c =="IFG" ~0,
        dm031c =="Untreated DIABETES" ~1,
        dm031c =="Treated DIABETES" ~1,
        TRUE ~ NA_real_)
    ) |>
    dplyr::select(idno, age, egfr, BMI, smoking, diabetes) |>
    dplyr::left_join(Covs, dplyr::join_by(idno)) |>
    dplyr::mutate(Exam = 1) |>
    dplyr::filter(idno %in% protein_and_MIND_ids$idno)
                       
  ##E5
  
  Traits_E5 <- E5_covs |>
    dplyr::mutate(
      age = age5c,
      egfr = cepgfr5c,
      BMI = bmi5c,
      smoking = dplyr::case_when(
        cig5c=="Never" ~ 0,
        cig5c=="Former" ~ 1,
        cig5c=="Current" ~ 2,
        TRUE ~ NA_real_
      ),
      diabetes = dplyr::case_when(
        dm035c =="NORMAL" ~0, 
        dm035c =="IMPAIRED FASTING GLUCOSE" ~0,
        dm035c =="UNTREATED DIABETES" ~1,
        dm035c =="TREATED DIABETES" ~1,
        TRUE ~ NA_real_)
    ) |>
    dplyr::select(idno, age, egfr, BMI, smoking, diabetes) |>
    dplyr::left_join(Covs, dplyr::join_by(idno)) |>
    dplyr::mutate(Exam = 5) |>
    dplyr::filter(idno %in% protein_and_MIND_ids$idno)
  
  ##E6
        
  Traits_E6 <- E6_covs |>
    dplyr::mutate(
      age = age6c,
      egfr = cepgfr6c,
      BMI = bmi6c,
      smoking = dplyr::case_when(
        cig6c=="NEVER" ~ 0,
        cig6c=="FORMER" ~ 1,
        cig6c=="CURRENT" ~ 2,
        TRUE ~ NA_real_
      ),
      diabetes = dplyr::case_when(
        dm036c =="NORMAL" ~0, 
        dm036c =="IMPAIRED FASTING GLUCOSE" ~0,
        dm036c =="UNTREATED DIABETES" ~1,
        dm036c =="TREATED DIABETES" ~1,
        TRUE ~ NA_real_)
    ) |>
    dplyr::select(idno, age, egfr, BMI, smoking, diabetes) |>
    dplyr::left_join(Covs, dplyr::join_by(idno)) |>
    dplyr::mutate(Exam = 6) |>
    dplyr::filter(idno %in% protein_and_MIND_ids$idno)
  
  Traits <- dplyr::bind_rows(Traits_E1, Traits_E5, Traits_E6)
  
  common_keys <- intersect(names(Traits), names(cleaned_proteins))
  
  protein_and_MIND_and_cov_ids <- Traits |>
    dplyr::full_join(cleaned_proteins, by = common_keys) |>
    dplyr::filter(idno %in% protein_ids$idno) |>
    dplyr::filter(!is.na(fa) | !is.na(wmh) | !is.na(epvs) | !is.na(mb_present)) |>
    dplyr::filter(!is.na(age) & !is.na(egfr) & !is.na(BMI) & !is.na(smoking) & !is.na(diabetes) & !is.na(edu) &
                    !is.na(gender) & !is.na(race) & !is.na(E4) & !is.na(site) & !is.na(ldl) & !is.na(AFprevalent) &
                    !is.na(MIprevalent) & !is.na(CHFprevalent) & !is.na(sbp)
    ) |>
    dplyr::select(idno, Exam)

      
  Traits <-  Traits |>
    dplyr::full_join(cleaned_proteins, by = common_keys) |>
    dplyr::filter(idno %in% protein_ids$idno) |>
    dplyr::filter(!is.na(fa) | !is.na(wmh) | !is.na(epvs) | !is.na(mb_present)) |>
    dplyr::filter(!is.na(age) & !is.na(egfr) & !is.na(BMI) & !is.na(smoking) & !is.na(diabetes) & !is.na(edu) &
                    !is.na(gender) & !is.na(race) & !is.na(E4) & !is.na(site) & !is.na(ldl) & !is.na(AFprevalent) &
                    !is.na(MIprevalent) & !is.na(CHFprevalent) & !is.na(sbp))

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
    Traits_table = Traits,
    protein_ids = protein_ids,
    protein_and_MIND_ids = protein_and_MIND_ids,
    protein_and_MIND_and_cov_ids = protein_and_MIND_and_cov_ids
    )
  
  
}