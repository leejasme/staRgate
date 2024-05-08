cleanMarkerNames = function(dat){
  #' Clean marker naming conventions
  #'
  #' Purpose: cleans column names to consistent formatting of
  #' PERCENT_<marker1>_POS/NEG_OF_<marker2>_POS/NEG
  #' to match the naming in the count matrix created from the pipeline
  #'
  #' Inputs:
  #' @param dat dataframe: markers in columns and samples in rows
  #' @return dataframe with updated column names
  #' @export


  # Check dat is of class data.frame
  if(!(inherits(dat, "data.frame"))){
    rlang::abort(message = "Error: `dat` must be of data.frame class")
  }

  # First grab the current colnames of dat
  current_cnames = colnames(dat)

  new_cnames =
    current_cnames %>%
    # Take care of random types
    stringr::str_replace_all(., c("\\% " = "PERCENT_",
                                  # 2022-12-27 Add Temra, Tn, Tcm, Tcm
                                  "TEMRA" = "CCR7_NEG_CD45RA_POS",
                                  "TCM" = "CCR7_POS_CD45RA_NEG",
                                  "TEM" = "CCR7_NEG_CD45RA_NEG",
                                  "TN" = "CCR7_POS_CD45RA_POS",
                                  # Correct for hyphens first bc they'll be seen as a "negative"
                                  "LAG-3" = "LAG3", "Ki-67" = "KI67", "KI-67" = "KI67",
                                  "PD-1" = "PD1", "Tim3" = "TIM3",
                                  "\\+" = "_POS_", "\\-" = "_NEG_", #"CD4_" = "CD4_POS", "CD8_" = "CD8_POS",
                                  # The CD4 and CD8 in denom need to be changed for indicator vars
                                  "CD4$" = "CD4_POS", "CD8$" = "CD8_POS", "CD3$" = "CD3_POS",
                                  "FoxP3" = "FOXP3",
                                  "FOX_P3" = "FOXP3",
                                  # Change Treg
                                  # 2023-04-10 Confirmed with MA: CD4 Treg is %CD4+FoxP3+ of CD3 Tcells
                                  "CD4 Treg" = "FOXP3_POS_CD4_POS_OF_CD3_POS"#,
                                  #"CD4_TREG" = "FOXP3_POS_OF_CD4_POS"
    )) %>%
    # More general janitor cleaning
    # If start with janitor, it seems to remove the +/- in the names
    janitor::make_clean_names(case = "all_caps")

  # Check if any are missing the _POS_ or _NEG_ in num
  # But this selects the non percentage cols as well
  # Dont want the MFI cols
  # exclude ratio for now
  aux_bad_num =
    !grepl("_POS_OF_|_NEG_OF_|MFI|RATIO", new_cnames)

  aux_bad_num_rename =
    new_cnames[aux_bad_num] %>%
    stringr::str_replace_all(., c("_OF_" = "_POS_OF_"))

  # replace the additional bad numerators
  new_cnames_renamed =
    replace(new_cnames,
            aux_bad_num,
            aux_bad_num_rename)

  new_dat =
    dat %>%
    dplyr::rename(stats::setNames(current_cnames,
                                  new_cnames_renamed))

  return(new_dat)

}
