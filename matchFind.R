

matchFind<-function (df1, df2, org1, org2, msa, ort = TRUE){
  
  warnOpt <- getOption("warn")
  options(warn = -1)
  
  # Find index of msa having protein ids which both have variants
  cat("\n Finding orthologous variants.. \n")
  df1_id<-msa$index[msa[[paste0(org1,"_ID")]] %in% unique(df1$Refseq_ID)]
  df2_id<-msa$index[msa[[paste0(org2,"_ID")]] %in% unique(df2$Refseq_ID)]
  common_id<-as.numeric(intersect(df1_id, df2_id))
  df1_df2_ort_list<-vector(mode = "list", length = length(common_id))
  
  j<-0
  #if(max(common_id) > min(common_id)){
  pb <- pbapply::timerProgressBar(max = length(common_id), style = 2)
  #}
  # Loop over only common alignments
  for (i in common_id){
    # aa positions of variants from selected protein
    df1_aapos<-df1$aapos[df1$Refseq_ID == msa[[paste0(org1,"_ID")]][msa$index == i]]
    df2_aapos<-df2$aapos[df2$Refseq_ID == msa[[paste0(org2,"_ID")]][msa$index == i]]
    
    # reference aa of variants
    df1_from<-df1$from[df1$Refseq_ID == msa[[paste0(org1,"_ID")]][msa$index == i]]
    df2_from<-df2$from[df2$Refseq_ID == msa[[paste0(org2,"_ID")]][msa$index == i]]
    
    # converted aa
    df1_to<-df1$to[df1$Refseq_ID == msa[[paste0(org1,"_ID")]][msa$index == i]]
    df2_to<-df2$to[df2$Refseq_ID == msa[[paste0(org2,"_ID")]][msa$index == i]]
    
    # tables for calculating aa conservation
    df1_conservation<-data.table::data.table(from = df1_from, to = df1_to, aapos = df1_aapos)
    df1_conservation<-dplyr::left_join(df1_conservation, aa_conservation, by = "from")
    df1_conservation<-dplyr::mutate(df1_conservation, conservation = ifelse(stringi::stri_detect(df1_conservation$to.y, fixed = stringr::fixed(df1_conservation$to.x)), 1, 0))
    
    df2_conservation<-data.table::data.table(from = df2_from, to = df2_to, aapos = df2_aapos)
    df2_conservation<-dplyr::left_join(df2_conservation, aa_conservation, by = "from")
    df2_conservation<-dplyr::mutate(df2_conservation, conservation = ifelse(stringi::stri_detect(df2_conservation$to.y, fixed = stringr::fixed(df2_conservation$to.x)), 1, 0))
    
    # get sequence of selected proteins
    df1_seq<-msa[[paste0(org1,"_seq")]][msa$index == i]
    df2_seq<-msa[[paste0(org2,"_seq")]][msa$index == i]
    
    # get real (trimmed) sequence of selected proteins
    df1_real_seq<-stringi::stri_replace_all_regex(df1_seq, "-", "")
    df2_real_seq<-stringi::stri_replace_all_regex(df2_seq, "-", "")
    
    # find index of aa and combine these
    df1_seq_ind<-data.table::as.data.table(stringi::stri_locate_all_regex(df1_seq, "[A-Z]")[[1]])[[1]]
    df2_seq_ind<-data.table::as.data.table(stringi::stri_locate_all_regex(df2_seq, "[A-Z]")[[1]])[[1]]
    df1_real_seq_ind<-data.table::as.data.table(stringi::stri_locate_all_regex(df1_real_seq, "[A-Z]")[[1]])[[1]]
    df2_real_seq_ind<-data.table::as.data.table(stringi::stri_locate_all_regex(df2_real_seq, "[A-Z]")[[1]])[[1]]
    
    df1_all_seq_ind<-data.table::data.table(seq_ind = df1_seq_ind,  aapos = df1_real_seq_ind) %>%
      dplyr::filter(df1_real_seq_ind %in% df1_aapos)
    
    df2_all_seq_ind<-data.table::data.table(seq_ind = df2_seq_ind, aapos = df2_real_seq_ind) %>%
      dplyr::filter(df2_real_seq_ind %in% df2_aapos)
    
    df1_all_seq_ind<-df1_all_seq_ind[df1_all_seq_ind$seq_ind %in% df2_all_seq_ind$seq_ind,]
    df2_all_seq_ind<-df2_all_seq_ind[df2_all_seq_ind$seq_ind %in% df1_all_seq_ind$seq_ind,]
    
    # merge by aa positions, get only aa having same reference aa
    if (length(df1_all_seq_ind[[1]]) != 0){
      
      df1_all_seq_ind<-dplyr::left_join(df1_all_seq_ind, df1_conservation, by = "aapos") %>%
        unique() %>%
        dplyr::select(!to.y)
      
      df2_all_seq_ind<-dplyr::left_join(df2_all_seq_ind, df2_conservation, by = "aapos") %>%
        unique() %>%
        dplyr::select(!to.y)
      
      all_seq_ind<-dplyr::left_join(df1_all_seq_ind, df2_all_seq_ind, by = "seq_ind") %>%
        dplyr::filter(from.x == from.y)
      # all_seq_ind$conservation_score<-NA
      # if(length(all_seq_ind$seq_ind)>0){
      #   cons_str_1<-substring(df1_seq, as.numeric(all_seq_ind$seq_ind)-5, as.numeric(all_seq_ind$seq_ind)+5)
      #   cons_str_2<-substring(df2_seq, as.numeric(all_seq_ind$seq_ind)-5, as.numeric(all_seq_ind$seq_ind)+5)
      #   for(k in 1:length(all_seq_ind$seq_ind)){
      #     cons_list<-stringi::stri_match(str_split(cons_str_1[k], "", simplify = TRUE), regex = str_split(cons_str_2[k], "", simplify = TRUE))
      #     all_seq_ind$conservation_score[k]<-length(na.omit(cons_list))/length(cons_list)
      #   }
      # }
      
      # should dismiss conserved aa changes
      # if(ort == TRUE){
      #   all_seq_ind<-dplyr::filter(all_seq_ind, conservation.x == conservation.y)
      # }
      
      # combine both organisms
      j<-j+1
      df1_df2_ort<-setNames(data.table::data.table(msa[[paste0(org1,"_ID")]][msa$index == i],
                                                   all_seq_ind$aapos.x,
                                                   all_seq_ind$from.x,
                                                   all_seq_ind$to.x.x,
                                                   msa[[paste0(org2,"_ID")]][msa$index == i],
                                                   all_seq_ind$aapos.y,
                                                   all_seq_ind$from.y,
                                                   all_seq_ind$to.x.y,
                                                   i
      ),
      c(paste0(org1,"_ID"), paste0(org1,"_aapos"), paste0(org1,"_from"), paste0(org1,"_to"),
        paste0(org2,"_ID"), paste0(org2,"_aapos"), paste0(org2,"_from"), paste0(org2,"_to"), "msa_id"
      ))
      
      df1_df2_ort<-df1_df2_ort[df1_df2_ort[[paste0(org1,"_aapos")]] != "",]
      df1_df2_ort_list[[j]]<-df1_df2_ort
    }
    #if(max(common_id) > min(common_id)){
    n<-which(common_id == i)
    pbapply::setTimerProgressBar(pb, n)
    #}
  }
  df1_df2_ort<-data.table::rbindlist(df1_df2_ort_list)
  df1_df2_ort<-unique(df1_df2_ort)
  on.exit(options(warn = warnOpt))
  return(df1_df2_ort)
}

# all_seq_ind$conservation_score
# conservation_score
