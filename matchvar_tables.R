
# This script prepare variant tables for matching variants

library(orthoVar)
library(tidyr)
library(data.table)
library(stringr)
library(seqinr)

source("orthoFind.R")

mutagen<-fread("mutagen_cron.txt", select = c(4,5,6))
mutagen$from<-str_split(mutagen$aa_change, " to ", simplify = TRUE)[,1]
mutagen$to<-str_split(mutagen$aa_change, " to ", simplify = TRUE)[,2]
mutagen<-mutagen[,-1]
colnames(mutagen)[c(1,2)]<-c("aapos","Refseq_ID")
mutagen<-unique(mutagen)
write.table(mutagen, "mutagen_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")

mouse<-fread("mouse_cron.txt", select = c(4,5,6))
mouse$from<-str_split(mouse$aa_change, " to ", simplify = TRUE)[,1]
mouse$to<-str_split(mouse$aa_change, " to ", simplify = TRUE)[,2]
mouse<-mouse[,-1]
colnames(mouse)[c(1,2)]<-c("aapos","Refseq_ID")
mouse<-unique(mouse)
write.table(mouse, "mouse_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")

celegans<-fread("celegans_cron.txt", select = c(2,3,4))
celegans$from<-str_split(celegans$aa_change, " to ", simplify = TRUE)[,1]
celegans$to<-str_split(celegans$aa_change, " to ", simplify = TRUE)[,2]
celegans<-celegans[,-3]
colnames(celegans)[c(1,2)]<-c("aapos","Refseq_ID")
celegans<-unique(celegans)
write.table(celegans, "celegans_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")

clinvar<-fread("clinvar_cron.txt", select = c(5,14,16))
clinvar$variation<-gsub("p\\.", "", clinvar$variation)
clinvar$from<-str_split(clinvar$variation, "\\d+", simplify = TRUE)[,1]
clinvar$to<-str_split(clinvar$variation, "\\d+", simplify = TRUE)[,2]
clinvar$from<-ifelse(clinvar$from == "Ter", "Stp", clinvar$from)
clinvar$to<-ifelse(clinvar$to == "Ter", "Stp", clinvar$to)
clinvar$to<-ifelse(clinvar$to == "=", clinvar$from, clinvar$to)

clinvar$from<-a(clinvar$from)
clinvar$to<-a(clinvar$to)
clinvar<-clinvar[!is.na(clinvar$to),]
clinvar<-clinvar[clinvar$from != "*",]
clinvar<-clinvar[,-1]
clinvar<-unique(clinvar)
colnames(clinvar)[c(1,2)]<-c("Refseq_ID","aapos")
write.table(clinvar, "clinvar_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")


topmed<-fread("topmed_cron.txt", select = c(3,5,6))
topmed<-unique(topmed)
topmed$from<-str_split(topmed$variation, "/", simplify = TRUE)[,1]
topmed$to<-str_split(topmed$variation, "/", simplify = TRUE)[,2]
topmed$to<-ifelse(topmed$to == "", topmed$from, topmed$to)
topmed<-topmed[,-3]

ensp_np<-fread("est_esp_np_table_grch37.txt")
colnames(ensp_np)[4]<-"ensm_transcript_id"
ensp_np<-ensp_np[,c(1,2,4)]
ensp_np$ensembl_peptide_id[ensp_np$ensembl_peptide_id == ""]<-NA
ensp_np<-ensp_np[!is.na(ensp_np$ensm_transcript_id),]

topmedx<-left_join(topmed, ensp_np, by = "ensm_transcript_id")
topmedx$Refseq_ID<-ifelse(is.na(topmedx$refseq), topmedx$ensembl_peptide_id, topmedx$refseq)
topmedx1<-unique(topmedx[,c(2,3,4,9)])
topmedx1<-topmedx1[!is.na(topmedx1$Refseq_ID),]
topmedx1<-topmedx1[topmedx1$from != "-",]
topmedx1<-topmedx1[topmedx1$to != "-",]
colnames(topmedx1)[1]<-"aapos"

write.table(topmedx1, "topmed_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")

library(RMySQL)
con <- dbConnect(RMySQL::MySQL(), host = 'localhost', user = 'user', password = 'password', dbname = 'convart', port = 3306)
cosmicx<-tbl(con, "CosmicMutantExport") %>%
  dplyr::select(c(3,20)) %>%
  collect()
cosmic$mutation_aa<-gsub("p\\.", "", cosmic$mutation_aa)
cosmic<-cosmic[cosmic$mutation_aa != "?",]
cosmic<-unique(cosmic)
cosmic$from<-str_split(cosmic$mutation_aa, "\\d+", simplify = TRUE)[,1]
cosmic$to<-str_split(cosmic$mutation_aa, "\\d+", simplify = TRUE)[,2]
cosmic$aapos<-as.numeric(str_extract(cosmic$mutation_aa, "\\d+"))
colnames(cosmic)[1]<-"ensm_transcript_id"
cosmic$ensm_transcript_id<-str_split(cosmic$ensm_transcript_id, "\\.", simplify = TRUE)[,1]
cosmic<-left_join(cosmic, ensp_np, by = "ensm_transcript_id")
cosmic$Refseq_ID<-ifelse(!is.na(cosmic$refseq), cosmic$refseq, 
                         ifelse(!is.na(cosmic$ensembl_peptide_id), cosmic$ensembl_peptide_id, cosmic$ensm_transcript_id))
cosmic<-cosmic %>% dplyr::select(Refseq_ID, aapos, from, to)
cosmic<-cosmic[nchar(cosmic$from) == 1,]
cosmic<-cosmic[cosmic$from != "(",]
cosmic<-cosmic[nchar(cosmic$to) == 1,]
cosmic$to[cosmic$to == "="]<-cosmic$from[cosmic$to == "="]
cosmic<-cosmic[cosmic$to != "_",]
cosmic<-cosmic[cosmic$to != "?",]
cosmic<-cosmic[cosmic$to != "X",]
cosmic<-cosmic[cosmic$to != "Z",]
cosmic<-cosmic[!is.na(cosmic$Refseq_ID),]
cosmic<-unique(cosmic)


con <- dbConnect(RMySQL::MySQL(), host = 'localhost', user = 'user', password = 'password', dbname = 'convart', port = 3306)
gnomad<-tbl(con, "gnomad") %>%
  dplyr::select(canonical_transcript, variation, position) %>%
  collect()
dbDisconnect(con)
gnomad$variation<-gsub("p\\.", "", gnomad$variation)
gnomad$from<-str_split(gnomad$variation, "\\d+", simplify = TRUE)[,1]
gnomad$to<-str_split(gnomad$variation, "\\d+", simplify = TRUE)[,2]
gnomad$from[gnomad$from == "Ter"]<-"Stp"
gnomad$to[gnomad$to == "Ter"]<-"Stp"
gnomad<-unique(gnomad)
gnomad$from<-a(gnomad$from)
gnomad$to<-a(gnomad$to)
gnomad<-gnomad[!is.na(gnomad$from),]
gnomad<-gnomad[!is.na(gnomad$to),]
gnomad<-gnomad[!is.na(gnomad$position),]
colnames(gnomad)[c(1,3)]<-c("ensm_transcript_id","aapos")
gnomad<-left_join(gnomad, ensp_np, by = "ensm_transcript_id")
gnomad$Refseq_ID<-ifelse(is.na(gnomad$refseq), gnomad$ensembl_peptide_id, gnomad$refseq)
gnomad<-gnomad %>% dplyr::select(Refseq_ID, aapos, from, to)
gnomad<-unique(gnomad)


mutagen_old<-fread("mutagen_analysis.txt")
mouse_old<-fread("mouse_analysis.txt")
celegans_old<-fread("celegans_analysis.txt")
clinvar_old<-fread("clinvar_analysis.txt")
topmed_old<-fread("topmed_analysis.txt")
cosmic_old<-fread("cosmic_analysis.txt")
gnomad_old<-fread("gnomad_analysis.txt")


mutagen_new<-anti_join(mutagen, mutagen_old)
mouse_new<-anti_join(mouse, mouse_old)
celegans_new<-anti_join(celegans, celegans_old)
clinvar_new<-anti_join(clinvar, clinvar_old)
topmed_new<-anti_join(topmedx1, topmed_old)
cosmic_new<-anti_join(cosmic, cosmic_old)
gnomad_new<-anti_join(gnomad, gnomad_old)


write.table(mutagen, "mutagen_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(mouse, "mouse_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(celegans, "celegans_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(clinvar, "clinvar_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(topmedx1, "topmed_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(cosmic, "cosmic_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gnomad, "gnomad_analysis.txt", row.names = FALSE, quote = FALSE, sep = "\t")


write.table(mutagen_new, "mutagen_analysis_new.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(mouse_new, "mouse_analysis_new.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(celegans_new, "celegans_analysis_new.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(clinvar_new, "clinvar_analysis_new.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(topmedx1_new, "topmed_analysis_new.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(cosmic_new, "cosmic_analysis_new.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gnomad_new, "gnomad_analysis_new.txt", row.names = FALSE, quote = FALSE, sep = "\t")

