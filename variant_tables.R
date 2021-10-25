

library(data.table)
library(stringr)
library(dplyr)
library(seqinr)
library(biomaRt)
library(Biostrings)
library(tidyr)
library(RMySQL)
library(ontologyIndex)


mart<-useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
mouse_ens<-getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
                 mart = mart, uniqueRows = TRUE)


mouse_ens_prot<-read.fasta(file = "Mus_musculus.GRCm39.pep.all.fa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
  unlist() %>%
  as.data.frame()
mouseseq<-read.fasta(file = "GCF_000001635.27_GRCm39_protein.faa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
  unlist() %>%
  as.data.frame()


mouse_ens_prot$ens<-rownames(mouse_ens_prot)
mouse_ens_prot1<-str_split(mouse_ens_prot$ens, "\\.", simplify = TRUE)
mouse_ens_prot$ens<-mouse_ens_prot1[,1]
colnames(mouse_ens_prot)<-c("seq","ensembl_peptide_id")

mouse_all<-merge(mouse_ens, mouse_ens_prot, by = "ensembl_peptide_id")

mouseseq$Refseq_ID<-rownames(mouseseq)
tempmouseseq<-str_split_fixed(mouseseq$Refseq_ID, "\\.", 2)
mouseseq$Refseq_ID<-tempmouseseq[,1]
colnames(mouseseq)[1]<-"seq"

mouse_all<-left_join(mouse_all, mouseseq, by = "seq")
colnames(mouse_all)[3]<-"transcript_id"


# Mutagenetix variants ----

datenow<-"2021-08-14"
idnow1<-1045
idnow2<-1046

numbOfDays<-as.numeric(difftime(Sys.Date(), as.Date(datenow)), units = "days")
mod<-numbOfDays %% 7
modid<-numbOfDays %/% 7

date1<-Sys.Date() - mod
date1<-gsub("-","",date1)

id1<-idnow1 + modid*3
id2<-idnow2 + modid*3

download.file(paste0("https://mutagenetix.utsouthwestern.edu/linkedfile/linked_file.cfm?id=", id1, "&fn=mutagenetix%5Fincidental%5Fmutations%5F", date1, "%2Etxt"), destfile = paste0(file.path(getwd()), "/mutagenetix_incidental.txt"),
              method = "curl", quiet = FALSE, extra='-k')

download.file(paste0("https://mutagenetix.utsouthwestern.edu/linkedfile/linked_file.cfm?id=", id2, "&fn=mutagenetix%5Fphenotypic%5Fmutations%5F", date1, "%2Etxt"), destfile = paste0(file.path(getwd()), "/mutagenetix_phenotypic.txt"),
              method = "curl", quiet = FALSE, extra='-k')


mut_inc<-fread("mutagenetix_incidental.txt")
mut_inc<-mut_inc[!is.na(mut_inc$amino_acid_position),]
mut_phe<-fread("mutagenetix_phenotypic.txt")
mut_phe<-mut_phe[!is.na(mut_phe$amino_acid_position),]

if(length(mut_inc[[1]]) > 0){
  mut_inc<-mut_inc %>%
    dplyr::select(mutagenetix_id, gene_symbol, transcript_id, amino_acid_position, ref_amino_acid, var_amino_acid, polyphen_score, predicted_effect,
                  mgi_accession_id, mutation_type)
  
  mut_phe<-mut_phe %>%
    dplyr::select(phenotypic_id, gene_symbol, transcript_id, amino_acid_position, ref_amino_acid, var_amino_acid, polyphen_score, predicted_effect,
                  mgi_accession_id, mutation_type, phenotypes)
  
  mut_inc$phenotypes<-NA
  mut_inc$variant_type<-"incidental"
  mut_phe$variant_type<-"phenotypic"
  
  colnames(mut_phe)[1]<-"mutagenetix_id"
  
  mutagen<-rbind(mut_inc, mut_phe)
  mutagen$aa_change<-paste(mutagen$ref_amino_acid, "to", mutagen$var_amino_acid)
  mutagen$variation<-paste0(mutagen$ref_amino_acid, mutagen$amino_acid_position, mutagen$var_amino_acid)
  
  mutagen<-left_join(mutagen, mouse_all[,c(3,6)], by = "transcript_id")
  mutagen$Refseq_ID<-ifelse(is.na(mutagen$Refseq_ID), mutagen$transcript_id, mutagen$Refseq_ID)
  
  colnames(mutagen)[c(1,2,4,11,15)]<-c("variation_id","gene_name","pos","phenotype","refseq_id")
  mutagen$source<-"Mutagenetix"
  
  mutagen<-mutagen %>%
    dplyr::select(variation_id, gene_name, phenotype, aa_change, pos, refseq_id, variant_type, variation, polyphen_score, 
                  predicted_effect, source)
  
  write.table(mutagen, "mutagen_cron.txt", quote = FALSE, row.names = FALSE, sep = "\t")
}

# APF ----

datenow<-"2021-08-15"
numbOfDays<-as.numeric(difftime(Sys.Date(), as.Date(datenow)), units = "days")
mod<-numbOfDays %% 7
modid<-numbOfDays %/% 7

date1<-Sys.Date() - mod
date1<-gsub("-","",date1)

apf<-fread("variant_files/apf1.txt")

download.file(paste0("https://databases.apf.edu.au/download/SNVs",date1,".txt"), destfile = paste0(file.path(getwd()), "/apf1.txt"),
              method = "curl", quiet = FALSE, extra='-k')


apf_phe<-apf[apf$Phenotype != "",]
apf_phe<-apf_phe[!is.na(apf_phe$`Amino Acid Position`),]

apf1<-fread("apf1.txt")
apf1_phe<-apf1[apf1$Phenotype != "",]
apf1_phe<-apf1_phe[!is.na(apf1_phe$`Amino Acid Position`),]

apf_phe_recent<-anti_join(apf1_phe, apf_phe, by = "Id")

library(data.table)
library(rvest)

if (length(apf_phe_recent[[1]]) != 0){
  genelist<-apf_phe_recent
  genelist$aa_position<-NA
  genelist$ccdsid<-NA
  txtt<-c()
  pb <- pbapply::timerProgressBar(min = 1, max = length(genelist[[1]]), style = 2)
  for (i in 1:length(genelist[[1]])){
    geneid<-genelist$Id[i]
    url<-paste0("https://databases.apf.edu.au/mutations/snpRow/show/", geneid, "?assembly=GRCm38")
    webpage <- read_html(url)
    ttll<-html_nodes(webpage, "span.property-value[aria-labelledby=\"vepAminoAcidPosition-label\"]")
    ccds<-html_nodes(webpage, "a.ext_link[target=\"_tab\"]")
    aaposition<-html_text(ttll)
    ccdsposition<-html_text(ccds)
    ccdsposition<-ccdsposition[grep("CCDS", ccdsposition)]
    if (length(aaposition != 0)){
      genelist$aa_position[which(genelist$Id == geneid)]<-paste(aaposition, collapse = ",")
      genelist$ccdsid[which(genelist$Id == geneid)]<-paste(ccdsposition, collapse = ",")
    }
    else {genelist$aa_position[which(genelist$Id == geneid)]<-0
    genelist$ccdsid[which(genelist$Id == geneid)]<-0}
    txtt<-c(txtt, i)
    pbapply::setTimerProgressBar(pb, i)
  }
  
  mouse<-genelist
  mouse1<-str_split_fixed(mouse$`Amino Acid Change`, "->", 2)
  mouse$from<-mouse1[,1]
  mouse$to<-mouse1[,2]
  
  mouse$to[mouse$to == "Stop"]<-"*"

  seq1<-read.fasta(file = "CCDS_protein.current.faa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
    unlist() %>%
    as.data.frame()
  seq1$ens<-rownames(seq1)
  ss<-str_split(seq1$ens, "\\|", simplify = TRUE)
  ss<-data.frame(ss, stringsAsFactors = FALSE)
  ss1<-ss[,1]
  ss1<-str_split(ss1, "\\.", simplify = TRUE)
  seq1$ens<-ss1[,1]
  
  ccdsids<-fread("CCDS.current.txt")
  ccdsids<-ccdsids[ccdsids$match_type == "Identical"]
  ccdsids1<-str_split(ccdsids$ccds_id, "\\.", simplify = TRUE)
  ccdsids$ccds_id<-ccdsids1[,1]
  
  mouseseq<-read.fasta(file = "GCF_000001635.27_GRCm39_protein.faa", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
    unlist() %>%
    as.data.frame()
  mouseseq$ens<-rownames(mouseseq)
  tempmouseseq<-str_split_fixed(mouseseq$ens, "\\.", 2)
  mouseseq$ens<-tempmouseseq[,1]
  
  mouseall<-merge(mouseseq, seq1, by = ".")
  mouseall<-mouseall[,2:3]
  colnames(mouseall)<-c("refseq_id","CCDS_ID")
  
  mouse<-mouse %>% dplyr::select(Id, `Gene Name`, Phenotype, `Amino Acid Position`, ccdsid, from, to, `Polyphen Score`, `Polyphen Prediction`)
  colnames(mouse)<-c("variation_id","gene_name","phenotype","pos","CCDS_ID","from","to","polyphen_score","predicted_effect")
  mouse$aa_change<-paste(mouse$from, "to", mouse$to)
  mouse$variation<-paste0(mouse$from, mouse$pos, mouse$to)
  mouse$source<-"APF"
  mouse$variant_type<-"phenotypic"
  
  mouse<-mouse %>% 
    mutate(CCDS_ID = strsplit(as.character(CCDS_ID), ",")) %>%
    unnest(CCDS_ID)
  tempmouse<-str_split_fixed(mouse$CCDS_ID, "\\.", 2)
  mouse$CCDS_ID<-tempmouse[,1]
  mouse<-mouse[mouse$CCDS_ID != "NO_CCDS",]
  mouse<-left_join(mouse, mouseall, by = "CCDS_ID")
  
  mouse<-mouse[!is.na(mouse$pos),]
  mouse<-unique(mouse)
  mouse$pos<-as.integer(mouse$pos)
  mouse<-mouse[mouse$to != "",]
  mouse<-mouse[mouse$from != "Disrupted splicing",]
  mouse<-mouse[mouse$from != "Stop",]
  
  mouse<-mouse[mouse$to != "",]
  mouse<-mouse[!is.na(mouse$refseq_id),]
  
  mouse<-mouse %>%
    dplyr::select(variation_id, gene_name, phenotype, aa_change, pos, refseq_id, variant_type, variation, polyphen_score, 
                  predicted_effect, source)
  
  mouse_old<-fread("mouse_cron.txt")
  
  mouse<-rbind(mouse_old, mouse)
  
  write.table(mouse, "mouse_cron.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  
}


# Wormbase ----

download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_protein.faa.gz",
              destfile = paste0(file.path(getwd()), "/GCF_000002985.6_WBcel235_protein.faa.gz"), 
              method = "curl", quiet = TRUE, extra = '-k')

download.file("ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/annotation/geneOtherIDs/c_elegans.PRJNA13758.current_development.geneOtherIDs.txt.gz",
              destfile = paste0(file.path(getwd()), "/c_elegans.PRJNA13758.current.geneOtherIDs.txt.gz"), 
              method = "curl", quiet = TRUE, extra = '-k')

download.file("ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/protein/c_elegans.PRJNA13758.current_development.protein.fa.gz",
              destfile = paste0(file.path(getwd()), "/c_elegans_protein_wb.fa.gz"), 
              method = "curl", quiet = TRUE, extra = '-k')

celegans_prot<-read.fasta(file = "c_elegans_protein_wb.fa.gz", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
  unlist() %>%
  as.data.frame()
celegans_prot$wb_id<-rownames(celegans_prot)
colnames(celegans_prot)[1]<-"seq"

celegans_np<-read.fasta(file = "GCF_000002985.6_WBcel235_protein.faa.gz", seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>%
  unlist() %>%
  as.data.frame()
celegans_np$Refseq_ID<-rownames(celegans_np)
tempceleseq<-str_split_fixed(celegans_np$Refseq_ID, "\\.", 2)
celegans_np$Refseq_ID<-tempceleseq[,1]
colnames(celegans_np)[1]<-"seq"

celegans_annot<-left_join(celegans_np, celegans_prot, by = "seq")
colnames(celegans_annot)[3]<-"Transcript_ID"


wormbase_aa<-fread("wormbase_variants_aachange.txt")

wormbase_aa$Change<-as.character(stringr::str_extract(wormbase_aa$V9, "(?<=aachange=).+?(?=;)"))
wormbase_aa$Transcript_ID<-as.character(stringr::str_extract(wormbase_aa$V9, "(?<=hgvsc=).+?(?=:)"))
wormbase_aa$Position<-as.character(stringr::str_extract(wormbase_aa$V9, "(?<=aa_position=).+?(?=;)"))
wormbase_aa$Wormbase_ID<-as.character(stringr::str_extract(wormbase_aa$V9, "(?<=variation=).+?(?=;)"))
wormbase_aa$SIFT<-as.character(stringr::str_extract(wormbase_aa$V9, "(?<=sift=).+?(?=;)"))
wormbase_aa$PolyPhen<-as.character(stringr::str_extract(wormbase_aa$V9, "(?<=polyphen=).+?(?=;)"))
wormbase_aa$Transcript_ID<-gsub("(.*)\\..*", "\\1", wormbase_aa$Transcript_ID)

wormbase_aa$Position<-as.numeric(wormbase_aa$Position)
wormbase_aa<-wormbase_aa[!is.na(wormbase_aa$Position),]
wormbase_aa<-wormbase_aa %>% filter(V3 == "point_mutation" | V3 == "SNP")
wormbase_aa<-wormbase_aa[nchar(wormbase_aa$Change)<4,]
wormbase_aa<-wormbase_aa[stri_sub(wormbase_aa$Change, 1, 1) != "*",]
colnames(wormbase_aa)[2:3]<-c("Source","Type")

wb_name<-fread("c_elegans.PRJNA13758.current.geneOtherIDs.txt.gz", header = FALSE)
wb_name$V4<-ifelse(wb_name$V4 == "", wb_name$V3, wb_name$V4)
wb_name<-wb_name[,c(1,3,4)]
colnames(wb_name)<-c("wb_gene_id","Transcript_ID","gene_name")
wb_name<-wb_name[wb_name$Transcript_ID != "",]

wormbase_aa$id<-wormbase_aa$Transcript_ID
wormbase_aa$Transcript_ID<-gsub("[^0-9]*$","",wormbase_aa$Transcript_ID)
wormbase_aa<-left_join(wormbase_aa, wb_name, by = "Transcript_ID")
wormbase_aa$Change<-gsub("/"," to ", wormbase_aa$Change)

wormbase<-wormbase_aa %>%
  dplyr::select(Wormbase_ID, id, Position, Change, SIFT, PolyPhen, Source, Type, gene_name) %>%
  distinct()
colnames(wormbase)[2]<-"Transcript_ID"
wormbase<-left_join(wormbase, celegans_annot[,2:3], by = "Transcript_ID")
wormbase$Refseq_ID<-ifelse(is.na(wormbase$Refseq_ID), wormbase$Transcript_ID, wormbase$Refseq_ID)
wormbase<-unique(wormbase)

download.file("ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/annotation/phenotype_association/c_elegans.PRJNA13758.current.phenotype_association.wb.gz",
              destfile = paste0(file.path(getwd()), "/c_elegans_phenotypes.wb.gz"), 
              method = "curl", quiet = TRUE, extra = '-k')
wb_phe<-fread("c_elegans_phenotypes.wb.gz")


download.file("https://raw.githubusercontent.com/obophenotype/c-elegans-phenotype-ontology/master/wbphenotype.obo",
              destfile = paste0(file.path(getwd()), "/wbphenotype.obo"), 
              method = "curl", quiet = TRUE, extra = '-k')


phenotypes<-get_ontology("wbphenotype.obo", extract_tags = "everything")

dfphe<-data.frame(id = phenotypes$id, def = phenotypes$def)

colnames(wb_phe)[5]<-"id"
wb_phe<-left_join(wb_phe, dfphe, by = "id")
wb_phe$V6<-ifelse(grepl("WB:WBVar", wb_phe$V6), wb_phe$V6, wb_phe$V8)

wb_phe1<-wb_phe %>% 
  mutate(V6 = strsplit(as.character(V6), "\\|")) %>%
  unnest(V6)
wb_phe2<-wb_phe1 %>%
  dplyr::select(V6, def) %>%
  distinct()

wb_phe2$V6<-gsub("WB:","",wb_phe2$V6)
colnames(wb_phe2)[1]<-"Wormbase_ID"
wormbase<-left_join(wormbase, wb_phe2, by = "Wormbase_ID")
wormbase$variant_type<-"Unknown"
wormbase$variant_type[!is.na(wormbase$def)]<-"Phenotypic"
wormbase$Source<-"Wormbase"
wormbase<-unique(wormbase)


con <- dbConnect(RMySQL::MySQL(), host = 'localhost', user = 'user', password = 'password', dbname = 'convart', port = 3306)
wormbase_old<-tbl(con, "celegans_variants") %>%
  collect()
dbDisconnect(con)

colnames(wormbase)[c(1,3,4,7,10,11)]<-c("ids","pos","aa_change","source","refseq_id","phenotype")
wormbase_old$phenotype[wormbase_old$phenotype == "NA"]<-NA
wormbase_old$phenotype[wormbase_old$phenotype == "NA\nNA"]<-NA

wormbase_old_new<-anti_join(wormbase_old, wormbase, by = c("ids","pos","refseq_id","aa_change","source"))
wormbase_old_new$phenotype[wormbase_old_new$phenotype == "No"]<-NA

wormbase<-wormbase %>% dplyr::select(ids, pos, refseq_id, aa_change, gene_name, variant_type, phenotype, SIFT, PolyPhen, source)
wormbase_old_new<-wormbase_old_new %>% dplyr::select(-primary_id)

wormbase_all<-rbind(wormbase, wormbase_old_new)

wormbase_all1<-as.data.table(wormbase_all)
wormbase_all1<-wormbase_all1[, list(aa_change, gene_name, variant_type, SIFT, PolyPhen, source, phenotype = paste0(phenotype, collapse = "\n")), 
                             by = c("ids","pos","refseq_id")]
wormbase_all1<-unique(wormbase_all1)
wormbase_all1$phenotype<-gsub("\"", "", wormbase_all1$phenotype)
wormbase_all1$primary_id<-as.numeric(rownames(wormbase_all1))

write.table(wormbase_all1, "celegans_cron.txt", row.names = FALSE, sep = "\t")


# Clinvar ----

download.file("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/hgvs4variation.txt.gz",
              destfile = paste0(file.path(getwd()), "/hgvs4variation.txt.gz"), 
              method = "curl", quiet = TRUE, extra = '-k')

clinvar<-fread("hgvs4variation.txt.gz") %>%
  filter(ProteinExpression != "-") %>%
  filter(across(ProteinExpression, ~ grepl("NP_", .))) %>%
  filter(UsedForNaming == "Yes")
clinvar<-clinvar[,c(1,2,3,4,7,8,9,10)]

download.file("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
              destfile = paste0(file.path(getwd()), "/variant_summary.txt.gz"), 
              method = "curl", quiet = TRUE, extra = '-k')

clinvarx<-fread("variant_summary.txt.gz")
clinvarx<-clinvarx[clinvarx$Type == "single nucleotide variant",]
clinvarx<-clinvarx %>%
  dplyr::select(Name, ClinicalSignificance, PhenotypeList, VariationID, Type, `RS# (dbSNP)`,
                RCVaccession, PhenotypeIDS, ReviewStatus, Assembly)

clinvarxy<-clinvarx[clinvarx$Assembly == "GRCh37",]

clinvara<-left_join(clinvar, clinvarxy, by = "VariationID")
clinvara<-unique(clinvara)

clinvarb<-clinvara[!is.na(clinvara$Name),]

clinvarb$np_id<-str_extract(clinvarb$ProteinExpression, ".+?(?=\\.)")
clinvarb$nm_id<-str_extract(clinvarb$NucleotideExpression, ".+?(?=\\.)")
clinvarb<-clinvarb[,c(-5,-6,-7,-17)]

colnames(clinvarb)<-c("symbol","gene_id","variation_id","allele_id","variation","name","clinical_significance","phenotypes",
                      "type","rs_number","rcv_accession","phenotype_ids","review_status","np_id","nm_id")
clinvarb$position<-str_extract(clinvarb$variation, "\\d+")
clinvarb$rs_number<-paste0("rs", clinvarb$rs_number)

clinvarb$clinvar_id<-as.numeric(rownames(clinvarb))
colnames(clinvarb)[which(colnames(clinvarb) == "type")]<-"variant_type"

write.table(clinvarb, "clinvar_cron.txt", quote = FALSE, row.names = FALSE, sep = "\t")


# TOPMed ----

topmed<-read.table("topmed/all_output_clean1.txt", header = FALSE, sep = " ")
colnames(topmed)<-c("rs_id","ensm_gene_id","ensm_transcript_id","variant_type","position","variation","ensm_protein_id",
                    "sift_score","polyphen_score")
topmed<-as.data.table(topmed)
topmed<-topmed[nchar(topmed$variation)<4,]
topmed$position<-as.numeric(topmed$position)
topmed<-topmed[!is.na(topmed$position),]
topmed<-topmed[stri_sub(topmed$variation, 1, 1) != "*",]
topmed<-unique(topmed)
topmed<-topmed[!duplicated(topmed[,c(1,3,6)])]
topmed$id<-as.numeric(rownames(topmed))

write.table(topmed, "topmed_cron.txt", quote = FALSE, row.names = FALSE, sep = "\t")


