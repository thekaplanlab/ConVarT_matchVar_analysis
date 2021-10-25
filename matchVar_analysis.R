

library(data.table)
library(RMySQL)
library(dplyr)
library(pbapply)
library(tidyr)

source("matchFind.R")

mutagen<-fread("mutagen_analysis.txt")
mouse<-fread("mouse_analysis.txt")
celegans<-fread("celegans_analysis.txt")
clinvar<-fread("clinvar_analysis.txt")
topmed<-fread("topmed_analysis.txt")
cosmic<-fread("cosmic_analysis.txt")
gnomad<-fread("gnomad_analysis.txt")

con <- dbConnect(RMySQL::MySQL(), host = 'localhost', user = 'user', password = 'password', dbname = 'convart', port = 3306)
msa<-tbl(con, "msa") %>%
  collect()

msa<-msa %>%
  filter(grepl('musculus', fasta)) %>%
  filter(grepl('sapiens', fasta)) %>%
  filter(grepl('elegans', fasta))

msa<-msa %>%
  separate(fasta, c("Human", "Mouse","C.elegans"), sep = "\n>")

# Remove "\n" to get single sequence text
msa<-as.data.frame(apply(msa[,2:4], 2, function(x) gsub("\n", "", x, perl = TRUE)))
msa<-as.data.frame(apply(msa, 2, function(x) gsub(">", "", x, perl = TRUE)))

# If some rows have mixed data (e.g Mouse column having C. elegans sequence data), swap them.
msa1<-msa
msa1$Mouse[grep("elegans", msa1$Mouse)]<-msa1$C.elegans[grep("elegans", msa1$Mouse)]
msa1$C.elegans[grep("musculus", msa$C.elegans)]<-msa$Mouse[grep("musculus", msa$C.elegans)]

# Two of "Human", "Mouse" and "C_elegans" names will be used in place of org1 and org2 in matchFinder function
msa1<-msa1 %>%
  separate(Human, c("Human_ID","Human_seq"), " \\[Homo sapiens\\]") %>%
  separate(Mouse, c("Mouse_ID","Mouse_seq"), " \\[Mus musculus\\]") %>%
  separate(C.elegans, c("C_elegans_ID","C_elegans_seq"), " \\[Caenorhabditis elegans\\]")

# set rownames as index numbers
msa1$index<-rownames(msa1)


celegans_topmed_ort<-matchFind(df1 = celegans, df2 = topmed, 
                               org1 = "C_elegans", org2 = "Human", 
                               msa = msa1)
celegans_topmed_ort<-unique(celegans_topmed_ort)
write.table(celegans_topmed_ort, "celegans_topmed_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


celegans_clinvar_ort<-matchFind(df1 = celegans, df2 = clinvar, 
                               org1 = "C_elegans", org2 = "Human", 
                               msa = msa1)
celegans_clinvar_ort<-unique(celegans_clinvar_ort)
write.table(celegans_clinvar_ort, "celegans_clinvar_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


celegans_cosmic_ort<-matchFind(df1 = celegans, df2 = cosmic, 
                                org1 = "C_elegans", org2 = "Human", 
                                msa = msa1)
celegans_cosmic_ort<-unique(celegans_cosmic_ort)
write.table(celegans_cosmic_ort, "celegans_cosmic_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


celegans_gnomad_ort<-matchFind(df1 = celegans, df2 = gnomad, 
                                org1 = "C_elegans", org2 = "Human", 
                                msa = msa1)
celegans_gnomad_ort<-unique(celegans_gnomad_ort)
write.table(celegans_gnomad_ort, "celegans_gnomad_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


mutagen_topmed_ort<-matchFind(df1 = mutagen, df2 = topmed, 
                              org1 = "Mouse", org2 = "Human", 
                              msa = msa1)
mutagen_topmed_ort<-unique(mutagen_topmed_ort)
write.table(mutagen_topmed_ort, "mutagen_topmed_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


mutagen_clinvar_ort<-matchFind(df1 = mutagen, df2 = clinvar, 
                              org1 = "Mouse", org2 = "Human", 
                              msa = msa1)
mutagen_clinvar_ort<-unique(mutagen_clinvar_ort)
write.table(mutagen_clinvar_ort, "mutagen_clinvar_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


mutagen_cosmic_ort<-matchFind(df1 = mutagen, df2 = cosmic, 
                              org1 = "Mouse", org2 = "Human", 
                              msa = msa1)
mutagen_cosmic_ort<-unique(mutagen_cosmic_ort)
write.table(mutagen_cosmic_ort, "mutagen_cosmic_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


mutagen_gnomad_ort<-matchFind(df1 = mutagen, df2 = gnomad, 
                              org1 = "Mouse", org2 = "Human", 
                              msa = msa1)
mutagen_gnomad_ort<-unique(mutagen_gnomad_ort)
write.table(mutagen_gnomad_ort, "mutagen_gnomad_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


mouse_clinvar_ort<-matchFind(df1 = mouse, df2 = clinvar, 
                            org1 = "Mouse", org2 = "Human", 
                            msa = msa1)
mouse_clinvar_ort<-unique(mouse_clinvar_ort)
write.table(mouse_clinvar_ort, "mouse_clinvar_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


mouse_cosmic_ort<-matchFind(df1 = mouse, df2 = cosmic, 
                            org1 = "Mouse", org2 = "Human", 
                            msa = msa1)
mouse_cosmic_ort<-unique(mouse_cosmic_ort)
write.table(mouse_cosmic_ort, "mouse_cosmic_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")


mouse_gnomad_ort<-matchFind(df1 = mouse, df2 = gnomad, 
                            org1 = "Mouse", org2 = "Human", 
                            msa = msa1)
mouse_gnomad_ort<-unique(mouse_gnomad_ort)
write.table(mouse_gnomad_ort, "mouse_gnomad_matchvar.txt", row.names = FALSE, quote = FALSE, sep = "\t")

