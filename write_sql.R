
# This script writes variant tables to MySQL database
# After generating initial tables, you need to uncomment to update existing tables.

library(data.table)
library(RMySQL)
library(dplyr)

mutagen<-fread("mutagen_cron.txt")
mouse<-fread("mouse_cron.txt")
celegans<-fread("celegans_cron.txt")
clinvar<-fread("clinvar_cron.txt")
topmed<-fread("topmed_cron.txt")

mouseall<-rbind(mutagen, mouse) %>%
  distinct
mouseall$id<-row.names(mouseall)
mouseall$id<-as.numeric(mouseall$id)

con <- dbConnect(RMySQL::MySQL(), host = 'localhost', user = 'user', password = 'password', dbname = 'convart', port = 3306)

# dbSendQuery(con, 'DROP TABLE `mouse_variants_old`')
# dbSendQuery(con, 'RENAME TABLE `convart`.`mouse_variants_new`
#                   TO `convart`.`mouse_variants_old`;')
dbSendQuery(con, 'CREATE TABLE `convart`.`mouse_variants_new` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `variation_id` bigint(20) DEFAULT NULL,
  `gene_name` varchar(255) DEFAULT NULL,
  `phenotype` varchar(255) DEFAULT NULL,
  `aa_change` varchar(255) DEFAULT NULL,
  `pos` bigint(20) DEFAULT NULL,
  `refseq_id` varchar(255) DEFAULT NULL,
  `variant_type` varchar(255) DEFAULT NULL,
  `variation` varchar(255) DEFAULT NULL,
  `polyphen_score` varchar(255) DEFAULT NULL,
  `predicted_effect` varchar(255) DEFAULT NULL,
  `source` varchar(255) DEFAULT NULL,
  PRIMARY KEY (id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;')

dbWriteTable(con, "mouse_variants_new", mouseall, row.names = FALSE, overwrite = FALSE, append = TRUE)


# C. elegans ----

celegans<-fread("celegans_cron.txt")

# dbSendQuery(con, 'DROP TABLE `celegans_variants_old`')
# dbSendQuery(con, 'RENAME TABLE `convart`.`celegans_variants`
#                   TO `convart`.`celegans_variants_old`;')

dbSendQuery(con, 'CREATE TABLE `convart`.`celegans_variants` (
  `primary_id` int(11) NOT NULL,
  `ids` varchar(255) DEFAULT NULL,
  `pos` bigint(20) DEFAULT NULL,
  `refseq_id` varchar(255) DEFAULT NULL,
  `aa_change` varchar(255) DEFAULT NULL,
  `gene_name` varchar(255) DEFAULT NULL,
  `variant_type` varchar(255) DEFAULT NULL,
  `phenotype` varchar(255) DEFAULT NULL,
  `SIFT` varchar(255) DEFAULT NULL,
  `PolyPhen` varchar(255) DEFAULT NULL,
  `source` varchar(255) DEFAULT NULL,
  PRIMARY KEY (primary_id),
  KEY `celegans_refseq_id_index` (`refseq_id`(191))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;')

dbWriteTable(con, "celegans_variants", celegans, row.names = FALSE, overwrite = FALSE, append = TRUE)


# Clinvar ----

clinvar<-fread("clinvar_cron.txt")

# dbSendQuery(con, 'DROP TABLE `clinvar_old`')
# dbSendQuery(con, 'RENAME TABLE `convart`.`clinvar`
#                   TO `convart`.`clinvar_old`;')

dbSendQuery(con, 'CREATE TABLE `convart`.`clinvar` (
  `clinvar_id` bigint(20) NOT NULL,
  `gene_id` bigint(20) DEFAULT NULL,
  `allele_id` bigint(20) DEFAULT NULL,
  `symbol` varchar(255) DEFAULT NULL,
  `rs_number` varchar(255) DEFAULT NULL,
  `rcv_accession` varchar(255) DEFAULT NULL,
  `variation_id` bigint(20) DEFAULT NULL,
  `variant_type` varchar(255) DEFAULT NULL,
  `phenotype_ids` varchar(255) DEFAULT NULL,
  `name` varchar(255) DEFAULT NULL,
  `clinical_significance` varchar(255) DEFAULT NULL,
  `review_status` varchar(255) DEFAULT NULL,
  `phenotypes` varchar(255) DEFAULT NULL,
  `nm_id` text,
  `variation` varchar(255) DEFAULT NULL,
  `position` bigint(20) DEFAULT NULL,
  `np_id` varchar(255) DEFAULT NULL,
  PRIMARY KEY (clinvar_id),
  UNIQUE KEY `clinvar_id` (`clinvar_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;')

dbWriteTable(con, "clinvar", clinvar, row.names = FALSE, overwrite = FALSE, append = TRUE)


# Topmed ----

topmed<-fread("topmed_cron.txt")

# dbSendQuery(con, 'DROP TABLE `dbsnp`')
# dbSendQuery(con, 'RENAME TABLE `convart`.`dbsnp_new`
#                   TO `convart`.`dbsnp`;')

dbSendQuery(con, 'CREATE TABLE `convart`.`dbsnp_new` (
  `id` int(20) NOT NULL AUTO_INCREMENT,
  `ensm_protein_id` varchar(255) DEFAULT NULL,
  `rs_id` varchar(50) DEFAULT NULL,
  `ensm_gene_id` varchar(255) DEFAULT NULL,
  `ensm_transcript_id` varchar(80) DEFAULT NULL,
  `variant_type` varchar(255) DEFAULT NULL,
  `position` int(20) DEFAULT NULL,
  `variation` varchar(50) DEFAULT NULL,
  `sift_score` varchar(255) DEFAULT NULL,
  `polyphen_score` varchar(255) DEFAULT NULL,
  PRIMARY KEY (id),
  UNIQUE KEY `canonical_transcript_2` (`rs_id`,`ensm_transcript_id`,`variation`) USING BTREE,
  KEY `ensm_transcript_id` (`ensm_transcript_id`),
  KEY `rs_id` (`rs_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8mb4;')

dbWriteTable(con, "dbsnp_new", topmed, row.names = FALSE, overwrite = FALSE, append = TRUE)

dbDisconnect(con)

