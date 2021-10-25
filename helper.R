
# Script for getting current Wormbase version number

library(RCurl)
library(XML)
library(stringr)

url<-"ftp://ftp.wormbase.org/pub/wormbase/releases/"
page<-getURL(url, verbose=TRUE, ftp.use.epsv=TRUE, dirlistonly = TRUE)
links<-getHTMLLinks(page)
allws<-str_split(page, "\n", simplify = TRUE)
allws<-allws[grep("WS", allws)]

ws<-allws[length(allws)]
write.table(ws, "wb_version", quote = FALSE, row.names = FALSE, col.names = FALSE)
