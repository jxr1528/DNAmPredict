# Downloading from TCGA CESC methylation and mutation data  ---------------
library(TCGA2STAT)
library(devtools)
library(CNTools)
library(BiocManager)

#modified merge function that returns transposed data (Rao et al. 2020)
#https://github.com/hhhxz305/DNAmPredict/blob/master/import_setupdata.0.R
OMICSBIND2<-function (dat1, dat2) 
{
    if (is.null(dat1) || is.null(dat2)) {
        message("Empy data object")
        return(NULL)
    }
    if ((class(dat1) == "data.frame" || class(dat1) == "matrix") & 
        (class(dat2) == "data.frame" || class(dat2) == "matrix")) {
        t1 <- sapply(colnames(dat1), function(s) nchar(s))
        t2 <- sapply(colnames(dat2), function(s) nchar(s))
        dat1.tumor <- dat1
        dat2.tumor <- dat2
        if (sum(t1 > 12) == ncol(dat1)) {
            dat1.type <- sapply(colnames(dat1), function(s) unlist(strsplit(s, 
                                                                            "-"))[4])
            dat1.tumor <- dat1[, grep("^01", dat1.type)]
            colnames(dat1.tumor) <- substr(colnames(dat1.tumor), 
                                           1, 12)
        }
        if (sum(t2 > 12) == ncol(dat2)) {
            dat2.type <- sapply(colnames(dat2), function(s) unlist(strsplit(s, 
                                                                            "-"))[4])
            dat2.tumor <- dat2[, grep("^01", dat2.type)]
            colnames(dat2.tumor) <- substr(colnames(dat2.tumor), 
                                           1, 12)
        }
        matching <- intersect(colnames(dat1.tumor), colnames(dat2.tumor))
        if (length(matching) == 0) {
            message("No matched samples")
            return(NULL)
        }
        dat1.good <- dat1.tumor[, matching]
        dat2.good <- dat2.tumor[, matching]
        rownames(dat1.good) <- paste("d1.", rownames(dat1.good), 
                                     sep = "")
        rownames(dat2.good) <- paste("d2.", rownames(dat2.good), 
                                     sep = "")
        mdata <- t(rbind(dat1.good, dat2.good))
        return(list(merged.data = t(mdata)))
    }
}
#Downloading somatic mutations (SNPs) and clinical data

cesc.mut <- getTCGA(disease="CESC", data.type="Mutation", type="somatic", 
                    clinical=TRUE)
luad.mut <- getTCGA(disease="LUAD", data.type="Mutation", type="somatic", 
                    clinical=TRUE)
coad.mut <- getTCGA(disease="COAD", data.type="Mutation", type="somatic", 
                    clinical=TRUE)
kirc.mut <- getTCGA(disease="KIRC", data.type="Mutation", type="somatic", 
                    clinical=TRUE)
lihc.mut <- getTCGA(disease="LIHC", data.type="Mutation", type="somatic", 
                    clinical=TRUE)
#Downloading methylation data
options(timeout=50000) # to increase timeout to 20,000 seconds 

cesc.met <- getTCGA(disease="CESC", data.type="Methylation", type="450K",
                    clinical=T)
luad.met <- getTCGA(disease="LUAD", data.type="Methylation", type="450K",
                    clinical=T)
coad.met <- getTCGA(disease="COAD", data.type="Methylation", type="450K",
                    clinical=T)
kirc.met <- getTCGA(disease="KIRC", data.type="Methylation", type="450K",
                    clinical=T)
lihc.met <- getTCGA(disease="LIHC", data.type="Methylation", type="450K",
                    clinical=T)

#Saving methylation for later use in Gini/Gap

save(cesc.met,file='cesc_met_raw.rdata')
save(luad.met,file='luad_met_raw.rdata')
save(coad.met,file='coad_met_raw.rdata')
save(kirc.met,file='kirc_met_raw.rdata')
save(lihc.met,file='lihc_met_raw.rdata')

#Merging methylation and mutation data and transposing the matrix

cesc.mut1 <- OMICSBIND2(dat1 =cesc.mut$dat, dat2 =t(cesc.mut$clinical))
#dimensions of cesc.mut1 : 15056 x 194 (15,006 genes / mutation data + 50 clinicals data) x 194 patients
luad.mut1 <- OMICSBIND2(dat1 =luad.mut$dat, dat2 =t(luad.mut$clinical))
coad.mut1 <- OMICSBIND2(dat1 =coad.mut$dat, dat2 =t(coad.mut$clinical))
kirc.mut1 <- OMICSBIND2(dat1 =kirc.mut$dat, dat2 =t(kirc.mut$clinical))
lihc.mut1 <- OMICSBIND2(dat1 =lihc.mut$dat, dat2 =t(lihc.mut$clinical))

CESC11 <- OMICSBIND2(dat1 =cesc.met$dat, dat2 =cesc.mut1$merged.data)
#dimensions of CESC11 : 49477 x 194 (15,006 genes / mutation data + 50 clinicals data + 482421 cpgs) x 194 patients
LUAD11 <- OMICSBIND2(dat1 =luad.met$dat, dat2 =luad.mut1$merged.data)
COAD11 <- OMICSBIND2(dat1 =coad.met$dat, dat2 =coad.mut1$merged.data)
KIRC11 <- OMICSBIND2(dat1 =kirc.met$dat, dat2 =kirc.mut1$merged.data)
LIHC11 <- OMICSBIND2(dat1 =lihc.met$dat, dat2 =lihc.mut1$merged.data)

# Writing combined data for second part of analysis

write.csv(CESC11$merged.data,'CESC.csv')
write.csv(LUAD11$merged.data,'LUAD.csv')
write.csv(COAD11$merged.data,'COAD.csv')
write.csv(KIRC11$merged.data,'KIRC.csv')
write.csv(LIHC11$merged.data,'LIHC.csv')
