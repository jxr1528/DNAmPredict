library(data.table)
setDTthreads(6)

#select cancer types to combine
cancers=c('LIHC','COAD','KIRC', 'CESC','LUAD')[c(1,4,5)]

CESC12=NULL
types=NULL
if(sum(cancers=='LIHC')==1){
    if(is.null(CESC12)){
        holder=as.matrix(fread('lihc.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1]
        #next need to filter race
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        CESC12=holder
        types=c(types,rep('LIHC',dim(CESC12)[1]))
    }else{
        holder=as.matrix(fread('lihc.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1] 
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        cols=intersect(colnames(CESC12),colnames(holder))
        indx1=match(cols,colnames(CESC12))
        indx2=match(cols,colnames(holder))
        CESC12=rbind(CESC12[,indx1],holder[,indx2])
        #print(paste('Combining',cancers[i]))
        types=c(types,rep('LIHC',dim(holder)[1]))
    }
    rm(holder)
    gc()
}


if(sum(cancers=='CESC')==1){
    if(is.null(CESC12)){
        holder=as.matrix(fread('cesc.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1]
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        CESC12=holder
        types=c(types,rep('CESC',dim(CESC12)[1]))
    }else{
        holder=as.matrix(fread('cesc.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1]
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        cols=intersect(colnames(CESC12),colnames(holder))
        indx1=match(cols,colnames(CESC12))
        indx2=match(cols,colnames(holder))
        CESC12=rbind(CESC12[,indx1],holder[,indx2])
        #print(paste('Combining',cancers[i]))
        types=c(types,rep('CESC',dim(holder)[1]))
    }
    rm(holder)
    gc()
}

if(sum(cancers=='COAD')==1){
    if(is.null(CESC12)){
        holder=as.matrix(fread('coad.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1]
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        CESC12=holder
        types=c(types,rep('COAD',dim(CESC12)[1]))
    }else{
        holder=as.matrix(fread('coad.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1] 
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        cols=intersect(colnames(CESC12),colnames(holder))
        indx1=match(cols,colnames(CESC12))
        indx2=match(cols,colnames(holder))
        CESC12=rbind(CESC12[,indx1],holder[,indx2])
        #print(paste('Combining',cancers[i]))
        types=c(types,rep('COAD',dim(holder)[1]))
    }
    rm(holder)
    gc()
}


if(sum(cancers=='KIRC')==1){
    if(is.null(CESC12)){
        holder=as.matrix(fread('kirc.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1]
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        CESC12=holder
        types=c(types,rep('KIRC',dim(CESC12)[1]))
    }else{
        holder=as.matrix(fread('kirc.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1] 
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        cols=intersect(colnames(CESC12),colnames(holder))
        indx1=match(cols,colnames(CESC12))
        indx2=match(cols,colnames(holder))
        CESC12=rbind(CESC12[,indx1],holder[,indx2])
        #print(paste('Combining',cancers[i]))
        types=c(types,rep('KIRC',dim(holder)[1]))
    }
    rm(holder)
    gc()
}


if(sum(cancers=='LUAD')==1){
    if(is.null(CESC12)){
        holder=as.matrix(fread('luad.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1]
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        CESC12=holder
        types=c(types,rep('LUAD',dim(CESC12)[1]))
    }else{
        holder=as.matrix(fread('luad.csv',stringsAsFactors = F))
        nam=holder[,1]
        nam=strsplit(nam,split='\\.')
        nam2=rep('',length(nam))
        for(i in 1:length(nam)){
            print(i)
            if(length(nam[[i]])>0){
                nam2[i]=nam[[i]][[length(nam[[i]])]]
            }
        }
        rownames(holder)=nam2
        rm(nam)
        rm(nam2)
        holder=t(holder[,-1])[,-1] 
        holder=holder[(holder[,'race']=='white'|holder[,'race']=='black or african american')&!is.na(holder[,'race']),]
        cols=intersect(colnames(CESC12),colnames(holder))
        indx1=match(cols,colnames(CESC12))
        indx2=match(cols,colnames(holder))
        CESC12=rbind(CESC12[,indx1],holder[,indx2])
        #print(paste('Combining',cancers[i]))
        types=c(types,rep('LUAD',dim(holder)[1]))
    }
    rm(holder)
    gc()
}

#cancer types
table(types)
## dim(CESC12) = 439 x 488800 ## CESC = 154 LUAD = 151 LIHC = 134
    ## MET= 482420 (1:482420) , SNP =6364 (482421 : 488784) and 
        ## Clinical = 16 (488785:488800)

## Separating genetic measures from clinical
CESC12_clin=CESC12[,-(1:488784)]
CESC12=CESC12[,(1:488784)]
storage.mode(CESC12)='numeric'

#set .7 as training and .3 as testing
set.seed(1314)
train.idx <- sample(1:dim(CESC12)[1],round(.7*dim(CESC12)[1]))
test.idx <- (1:dim(CESC12)[1])[-train.idx]




    #Pre-processing data and filtering data

    ###Based on Training dataset
all.var <- apply(CESC12[train.idx,1:488784], 2, var)

######CESC:482420 MET, 6364 SNP; LUAD:482420 MET, 6364 SNP ; LIHC:482420 MET, 6364 SNP

met.var <- all.var[1:482420]
snp.var <- all.var[482421:488784]

hist(met.var,breaks=50)
hist(met.var,breaks=50,ylim=c(0,50))
abline(v=quantile(met.var, 0.9999, na.rm=T),col='red')

hist(snp.var,breaks=50)
hist(snp.var,breaks=50,ylim=c(0,50))
abline(v=quantile(snp.var, 0.995, na.rm=T),col='red')


met.data <- CESC12[,1:482420][, met.var >= quantile(met.var, 0.9999, na.rm=T) & !is.na(met.var)]
snp.data <- CESC12[,482421:488784][, snp.var >= quantile(snp.var, 0.995, na.rm=T) & !is.na(snp.var)]
race.data <-CESC12_clin[,'race']
type.data <-types
rm(CESC12)
gc()

#Saving final combined dataset under: LIHC_CESC_LUAD_combined.rdata
save.image(paste(paste(cancers,collapse = '_'),'_combined.rdata',sep=''))


