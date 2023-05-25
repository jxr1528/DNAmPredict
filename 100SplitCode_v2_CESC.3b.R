
load("~/CESC_LUAD_combined.rdata")
options(error=traceback)
library(glmnet)
library(cluster)
library(cmmp)
library(nlme)
library(xtable)
library(gplots)
library(corrplot)
library(randomForest)
dim(met.data)
colnames(met.data)
colnames(snp.data)
split.num=100
seed=1516
set.seed(seed)
collector=list()

## Subset to only CESC

met.data=met.data[types=='CESC',]
snp.data=snp.data[types=='CESC',]
race.data=race.data[types=='CESC']
type.data=type.data[types=='CESC']

#100 split runs for CESC only analysis

for (z in 1:split.num) {
    
    #Set Training dataset index
    
    train.idx=sample(1:dim(met.data)[1],round(.7*dim(met.data)[1]))
    
    snp.data_filt<-snp.data
    dim(snp.data_filt) #154  56
    #To eliminate CpG (columns with NA values)
    which(is.na(met.data),arr.ind = T)
    if(length(which(is.na(met.data),arr.ind = T)[,2])!=0){
        met.data2=met.data[,-which(is.na(met.data),arr.ind = T)[,2]]} else{met.data2=met.data}
    dim(met.data2)  #154 37
    
    #Gap clustering                                                                                                    
    pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
    #set.seed(1315)
    gap2 <- clusGap(met.data2[train.idx,],pam1,10)
    gap2 #suggesting K=5, may want to make this automatic
    plot(gap2)
    
    #set.seed(1111)
    nc <- maxSE(f = gap2$Tab[, "gap"], SE.f = gap2$Tab[, "SE.sim"], method = "firstSEmax", SE.factor = 1)
    groups=kmeans(met.data2,nc)  #Set clustering groups
    groups$size # 3 groups clustered based on the profile of 37 CpG (Pam + Gap): 1=25 2=94 3=35
    ## total =154
    
    #Set Training dataset
    train2 <- data.frame(met.data2[train.idx,],snp.data_filt[train.idx,],race=race.data[train.idx],type=type.data[train.idx],groupid=groups$cluster[train.idx])
    test2 <- data.frame(met.data2[-train.idx,],snp.data_filt[-train.idx,],race=race.data[-train.idx],type=type.data[-train.idx],groupid=groups$cluster[-train.idx])
    
    #Set Clinical variables
    CESC12_clin2=as.data.frame(CESC12_clin)
    CESC12_clin2$Age=as.numeric(CESC12_clin2$yearstobirth)
    CESC12_clin2$Gender=as.factor(CESC12_clin2$gender)
    CESC12_clin2$pathologicstage=as.factor(CESC12_clin2$pathologicstage)
    CESC12_clin2$pathologyTstage=as.factor(CESC12_clin2$pathologyTstage)
    CESC12_clin2$pathologyNstage=as.factor(CESC12_clin2$pathologyNstage)
    CESC12_clin2$pathologyMstage=as.factor(CESC12_clin2$pathologyMstage)
    sum(!complete.cases(CESC12_clin2[,8:10])) #lost 40
    
    #Set staging variables
    table(CESC12_clin2$pathologyTstage,useNA = 'always')
    PathoT=rep(NA,dim(train2)[1])
    PathoT[CESC12_clin2$pathologyTstage=='t1'|CESC12_clin2$pathologyTstage=='t1a'|CESC12_clin2$pathologyTstage=='t1a1'|CESC12_clin2$pathologyTstage=='t1b'|CESC12_clin2$pathologyTstage=='t1b1'|CESC12_clin2$pathologyTstage=='t1b2'|CESC12_clin2$pathologyTstage=='tis']='t1'
    PathoT[CESC12_clin2$pathologyTstage=='t2'|CESC12_clin2$pathologyTstage=='t2a'|CESC12_clin2$pathologyTstage=='t2a1'|CESC12_clin2$pathologyTstage=='t2a2'|CESC12_clin2$pathologyTstage=='t2b']='t2'
    PathoT[CESC12_clin2$pathologyTstage=='t3'|CESC12_clin2$pathologyTstage=='t3a'|CESC12_clin2$pathologyTstage=='t3b']='t3'
    PathoT[CESC12_clin2$pathologyTstage=='t4']='t4'
    PathoT[CESC12_clin2$pathologyTstage=='tx']='tx'
    
    PathoM=rep(NA,dim(train2)[1])
    PathoM[CESC12_clin2$pathologyMstage=='m0']='m0'
    PathoM[CESC12_clin2$pathologyMstage=='m1'|CESC12_clin2$pathologyMstage=='m1a'|CESC12_clin2$pathologyMstage=='m1b']='m0'
    PathoM[CESC12_clin2$pathologyMstage=='mx']='mx'
    
    PathoN=CESC12_clin2$pathologyNstage
    PathoN[PathoN=='n3']=NA #exlucding single person category
    PathoN=droplevels(PathoN)
    
    CESC12_clin3=data.frame(Age=CESC12_clin2$Age,Gender=CESC12_clin2$Gender,PathoT=PathoT,PathoN=PathoN,PathoM=PathoM)
    
    comptrain=complete.cases(train2)
    comptest=complete.cases(test2)
    
    dim(CESC12_clin2)
    dim(all)
    table(types)
    table(CESC12_clin2$pathologicstage)
    table(types,CESC12_clin2$pathologyTstage)
    table(types,CESC12_clin2$pathologyNstage)
    table(types,CESC12_clin2$pathologyMstage)
    table(CESC12_clin2$pathologyTstage,CESC12_clin2$pathologyMstage)
    newstage=rep(NA,dim(train2)[1])
    str(CESC12_clin2$pathologicstage[types=='LUAD'])
    newstage[types=='LUAD']=as.character(CESC12_clin2$pathologicstage[types=='LUAD'])
    newstage[types=='CESC'&(CESC12_clin2$pathologyMstage=='mx'|
                                CESC12_clin2$pathologyNstage=='nx'|
                                CESC12_clin2$pathologyTstage=='tx'|
                                is.na(CESC12_clin2$pathologyMstage)|
                                is.na(CESC12_clin2$pathologyTstage))]='sx'
    newstage[types=='CESC'&is.na(newstage)&(CESC12_clin2$pathologyMstage=='m1'|
                                                CESC12_clin2$pathologyMstage=='m1a'|
                                                CESC12_clin2$pathologyMstage=='m1b'|
                                                CESC12_clin2$pathologyTstage=='t4')]='s4'
    newstage[types=='CESC'&is.na(newstage)&(CESC12_clin2$pathologyTstage=='t3'|
                                                CESC12_clin2$pathologyTstage=='t3a'|
                                                CESC12_clin2$pathologyTstage=='t3b'|
                                                CESC12_clin2$pathologyNstage=='n1')]='s3'
    newstage[types=='CESC'&is.na(newstage)&(CESC12_clin2$pathologyTstage=='t2'|
                                                CESC12_clin2$pathologyTstage=='t2a'|
                                                CESC12_clin2$pathologyTstage=='t2a1'|
                                                CESC12_clin2$pathologyTstage=='t2a2'|
                                                CESC12_clin2$pathologyTstage=='t2b')]='s2'
    newstage[types=='CESC'&is.na(newstage)&(CESC12_clin2$pathologyTstage=='t1'|
                                                CESC12_clin2$pathologyTstage=='t1a'|
                                                CESC12_clin2$pathologyTstage=='t1a1'|
                                                CESC12_clin2$pathologyTstage=='t1b'|
                                                CESC12_clin2$pathologyTstage=='t1b1'|
                                                CESC12_clin2$pathologyTstage=='t1b2')]='s1'
    New.Stage=newstage
    New.Stage[newstage=='s1'|newstage=='stage i'|newstage=='stage ia'|newstage=='stage ib']='Stage.1'
    New.Stage[newstage=='s2'|newstage=='stage ii'|newstage=='stage iia'|newstage=='stage iib']='Stage.2'
    New.Stage[newstage=='s3'|newstage=='stage iiia'|newstage=='stage iiib']='Stage.3'
    New.Stage[newstage=='s4'|newstage=='stage iv']='Stage.4'
    New.Stage[newstage=='sx'|is.na(newstage)]='Stage.NA'
    
    #Add clinical and stage variables
    
    CESC12_clin3=data.frame(Age=CESC12_clin2$Age,Gender=CESC12_clin2$Gender,Stage=New.Stage)[CESC12_clin2$tumortissuesite=='cervical',]
    
    train.clic=CESC12_clin3[train.idx[comptrain],]
    train2=data.frame(train2,train.clic)
    train2=train2[complete.cases(train2),]
    test2.clic=CESC12_clin3[-train.idx[comptest],]
    test2=data.frame(test2,test2.clic)
    test3=test2[test2$type=='CESC',]
    test2=test2[complete.cases(test2),]
    test3=test3[complete.cases(test3),]
    testrace=test3$race
    races=unique(train2$race)
    races #races[1]=white,races[2]=black
    
    intcol=list()
    coefcol=list()
    coefcol_m=list()
    coefcol_rau=list()
    
    #Set Beta, M and RAU variables
    
    #note matrix is overall, white, black
    cmmpcol=matrix(NA,nrow=37,ncol=3)
    lmcol=matrix(NA,nrow=37,ncol=3)
    lascol=matrix(NA,nrow=37,ncol=3)
    las_cmmpcol=matrix(NA,nrow=37,ncol=3)
    rfcol=matrix(NA,nrow=37,ncol=3)
    
    cmmpcol_m=matrix(NA,nrow=37,ncol=3)
    lmcol_m=matrix(NA,nrow=37,ncol=3)
    lascol_m=matrix(NA,nrow=37,ncol=3)
    las_cmmpcol_m=matrix(NA,nrow=37,ncol=3)
    rfcol_m=matrix(NA,nrow=37,ncol=3)
    
    cmmpcol_rau=matrix(NA,nrow=37,ncol=3)
    lmcol_rau=matrix(NA,nrow=37,ncol=3)
    lascol_rau=matrix(NA,nrow=37,ncol=3)
    las_cmmpcol_rau=matrix(NA,nrow=37,ncol=3)
    rfcol_rau=matrix(NA,nrow=37,ncol=3)
    
    colnames(cmmpcol)=colnames(lmcol)=colnames(lascol)=colnames(las_cmmpcol)=colnames(rfcol)=
        colnames(cmmpcol_m)=colnames(lmcol_m)=colnames(lascol_m)=colnames(las_cmmpcol_m)=colnames(rfcol_m)=
        colnames(cmmpcol_rau)=colnames(lmcol_rau)=colnames(lascol_rau)=colnames(las_cmmpcol_rau)=colnames(rfcol_rau)=c('all','white','black')
    
    convert_m=function(p){
        p2=log(p/(1-p),2)
        p2[p2==Inf|p2==-Inf]=NA
        p2
    }
    
    convert_rau=function(p){
        t = 2 * asin(sqrt(p))
        (46.47324337*t) - 23
    }
    
    
    #Final analysis for 37 filtered methylations
    
    for(outcome in 1:37){
        #set dataset and formula
        
        train6=train2
        train6=train6[,-which(names(train6)=='groupid')]
        train6=train6[,-which(names(train6)=='type')]
        train6=train6[,-which(names(train6)=='Gender')]
        test6=test3
        ##f.in <- as.formula(paste(colnames(train2)[outcome] ,'~ -1+',paste(c(colnames(train6)[-grep('cg',colnames(train6))],'Stage*type'),collapse='+')))
        f.in <- as.formula(paste(colnames(train2)[outcome] ,'~ -1+',paste(c(colnames(train6)[-grep('cg',colnames(train6))]),collapse='+')))
        train3=data.frame(train2[,1:37],model.matrix(f.in,data=train2),groupid=train2$groupid)
        
        train3_m=train3
        train3_m[,outcome]=convert_m(train3_m[,outcome])
        train3_m=train3_m[complete.cases(train3_m),]
        
        train3_rau=train3
        train3_rau[,outcome]=convert_rau(train3_rau[,outcome])
        train3_rau=train3_rau[complete.cases(train3_rau),]
        

        
        test4=data.frame(test3[,1:37],model.matrix(f.in,data=test3)) ##,groupid=test3$groupid)
        if ("StageStage.4" %in% colnames (test4)){test4=test4} else {test4 <- cbind(test4, StageStage.4=0)}
        if ("StageStage.3" %in% colnames (test4)){test4=test4} else {test4 <- cbind(test4, StageStage.3=0)}
        if ("StageStage.2" %in% colnames (test4)){test4=test4} else {test4 <- cbind(test4, StageStage.2=0)}
        if ("StageStage.NA" %in% colnames (test4)){test4=test4} else {test4 <- cbind(test4, StageStage.NA=0)}
        
        ##test4b=data.frame(test2[,1:37],model.matrix(f.in,data=test2),groupid=test2$groupid)
        
        
        
        
        test4_m=test4
        test4_m[,outcome]=convert_m(test4_m[,outcome])
        ## test4_m=test4_m[complete.cases(test4_m),]  
        
        test4_rau=test4
        test4_rau[,outcome]=convert_rau(test4_rau[,outcome])
        ## test4_rau=test4_rau[complete.cases(test4_rau),]
        
        #test4b_m=test4b
        #test4b_m[,outcome]=convert_m(test4b_m[,outcome])
        ##test4b_m=test4b_m[complete.cases(test4b_m),]  
        
        #test4b_rau=test4b
        #test4b_rau[,outcome]=convert_rau(test4b_rau[,outcome])
        ##test4b_rau=test4b_rau[complete.cases(test4b_rau),]
        
        ## CMMP
        
        r.in <- as.formula(~1|groupid) 
        
        test.asian=NULL
        f.in=as.formula(paste(colnames(train3)[outcome] ,'~ -1+',paste(c(colnames(train3)[c(38:92,95:(dim(train3)[2]-1))]),collapse='+')))
        # in above funcion we are leaving out column 93 (KEAP1) due to collinearity
        
        try(test.asian <- cmmp(f.in, r.in, train = train3[,c(outcome,38:dim(train3)[2])], x.new = test4[,c(38:dim(test4)[2])], y.new = test4[,outcome], x.fut = NULL, match.train = NULL, a1 = NULL, n.new = NULL, interval = TRUE))
        
        test.asian_m=NULL
        try(test.asian_m <- cmmp(f.in, r.in, train = train3_m[,c(outcome,38:dim(train3)[2])], x.new = test4_m[,c(38:dim(test4)[2])], y.new = test4_m[,outcome], x.fut = NULL, match.train = NULL, a1 = NULL, n.new = NULL, interval = TRUE))
        
        test.asian_rau=NULL
        try(test.asian_rau <- cmmp(f.in, r.in, train = train3_rau[,c(outcome,38:dim(train3)[2])], x.new = test4_rau[,c(38:dim(test4)[2])], y.new = test4_rau[,outcome], x.fut = NULL, match.train = NULL, a1 = NULL, n.new = NULL, interval = TRUE))
        
      
        
        if(is.null(test.asian)|is.null(test.asian_m)|is.null(test.asian_rau)){
            break()}
        
        testext=test4[-(1:dim(test4)[1]),]
        testext[1,]=rep(0,dim(testext)[2])
        testext=testext[rep(1,length(unique((train2$groupid)))),]
        testext$groupid=unique((train2$groupid))
        
        ##f.in_int <- as.formula(paste(colnames(train6)[outcome] ,'~ -1+',paste(c(colnames(train6)[c(38:dim(train6)[2])],paste('race*',colnames(train6)[c(38:93)],sep=''),'type*Stage'),collapse='+')))
        f.in_int <- as.formula(paste(colnames(train6)[outcome] ,'~ -1+',paste(c(colnames(train6)[c(38:92, 94:dim(train6)[2])],paste('race*',colnames(train6)[c(38:92)],sep='')),collapse='+')))
        train5=model.matrix(f.in_int,data=train6)
        if ("StageStage.4" %in% colnames (train5)){train5=train5} else {train5 <- cbind(train5, StageStage.4=0)}
        if ("StageStage.3" %in% colnames (train5)){train5=train5} else {train5 <- cbind(train5, StageStage.3=0)}
        if ("StageStage.2" %in% colnames (train5)){train5=train5} else {train5 <- cbind(train5, StageStage.2=0)}
        if ("StageStage.NA" %in% colnames (train5)){train5=train5} else {train5 <- cbind(train5, StageStage.NA=0)}
        test5=model.matrix(f.in_int,data=test2)
        if ("StageStage.4" %in% colnames (test5)){test5=test5} else {test5 <- cbind(test5, StageStage.4=0)}
        if ("StageStage.3" %in% colnames (test5)){test5=test5} else {test5 <- cbind(test5, StageStage.3=0)}
        if ("StageStage.2" %in% colnames (test5)){test5=test5} else {test5 <- cbind(test5, StageStage.2=0)}
        if ("StageStage.NA" %in% colnames (test5)){test5=test5} else {test5 <- cbind(test5, StageStage.NA=0)}
        storage.mode(train5)='numeric'
        storage.mode(test5)='numeric'
        
        ##ENET
        
        cvfit=cv.glmnet(train5,train2[,outcome],intercept=T)
        elaspred=predict(cvfit,newx=test5,s='lambda.min')
        
        
        comp_m=complete.cases(convert_m(train2[,outcome]))
        cvfit_m=cv.glmnet(train5[comp_m,],convert_m(train2[,outcome])[comp_m],intercept=T)
        elaspred_m=predict(cvfit_m,newx=test5,s='lambda.min')
        
        comp_rau=complete.cases(convert_rau(train2[,outcome]))
        cvfit_rau=cv.glmnet(train5[comp_rau,],convert_rau(train2[,outcome])[comp_rau],intercept=T)
        elaspred_rau=predict(cvfit_rau,newx=test5,s='lambda.min')
        
        
        ### LR (LM)
        
        lm_mod=glm(f.in,data=train3[,c(outcome,38:dim(train3)[2])])
        lm_modpred=predict(lm_mod,newdata = test4[,c(38:dim(test4)[2])])
        
        lm_mod_m=glm(f.in,data=train3_m[,c(outcome,38:dim(train3)[2])])
        lm_modpred_m=predict(lm_mod_m,newdata = test4[,c(38:dim(test4)[2])])
        
        lm_mod_rau=glm(f.in,data=train3_rau[,c(outcome,38:dim(train3)[2])])
        lm_modpred_rau=predict(lm_mod_rau,newdata = test4[,c(38:dim(test4)[2])])
        
        ### Random Forest (RF)
        
        rf_mod=randomForest(f.in,data=train3[,c(outcome,38:dim(train3)[2])])
        rf_modpred=predict(rf_mod,newdata = test4[,c(38:dim(test4)[2])])
        
        rf_mod_m=randomForest(f.in,data=train3_m[,c(outcome,38:dim(train3)[2])])
        rf_modpred_m=predict(rf_mod_m,newdata = test4[,c(38:dim(test4)[2])])
        
        rf_mod_rau=randomForest(f.in,data=train3_rau[,c(outcome,38:dim(train3)[2])])
        rf_modpred_rau=predict(rf_mod_rau,newdata = test4[,c(38:dim(test4)[2])])
        
        
        

        
        ## Taking residual from elastic net
        
        resid=train3[,outcome]-predict(cvfit,newx=train5,s='lambda.min')
        resid_m=train3_m[,outcome]-predict(cvfit_m,newx=train5,s='lambda.min')
        resid_rau=train3_rau[,outcome]-predict(cvfit_rau,newx=train5,s='lambda.min')
        
        ### Looking at resid, the column name is "lambda.min"
        ### Here we change the column name to X1 ; to match the new f.in formula f.in<-X1~-1
        
        colnames(resid)='X1'
        colnames(resid_m)='X1'
        colnames(resid_rau)='X1'
        
        
        f.in<-X1~-1
        combocmmp <- cmmp(f.in, r.in, train = data.frame(resid,groupid=train2[,'groupid']), x.new = test4[,c(38:100)], y.new = resid, x.fut = NULL, match.train = NULL, a1 = NULL, n.new = NULL, interval = TRUE)
        combocmmp_m <- cmmp(f.in, r.in, train = data.frame(resid_m,groupid=train2[,'groupid']), x.new = test4[,c(38:100)], y.new = resid_m, x.fut = NULL, match.train = NULL, a1 = NULL, n.new = NULL, interval = TRUE)
        
        combocmmp_rau=NULL
        try(combocmmp_rau <- cmmp(f.in, r.in, train = data.frame(resid_rau,groupid=train2[,'groupid']), x.new = test4[,c(38:100)], y.new = resid_rau, x.fut = NULL, match.train = NULL, a1 = NULL, n.new = NULL, interval = TRUE))
        if(is.null(combocmmp_rau)){
            try(combocmmp_rau <- cmmp2(f.in, r.in, train = data.frame(resid_rau,groupid=train2[,'groupid']), x.new = test4[,c(38:100)], y.new = resid_rau, x.fut = NULL, match.train = NULL, a1 = NULL, n.new = NULL, interval = TRUE))
        }
        
        
        ### Collecting MSPE to compare the different methods
        ### a) for beta values w/o transforming into M or Rau and for LR,  CMMP, ENET and CMMP + ENET
        
        ## LR
        
        lmcol[outcome,]=c(mean((lm_modpred-test4[,outcome])^2),
                          mean((lm_modpred-test4[,outcome])[testrace==races[1]]^2), #white
                          mean((lm_modpred-test4[,outcome])[testrace==races[2]]^2)) #black
        
        ## RF
        
        rfcol[outcome,]=c(mean((rf_modpred-test4[,outcome])^2),
                          mean((rf_modpred-test4[,outcome])[testrace==races[1]]^2), #white
                          mean((rf_modpred-test4[,outcome])[testrace==races[2]]^2)) #black
        
        
        ## CMMP
        cmmpcol[outcome,]=c( mean((test.asian$mixed.pred-test4[,outcome])^2),
                             mean((test.asian$mixed.pred-test4[,outcome])[testrace==races[1]]^2),
                             mean((test.asian$mixed.pred-test4[,outcome])[testrace==races[2]]^2))
        
        ## ENET
        
        lascol[outcome,]=c(mean((elaspred-test4[,outcome])^2),
                           mean((elaspred-test4[,outcome])[testrace==races[1]]^2),
                           mean((elaspred-test4[,outcome])[testrace==races[2]]^2))
        
        ## CMMP + ENET
        
        las_cmmpcol[outcome,]=c(mean((elaspred+combocmmp$mixed.pred-test4[,outcome])^2),
                                mean((elaspred+combocmmp$mixed.pred-test4[,outcome])[testrace==races[1]]^2),
                                mean((elaspred+combocmmp$mixed.pred-test4[,outcome])[testrace==races[2]]^2))
        
        
        ### b) for m values
        
        ##LR _m
        
        lmcol_m[outcome,]=c(mean((lm_modpred_m-test4_m[,outcome])^2),
                            mean((lm_modpred_m-test4_m[,outcome])[testrace==races[1]]^2), #white
                            mean((lm_modpred_m-test4_m[,outcome])[testrace==races[2]]^2)) #black
        
        ##RF_m
        
        rfcol_m[outcome,]=c(mean((rf_modpred_m-test4_m[,outcome])^2),
                          mean((rf_modpred_m-test4_m[,outcome])[testrace==races[1]]^2), #white
                          mean((rf_modpred_m-test4_m[,outcome])[testrace==races[2]]^2)) #black
        
        
        
        ##CMMP _m
        
        cmmpcol_m[outcome,]=c( mean((test.asian_m$mixed.pred-test4_m[,outcome])^2),
                               mean((test.asian_m$mixed.pred-test4_m[,outcome])[testrace==races[1]]^2),
                               mean((test.asian_m$mixed.pred-test4_m[,outcome])[testrace==races[2]]^2))
        
        ## ENET _m
        
        lascol_m[outcome,]=c(mean((elaspred_m-test4_m[,outcome])^2),
                             mean((elaspred_m-test4_m[,outcome])[testrace==races[1]]^2),
                             mean((elaspred_m-test4_m[,outcome])[testrace==races[2]]^2))
        
        ## ENET + CMMP
        
        las_cmmpcol_m[outcome,]=c(mean((elaspred_m+combocmmp_m$mixed.pred-test4_m[,outcome])^2),
                                  mean((elaspred_m+combocmmp_m$mixed.pred-test4_m[,outcome])[testrace==races[1]]^2),
                                  mean((elaspred_m+combocmmp_m$mixed.pred-test4_m[,outcome])[testrace==races[2]]^2))
        
        ### c) for rau values
        
        ##LR _rau
        
        lmcol_rau[outcome,]=c(mean((lm_modpred_rau-test4_rau[,outcome])^2),
                              mean((lm_modpred_rau-test4_rau[,outcome])[testrace==races[1]]^2), #white
                              mean((lm_modpred_rau-test4_rau[,outcome])[testrace==races[2]]^2)) #black
        
        ##RF_rau
        
        rfcol_rau[outcome,]=c(mean((rf_modpred_rau-test4_rau[,outcome])^2),
                          mean((rf_modpred_rau-test4_rau[,outcome])[testrace==races[1]]^2), #white
                          mean((rf_modpred_rau-test4_rau[,outcome])[testrace==races[2]]^2)) #black
        
        ##CMMP _rau
        
        cmmpcol_rau[outcome,]=c( mean((test.asian_rau$mixed.pred-test4_rau[,outcome])^2),
                                 mean((test.asian_rau$mixed.pred-test4_rau[,outcome])[testrace==races[1]]^2),
                                 mean((test.asian_rau$mixed.pred-test4_rau[,outcome])[testrace==races[2]]^2))
        
        ## ENET _m
        
        lascol_rau[outcome,]=c(mean((elaspred_rau-test4_rau[,outcome])^2),
                               mean((elaspred_rau-test4_rau[,outcome])[testrace==races[1]]^2),
                               mean((elaspred_rau-test4_rau[,outcome])[testrace==races[2]]^2))
        
        ## ENET + CMMP
        
        las_cmmpcol_rau[outcome,]=c(mean((elaspred_rau+combocmmp_rau$mixed.pred-test4_rau[,outcome])^2),
                                    mean((elaspred_rau+combocmmp_rau$mixed.pred-test4_rau[,outcome])[testrace==races[1]]^2),
                                    mean((elaspred_rau+combocmmp_rau$mixed.pred-test4_rau[,outcome])[testrace==races[2]]^2))
        
        coefcol[[outcome]]=coef(cvfit,s='lambda.min')
        coefcol_m[[outcome]]=coef(cvfit_m,s='lambda.min')
        coefcol_rau[[outcome]]=coef(cvfit_rau,s='lambda.min')
        
    }
  
    if(length(coefcol)==0){
        z=z-1
        next()
    }
    
    print(z)
    
    collector[[z]]=list(coefcol=list(beta=coefcol,m=coefcol_m,rau=coefcol_rau),errors=list(cmmpcol=cmmpcol,lmcol=lmcol,lascol=lascol,las_cmmpcol=las_cmmpcol,rfcol=rfcol,
                                                                                           cmmpcol_m=cmmpcol_m,lmcol_m=lmcol_m,lascol_m=lascol_m,las_cmmpcol_m=las_cmmpcol_m,rfcol_m=rfcol_m,
                                                                                           cmmpcol_rau=cmmpcol_rau,lmcol_rau=lmcol_rau,lascol_rau=lascol_rau,las_cmmpcol_rau=las_cmmpcol_rau,rfcol_rau=rfcol_rau))
    
}


## Need name all snps used
### To change to i in 1:100 depending on how many runs

nam=rownames(collector[[1]]$coefcol$beta[[1]])
for(i in 1:100){
    for(j in 1:37){
        nam=union(nam,rownames(collector[[i]]$coefcol$beta[[j]]))
    }
}


#coeficient plot
coefmats=list()
for( i in 1:37){
    coefmats[[i]]=matrix(NA,nrow=length(nam),ncol=100)
    rownames(coefmats[[i]])=nam
}

for(j in 1:37){
    ##for(i in c(1:100)){
    for(i in c(1:8,10:31,33:53,55:93,95:100)){
        ##coefmats[[j]][,i]=matrix(collector[[i]]$coefcol[[j]])[match(rownames(coefmats[[j]]),rownames(collector[[i]]$coefcol$beta[[j]])),1]
        coefmats[[j]][,i]=matrix(collector[[i]]$coefcol$beta[[j]])
    }
    
}

se <- function(x) sqrt(var(x,na.rm=T)/sum(!is.na(x)))

#collecting mean and se 
meanholder=matrix(NA,nrow=length(nam),ncol=37)
#seholder=matrix(NA,nrow=length(nam),nrow=37)
rownames(meanholder)=nam
#rownames(seholder)=nam

for(i in 1:length(coefmats)){
    meanholder[,i]=rowMeans(coefmats[[i]],na.rm=T)
    #seholder[,i]=apply(coefmats[[i]],1)
}

meanholder2=meanholder[c(1:118),]

vars=apply(meanholder2[-1,],1,var)
mvar=which(vars>quantile(vars,.9))

par(mar=c(6.5,4.1,2.5,1))
plot(1:(dim(meanholder2)[1]-1),meanholder2[-1,1],col=c(rep(2,56),rep(4,7),rep(3,55)),ylim=range(meanholder2),
     xlab='',ylab='Coefficent',xaxt='n',pch=20,cex=.7,main='Elastic Net Coefficients, mean over 100 samples, shaded=interactions-CESC only')
for(i in 2:37){
    points(1:(dim(meanholder2)[1]-1),meanholder2[-1,i],col=c(rep(2,56),rep(4,7),rep(3,55)),pch=20,cex=.7)
}
rect(64, -5, 150, 5, density = 8)
legend('topleft',legend=c('SNP','Clinical','SNP*Race'),lty=1,col=c(2,4,3),bg='transparent')
axis(1,at=mvar,rownames(meanholder2)[-1][mvar],las=2,cex.axis=.46)

##_______________________________________________________________________________
##Optional to rotate the Elastic net coefficient plot

par(mar=c(4,6,4,2))
plot(meanholder2[-1,1],1:(dim(meanholder2)[1]-1),col=c(rep(2,56),rep(4,7),rep(3,55)), xlim= range(-0.5:0.5),
     ylab='Variables',xlab='Coefficent',yaxt='n',pch=20,cex=.7,main='Elastic Net Coefficients, mean over 100 samples, shaded are interactions')
for(i in 2:37){
    points(meanholder2[-1,i], 1:(dim(meanholder2)[1]-1),col=c(rep(2,56),rep(4,7),rep(3,55)),pch=20,cex=.7)
}
rect(ybottom=64, ytop =118 , xright=-0.7, 150, 5, density = 8)
legend('bottomright',legend=c('SNP','Clinical','SNP*Race'),lty=1,col=c(2,4,3),bg='transparent')
axis(2,at=mvar,rownames(meanholder2)[-1][mvar],las=2,cex.axis=.4)

##------------------------------------------------------------------------------

#error plot

errcol=list()
errcol$cmmpcol=matrix(NA,ncol=100,nrow=37)
errcol$lmcol=matrix(NA,ncol=100,nrow=37)
errcol$lascol=matrix(NA,ncol=100,nrow=37)
errcol$las_cmmpcol=matrix(NA,ncol=100,nrow=37)
errcol$rfcol=matrix(NA,ncol=100,nrow=37)

errcol$cmmpcolw=matrix(NA,ncol=100,nrow=37)
errcol$lmcolw=matrix(NA,ncol=100,nrow=37)
errcol$lascolw=matrix(NA,ncol=100,nrow=37)
errcol$las_cmmpcolw=matrix(NA,ncol=100,nrow=37)
errcol$rfcolw=matrix(NA,ncol=100,nrow=37)

errcol$cmmpcolb=matrix(NA,ncol=100,nrow=37)
errcol$lmcolb=matrix(NA,ncol=100,nrow=37)
errcol$lascolb=matrix(NA,ncol=100,nrow=37)
errcol$las_cmmpcolb=matrix(NA,ncol=100,nrow=37)
errcol$rfcolb=matrix(NA,ncol=100,nrow=37)


errcol$cmmpcol_m=matrix(NA,ncol=100,nrow=37)
errcol$lmcol_m=matrix(NA,ncol=100,nrow=37)
errcol$lascol_m=matrix(NA,ncol=100,nrow=37)
errcol$las_cmmpcol_m=matrix(NA,ncol=100,nrow=37)
errcol$rfcol_m=matrix(NA,ncol=100,nrow=37)

errcol$cmmpcolw_m=matrix(NA,ncol=100,nrow=37)
errcol$lmcolw_m=matrix(NA,ncol=100,nrow=37)
errcol$lascolw_m=matrix(NA,ncol=100,nrow=37)
errcol$las_cmmpcolw_m=matrix(NA,ncol=100,nrow=37)
errcol$rfcolw_m=matrix(NA,ncol=100,nrow=37)

errcol$cmmpcolb_m=matrix(NA,ncol=100,nrow=37)
errcol$lmcolb_m=matrix(NA,ncol=100,nrow=37)
errcol$lascolb_m=matrix(NA,ncol=100,nrow=37)
errcol$las_cmmpcolb_m=matrix(NA,ncol=100,nrow=37)
errcol$rfcolb_m=matrix(NA,ncol=100,nrow=37)


errcol$cmmpcol_rau=matrix(NA,ncol=100,nrow=37)
errcol$lmcol_rau=matrix(NA,ncol=100,nrow=37)
errcol$lascol_rau=matrix(NA,ncol=100,nrow=37)
errcol$las_cmmpcol_rau=matrix(NA,ncol=100,nrow=37)
errcol$rfcol_rau=matrix(NA,ncol=100,nrow=37)

errcol$cmmpcolw_rau=matrix(NA,ncol=100,nrow=37)
errcol$lmcolw_rau=matrix(NA,ncol=100,nrow=37)
errcol$lascolw_rau=matrix(NA,ncol=100,nrow=37)
errcol$las_cmmpcolw_rau=matrix(NA,ncol=100,nrow=37)
errcol$rfcolw_rau=matrix(NA,ncol=100,nrow=37)

errcol$cmmpcolb_rau=matrix(NA,ncol=100,nrow=37)
errcol$lmcolb_rau=matrix(NA,ncol=100,nrow=37)
errcol$lascolb_rau=matrix(NA,ncol=100,nrow=37)
errcol$las_cmmpcolb_rau=matrix(NA,ncol=100,nrow=37)
errcol$rfcolb_rau=matrix(NA,ncol=100,nrow=37)

for(i in c(1:8,10:31,33:53,55:93,95:100)){
    ##for(i in c(1:17,19:62,64:72,74:84,86:89)){
    errcol$cmmpcol[,i]=collector[[i]]$errors$cmmpcol[,1]
    errcol$cmmpcolw[,i]=collector[[i]]$errors$cmmpcol[,2]
    errcol$cmmpcolb[,i]=collector[[i]]$errors$cmmpcol[,3]
    
    errcol$lmcol[,i]=collector[[i]]$errors$lmcol[,1]
    errcol$lmcolw[,i]=collector[[i]]$errors$lmcol[,2]
    errcol$lmcolb[,i]=collector[[i]]$errors$lmcol[,3]
    
    errcol$lascol[,i]=collector[[i]]$errors$lascol[,1]
    errcol$lascolw[,i]=collector[[i]]$errors$lascol[,2]
    errcol$lascolb[,i]=collector[[i]]$errors$lascol[,3]
    
    errcol$las_cmmpcol[,i]=collector[[i]]$errors$las_cmmpcol[,1]
    errcol$las_cmmpcolw[,i]=collector[[i]]$errors$las_cmmpcol[,2]
    errcol$las_cmmpcolb[,i]=collector[[i]]$errors$las_cmmpcol[,3]
    
    errcol$rfcol[,i]=collector[[i]]$errors$rfcol[,1]
    errcol$rfcolw[,i]=collector[[i]]$errors$rfcol[,2]
    errcol$rfcolb[,i]=collector[[i]]$errors$rfcol[,3]
    
    
    errcol$cmmpcol_m[,i]=collector[[i]]$errors$cmmpcol_m[,1]
    errcol$cmmpcolw_m[,i]=collector[[i]]$errors$cmmpcol_m[,2]
    errcol$cmmpcolb_m[,i]=collector[[i]]$errors$cmmpcol_m[,3]
    
    errcol$lmcol_m[,i]=collector[[i]]$errors$lmcol_m[,1]
    errcol$lmcolw_m[,i]=collector[[i]]$errors$lmcol_m[,2]
    errcol$lmcolb_m[,i]=collector[[i]]$errors$lmcol_m[,3]
    
    errcol$lascol_m[,i]=collector[[i]]$errors$lascol_m[,1]
    errcol$lascolw_m[,i]=collector[[i]]$errors$lascol_m[,2]
    errcol$lascolb_m[,i]=collector[[i]]$errors$lascol_m[,3]
    
    errcol$las_cmmpcol_m[,i]=collector[[i]]$errors$las_cmmpcol_m[,1]
    errcol$las_cmmpcolw_m[,i]=collector[[i]]$errors$las_cmmpcol_m[,2]
    errcol$las_cmmpcolb_m[,i]=collector[[i]]$errors$las_cmmpcol_m[,3]
    
    errcol$rfcol_m[,i]=collector[[i]]$errors$rfcol_m[,1]
    errcol$rfcolw_m[,i]=collector[[i]]$errors$rfcol_m[,2]
    errcol$rfcolb_m[,i]=collector[[i]]$errors$rfcol_m[,3]
    
    
    errcol$cmmpcol_rau[,i]=collector[[i]]$errors$cmmpcol_rau[,1]
    errcol$cmmpcolw_rau[,i]=collector[[i]]$errors$cmmpcol_rau[,2]
    errcol$cmmpcolb_rau[,i]=collector[[i]]$errors$cmmpcol_rau[,3]
    
    errcol$lmcol_rau[,i]=collector[[i]]$errors$lmcol_rau[,1]
    errcol$lmcolw_rau[,i]=collector[[i]]$errors$lmcol_rau[,2]
    errcol$lmcolb_rau[,i]=collector[[i]]$errors$lmcol_rau[,3]
    
    errcol$lascol_rau[,i]=collector[[i]]$errors$lascol_rau[,1]
    errcol$lascolw_rau[,i]=collector[[i]]$errors$lascol_rau[,2]
    errcol$lascolb_rau[,i]=collector[[i]]$errors$lascol_rau[,3]
    
    errcol$las_cmmpcol_rau[,i]=collector[[i]]$errors$las_cmmpcol_rau[,1]
    errcol$las_cmmpcolw_rau[,i]=collector[[i]]$errors$las_cmmpcol_rau[,2]
    errcol$las_cmmpcolb_rau[,i]=collector[[i]]$errors$las_cmmpcol_rau[,3]
    
    errcol$rfcol_rau[,i]=collector[[i]]$errors$rfcol_rau[,1]
    errcol$rfcolw_rau[,i]=collector[[i]]$errors$rfcol_rau[,2]
    errcol$rfcolb_rau[,i]=collector[[i]]$errors$rfcol_rau[,3]
    
}

error=list()
error$cmmp=rowMeans(errcol$cmmpcol,na.rm=T)
error$cmmpw=rowMeans(errcol$cmmpcolw,na.rm=T)
error$cmmpb=rowMeans(errcol$cmmpcolb,na.rm=T)

error$lm=rowMeans(errcol$lmcol,na.rm=T)
error$lmw=rowMeans(errcol$lmcolw,na.rm=T)
error$lmb=rowMeans(errcol$lmcolb,na.rm=T)

error$las=rowMeans(errcol$lascol,na.rm=T)
error$lasw=rowMeans(errcol$lascolw,na.rm=T)
error$lasb=rowMeans(errcol$lascolb,na.rm=T)

error$las_cmmp=rowMeans(errcol$las_cmmpcol,na.rm=T)
error$las_cmmpw=rowMeans(errcol$las_cmmpcolw,na.rm=T)
error$las_cmmpb=rowMeans(errcol$las_cmmpcolb,na.rm=T)

error$rf=rowMeans(errcol$rfcol,na.rm=T)
error$rfw=rowMeans(errcol$rfcolw,na.rm=T)
error$rfb=rowMeans(errcol$rfcolb,na.rm=T)

error$cmmp_m=rowMeans(errcol$cmmpcol_m,na.rm=T)
error$cmmpw_m=rowMeans(errcol$cmmpcolw_m,na.rm=T)
error$cmmpb_m=rowMeans(errcol$cmmpcolb_m,na.rm=T)

error$lm_m=rowMeans(errcol$lmcol_m,na.rm=T)
error$lmw_m=rowMeans(errcol$lmcolw_m,na.rm=T)
error$lmb_m=rowMeans(errcol$lmcolb_m,na.rm=T)

error$las_m=rowMeans(errcol$lascol_m,na.rm=T)
error$lasw_m=rowMeans(errcol$lascolw_m,na.rm=T)
error$lasb_m=rowMeans(errcol$lascolb_m,na.rm=T)

error$las_cmmp_m=rowMeans(errcol$las_cmmpcol_m,na.rm=T)
error$las_cmmpw_m=rowMeans(errcol$las_cmmpcolw_m,na.rm=T)
error$las_cmmpb_m=rowMeans(errcol$las_cmmpcolb_m,na.rm=T)

error$rf_m=rowMeans(errcol$rfcol_m,na.rm=T)
error$rfw_m=rowMeans(errcol$rfcolw_m,na.rm=T)
error$rfb_m=rowMeans(errcol$rfcolb_m,na.rm=T)

error$cmmp_rau=rowMeans(errcol$cmmpcol_rau,na.rm=T)
error$cmmpw_rau=rowMeans(errcol$cmmpcolw_rau,na.rm=T)
error$cmmpb_rau=rowMeans(errcol$cmmpcolb_rau,na.rm=T)

error$lm_rau=rowMeans(errcol$lmcol_rau,na.rm=T)
error$lmw_rau=rowMeans(errcol$lmcolw_rau,na.rm=T)
error$lmb_rau=rowMeans(errcol$lmcolb_rau,na.rm=T)

error$las_rau=rowMeans(errcol$lascol_rau,na.rm=T)
error$lasw_rau=rowMeans(errcol$lascolw_rau,na.rm=T)
error$lasb_rau=rowMeans(errcol$lascolb_rau,na.rm=T)

error$las_cmmp_rau=rowMeans(errcol$las_cmmpcol_rau,na.rm=T)
error$las_cmmpw_rau=rowMeans(errcol$las_cmmpcolw_rau,na.rm=T)
error$las_cmmpb_rau=rowMeans(errcol$las_cmmpcolb_rau,na.rm=T)

error$rf_rau=rowMeans(errcol$rfcol_rau,na.rm=T)
error$rfw_rau=rowMeans(errcol$rfcolw_rau,na.rm=T)
error$rfb_rau=rowMeans(errcol$rfcolb_rau,na.rm=T)

errors=list()
errors$cmmp=apply(errcol$cmmpcol,1,se)
errors$cmmpw=apply(errcol$cmmpcolw,1,se)
errors$cmmpb=apply(errcol$cmmpcolb,1,se)

errors$lm=apply(errcol$lmcol,1,se)
errors$lmw=apply(errcol$lmcolw,1,se)
errors$lmb=apply(errcol$lmcolb,1,se)

errors$las=apply(errcol$lascol,1,se)
errors$lasw=apply(errcol$lascolw,1,se)
errors$lasb=apply(errcol$lascolb,1,se)

errors$las_cmmp=apply(errcol$las_cmmpcol,1,se)
errors$las_cmmpw=apply(errcol$las_cmmpcolw,1,se)
errors$las_cmmpb=apply(errcol$las_cmmpcolb,1,se)

errors$rf=apply(errcol$rfcol,1,se)
errors$rfw=apply(errcol$rfcolw,1,se)
errors$rfb=apply(errcol$rfcolb,1,se)

errors$cmmp_m=apply(errcol$cmmpcol_m,1,se)
errors$cmmpw_m=apply(errcol$cmmpcolw_m,1,se)
errors$cmmpb_m=apply(errcol$cmmpcolb_m,1,se)

errors$lm_m=apply(errcol$lmcol_m,1,se)
errors$lmw_m=apply(errcol$lmcolw_m,1,se)
errors$lmb_m=apply(errcol$lmcolb_m,1,se)

errors$las_m=apply(errcol$lascol_m,1,se)
errors$lasw_m=apply(errcol$lascolw_m,1,se)
errors$lasb_m=apply(errcol$lascolb_m,1,se)

errors$las_cmmp_m=apply(errcol$las_cmmpcol_m,1,se)
errors$las_cmmpw_m=apply(errcol$las_cmmpcolw_m,1,se)
errors$las_cmmpb_m=apply(errcol$las_cmmpcolb_m,1,se)

errors$rf_m=apply(errcol$rfcol_m,1,se)
errors$rfw_m=apply(errcol$rfcolw_m,1,se)
errors$rfb_m=apply(errcol$rfcolb_m,1,se)

errors$cmmp_rau=apply(errcol$cmmpcol_rau,1,se)
errors$cmmpw_rau=apply(errcol$cmmpcolw_rau,1,se)
errors$cmmpb_rau=apply(errcol$cmmpcolb_rau,1,se)

errors$lm_rau=apply(errcol$lmcol_rau,1,se)
errors$lmw_rau=apply(errcol$lmcolw_rau,1,se)
errors$lmb_rau=apply(errcol$lmcolb_rau,1,se)

errors$las_rau=apply(errcol$lascol_rau,1,se)
errors$lasw_rau=apply(errcol$lascolw_rau,1,se)
errors$lasb_rau=apply(errcol$lascolb_rau,1,se)

errors$las_cmmp_rau=apply(errcol$las_cmmpcol_rau,1,se)
errors$las_cmmpw_rau=apply(errcol$las_cmmpcolw_rau,1,se)
errors$las_cmmpb_rau=apply(errcol$las_cmmpcolb_rau,1,se)

errors$rf_rau=apply(errcol$rfcol_rau,1,se)
errors$rfw_rau=apply(errcol$rfcolw_rau,1,se)
errors$rfb_rau=apply(errcol$rfcolb_rau,1,se)

#Vioplot

library(vioplot)
par(mar=c(4.5, 4.1, 6.5, 1.0))
labs=apply(expand.grid(c('CMMP','LM','ENET','CMMP+ENET','RF'),c('All','White','Black')),1,paste,collapse='\n')
labs=labs[c(1,6,11,2,7,12,3,8,13,4,9,14,5,10,15)]
vioplot(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
        error$las_cmmp,error$las_cmmpw,error$las_cmmpb,error$rf,error$rfw,error$rfb,names=labs,ylab='MSPE',main='MSPE averaged over 100 samples, Beta values',cex.axis=.8 ,col=rep(c('coral3','lightblue','orchid'),5))

#overlay points
holder=c(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
         error$las_cmmp,error$las_cmmpw,error$las_cmmpb, error$rf,error$rfw,error$rfb)
holders=c(errors$cmmp,errors$cmmpw,errors$cmmpb,errors$lm,errors$lmw,errors$lmb,errors$las,errors$lasw,errors$lasb,
          errors$las_cmmp,errors$las_cmmpw,errors$las_cmmpb,errors$rf,errors$rfw,errors$rfb)

holders_scaled=holders/max(holders)*.25

for(i in 1:(15*37)){
    arrows(floor(i/37)+1-holders_scaled[i],holder[i],floor(i/37)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
    pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
    sca=c(sca,0+c(-max(holders),-.5*max(holders),.5*max(holders),max(holders)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,3),cex.axis=.6,las=2)




vioplot(error$cmmp_m,error$cmmpw_m,error$cmmpb_m,error$lm_m,error$lmw_m,error$lmb_m,error$las_m,error$lasw_m,error$lasb_m,
        error$las_cmmp_m,error$las_cmmpw_m,error$las_cmmpb_m,error$rf_m,error$rfw_m,error$rfb_m, names=labs,ylab='MSPE',main='MSPE averaged over 100 samples, M values',cex.axis=.8,col=rep(c('coral3','lightblue','orchid'),5))

#overlay points
holder=c(error$cmmp_m,error$cmmpw_m,error$cmmpb_m,error$lm_m,error$lmw_m,error$lmb_m,error$las_m,error$lasw_m,error$lasb_m,
         error$las_cmmp_m,error$las_cmmpw_m,error$las_cmmpb_m,error$rf_m,error$rfw_m,error$rfb_m)
holders=c(errors$cmmp_m,errors$cmmpw_m,errors$cmmpb_m,errors$lm_m,errors$lmw_m,errors$lmb_m,errors$las_m,errors$lasw_m,errors$lasb_m,
          errors$las_cmmp_m,errors$las_cmmpw_m,errors$las_cmmpb_m,errors$rf_m,errors$rfw_m,errors$rfb_m )
holders_scaled=holders/max(holders)*.25

for(i in 1:(15*37)){
    arrows(floor(i/37)+1-holders_scaled[i],holder[i],floor(i/37)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
    pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
    sca=c(sca,0+c(-max(holders),-.5*max(holders),.5*max(holders),max(holders)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,3),cex.axis=.6,las=2)



vioplot(error$cmmp_rau,error$cmmpw_rau,error$cmmpb_rau,error$lm_rau,error$lmw_rau,error$lmb_rau,error$las_rau,error$lasw_rau,error$lasb_rau,
        error$las_cmmp_rau,error$las_cmmpw_rau,error$las_cmmpb_rau,error$rf_rau,error$rfw_rau,error$rfb_rau,names=labs,ylab='MSPE',main='MSPE averaged over 100 samples, RAU values',cex.axis=.8,col=rep(c('coral3','lightblue','orchid'),5))

#overlay points
holder=c(error$cmmp_rau,error$cmmpw_rau,error$cmmpb_rau,error$lm_rau,error$lmw_rau,error$lmb_rau,error$las_rau,error$lasw_rau,error$lasb_rau,
         error$las_cmmp_rau,error$las_cmmpw_rau,error$las_cmmpb_rau,error$rf_rau,error$rfw_rau,error$rfb_rau)
holders=c(errors$cmmp_rau,errors$cmmpw_rau,errors$cmmpb_rau,errors$lm_rau,errors$lmw_rau,errors$lmb_rau,errors$las_rau,errors$lasw_rau,errors$lasb_rau,
          errors$las_cmmp_rau,errors$las_cmmpw_rau,errors$las_cmmpb_rau,errors$rf_rau,errors$rfw_rau,errors$rfb_rau )

holders_scaled=holders/max(holders)*.25

for(i in 1:(15*37)){
    arrows(floor(i/37)+1-holders_scaled[i],holder[i],floor(i/37)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
    pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
    sca=c(sca,0+c(-max(holders),-.5*max(holders),.5*max(holders),max(holders)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,3),cex.axis=.6,las=2)


vioplot(error$cmmp_rau/100^2,error$cmmpw_rau/100^2,error$cmmpb_rau/100^2,error$lm_rau/100^2,error$lmw_rau/100^2,error$lmb_rau/100^2,error$las_rau/100^2,error$lasw_rau/100^2,error$lasb_rau/100^2,
        error$las_cmmp_rau/100^2,error$las_cmmpw_rau/100^2,error$las_cmmpb_rau/100^2, error$rf_rau/100^2,error$rfw_rau/100^2,error$rfb_rau/100^2,names=labs,ylab='MSPE',main='MSPE averaged over 100 samples, RAU values/10,000',cex.axis=.8,col=rep(c('coral3','lightblue','orchid'),5))

#overlay points
holder=c(error$cmmp_rau,error$cmmpw_rau,error$cmmpb_rau,error$lm_rau,error$lmw_rau,error$lmb_rau,error$las_rau,error$lasw_rau,error$lasb_rau,
         error$las_cmmp_rau,error$las_cmmpw_rau,error$las_cmmpb_rau,error$rf_rau,error$rfw_rau,error$rfb_rau )
holders=c(errors$cmmp_rau,errors$cmmpw_rau,errors$cmmpb_rau,errors$lm_rau,errors$lmw_rau,errors$lmb_rau,errors$las_rau,errors$lasw_rau,errors$lasb_rau,
          errors$las_cmmp_rau,errors$las_cmmpw_rau,errors$las_cmmpb_rau,errors$rf_rau,errors$rfw_rau,errors$rfb_rau )

holders_scaled=holders/max(holders)*.25

for(i in 1:(15*37)){
    arrows(floor(i/37)+1-holders_scaled[i],holder[i],floor(i/37)+1+holders_scaled[i],holder[i],code=1,angle=90,length=.0,col='grey50',lwd=.5)
}

pos=NULL
sca=NULL
for(i in 1:15){
    pos=c(pos,((1:15)[i]+c(-.25,-.125,.125,.25)))
    sca=c(sca,0+c(-max(holders),-.5*max(holders),.5*max(holders),max(holders)))
}
mtext('SE',at=c(0.1))
axis(3,at=pos,round(sca,1)/100,cex.axis=.6,las=2)



#points(c(rep(1,37),rep(2,37),rep(3,37),rep(4,37),rep(5,37),rep(6,37),rep(7,37),rep(8,37),rep(9,37),rep(10,37),rep(11,37),rep(12,37))[i],c(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
#error$las_cmmp,error$las_cmmpw,error$las_cmmpb)[i],col='black',pch=20,cex=.8)

#for(i in 1:(12*37)){
#points(c(rep(1,37),rep(2,37),rep(3,37),rep(4,37),rep(5,37),rep(6,37),rep(7,37),rep(8,37),rep(9,37),rep(10,37),rep(11,37),rep(12,37))[i],c(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
#error$las_cmmp,error$las_cmmpw,error$las_cmmpb)[i],col='black',pch=20,cex=(holders/median(holders))[i])
#}







#cut(sds[,1],breaks=seq(min(sds),max(sds),(max(sds)-min(sds))/4))



sds=data.frame(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,error$las_cmmp,error$las_cmmpw,error$las_cmmpb,error$rf,error$rfw,error$rfb)

##cols=NULL
##for(i in 1:dim(sds)[2]){
##    cols=c(cols,cut(sds[,i],breaks=seq(min(sds),max(sds),(max(sds)-min(sds))/4)))
##}

##points(c(rep(1,37),rep(2,37),rep(3,37),rep(4,37),rep(5,37),rep(6,37),rep(7,37),rep(8,37),rep(9,37),rep(10,37),rep(11,37),rep(12,37),rep(13,37),rep(14,37),rep(15,37))[i],c(error$cmmp,error$cmmpw,error$cmmpb,error$lm,error$lmw,error$lmb,error$las,error$lasw,error$lasb,
  ##                                                                                                                                        error$las_cmmp,error$las_cmmpw,error$las_cmmpb,error$rf,error$rfw,error$rfb)[i],pch=20,cex=.8,col=cols)

##names(sds)=labs
library(xtable)
xtable(sds)

