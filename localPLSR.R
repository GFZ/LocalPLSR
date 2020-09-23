#############################################################################################################
###                 localPLSR                 ###############################################################
#############################################################################################################
##
## Summary: Using spectral neighbours in a large scale soil spectral library to estimate soil properties in a dataset.
##
## Details on this method can be found in:
## Ward, K.J.; Chabrillat, S.; Neumann, C.; Foerster, S. A remote sensing adapted approach for soil organic 
## carbon prediction based on the spectrally clustered lucas soil database. Geoderma 2019, 353, 297-307 



###############################################################
## Author: Kathrin J. Ward
## GFZ Potsdam, section 1.4 Remote sensing and Geoinformatics
## kathrin.ward@gfz-potsdam.de
## 15th July 2020
###############################################################



dir <- "/.../" ## change directory

library(signal) ## version 0.7-6
library(pls) ## version 2.6-0
library(resemble) # version 1.2.2
library(e1071) ## version 1.6-8


################################################################################################################################
## 1. Prepare data including calibration/validation samples ####################################################################

## input data (here called input) with following columns: ID, OC, cal, spc
## input data should include both the spectral library as calibration samples and the validation samples (e.g. field samples)
## OC = numeric; organic carbon
## cal: binary; 1=sample belongs to calibration, 0=sample belongs to validation
## spc = numeric; spectra


## pre-processing: example for Savitzky-Golay smoothing and first derivative
rownames(input) <- input$ID
spc <- as.matrix(input[,4:ncol(input)]) # extract spectra: lines = samples; columns = spectral bands
d1_2 <- spc
for (i in 1:nrow(spc)) { d1_2[i,] <- sgolayfilt(spc[i,], p=2,n=11,m=1,ts=1) }  # p=polynom,n=window size,m=order derivative; adapt to level of noise! (package signal)
d1_2 <- cbind("ID"=row.names(d1_2),d1_2) ## add column with ID
d1 <- as.data.frame(d1_2)
d1[2:ncol(d1)] <- lapply(d1[2:ncol(d1)], function(x) as.numeric(as.character(x))) ## factor to numeric of spectra
d1[,1] <- as.character(d1[,1]) ## factor to character for ID
d1 <- merge(input[,c(1:3)],d1,by="ID") ## merge with OC and cal

matplot(colnames(d1[,4:ncol(d1)]),t(d1[1:20,4:ncol(d1)]),type="l",xlab="Wavelength [nm]",ylab="1st derivative of absorbance") ## test plot
abline(v=seq(500,2500,100),lty=2,col="lightgrey")

## subset for cal and val and save
d1_cal <- subset(d1,d1$cal==1)
d1_cal <- d1_cal[,-3] ## remove column cal
d1_val <- subset(d1,d1$cal==0)
d1_val <- d1_val[,-3] # remove column cal

setwd(dir)
save(d1_cal,file="CAL_SG1D_n11.Rdata")
save(d1_val,file="VAL_SG1D_n11.Rdata")

## calculate Mahalanobis distance based on first principal components of PCA
pca1d <- princomp(input[,3:ncol(input)]) ## PCA of spectra; hint: if there are more bands than samples, duplicate samples
cums <- cumsum(pca1d$sdev^2 / sum(pca1d$sdev^2)) 
cums[1:40] ## check number of PC to be used (e.g. which explain > 99.5% of variance)
pca1dn <- pca1d$scores[,1:7] ## select number of PC to be used
mdist <- fDiss(pca1dn,method="mahalanobis",center=T,scaled=T) ## calculate distance matrix, e.g. Mahalanobis (package resemble)
rownames(mdist) <- rownames(pca1dn) ## add rownames
colnames(mdist) <- rownames(pca1dn) ## add columnnames
mdist <- as.data.frame(mdist)

setwd(dir)
save(mdist,file="MDist_PCA_calval.Rdata")


###############################################################################################################################
##  2. Define localPLSR function #############################################################################################
LOCAL <- function(cal,val,diss,K,ncom,direc,meth,out,fixed,gt,minsamp,repe){
  
  # cal: potential calibration data set (e.g. LUCAS); ID, SOC, spectra; include rownames
  # val: validation samples; ID, SOC, spectra, include rownames
  # diss: distance matrix with cal and val samples included; include column and rownames (e.g. Mahalanobis distance matrix)
  # K: number of most similar samples to use if fixed=T OR threshold in dissimilarity measure if fixed=F; can be sequence of values
  #    e.g. seq(300,600,50) if fixed=T or seq(0.17,0.2,0.01) if fixed=F
  # ncom: maximum number of PLSR components to be used
  # meth: character; name of method, to be included in output filenames, choose whatever makes sense for you
  # out: out=T: plot graphs showing number of PLSR comp; out=F: don't plot
  # fixed: fixed=T: fixed number of samples, fixed=F: threshold in dissimilarity measure
  # gt: only applied if fixed=F; gt=T: greater than threshold, gt=F: lower than threshold
  # minsamp: only applied if fixed=F; minimum number of similar samples to be used, meaning a threshold in the dissimilarity 
  #          measure is used
  # repe: repeat model calibration x-times and average SOC predictions (ensemble modelling) -> more stable results for smaller 
  #       number of validation samples
  
  
  j=1 ## extra iterator connected to i but with intervals of 1 --> iterates different thresholds/ nb of nearest samples
  
  for(i in K){ ## thresholds in distance measure resp. number of nearest samples; i will equal all values of K one after the other
    
    par(mfrow=c(1,3),mar=c(30, 4, 1, 1) + 0.1) ## prepare plot parameters
    
    p <- as.numeric() ## container for variables
    res <- as.numeric()
    se <- as.numeric()
    res2 <- as.numeric()
    se2 <- as.numeric()
    st <- as.numeric()
    names <- list()
    ncoma <- as.numeric()
    ocp_bt <- as.numeric()
    oco_bt <- as.numeric()
    
    for(n in 1:nrow(val)){  ## do for each validation sample , nrow(val)
      
      if(n %% 50==0){cat(paste(n," "))} ## plot progress
      
      
      ## extract distance measure for validation sample #####
      #######################################################
      v1 <- row.names(val[n,])  # name of validation sample
      dist2 <- diss[which(row.names(diss)==v1),] ## extract distance value for validation sample
      dist <- t(dist2[,colnames(dist2) %in% rownames(cal)]) ## reduce to calibration subset
      
      ## find most similar K samples ########
      #######################################
      
      ## create vector for thresholds
      if(n==1){thresh <- as.numeric(rep(0,nrow(val)))
      names(thresh) <- row.names(val)}
      
      ## fixed number, like 300 samples
      if(fixed==T){ 
        simsam <- order(dist[,1])[1:i] ## positions of most similar K samples
        s1 <- cal[simsam,]  ## subset of most similar spectra
      }
      
      ## threshold 1
      if(fixed==F & gt==T){ ## gt threshold (greater than)
        simsam <- which(dist[,1]>i) ## positions of most similar K samples using threshold i from K
        s1 <- cal[simsam,]  ## subset of most similar spectra
        thresh[n] <- i
        
        ## if cal dataset too small use fixed number of cal samples
        if(nrow(s1)<minsamp){ 
          simsam <- order(dist[,1])[1:minsamp] ## positions of most similar K samples using the fixed number
          s1 <- cal[simsam,]  ## subset of most similar spectra
          thresh[n] <- threshold
        }
      }
      
      ## threshold 2 (used for Mahalanobis)
      if(fixed==F & gt==F){ ## lt threshold (lower than)
        simsam <- which(dist[,1]<i) ## positions of most similar K samples using threshold i from K
        s1 <- cal[simsam,]  ## subset of most similar spectra
        thresh[n] <- i
        
        ## if cal dataset too small use fixed number of cal samples
        if(nrow(s1)<minsamp){
          simsam <- order(dist[,1])[1:minsamp] ## positions of most similar K samples using fixed number
          s1 <- cal[simsam,]  ## subset of most similar spectra
          thresh[n] <- minsamp
        }
      }
      
      
      ## repeat model calibration x times and average SOC result to get stable results
      for(rr in 1:repe){
        
        if(rr==1){
          psoc <- as.numeric()
          ncomv <- as.numeric()}
        
        ## calibrate PLSR model ####
        ############################
        train2 <- data.frame(OC=s1[,2],spct=I(as.matrix(s1[,3:ncol(s1)])))  ## bring data to correct format for PLSR
        
        if(nrow(train2)<ncom+1){ ## if calibration subset is smaller than maximum number of latent variables
          ## started to run 10 CV in parallel 
          ## ! Check if 10 CPUs are available, otherwise reduce number !
          pls.options(parallel = 10) ## doesn't need to be stoped, automatically destroyed after usage
          m <- plsr(OC ~ spct, data=train2,ncomp=nrow(train2)-1,validation="CV")  ## CV with 10 segments (pls package)
        } else {
          pls.options(parallel = 10) 
          m <- plsr(OC ~ spct, data=train2,ncomp=ncom,validation="CV") ## CV with 10 segments
        }
        
        ## best number of model components #### -> dependent on CV segments -> random -> different number of model comp.
        #######################################
        # maximum of adjusted r2 #
        r2 <- R2(m,estimate="CV")
        r22 <- r2$val[,,1:(m$ncomp+1)]
        adjr2 <- as.numeric()
        for(v in 1:length(r22)){
          r23 <- r22[1:v]
          adjr2[v] <- 1-(1-r23[v])*((nrow(s1)-1)/(nrow(s1)-(length(r23)-1)-1))  ## (length(r23)-1) because first entry is intercept
        }
        ncom_adjR2 <- which.max(adjr2)
        
        # adj. Wold's R #  R = PRESS(m + 1)/PRESS(m); m=number of included model components #
        wr <- as.numeric()
        for(w in 1:(m$ncomp-1)){
          wr[w] <- m$validation$PRESS[,w+1]/m$validation$PRESS[,w]
        }
        ncom_wr <- which(wr>0.95)[1]
        
        # minimum of RMSEP #
        rmse1 <- RMSEP(m, estimate="CV")  ## calculate RMSE (pls package)
        rmse2 <- rmse1$val[,,1:m$ncomp+1]
        ncom3 <- match(min(rmse2),rmse2)  ## calculate minimum RMSE
        rmse3 <- data.frame(rmse=rmse2[1:ncom3],comp=c(1:ncom3)) ## cut RMSE vector until minimum and add column with comp. nb
        if(nrow(rmse3)>1){
          subrmse <- subset(rmse3, rmse3[,1]<min(rmse3[,1])+(sd(rmse3[,1])))  ##  /2; subset cut vector for all values < min+sd or min+sd/2
        } else {
          subrmse <- rmse3
        }
        ncom_RMSE <- min(subrmse[,2]) ## select minimum number of components of remaining vector 
        
        ## average outcomes
        ncom2 <- round((ncom_RMSE+ncom_wr+ncom_adjR2)/3,0)
        ncomv[rr] <- ncom2
        
        ## plot best number of model comp. ####
        if(out==T){
          
          plot(c(0:(length(adjr2)-1)),adjr2,type="l",xlab="Number of PLSR components",ylab="adj. R2", 
               main=paste("K = ",i,"; ",meth,sep=""))
          abline(v=ncom_adjR2-1)
          abline(v=ncom2,col="red") ## counts from 0 on
          
          #plot(rmse3$comp,rmse3$rmse,type="l",xlab="Number of PLSR components",ylab="RMSE")
          plot(c(1:length(rmse2)),rmse2,type="l",xlab="Number of PLSR components",ylab="RMSE")
          abline(v=ncom_RMSE)
          abline(v=ncom2,col="red")
          
          plot(c(1:length(wr)),wr,type="l",xlab="Number of PLSR components",ylab="adj. Wold's R (0.95)")
          abline(v=ncom_wr)
          abline(v=ncom2,col="red")
        }  
        
        ## predict validation data ####
        ###############################
        vv1 <- val[row.names(val)==v1,]
        pred2 <- data.frame(OC=vv1[,2],spct=I(as.matrix(vv1[,3:ncol(vv1)])))
        psoc[rr] <- predict(m,ncomp=ncom2,newdata=pred2)  ## vector with predicted SOC values
        
      } ## end repeat model calibration x times
      
     
      p[n] <- mean(psoc) ## mean of predicted SOC of all repetitions
      
      ## back-transformation of normalized SOC
      ## ! might need to be adapted !!
      ocp_bt[n] <- exp(p[n])-1  ## back transformation predicted (normalization was log(x+1))
      oco_bt[n] <- exp(pred2$OC)-1 ## back transformation observed
      
      ## back-transformed for measures with units
      res[n] <- apply(rbind(oco_bt[n],ocp_bt[n]),2,diff) ## residuals per sample: predicted - measured (underestimation neg number)
      se[n] <- res[n]^2 ## squared residuals
      ## not back-transformed for unitless measures
      res2[n] <- apply(rbind(pred2$OC,p[n]),2,diff) ## residuals per sample: predicted - measured (underestimation neg number)
      se2[n] <- res2[n]^2 ## squared residuals
      st[n] <- (pred2$OC-mean(val[,2]))^2  ## total sum of squares not back-transformed (observed - mean(observed))
      
      names[[n]] <- row.names(s1) ## names of samples used for calibration
      ncoma[n] <- mean(ncomv) ## mean of latent variables of all repetitions
      
      
    } ## end for(n in 1:nrow(val))
    
    ## prediction results ####
    ##########################
    
    if(j==1){
      result <- as.data.frame(matrix(NA,nrow=length(K),ncol=8,dimnames=list(c(K),
                       c("RMSEP","nRMSEP","rRMSEP","R2","RPD","RPIQ","Bias","NbComponents")))) ## create empty data frame
      reslist <- list()
    }
    
    reslist[[j]] <- res
    
    dat <- data.frame(se,se2,st)
    dat <- na.omit(dat)
    res <- na.omit(res)
    ncoma <- na.omit(ncoma)
    rmsep_norm <- sqrt(sum(dat$se2)/length(dat$se2))
    
    ## R2, RPD and RPIQ are calculated based on normalized data; for others data is back-transformed
    result[j,1] <- round(sqrt(sum(dat$se)/length(dat$se)),4) ## RMSEP
    result[j,2] <- round((result[j,1]/(max(oco_bt)-min(oco_bt)))*100,4)  ## nRMSEP
    result[j,3] <- round((result[j,1]/mean(oco_bt))*100,4)  ## rRMSEP
    result[j,4] <- round(1-(sum(dat$se2)/sum(dat$st)),4) ## R2 =1-(sum of squared errors, fitted)/(sum of squared, target)
    result[j,5] <- round(sd(val[,2])/rmsep_norm,4) ## RPD
    result[j,6] <- round(IQR(val[,2])/rmsep_norm,4) ## RPIQ
    result[j,7] <- round(sum(res)/length(res),4) ## bias 
    result[j,8] <- round(mean(ncoma),1) ## Number of model components
    
    cat(paste("\n ******************************************** \n",
              "\n","iteration",j,"/",length(K),"\n",
              meth,": R2:",round(result[j,4],2)," *** RMSE:",round(result[j,1],2)," *** RPD:",round(result[j,5],2),
              "\n","no. val samples",length(se),"\n"))
    
    j <- j+1
    
    ## save results
    setwd(direc)
    save(reslist,file=paste(Sys.time(),"_",meth,"_K",i,"th",minsamp,"residuals.Rdata",sep="")) ## residuals
    save(names,file=paste(Sys.time(),"_",meth,"_K",i,"th",minsamp,"names.Rdata",sep="")) ## names of calibration samples
    save(thresh,file=paste(Sys.time(),"_",meth,"_K",i,"th",minsamp,"thresholds.Rdata",sep="")) ## thresholds
    write.csv(result,file=paste(Sys.time(),"_",meth,"_K",i,"th",minsamp,".csv",sep="")) ## validation results
    save(ocp_bt,file=paste(Sys.time(),"_",meth,"_K",i,"th",minsamp,"predicted_bt.Rdata",sep="")) ## predicted SOC contents back-transformed
    
  }
  return(result)
  
}

###############################################################################################################################
## 3. Apply localPLSR to predict SOC of validation soil samples ###############################################################

## load prepared data
setwd(dir)
load("CAL_SG1D_n11.Rdata")  ## d1_cal; calibration set, e.g. LUCAS
rownames(d1_cal) <- d1_cal$ID
d1_cal2 <- d1_cal
skewness(d1_cal$OC) ## estimate skewness and find transformation towards normal distribution, here e.g. x=log(SOC+1) (e1071 package)
## ! back-transformation has to be adapted in local function, line 252 ff ! here e.g. SOC=exp(x)-1
OCnew <- d1_cal2$OC +1 ## add one to all SOC values as log(0)=Inf and log 0-1 is negative
d1_cal2$OC <- log(OCnew)
skewness(d1_cal2$OC)

load("VAL_SG1D_n11.Rdata") ## d1_val; validation set, e.g. field samples for which SOC content should be estimated
rownames(d1_val) <- d1_val$ID
d1_val2 <- d1_val
OCnew2 <- d1_val2$OC +1 
d1_val2$OC <- log(OCnew2)

load("MDist_PCA_calval.Rdata") ## mdist; data frame with Mahalanobis distance including cal and val


## apply local PLSR function ###
## example with a threshold in Mahalanobis distance of < 0.19 and a minumum of 200 samples used for calibration
t1 <- Sys.time() ## takes a while depending on size of dataset ;-)
result <- LOCAL(cal=d1_cal2,val=d1_val2,diss=mdist,K=rep(seq(0.19,0.19,0.01),1),ncom=40,direc=dir06,meth="MahalDist",out=F, 
                fixed=F,gt=F,minsamp=200,repe=100) 
Sys.time()-t1



