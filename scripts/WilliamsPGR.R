## recalculate Sahul-wide human population growth values from:
## Williams, A.N., 2013. A new population curve for prehistoric Australia. Proc. R. Soc. Lond. B 280, 20130486
## https://royalsocietypublishing.org/doi/full/10.1098/rspb.2013.0486
## approach based on number of dates per temporal bin
## dates calibrated and taphonomically corrected

# Corey Bradshaw
# October 2022

# libraries
library(ggplot2)
library(stratigraphr)
library(DescTools)
library(rcarbon)
library(stringr)

# import raw radiometric dates
RCdates <- read.csv("radiometricdates.csv")
head(RCdates)

# test if calculated SDs are correct # expect similar CVs
datSDraw <- subset(RCdates, calc==0)
datSDcalc <- subset(RCdates, calc==1)

datSDraw$CV <- 100*datSDraw$RCerr/datSDraw$RCage
datSDcalc$CV <- 100*datSDcalc$RCerr/datSDcalc$RCage
par(mfrow=c(1,2))
hist(log(datSDraw$CV), main="", xlab="log raw CV")
mean(datSDraw$CV)
hist(log(datSDcalc$CV), main="", xlab="log calculated CV")
mean(datSDcalc$CV)
par(mfrow=c(1,1))

# sort from oldest to youngest
RCsort <- RCdates[order(RCdates$RCage, decreasing=T),]
head(RCsort)

# calendar year calibration
# test
firstdate <- calibrate(x=RCsort$RCage[1], errors=RCsort$RCerr[1], ids=RCsort$labcode[1], calCurves="shcal20")
plot(firstdate, HPD=T, credMass=0.95)
summary(firstdate)

table(RCsort$sitetype)
table(RCsort$material)

# re-classify material by major class, type, species, specs
RCsort$matClass <- NA
RCsort$matClass[grep("marine", RCsort$material)] <- "marine"
RCsort$matClass[grep("Esuarine", RCsort$material)] <- "marine"
RCsort$matClass[grep("Anadara", RCsort$material)] <- "marine"
RCsort$matClass[grep("charcoal", RCsort$material)] <- "charcoal"
RCsort$matClass[grep("wood", RCsort$material)] <- "wood"
RCsort$matClass[grep("bark", RCsort$material)] <- "wood"
RCsort$matClass[grep("soil", RCsort$material)] <- "soil"
RCsort$matClass[grep("organic", RCsort$material)] <- "organic"
RCsort$matClass[grep("bone", RCsort$material)] <- "bone"
RCsort$matClass[grep("humic", RCsort$material)] <- "soil"
RCsort$matClass[grep("freshwater", RCsort$material)] <- "freshwater"
RCsort$matClass[grep("eggshell", RCsort$material)] <- "eggshell"
RCsort$matClass[grep("nodule", RCsort$material)] <- "nodule"
RCsort$matClass[grep("protein", RCsort$material)] <- "protein"
RCsort$matClass[grep("plant", RCsort$material)] <- "plant"
RCsort$matClass[grep("Peat", RCsort$material)] <- "plant"
RCsort$matClass[grep("otolith", RCsort$material)] <- "otolith"
RCsort$matClass[grep("carbon", RCsort$material)] <- "carbon"
RCsort$matClass[grep("ash", RCsort$material)] <- "ash"
RCsort$matClass <- ifelse(RCsort$material == "carbonised Eucalypt sp. ", "wood", RCsort$matClass)
RCsort$matClass <- ifelse(RCsort$material == "wood charcoal", "charcoal", RCsort$matClass)
RCsort$matClass <- ifelse(RCsort$material == "soil organic matter", "soil", RCsort$matClass)
RCsort$matClass <- ifelse(RCsort$material == "organic carbon (from animal bone)", "bone", RCsort$matClass)


RCsort$matType <- NA
RCsort$matType[grep("shell", RCsort$material)] <- "shell"
RCsort$matType[grep("carbonate", RCsort$material)] <- "carbonate"
RCsort$matType[grep("humic", RCsort$material)] <- "humics"
RCsort$matType[grep("blood", RCsort$material)] <- "blood"
RCsort$matType[grep("bark", RCsort$material)] <- "bark"
RCsort$matType[grep("carbonised", RCsort$material)] <- "carbonised"
RCsort$matType[grep("apatite", RCsort$material)] <- "apatite"
RCsort$matType[grep("residue", RCsort$material)] <- "residue"
RCsort$matType[grep("collagen", RCsort$material)] <- "collagen"
RCsort$matType[grep("gelatin", RCsort$material)] <- "gelatin"
RCsort$matType[grep("charred", RCsort$material)] <- "charred"
RCsort$matType[grep("Peat", RCsort$material)] <- "peat"
RCsort$matType[grep("Cellana", RCsort$material)] <- "shell"
RCsort$matType <- ifelse(RCsort$material == "freshwater mussel", "shell", RCsort$matType)
RCsort$matType <- ifelse(RCsort$material == "freshwater mussel (Velesunio sp)", "shell", RCsort$matType)
RCsort$matType <- ifelse(RCsort$material == "wood charcoal", "wood", RCsort$matType)
RCsort$matType <- ifelse(RCsort$material == "soil organic matter", "organic", RCsort$matType)
RCsort$matType <- ifelse(RCsort$material == "organic carbon (from animal bone)", "organic", RCsort$matType)


RCsort$sp <- NA
RCsort$sp[grep("Subinella", RCsort$material)] <- "Subinella"
RCsort$sp[grep("Subninella", RCsort$material)] <- "Subinella"
RCsort$sp[grep("Turbo", RCsort$material)] <- "Turbo"
RCsort$sp[grep("Saccostrea", RCsort$material)] <- "Saccostrea"
RCsort$sp[grep("Pyrazus", RCsort$material)] <- "Pyrazus"
RCsort$sp[grep("Plebidonax", RCsort$material)] <- "Plebidonax"
RCsort$sp[grep("Patella", RCsort$material)] <- "Patella"
RCsort$sp[grep("Paphies", RCsort$material)] <- "Paphies"
RCsort$sp[grep("Ostrea", RCsort$material)] <- "Ostrea"
RCsort$sp[grep("oyster", RCsort$material)] <- "oyster"
RCsort$sp[grep("Notohaliotus", RCsort$material)] <- "Notohaliotus"
RCsort$sp[grep("Nerita", RCsort$material)] <- "Nerita"
RCsort$sp[grep("Mylitus", RCsort$material)] <- "Mylitus"
RCsort$sp[grep("Haliotis", RCsort$material)] <- "Haliotis"
RCsort$sp[grep("Donax", RCsort$material)] <- "Donax"
RCsort$sp[grep("Dicathais", RCsort$material)] <- "Dicathais"
RCsort$sp[grep("Dicatgaus", RCsort$material)] <- "Dicathais"
RCsort$sp[grep("catrut", RCsort$material)] <- "Dicathais"
RCsort$sp[grep("Crassostera", RCsort$material)] <- "Crassostera"
RCsort$sp[grep("Cellana", RCsort$material)] <- "Cellana"
RCsort$sp[grep("Cabesatana", RCsort$material)] <- "Cabesatana"
RCsort$sp[grep("cabesatana", RCsort$material)] <- "Cabesatana"
RCsort$sp[grep("Brachidontes", RCsort$material)] <- "Brachidontes"
RCsort$sp[grep("Austomylitus", RCsort$material)] <- "Austomylitus"
RCsort$sp[grep("Anadara", RCsort$material)] <- "Anadara"
RCsort$sp[grep("Brachidontes", RCsort$material)] <- "Brachidontes"
RCsort$sp[grep("mussel", RCsort$material)] <- "mussel"
RCsort$sp[grep("Velesunio", RCsort$material)] <- "Velesunio"
RCsort$sp[grep("Unio", RCsort$material)] <- "Unio"
RCsort$sp[grep("Eucalypt", RCsort$material)] <- "Eucalyptus"
RCsort$sp[grep("Acacia", RCsort$material)] <- "Acacia"
RCsort$sp[grep("Casuarina", RCsort$material)] <- "Casuarina"
RCsort$sp[grep("Melaleuca", RCsort$material)] <- "Melaleuca"
RCsort$sp[grep("Dugong", RCsort$material)] <- "Dugong"
RCsort$sp[grep("human", RCsort$material)] <- "Homo"
RCsort$sp[grep("Vombatus", RCsort$material)] <- "Vombatus"
RCsort$sp[grep("Macropus", RCsort$material)] <- "Macropus"
RCsort$sp[grep("emu", RCsort$material)] <- "Dromaius"

head(RCsort)
table(RCsort$matClass)
RCsort[which(is.na(RCsort$matClass)==T),]

table(RCsort$matType)
RCsort[which(is.na(RCsort$matType)==T & RCsort$material != RCsort$matClass & RCsort$matClass != "wood" & RCsort$matClass != "charcoal"),]

table(RCsort$sp)
plot(RCsort$RCage, RCsort$RCerr, pch=19, cex=0.7, xlab="radiocarbon years BP", ylab="SD")

# make ids unique
length(RCsort$labcode)
length(unique(RCsort$labcode))
sort(table(RCsort$labcode), decreasing=T)[1:20]

RCsort$labcode[which(RCsort$labcode == "Not given")] <- paste("not given ", seq(1,34,1), sep="")
RCsort$labcode[which(RCsort$labcode == "V-5")] <- paste("V-5", letters[1:4], sep="")
RCsort$labcode[which(RCsort$labcode == "ANU 404A")] <- paste("ANU 404A", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "ANU 404B")] <- paste("ANU 404B", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "ANU 405")] <- paste("ANU 405", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "ANU 422")] <- paste("ANU 422", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "ANU-1038")] <- paste("ANU-1038", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "ANU-1242")] <- paste("ANU-1242", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "ANU-648")] <- paste("ANU-648", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "SUA-1109")] <- paste("SUA-1109", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "SUA-1110")] <- paste("SUA-1110", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "SUA-1166")] <- paste("SUA-1166", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "SUA-1167")] <- paste("SUA-1167", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "SUA-2173")] <- paste("SUA-2173", letters[1:2], sep="")
RCsort$labcode[which(RCsort$labcode == "Wk-16145")] <- paste("Wk-16145", letters[1:2], sep="")

sort(table(RCsort$labcode), decreasing=T)[1:20]

## only dates in 4 ka to 40 ka range
RCsortlim <- subset(RCsort, RCage > 4000)
RCsortlim$CV <- 100*(RCsortlim$RCerr/RCsortlim$RCage)
range(RCsortlim$CV)
hist(RCsortlim$CV)
dim(RCsortlim)
head(RCsortlim)

## calibrate
# CALdates <- oxcalCalibrate(RCsort$RCage, RCsort$RCerr, RCsort$labcode)
RCnotmarine <- RCsortlim[-which(RCsortlim$matClass == "marine" | RCsortlim$matClass == "otolith"),]
factor(RCnotmarine$matClass)
dim(RCnotmarine)
table(RCnotmarine$matClass)
CALnotmarine <- calibrate(x=RCnotmarine$RCage, errors=RCnotmarine$RCerr, ids=RCnotmarine$labcode, calCurves="shcal20")
CALnotmarine.summ <- summary(CALnotmarine)

RCmarine <- subset(RCsortlim, matClass=="marine" | matClass=='otolith')
dim(RCmarine)
table(RCmarine$matClass)
CALmarine <- calibrate(x=RCmarine$RCage, errors=RCmarine$RCerr, ids=RCmarine$labcode, calCurves="marine20", verbose=T)
CALmarine.summ <- summary(CALmarine)
dim(CALmarine.summ)
head(CALnotmarine.summ)
head(CALmarine.summ)

CALdates <- rbind(CALnotmarine.summ[,c(1:4,9)], CALmarine.summ)
head(CALdates)

CALdates.sort <- CALdates[order(CALdates$MedianBP, decreasing=T),]
head(CALdates.sort)

CALdates.sort[which(CALdates.sort$DateID=="V-64"),]

CALdatMD <- CALdatUP <- CALdatLO <- CALdatSD <- rep(NA, length(CALdates.sort$OneSigma_BP_1))
for (i in 1:length(CALdates.sort$OneSigma_BP_1)) {
  CALdatMD[i] <- CALdates.sort$MedianBP[i]
  CALdatUP[i] <- as.numeric(str_split(CALdates.sort$OneSigma_BP_1[i], " to ", simplify=T))[1]
  CALdatLO[i] <- as.numeric(str_split(CALdates.sort$OneSigma_BP_1[i], " to ", simplify=T))[2]
  CALdatSD[i] <- (mean(c((CALdatUP[i] - CALdatMD[i]), (CALdatMD[i] - CALdatLO[i]))))
  
  if (100*(((CALdatUP[i] - CALdatLO[i])/2)/CALdatMD[i]) < quantile(RCsortlim$CV, probs=0.125, na.rm=T)) {
    if (is.na(as.numeric(str_split(CALdates.sort$OneSigma_BP_2[i], " to ", simplify=T))[1]) == F) {
      CALdatUP[i] <- as.numeric(str_split(CALdates.sort$OneSigma_BP_2[i], " to ", simplify=T))[1]
      CALdatLO[i] <- as.numeric(str_split(CALdates.sort$OneSigma_BP_2[i], " to ", simplify=T))[2]
      CALdatSD[i] <- (mean(c((CALdatUP[i] - CALdatMD[i]), (CALdatMD[i] - CALdatLO[i]))))
    }
  }
  
  if (100*(((CALdatUP[i] - CALdatLO[i])/2)/CALdatMD[i]) < quantile(RCsortlim$CV, probs=0.125, na.rm=T)) {
    if (is.na(as.numeric(str_split(CALdates.sort$TwoSigma_BP_1[i], " to ", simplify=T))[1]) == F) {
      CALdatUP[i] <- as.numeric(str_split(CALdates.sort$TwoSigma_BP_1[i], " to ", simplify=T))[1]
      CALdatLO[i] <- as.numeric(str_split(CALdates.sort$TwoSigma_BP_1[i], " to ", simplify=T))[2]
      CALdatSD[i] <- (mean(c((CALdatUP[i] - CALdatMD[i]), (CALdatMD[i] - CALdatLO[i]))))/1.96
    }
  }
  
}

CAL.out <- data.frame(CALdates.sort$DateID, CALdatMD, CALdatSD, CALdatLO, CALdatUP)
colnames(CAL.out)[1] <- "labcode"
head(CAL.out)
head(CAL.out[order(CAL.out$CALdatSD, decreasing=F),])

CALsort <- merge(RCsortlim, CAL.out, by="labcode")
head(CALsort)
tail(CALsort)
hist(CALsort$CALdatSD)

plot((CALsort$CALdatMD), log10(CALsort$CALdatSD), pch=19, cex=0.7, xlab="calendar years BP", ylab="log10 SD", xlim=c(5000,40000), ylim=c(0,4))
fitCAL <- lm(log10(CALsort$CALdatSD[which(CALsort$CALdatMD >= 5000 & CALsort$CALdatMD <= 40000)]) ~
               CALsort$CALdatMD[which(CALsort$CALdatMD >= 5000 & CALsort$CALdatMD <= 40000)])
summary(fitCAL)
abline(fitCAL, lty=2, lwd=2, col="red")
10^as.numeric(1000*coef(fitCAL)[2]) # for every 1000 years in the past, SD increases by this many times

  
# set up x-year intervals
intSTvec <- seq(0,39800,200)
intENvec <- seq(200,40000,200)
intMDvec <- (intSTvec + intENvec)/2

# taphonomic correction (Williams 2012)
# nt = 2.107e7 * (t + 2754)^-1.526; where nt = number of radiocarbon dates surviving; t = time
taph.corr <- read.csv("corrWilliams.csv")

corr.intp <- approx(x=taph.corr$CALt, y=taph.corr$corr, xout=intMDvec) 
corr.UPintp <- approx(x=taph.corr$CALt, y=taph.corr$corrUP, xout=intMDvec) 
corr.LOintp <- approx(x=taph.corr$CALt, y=taph.corr$corrLO, xout=intMDvec) 
taph.corr.interp <- data.frame(intMDvec, corr.intp$y, corr.UPintp$y, corr.LOintp$y)
colnames(taph.corr.interp) <- c("t","corr","corrUP","corrLO")
taph.corr.interp$sd <- (mean(c((taph.corr.interp$corrUP - taph.corr.interp$corr), (taph.corr.interp$corr - taph.corr.interp$corrLO))))/1.96
head(taph.corr.interp)

iter <- 1000
GRann.mat <- r.mat <- nDatesCorr.mat <- matrix(data=NA, ncol=length(intMDvec), nrow=iter)
cor.it <- rep(NA,iter)

for (i in 1:iter) {
  dates.it <- rnorm(dim(CALsort)[1], CALsort$CALdatMD, CALsort$CALdatSD)
  
  num.dates <- rep(0,length(intSTvec))
  for (s in 1:length(intSTvec)) {
    num.dates[s] <- length(which(dates.it >= intSTvec[s] & dates.it < intENvec[s]))
  } # end s
  
  # Williams' approach
  # apply df=25 smoothing spline to time series
  num.dates.ss <- smooth.spline(x=intMDvec, y=num.dates, df=25)
  
  # taphonomic correction
  #num.dates.taphc <- num.dates.ss$y / rnorm(length(taph.corr.interp$t), taph.corr.interp$corr, taph.corr.interp$sd)
  nDatesCorr.mat[i,] <- num.dates.ss$y / taph.corr.interp$corr 
    
  # mean annual growth rate
  GRann.mat[i,2:(length(intMDvec))] <- (diff(nDatesCorr.mat[i,])*0.5)/nDatesCorr.mat[i,][1:(length(nDatesCorr.mat[i,])-1)]
  
  # instantaneous exponential rate of increase
  r.mat[i,2:(length(intMDvec))] <- log(nDatesCorr.mat[i,][2:length(nDatesCorr.mat[i,])] / nDatesCorr.mat[i,][1:(length(nDatesCorr.mat[i,])-1)])
  
  # correlation between GRann & r
  cor.it[i] <- cor(GRann.mat[i,], r.mat[i,], method="pearson", use="na.or.complete")
  
} # end i

hist(cor.it, main="", xlab="Pearson's correlation GRann vs. r")
hist(Logit(cor.it), main="", xlab="logit Pearson's correlation GRann vs. r")
median(cor.it, na.rm=T)
quantile(cor.it, probs=0.025, na.rm=T)
quantile(cor.it, probs=0.975, na.rm=T)


GRannMD <- apply(GRann.mat, MARGIN=2, median, na.rm=T)
GRannLO <- apply(GRann.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
GRannUP <- apply(GRann.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
rMD <- apply(r.mat, MARGIN=2, median, na.rm=T)
rLO <- apply(r.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
rUP <- apply(r.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
nDatesMD <- apply(nDatesCorr.mat, MARGIN=2, median, na.rm=T)
nDatesLO <- apply(nDatesCorr.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
nDatesUP <- apply(nDatesCorr.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)

out.dat <- data.frame(intMDvec, GRannMD, GRannUP, GRannLO, rMD, rUP, rLO, nDatesMD, nDatesUP, nDatesLO)
head(out.dat)

# GRann
plot(out.dat$intMDvec, out.dat$GRannMD, xlab="calendar years BP", ylab="rate of change (GRann)", type="l", ylim=c(-0.2,0.25), xlim=c(5000,40000))
abline(h=0, lwd=2, lty=3)
lines(out.dat$intMDvec, out.dat$GRannLO, lty=2, col="grey")
lines(out.dat$intMDvec, out.dat$GRannUP, lty=2, col="grey")

# r
plot(out.dat$intMDvec, out.dat$rMD, xlab="calendar years BP", ylab="rate of change (r)", type="l", ylim=c(-0.2,0.4), xlim=c(5000,40000))
lines(out.dat$intMDvec, out.dat$rLO, lty=2, col="pink")
lines(out.dat$intMDvec, out.dat$rUP, lty=2, col="pink")
abline(h=0, lwd=2, lty=3)

# number of dates
plot(out.dat$intMDvec, out.dat$nDatesMD, xlab="calendar years BP", ylab="taphonomically corrected number of dates", type="l", ylim=c(0,150), xlim=c(5000,40000))
lines(out.dat$intMDvec, out.dat$nDatesLO, lty=2, col="grey")
lines(out.dat$intMDvec, out.dat$nDatesUP, lty=2, col="grey")

