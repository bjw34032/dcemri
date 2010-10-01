###################################################
### chunk number 1: preliminaries
###################################################
library("oro.dicom")
library("XML")
library("oro.nifti")
library("dcemriS4")
library("bitops")
library("minpack.lm")
library("splines")
## options(niftiAuditTrail=FALSE)
options(prompt="R> ")
options(width=76)


###################################################
### chunk number 2: doubleanglemethod
###################################################
f60 <- file.path("nifti", "SDAM_ep2d_60deg_26slc.nii.gz")
sdam60 <- readNIfTI(system.file(f60, package="dcemriS4"))
f120 <- file.path("nifti", "SDAM_ep2d_120deg_26slc.nii.gz")
sdam120 <- readNIfTI(system.file(f120, package="dcemriS4"))
sdam.image <- rowMeans(dam(sdam60, sdam120, 60), dims=3)
mask <- (rowSums(sdam60, dims=3) > 500) # nifti(rowSums(sdam60, dims=3) > 500)


###################################################
### chunk number 3: figure5-png
###################################################
png("sdam.png", width=400, height=400)


###################################################
### chunk number 4: figure5-code
###################################################
# A smooth version of "sdam.image"
fsmooth <- file.path("nifti", "SDAM_smooth.nii.gz")
SDAM <- readNIfTI(system.file(fsmooth, package="dcemriS4"))
overlay(sdam120, ifelse(mask, SDAM, NA), z=13, zlim.x=range(sdam120), 
        zlim.y=c(0.5,1.5), plot.type="single")


###################################################
### chunk number 5: figure5-dev.off
###################################################
dev.off()


###################################################
### chunk number 6: t1estimation
###################################################
alpha <- c(5,10,20,25,15)
nangles <- length(alpha)
fnames <- file.path("nifti", paste("fl3d_vibe-", alpha, "deg.nii.gz", sep=""))
X <- Y <- 64
Z <- 36
flip <- fangles <- array(0, c(X,Y,Z,nangles))
for (w in 1:nangles) {
  vibe <- readNIfTI(system.file(fnames[w], package="dcemriS4"))
  flip[,,1:nsli(vibe),w] <- vibe
  fangles[,,,w] <- array(alpha[w], c(X,Y,Z))
}
TR <- 4.22 / 1000 # in seconds
fanglesB1 <- fangles * array(SDAM, c(X,Y,Z,nangles))
zi <- 10:13
maskzi <- mask
maskzi[,,(! 1:Z %in% zi)] <- FALSE
R1 <- R1.fast(flip, maskzi, fanglesB1, TR, verbose=TRUE)


###################################################
### chunk number 7: figure6-png
###################################################
png("t1_phantom.png", width=400, height=400)


###################################################
### chunk number 8: figure6-code
###################################################
overlay(vibe, 1/R1$R10[,,1:nsli(vibe)], z=13, zlim.x=c(0,1024), 
        zlim.y=c(0,2.5), plot.type="single")


###################################################
### chunk number 9: figure6-dev.off
###################################################
dev.off()


###################################################
### chunk number 10: FSLmask
###################################################
t1pmask <- readNIfTI(system.file("nifti/t1_phantom_mask.nii.gz", 
                                 package="dcemriS4"))
pmask <- nifti(array(t1pmask[,,25], dim(t1pmask))) # repeat slice #25


###################################################
### chunk number 11: figure7
###################################################
T1 <- c(.484,.350,1.07,.648,.456,1.07,.660,1.543,1.543,.353)
par(mfrow=c(1,1), mar=c(5,4,4,2)+.1)
boxplot(split(1/drop(R1$R10), as.factor(drop(pmask)))[-1], 
        ylim=c(0,2.5), xlab="Region of Interest", ylab="T1 (seconds)")
points(1:10, T1, col=rainbow(10), pch=16, cex=2)


###################################################
### chunk number 12: buckley.aif
###################################################
data("buckley")
aifparams <- with(buckley, orton.exp.lm(time.min, input))
fit.aif <- with(aifparams, aif.orton.exp(buckley$time.min, AB, muB, AG, muG))


###################################################
### chunk number 13: figure8
###################################################
with(buckley, plot(time.min, input, type="l", lwd=2, xlab="Time (minutes)", 
                   ylab=""))
with(buckley, lines(time.min, fit.aif, lwd=2, col=2))
legend("topright", c("Simulated AIF", "Estimated AIF"), lwd=2, col=1:2, 
       bty="n")


###################################################
### chunk number 14: buckley.kinetic
###################################################
xi <- seq(5, 300, by=5)
img <- array(t(breast$data)[,xi], c(13,1,1,60))
time <- buckley$time.min[xi]
aif <- buckley$input[xi]
mask <- array(TRUE, dim(img)[1:3])
aifparams <- orton.exp.lm(time, aif)
fit <- dcemri.lm(img, time, mask, model="orton.exp",
                 aif="user", user=aifparams)


###################################################
### chunk number 15: figure9
###################################################
par(mfrow=c(4,4), mar=c(5,4,4,2)/1.25, mex=0.8)
for (x in 1:nrow(img)) {
  plot(time, img[x,1,1,], ylim=range(img), xlab="Time (minutes)",
       ylab="", main=paste("Series", x))
  kinparams <- with(fit, c(vp[x,1,1], ktrans[x,1,1], kep[x,1,1]))
  lines(time, model.orton.exp(time, aifparams[1:4], kinparams), 
        lwd=1.5, col=2)
}


###################################################
### chunk number 16: RIDER_Neuro_MRI+pre1
###################################################
perf <- paste("281949", "19040721", "perfusion", sep="_")
mask <- readANALYZE(paste(perf, "mask", sep="_"))
mask <- ifelse(mask > 0, TRUE, FALSE)
dynamic <- readNIfTI(perf)
TR <- 4.43 / 1000 # taken from CSV file
dangle <- 25      # taken from CSV file
fflip <- list.files(pattern="ax[0-9]")
fangles <- as.numeric(sub(".*ax([0-9]+).*", "\\1", fflip))
flip <- array(NA, c(dim(dynamic)[1:3], length(fangles)))
for (fa in 1:length(fangles)) {
  flip[,,,fa] <- readNIfTI(fflip[fa])
}


###################################################
### chunk number 17: RIDER_Neuro_MRI+pre2 eval=FALSE
###################################################
## ca <- CA.fast(dynamic, mask, dangle, flip, fangles, TR, verbose=TRUE)
## writeNIfTI(ca$M0, paste(perf, "m0", sep="_"))
## writeNIfTI(ca$R10, paste(perf, "r10", sep="_"))
## writeNIfTI(ca$conc, paste(perf, "gdconc", sep="_"))


###################################################
### chunk number 18: RIDER_Neuro_MRI+lm0
###################################################
acqtimes <- str2time(unique(sort(scan("rawtimes.txt", quiet=TRUE))))$time
acqtimes <- (acqtimes - acqtimes[9]) / 60 # minutes
conc <- readNIfTI(paste(perf, "gdconc", sep="_"))


###################################################
### chunk number 19: RIDER_Neuro_MRI+lm1 eval=FALSE
###################################################
## fit.lm <- dcemri.lm(conc, acqtimes, mask, model="extended", 
##                     aif="fritz.hansen", verbose=TRUE)
## writeNIfTI(fit.lm$ktrans, paste(perf, "ktrans", sep="_"))


###################################################
### chunk number 20: RIDER_Neuro_MRI+lm2 eval=FALSE
###################################################
## writeNIfTI(fit.lm$kep, paste(perf, "kep", sep="_"))
## writeNIfTI(fit.lm$vp, paste(perf, "vp", sep="_"))
## writeNIfTI(fit.lm$ve, paste(perf, "ve", sep="_"))
## writeNIfTI(fit.lm$sse, paste(perf, "sse", sep="_"))
## rm(fit.lm)


###################################################
### chunk number 21: RIDER_Neuro_MRI+lm3
###################################################
png(file=paste(paste(perf, "ktrans", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)


###################################################
### chunk number 22: RIDER_Neuro_MRI+lm4
###################################################
fit.lm <- list(ktrans=readNIfTI(paste(perf, "ktrans", sep="_")))


###################################################
### chunk number 23: RIDER_Neuro_MRI+lm5
###################################################
overlay(dynamic, ifelse(fit.lm$ktrans < 0.1, fit.lm$ktrans, NA),
        w=11, zlim.x=c(32,512), zlim.y=c(0,0.1))


###################################################
### chunk number 24: RIDER_Neuro_MRI+lm6
###################################################
dev.off()
fit.lm$kep <- readNIfTI(paste(perf, "kep", sep="_"))
fit.lm$vp <- readNIfTI(paste(perf, "vp", sep="_"))
fit.lm$ve <- readNIfTI(paste(perf, "ve", sep="_"))
fit.lm$sse <- readNIfTI(paste(perf, "sse", sep="_"))
xx <- 41:220
yy <- 21:220
png(file=paste(paste(perf, "kep", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
overlay(as(dynamic[xx,yy,,], "nifti"), 
        ifelse(fit.lm$kep[xx,yy,] < 1.25, fit.lm$kep[xx,yy,], NA),
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,1.25), plot.type="single")
dev.off()
png(file=paste(paste(perf, "vp", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
overlay(as(dynamic[xx,yy,,], "nifti"), 
        ifelse(fit.lm$vp[xx,yy,] < 0.3, fit.lm$vp[xx,yy,], NA),
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,0.03), plot.type="single")
dev.off()
png(file=paste(paste(perf, "ve", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
overlay(as(dynamic[xx,yy,,], "nifti"), 
        ifelse(fit.lm$ve[xx,yy,] < 0.3, fit.lm$ve[xx,yy,], NA),
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,0.3), plot.type="single")
dev.off()
png(file=paste(paste(perf, "sse", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
overlay(as(dynamic[xx,yy,,], "nifti"), 
        ifelse(fit.lm$sse[xx,yy,] < 0.05, fit.lm$sse[xx,yy,], NA),
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,0.05), plot.type="single")
dev.off()


###################################################
### chunk number 25: RiderNeuroMRI+map1 eval=FALSE
###################################################
## fit.map <- dcemri.map(conc, acqtimes, mask, model="extended", 
##                       aif="fritz.hansen", ab.ktrans=c(log(0.05),1),
##                       ab.kep=c(log(0.7),1), ab.vp=c(1,19),
##                       multicore=TRUE)
## writeNIfTI(fit.map$ktrans, paste(perf, "ktrans", "map", sep="_"))


###################################################
### chunk number 26: RiderNeuroMRI+map2
###################################################
fit.map <- list(ktrans=readNIfTI(paste(perf, "ktrans", "map", sep="_")))
png(file=paste(paste(perf, "ktrans", "map", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
overlay(as(dynamic[xx,yy,,], "nifti"), 
        ifelse(fit.map$ktrans[xx,yy,] < 0.1, fit.map$ktrans[xx,yy,], NA),
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,0.1), plot.type="single")
dev.off()
png(file=paste(paste(perf, "ktrans", "compare", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
plot(fit.lm$ktrans, fit.map$ktrans, xlim=c(0,0.3), ylim=c(0,0.3),
     xlab=expression(paste(K^{trans}, " (Levenberg-Marquardt)")), 
     ylab=expression(paste(K^{trans}, " (Maximum A Posteriori)")), 
     pch=19)
abline(0, 1, col="red", lwd=2)
dev.off()


###################################################
### chunk number 27: RiderNeuroMRI+map3
###################################################
x <- apply(mask, 1, sum) != 0
y <- apply(mask, 2, sum) != 0
print(c(sum(mask[x,y,]), 
        sum(is.na(fit.lm$ktrans[x,y,])),
        sum(is.na(fit.map$ktrans[x,y,]))))


###################################################
### chunk number 28: RiderNeuroMRI+Bayes0 eval=FALSE
###################################################
## fit.bayes <- dcemri.bayes(conc, acqtimes, mask, model="extended", 
##                           aif="fritz.hansen", ab.ktrans=c(log(0.05),1), 
##                           ab.kep=c(log(0.7),1), ab.vp=c(1,19),
##                           multicore=TRUE)
## writeNIfTI(fit.bayes$ktrans, paste(perf, "ktrans", "bayes", sep="_"))


###################################################
### chunk number 29: RiderNeuroMRI+Bayes1 eval=FALSE
###################################################
## writeNIfTI(fit.bayes$ktranserror, paste(perf, "ktrans", "bayes", "sd", sep="_"))
## writeNIfTI(fit.bayes$kep, paste(perf, "kep", "bayes", sep="_"))
## writeNIfTI(fit.bayes$keperror, paste(perf, "kep","bayes", "sd", sep="_"))
## rm(fit.bayes)


###################################################
### chunk number 30: RiderNeuroMRI+Bayes2
###################################################
fit.bayes <- list(ktrans=readNIfTI(paste(perf, "ktrans","bayes", sep="_")))
png(file=paste(paste(perf, "ktrans", "bayes", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
overlay(as(dynamic[xx,yy,,], "nifti"), 
        ifelse(fit.bayes$ktrans[xx,yy,] < 0.1, fit.bayes$ktrans[xx,yy,], NA), 
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,0.1), plot.type="single")
dev.off()
fit.bayes$ktranserror <- readNIfTI(paste(perf, "ktrans","bayes","sd", sep="_"))
png(file=paste(paste(perf, "ktrans", "bayes", "sd", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
overlay(dynamic[xx,yy,,], 
        ifelse(fit.bayes$ktranserror[xx,yy,] < 0.0075, fit.bayes$ktranserror[xx,yy,], NA),
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,0.0075), plot.type="single")
dev.off()
png(file=paste(paste(perf, "ktrans", "bayes", "cv", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
overlay(as(dynamic[xx,yy,,], "nifti"), 
        ifelse(fit.bayes$ktranserror[xx,yy,] / fit.bayes$ktrans[xx,yy,] < 0.2, 
               fit.bayes$ktranserror[xx,yy,] / fit.bayes$ktrans[xx,yy,], NA),
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,0.2), plot.type="single")
dev.off()


###################################################
### chunk number 31: RiderNeuroMRI+Spline
###################################################
mask.spline <- array(FALSE, dim(mask))
z <- 7
mask.spline[,,z] <- mask[,,z]
fit.spline <- dcemri.spline(conc[,,,-(1:8)], acqtimes[-(1:8)], mask.spline,
                            model="weinmann", aif="fritz.hansen", 
                            ab.tauepsilon=c(1/10,1/1000), 
                            ab.hyper=c(0.01,0.01), multicore=TRUE, 
                            nlr=TRUE)
writeNIfTI(fit.spline$ktrans, paste(perf, "ktrans","spline", sep="_"))
writeNIfTI(fit.spline$Fp, paste(perf, "Fp","spline", sep="_"))


###################################################
### chunk number 32: RiderNeuroMRI+Spline2
###################################################
fit.spline <- list(ktrans=readNIfTI(paste(perf, "ktrans","spline", sep="_")))
png(file=paste(paste(perf, "ktrans", "spline", sep="_"), "png", sep="."), 
    width=2*480, height=2*480)
overlay(as(dynamic[xx,yy,,], "nifti"), 
        ifelse(fit.spline$ktrans[xx,yy,] < 0.1, fit.spline$ktrans[xx,yy,], NA),
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,0.1), plot.type="single")
dev.off()
fit.spline$Fp <- readNIfTI(paste(perf, "Fp","spline", sep="_"))
overlay(as(dynamic[xx,yy,,], "nifti"), 
        ifelse(fit.spline$Fp[xx,yy,] < 0.15, fit.spline$Fp[xx,yy,], NA),
        z=7, w=11, zlim.x=c(32,512), zlim.y=c(0,0.15), plot.type="single")
dev.off()


###################################################
### chunk number 33: AuditTrail
###################################################
print(audit.trail(fit$ktrans))


###################################################
### chunk number 34: end eval=FALSE
###################################################
## options(prompt="R> ")


