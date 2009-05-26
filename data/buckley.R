buckley <- read.table("buckley.txt", header=TRUE)
breast <- list(data = buckley[,4:16])
breast$Fp <- c(.17, .37, .57, .77, .97, rep(.57, 8))
breast$PS <- c(rep(.33, 5), .01, .17, .49, .65, rep(.33, 4))
breast$Vp <- c(rep(.06, 9), .0001, .03, .09, .12)
breast$Ve <- rep(.45, 13)
breast$Ktrans <- (1 - exp(-breast$PS/breast$Fp)) * breast$Fp
breast$Kep <- breast$Ktrans / breast$Ve
meningioma <- list(data = buckley[,17:29])
meningioma$Fp <- c(.4, .8, 1.2, 1.6, 2.0, rep(1.2, 8))
meningioma$PS <- c(rep(.34, 5), 0, .17, .51, .68, rep(.34, 4))
meningioma$Vp <- c(rep(.08, 9), .0001, .04, .12, .16)
meningioma$Ve <- rep(.4, 13)
meningioma$Ktrans <- (1 - exp(-meningioma$PS/meningioma$Fp)) * meningioma$Fp
meningioma$Kep <- meningioma$Ktrans / meningioma$Ve
