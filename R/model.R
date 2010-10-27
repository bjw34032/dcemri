
aifParameters <- function(type, user=NULL) {
  switch(type,
         tofts.kermode = list(D=0.1, a1=3.99, a2=4.78, m1=0.144, m2=0.0111),
         fritz.hansen = list(D=1.0, a1=2.4, a2=0.62, m1=3.0, m2=0.016),
         orton.exp = list(D=1.0, AB=323, muB=20.2, AG=1.07, muG=0.172),
         orton.cos = list(D=1.0, AB=2.84, muB=22.8, AG=1.36, muG=0.171),
         user = user,
         empirical = user,
         stop("AIF parameters must be specified!"))
}

compartmentalModel <- function(type) {
  switch(type,
         weinmann =
         function(time, th1, th3, p) {
           ## Convolution of Tofts & Kermode AIF with single-compartment model
           erg <- p$D * exp(th1) * ((p$a1 / (p$m1 - exp(th3))) * (exp(-(time * exp(th3))) - exp(-(time * p$m1))) + (p$a2 / (p$m2 - exp(th3))) * (exp(-(time * exp(th3))) - exp(-(time * p$m2))))
           erg[time <= 0] <- 0
           return(erg)
         },
         extended =
         function(time, th0, th1, th3, p) {
           ## Extended Tofts & Kermode model including the concentration of
           ## contrast agent in the blood plasma (vp)
           Cp <- function(tt, D, a1, a2, m1, m2) {
             D * (a1 * exp(-m1 * tt) + a2 * exp(-m2 * tt))
           }
           erg <- exp(th0) * Cp(time, p$D, p$a1, p$a2, p$m1, p$m2) +
             p$D * exp(th1) * ((p$a1 / (p$m1 - exp(th3))) * (exp(-(time * exp(th3))) - exp(-(time * p$m1))) + (p$a2 / (p$m2 - exp(th3))) * (exp(-(time * exp(th3))) - exp(-(time * p$m2))))
           erg[time <= 0] <- 0
           return(erg)
         },
         kety.orton.exp =
         function(time, th1, th3, p) {
           ## Kety model using the exponential AIF from Matthew Orton (ICR)
           ktrans <- exp(th1)
           kep <- exp(th3)
           T1 <- p$AB * kep / (kep - p$muB)
           T2 <- time * exp(-p$muB * time) -
             (exp(-p$muB * time) - exp(-kep * time)) / (kep - p$muB)
           T3 <- p$AG * kep
           T4 <- (exp(-p$muG * time) - exp(-kep * time)) / (kep - p$muG) -
             (exp(-p$muB * time) - exp(-kep * time)) / (kep - p$muB)
           erg <- ktrans * (T1 * T2 + T3 * T4)
           erg[time <= 0] <- 0
           return(erg)
         },
         orton.exp =
         function(time, th0, th1, th3, p) {
           ## Extended model using the exponential AIF from Matthew Orton (ICR)
           Cp <- function(tt, p) {
             p$AB * tt * exp(-p$muB * tt) + p$AG *
               (exp(-p$muG * tt) - exp(-p$muB * tt))
           }
           vp <- exp(th0)
           ktrans <- exp(th1)
           kep <- exp(th3)
           T1 <- p$AB * kep / (kep - p$muB)
           T2 <- time * exp(-p$muB * time) -
             (exp(-p$muB * time) - exp(-kep * time)) / (kep - p$muB)
           T3 <- p$AG * kep
           T4 <- (exp(-p$muG * time) - exp(-kep * time)) /
             (kep - p$muG) - (exp(-p$muB * time) - exp(-kep * time)) /
               (kep - p$muB)
           erg <- vp * Cp(time, p) + ktrans * (T1 * T2 + T3 * T4)
           erg[time <= 0] <- 0
           return(erg)
         },
         kety.orton.cos =
         function(time, th1, th3, p) {
           ## Extended model with the raised cosine AIF from Matthew Orton (ICR)
           A2 <- function(time, alpha, p) {
             (1 - exp(-alpha * time)) / alpha -
               (alpha * cos(p$muB * time) + p$muB * sin(p$muB * time) -
                alpha * exp(-alpha * time)) / (alpha^2 + p$muB^2)
           }
           ktrans <- exp(th1)
           kep <- exp(th3)
           tB <- 2 * pi / p$muB
           cp <- ifelse(time <= tB,
                        p$aB * (1 - cos(p$muB*time)) + p$aB * p$aG * A2(time, p$muG, p),
                        p$aB * p$aG * A2(time, p$muG, p) * exp(-p$muB * (time - tB)))
           erg <- ifelse(time <= tB,
                         p$aB * p$aG * ktrans / (kep - p$muG) * ((A2(time, p$muG, p) + (kep - p$muG) / p$aG - 1) * A2(time, kep, p)),
                         p$aB * p$aG * ktrans / (kep - p$muG) * (A2(tB, p$muG, p) * exp(-p$muB * (time - tB)) + ((kep - p$muG) / p$aG - 1) * A2(tB, kep, p) * exp(-kep * (time - tB))))
           erg[time <= 0] <- 0
           return(erg)
         },
         orton.cos =
         function(time, th0, th1, th3, p) {
           ## Extended model with the raised cosine AIF from Matthew Orton (ICR)
           A2 <- function(time, alpha, p) {
             (1 - exp(-alpha * time)) / alpha -
               (alpha * cos(p$muB * time) + p$muB * sin(p$muB * time) -
                alpha * exp(-alpha * time)) / (alpha^2 + p$muB^2)
           }
           vp <- exp(th0)
           ktrans <- exp(th1)
           kep <- exp(th3)
           tB <- 2 * pi / p$muB
           cp <- ifelse(time <= tB,
                        p$aB * (1 - cos(p$muB * time)) + p$aB * p$aG * A2(time, p$muG, p),
                        p$aB * p$aG * A2(time, p$muG, p) * exp(-p$muB * (time - tB)))
           erg <- ifelse(time <= tB,
                         vp * cp + p$aB * p$aG * ktrans / (kep - p$muG) * ((A2(time, p$muG, p) + (kep - p$muG) / p$aG - 1) * A2(time, kep, p)),
                         vp * cp + p$aB * p$aG * ktrans / (kep - p$muG) * (A2(tB, p$muG, p) * exp(-p$muB * (time - tB)) + ((kep - p$muG) / p$aG - 1) * A2(tB, kep, p) * exp(-kep * (time - tB))))
           erg[time <= 0] <- 0
           return(erg)
         },
         weinmann.empirical =
         function(time, th1, th3, aif) {
           tsec <- seq(min(time*60), ceiling(max(time*60)), by=1)
           aif.new <- approx(time*60, aif, tsec)$y
           erg <- approx(tsec, expConv(aif.new, exp(th1)/60, exp(th3)/60), time*60)$y
           erg[time <= 0] <- 0
           return(erg)
         },
         extended.empirical =
         function(time, th0, th1, th3, aif) {
           tsec <- seq(min(time*60), ceiling(max(time*60)), by=1)
           aif.new <- approx(time*60, aif, tsec)$y
           erg <- approx(tsec, exp(th0) * aif.new + expConv(aif.new, exp(th1)/60, exp(th3)/60), time*60)$y
           erg[time <= 0] <- 0
           return(erg)
         })
}

expConv <- function(input, k1, k2) {
  n <- length(input)
  convolution <- rep(0, n)
  ek2 <- exp(-k2)
  if (ek2 == 1) {
    convolution <- k1 * cumsum(input)
  } else {
    prev <- 0
    for (i in 1:n) {
      prev <- prev * ek2 + k1 * input[i] * (1 - ek2) / k2
      convolution[i] <- prev
    }
  }
  return(convolution)
}
  
