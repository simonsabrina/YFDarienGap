library(epiR)
library(deSolve)

setwd("C:/Users/simon/Documents/Doutorado/1_Tese/4darien")
source("4PRCCA_model.R")

# Assuming uniform distribution
set.seed(1234)

a <- runif(n=100, min=2/30, max=10/30)
b1 <- runif(n=100, min=0.2, max=0.6)
b2 <- runif(n=100, min=0.2, max=0.6)
c <- runif(n=100, min=0.2, max=0.6)
mum <- runif(n=100, min=0.001/30, max=0.01/30)
muv <- runif(n=100, min=0.4/30, max=1/30)
muh <- runif(n=100, min=3.04414e-05, max=5.479452e-05) #1/y/365: 50y to 90y
epsilonm <- runif(n=100, min=1/10, max=1/4)
epsilonv <- runif(n=100, min=1/10, max=1/4)
epsilonh <- runif(n=100, min=1/10, max=1/4)
gamam <- runif(n=100, min=1/12, max=1/5)
gamah <- runif(n=100, min=1/12, max=1/5)
alfam <- runif(n=100, min=0.15/30, max=0.5/30) #0.25/30
alfah <- runif(n=100, min=0.01/30, max=0.05/30) #0.0244/30
chi <- runif(n=100, min=1/10, max=1/2)
deltah <- runif(n=100, min=0.5, max=1)

#### Defining a vector for Infh #----
Infh <- numeric(length(deltah))

for (i in 1:length(deltah)){
  mod <- darien(deltah[i], a[i],b1[i],b2[i],c[i],muv[i],mum[i],muh[i],
                epsilonv[i],epsilonm[i],epsilonh[i],gamam[i],gamah[i],
                alfam[i],alfah[i],chi[i])
                           
  Infh[i] <- mod$Infh
}

# Sensitivity Analysis PRCC for Infh ----
dat <- data.frame(cbind(X1 = deltah, X2 = a, X3 = b1, X4 = b2, X5 = c,
                        X6 = muv, X7 = mum, X8 = muh , X9 = epsilonv,
                        X10 = epsilonm, X11 = epsilonh, X12 = gamam,
                        X13 = gamah, X14 = alfam, X15 = alfah, X16 = chi,
                        Y = Infh))
colnames(dat) <- c("migration rate", "biting rate", 
                   "Transmission Haem -> Humans", "Transmission Haem -> Monkeys",
                   "mosquito susceptibility", "mortality mosquitoes", 
                   "mortality monkeys", "mortality humans", "latency in mosquitoes",
                   "latency in monkeys", "latency in humans", "recovery monkeys",
                   "recovery humans", "YF death monkeys", "YF death humans",
                   "pace", "cases")
est <- epi.prcc(dat, sided.test = 2, conf.level = 0.95); est

# PLOTTING

par(mfrow=c(2,1))
bp <- barplot(est$est, ylim=c(-1,1), main="A) Correlation coefficient between inputs and number of cases",
            # cex.main = 0.95, cex.axis = 0.9, cex.names = 0.9, cex.lab = 0.9,
              names.arg=c(expression(delta[h]),"a", "b1", "b2", "c", expression(mu[v]), 
                          expression(mu[m]),expression(mu[h]),expression(epsilon[v]),
                          expression(epsilon[m]),expression(epsilon[h]),
                          expression(gamma[m]),expression(gamma[h]),
                          expression(alpha[m]),expression(alpha[h]),
                          expression(chi[h])),

              ylab='Coefficient', axes=FALSE)
axis(2)
mtext(text='Parameters', side=1, line=2.5)#, cex=0.9
box(); abline(h=0)
text(bp, -0.9, paste("r = ", round(est$est,2), sep=""), cex=0.9,pos=3) #cex=0.95,
grid(nx = 17, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.8)      # Grid line width


### # Defining a vector for Deaths #############################################
source("4PRCCA_model.R")

Death <- numeric(length(deltah))

for (i in 1:length(deltah)){
  mod <- darien(deltah[i], a[i],b1[i],b2[i],c[i],muv[i],mum[i],muh[i],
                epsilonv[i],epsilonm[i],epsilonh[i],gamam[i],gamah[i],
                alfam[i],alfah[i],chi[i])
  
  Death[i] <- mod$Death
}

# Sensitivity Analysis PRCC for Deaths ----
dat.death <- data.frame(cbind(X1 = deltah, X2 = a, X3 = b1, X4 = b2, X5 = c,
                        X6 = muv, X7 = mum, X8 = muh , X9 = epsilonv,
                        X10 = epsilonm, X11 = epsilonh, X12 = gamam,
                        X13 = gamah, X14 = alfam, X15 = alfah, X16 = chi,
                        Y = Death))
colnames(dat.death) <- c("migration rate", "biting rate", 
                   "Transmission Haem -> Humans", "Transmission Haem -> Monkeys",
                   "mosquito susceptibility", "mortality mosquitoes", 
                   "mortality monkeys", "mortality humans", "latency in mosquitoes",
                   "latency in monkeys", "latency in humans", "recovery monkeys",
                   "recovery humans", "YF death monkeys", "YF death humans",
                   "pace", "deaths")
est.death <- epi.prcc(dat.death, sided.test = 2, conf.level = 0.95); est.death

# PLOTTING

#par(mfrow=c(1,1))
bp2 <- barplot(est.death$est, ylim=c(-1,1), main="B) Correlation coefficient between inputs and number of deaths",
               #cex.main = 0.95, cex.axis = 0.9, cex.names = 0.9, cex.lab = 0.9,
              names.arg=c(expression(delta[h]),"a", "b1", "b2", "c", expression(mu[v]), 
                          expression(mu[m]),expression(mu[h]),expression(epsilon[v]),
                          expression(epsilon[m]),expression(epsilon[h]),
                          expression(gamma[m]),expression(gamma[h]),
                          expression(alpha[m]),expression(alpha[h]),
                          expression(chi[h])),
              
              ylab='Coefficient', axes=FALSE)
axis(2)
mtext(text='Parameters', side=1, line=2.5)#, cex = 0.9
box(); abline(h=0)
text(bp2, -0.9, paste("r = ", round(est.death$est,2), sep=""),cex=0.9,pos=3) #cex=0.95, 
grid(nx = 17, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.8)      # Grid line width


#### UNCERTANTY ANALYSIS OF b (MONTE CARLO SIMULATION) ----
source("PRCCA_model.R")

b <- runif(n=100, min=0.2, max=0.5)
Infh2 <- numeric(length(b))

for (i in 1:length(b)){
  Infh2[i] <- inf.bites(b[i])$Infh2
}

par(mfrow=c(1,2))
plot(b, Infh2,
        xlab= "parameter b",
        ylab= "Cases in host community")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.5)      # Grid line width


# Probability interval (95% PI) 
quantile(Infh2, probs=c(0.025,0.975))


#### UNCERTANTY ANALYSIS OF MIGRATION RATE deltah_12 (MONTE CARLO SIMULATION) ----
source("PRCCA_model.R")

deltah_12 <- runif(n=100, min=0.1/365, max=0.9/365)
Infh2 <- numeric(length(deltah_12))

for (i in 1:length(deltah_12)){
  Infh2[i] <- delta(deltah_12[i])$Infh2
}

plot(deltah_12, Infh2,
     xlab= "migration rate",
     ylab= "Cases in host community")
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.5)      # Grid line width

# Probability interval (95% PI) 
quantile(Infh2, probs=c(0.025,0.975))