# DARIEN GAP: SEIR-SEI Multihost Spatiotemporal

# Loading packages
library(deSolve)

# 1. PARAMETERS in daily rates From Esteva et al, 2019: https://onlinelibrary.wiley.com/doi/full/10.1002/cmm4.1059

# Vectors' parameters
a	<- 6/30  # Average biting rate	(Haemagogus)
b1 <- 0.25 # Transmission coef. Haem -> Humans (Fraction of bites actually infective)
b2 <- 0.4 # Transmission coef. Haem -> Monkeys
c1 <- 0.4 # Transm. coef. Humans -> Haem. (Vector susceptibility to the virus)
c2 <- 0.4 # Transm. coef. Monkeys -> Haem. 
muv <- 0.46/30 # Natural mortality rate of Haemagogu
epsilonv <- 4/30 #	Latency rate in vectors	(4.0 month-1 ~7.5 days)

# Monkeys' parameters
mum	<- 0.0048/30 # Alouatta natural mortality rate (Esteva, 2019)
gamam <- 1/10 # Alouatta recovery rate (3 month-1, Esteva, 2019)
epsilonm <- 1/6 # YF latency rate in monkeys (= humans)
alfam <- 0.25/30  # Average Disease-induced mortality rate among monkey species

# Humans' parameters
muh <- 0.00119/30 #(3.966667e-05 day-1) NATURAL MORTALITY
gamah <- 4/30 # Human recovery rate	(4.0 month-1 ~7 days)
epsilonh <- 1/6 # YF latency rate for humans (~6 days)
alfah <- 0.0244/30 #	Disease-induced mortality rate (0.0244 month-1)

# Time of exposure:
Sh0 <- 1000
pace <- 2 #days
chi <- (1/pace)#/Sh0 # have to transform to N humans
deltah <- (180000/365)/Sh0 # new people/day = 500 = 0.2 of the Sh0 inside the forest


par.darien <- c(a=a, b1=b1, b2=b2, c1=c1, c2=c2, muv=muv, epsilonv=epsilonv,
             epsilonh=epsilonh, epsilonm=epsilonm, gamam=gamam, alfam=alfam,
             deltah=deltah, muh=muh, gamah=gamah, alfah=alfah, mum=mum, chi=chi)

# 2. INITIAL CONDITIONS

# Population - Darien Monkeys
Sm0 <- 1000; Sv0 <- Sm0*1.5
m2 <- Sv0/Sm0
R0_monkeys <- a*a*b2*c2*m2*epsilonv/(muv*(muv+epsilonv)*(mum+alfam+gamam))
R0_monkeys

# Population - Migrant Humans
Sh0 <- 1000; # There are 2,500 people inside the forest every day
m1 <- Sv0/Sh0
R0_humans <- a*a*b1*c1*m1*epsilonm/(muv*(muv+epsilonv)*(muh+alfah+gamah))
R0_humans


state.darien <- c(Sh=Sh0,Eh=0,Ih=1,Rh=0,Infh=0,Risk=0,Death=0,Expos=0,dTotal=0,
               Sm=Sm0, Em=0, Im=0, Rm=0, Infm=0, Deathm=0,
               Sv=Sv0, Ev=0, Iv=0)

# Time of simulation
tsim <- 1*365
Dt <- 1

# 3. DETERMINISTIC ROSS-MACDONALD MODEL #

# Frequency-dependent transmission
model.darien <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # Human hosts
    dSh <- -a*b1*Iv*(Sh/(Sh+Eh+Ih+Rh)) - muh*Sh + deltah*Sh - chi*Sh
    dEh <- a*b1*Iv*(Sh/(Sh+Eh+Ih+Rh)) - epsilonh*Eh - muh*Eh - chi*Eh
    dIh <- epsilonh*Eh - gamah*Ih - muh*Ih - alfah*Ih - chi*Ih
    dRh <- gamah*Ih - muh*Rh - chi*Rh
    dInfh <- a*b1*Iv*(Sh/(Sh+Eh+Ih+Rh)) # total cases inside and outside forest
    dRisk <-  a*b1*Iv*(Sh/(Sh+Eh+Ih+Rh))/chi # = beta*pace = casos x 10
    dDeath <- alfah*Ih
    dExpos <- chi*(Eh + Ih)#Eh + Ih leaving the forest
    dTotal <- deltah*Sh #have to equal 180,000 in 1 year
    
    # Monkeys hosts   
    dSm <- -a*b2*Iv*(Sm/(Sm+Em+Im+Rm)) + mum*(Em+Im+Rm) + alfam*Im
    dEm <- a*b2*Iv*(Sm/(Sm+Em+Im+Rm)) - epsilonm*Em - muh*Em
    dIm <- epsilonm*Em - gamam*Im - mum*Im - alfam*Im
    dRm <- gamam*Im - mum*Rm
    dInfm <- a*b2*Iv*(Sm/(Sm+Em+Im+Rm))
    dDeathm <- alfam*Im
    
    # Mosquitoes vectors
    dSv <- -(a*c1*Sv*(Ih/Sh+Eh+Ih+Rh)+a*c2*Sv*(Im/Sm+Em+Im+Rm)) + muv*(Ev + Iv)
    dEv <- (a*c1*Sv*(Ih/Sh+Eh+Ih+Rh)+a*c2*Sv*(Im/Sm+Em+Im+Rm)) - (epsilonv + muv)*Ev
    dIv <- epsilonv*Ev - muv*Iv
    
    # return the output of the model
    return(list(c(dSh, dEh, dIh, dRh, dInfh, dRisk, dDeath, dExpos, dTotal,
                  dSm, dEm, dIm, dRm, dInfm, dDeathm,
                  dSv, dEv, dIv)))
    
  })
}

tempo <- seq(from=0,to=tsim,by=Dt)

mod <- ode(y = state.darien, times = tempo, func = model.darien, 
           parms = par.darien, method = "lsode")

mod <- as.data.frame(mod)

names(mod) <- c("t","Sh","Eh","Ih","Rh","Infh","Risk", "Death", "Expos", "Total",
                "Sm","Em","Im","Rm","Infm","Deathm",
                "Sv","Ev","Iv")

# 4. PLOTTING

paste("Pace=",pace,", Cases=",round(max(mod$Infh)),
      "Deaths=",round(max(mod$Death)),
      "Number of exposed leaving=",round(max(mod$Expos)))

par(mfrow=c(1,2))

# Human cases
matplot(mod$Ih, type = "l", lwd=2, col = "gray", # 10 days
        main=paste("(A) Infected humans"),cex.main=0.8,
        #sub=paste("Cases=",round(max(mod$Infh))),cex.sub=0.8,
        xlab= "Time (days)",cex.lab = 0.8,
        ylab= "N. cases", cex.axis = 0.8)
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.1)      # Grid line width

lines(mod$Ih, type = "l", lwd=2, lty=2, col = "black") # 7 days
lines(mod$Ih, type = "l", lwd=2, lty=3, col = "black") # 5 days

# Human deaths

matplot(mod$Death, type = "l", lwd=2, col = "gray", #lwd=1.5, #infected humans
        main=paste("(B) Human deaths"),cex.main=0.8,
        #sub=paste("Deaths=",round(max(mod$Death))),cex.sub=0.8,
        xlab= "Time (days)",cex.lab = 0.8,
        ylab= "N. deaths", cex.axis = 0.8)
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.1)      # Grid line width

lines(mod$Death, type = "l", lwd=2, lty=2, col = "black") # 7 days
lines(mod$Death, type = "l", lwd=2, lty=3, col = "black") # 5 days


# Monkeys 
matplot(mod$Im, type = "l", lwd=2, col = "gray", #lwd=1.5, #infected monkeys
        main="(B) Infected Monkeys",cex.main=0.8,
        sub=paste("Cases=",round(max(mod$Infm)),
                  "Deaths=",round(max(mod$Deathm))),cex.sub=0.8,
        xlab= "Time (days)",cex.lab = 0.8,
        ylab= "N. cases", cex.axis = 0.8)
grid(nx = NULL, ny = NULL,
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.1)      # Grid line width

paste("Cases=",round(max(mod$Infm)),
      "Deaths=",round(max(mod$Deathm)))
