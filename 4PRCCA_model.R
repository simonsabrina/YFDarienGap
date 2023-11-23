# Model for PRCCA

darien <- function(a,b1,b2,c,muv,mum,muh,epsilonv,epsilonm,epsilonh,
                              gamam,gamah,alfam,alfah,chi,deltah){
  
  
  # Vectors' parameters
#  a	<- 6/30  # Average biting rate	(Haemagogus)
#  b1 <- 0.25 # Transmission coef. Haem -> Humans (Fraction of bites actually infective)
#  b2 <- 0.4 # Transmission coef. Haem -> Monkeys
#  c <- 0.4 # Transm. coef. Humans/Monkeys -> Haem. (Vector susceptibility to the virus)
#  muv <- 0.46/30 # Natural mortality rate of Haemagogu
#  epsilonv <- 4/30 #	Latency rate in vectors	(4.0 month-1 ~7.5 days)
  
  # Monkeys' parameters
#  mum	<- 0.0048/30 # Alouatta natural mortality rate (Esteva, 2019)
#  gamam <- 1/10 # Alouatta recovery rate (3 month-1, Esteva, 2019)
#  epsilonm <- 1/6 # YF latency rate in monkeys (= humans)
#  alfam <- 0.25/30  # Average Disease-induced mortality rate among monkey species
  
  # Humans' parameters
#  muh <- 0.00119/30 #(3.966667e-05 day-1) NATURAL MORTALITY
#  gamah <- 4/30 # Human recovery rate	(4.0 month-1 ~7 days)
#  epsilonh <- 1/6 # YF latency rate for humans (~6 days)
#  alfah <- 0.0244/30 #	Disease-induced mortality rate (0.0244 month-1)
  
  # Time of exposure:
#  pace <- 10 #days
#  chi <- 1/pace
#  Sh0 <- 1000
#  deltah <- 0.2 #have to be bigger than chi
  #(180000/365)/Sh0 # new people/day = 500 = 0.2 of the Sh0 inside the forest
  
  par.darien <- c(a=a, b1=b1, b2=b2, c=c, c=c, muv=muv, epsilonv=epsilonv,
                  epsilonh=epsilonh, epsilonm=epsilonm, gamam=gamam, alfam=alfam,
                  deltah=deltah, muh=muh, gamah=gamah, alfah=alfah, mum=mum, chi=chi)
  
  # 2. INITIAL CONDITIONS
  
  # Population - Darien Monkeys
  Sm0 <- 1000; Sv0 <- Sm0*1.5
  m2 <- Sv0/Sm0
  R0_monkeys <- a*a*b2*c*m2*epsilonv/(muv*(muv+epsilonv)*(mum+alfam+gamam))
  R0_monkeys
  
  # Population - Migrant Humans
  Sh0 <- 1000; 
  m1 <- Sv0/Sh0
  R0_humans <- a*a*b1*c*m1*epsilonm/(muv*(muv+epsilonv)*(muh+alfah+gamah))
  R0_humans
  
  
  state.darien <- c(Sh=Sh0,Eh=0, Ih=0, Rh=0, Infh=0, Death=0,
                    Sm=Sm0, Em=0, Im=0, Rm=0, Infm=0, Deathm=0,
                    Sv=Sv0, Ev=0, Iv=1)
  
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
      dInfh <- a*b1*Iv*(Sh/(Sh+Eh+Ih+Rh))
      dDeath <- alfah*Ih
      
      # Monkeys hosts   
      dSm <- -a*b2*Iv*(Sm/(Sm+Em+Im+Rm)) + mum*(Em+Im+Rm) + alfam*Im
      dEm <- a*b2*Iv*(Sm/(Sm+Em+Im+Rm)) - epsilonm*Em - muh*Em
      dIm <- epsilonm*Em - gamam*Im - mum*Im - alfam*Im
      dRm <- gamam*Im - mum*Rm
      dInfm <- a*b2*Iv*(Sm/(Sm+Em+Im+Rm))
      dDeathm <- alfam*Im
      
      # Mosquitoes vectors
      dSv <- -(a*c*Sv*(Ih/Sh+Eh+Ih+Rh)+a*c*Sv*(Im/Sm+Em+Im+Rm)) + muv*(Ev + Iv)
      dEv <- (a*c*Sv*(Ih/Sh+Eh+Ih+Rh)+a*c*Sv*(Im/Sm+Em+Im+Rm)) - (epsilonv + muv)*Ev
      dIv <- epsilonv*Ev - muv*Iv
      
      # return the output of the model
      return(list(c(dSh, dEh, dIh, dRh, dInfh, dDeath,
                    dSm, dEm, dIm, dRm, dInfm, dDeathm,
                    dSv, dEv, dIv)))
      
    })
  }
  
  tempo <- seq(from=0,to=tsim,by=Dt)
  
  mod <- ode(y = state.darien, times = tempo, func = model.darien, 
             parms = par.darien, method = "lsode")
  
  mod <- as.data.frame(mod)
  
  names(mod) <- c("t","Sh","Eh","Ih","Rh","Infh", "Death",
                  "Sm","Em","Im","Rm","Infm", "Deathm",
                  "Sv","Ev","Iv")
  
#  return(data.frame(Infh = mod$Infh[length(mod$Infh)])) #switch between both
  return(data.frame(Death = mod$Death[length(mod$Death)]))
  
}
