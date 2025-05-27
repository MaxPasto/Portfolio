# Data setting ------------------------------------------------------------

a.min <- 0 # età minima dell'intervallo di fit
a.max <- 100 # età massima dell'intervallo di fit
omega <-121 # età estrema delle tavole
A.fit <-c(a.min:a.max) 
y.fit.min <- 1969 # anno minimo dell'intervallo di fit
y.fit.max <- 2019 # anno massimo dell'intervallo di fit
Y.fit <- c(y.fit.min:y.fit.max)
y.pred <- 50 # orizzonte di proiezione
Y.pred <- c((y.fit.max+1):(y.fit.max+y.pred))
n.sim <- 1000 # numero di simulazioni
n.boot <- 20 # numero di simulazioni per il bootstrap
int_vect <- c(rep(0.02,omega-a.min)) # vettore di tassi d'interesse per scadenza
v_vect <- cumprod((1+int_vect)^-1) # vettore di fattori di sconto per scadenza
#################################### fine data setting


# Importazione dati demografici HMD ---------------------------------------

library(demography)

Francedata <- read.demogdata(file = "France_Mx_1x1.txt",
                          popfile = "Exposures_1x1.txt",
                          type = "mortality",
                          label = "France")
names(Francedata)


# Il fit dei modelli ------------------------------------------------------

library("StMoMo")
source("Functions.r") # caricamento di funzioni da un altro script

# Trasfomazione dei dati e definizione dell'intervallo di stima

FrancemStMoMoC <- StMoMoData(Francedata, series = "male")     # esposizioni centrali
FrancefStMoMoC <- StMoMoData(Francedata, series = "female")   
FRANCEmStMoMoI <- central2initial(FrancemStMoMoC)             # esposizioni iniziali
FRANCEfStMoMoI <- central2initial(FrancefStMoMoC )

wxt <- genWeightMat(ages = A.fit, years = Y.fit, clip = 4) # pesi
FRANCEmRates <- FrancemStMoMoC$Dxt/FrancemStMoMoC$Ext
FRANCEmRates <- FRANCEmRates[A.fit+1,tail(FrancemStMoMoC$years+1,
                                    length(Y.fit))-FrancemStMoMoC$years[1]]
FRANCEfRates <- FrancefStMoMoC$Dxt/FrancefStMoMoC$Ext
FRANCEfRates <- FRANCEfRates[A.fit+1,tail(FrancefStMoMoC$years+1,
                                    length(Y.fit))-FrancefStMoMoC$years[1]]

################### Models fitting

LCfit_FRANCEm <- fit(lc(link = "logit"),
                  data = FRANCEmStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt)


LCfit_FRANCEf <- fit(lc(link = "logit"),
                  data = FRANCEfStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt)




RHfit_FRANCEm <- fit(rh(link = "logit", cohortAgeFun="1"),
                  data = FRANCEmStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt,
                  start.ax = LCfit_FRANCEm$ax,
                  start.bx = LCfit_FRANCEm$bx,
                  start.kt = LCfit_FRANCEm$kt)

RHfit_FRANCEf <- fit(rh(link = "logit", cohortAgeFun="1"),
                  data = FRANCEfStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt,
                  start.ax = LCfit_FRANCEf$ax,
                  start.bx = LCfit_FRANCEf$bx,
                  start.kt = LCfit_FRANCEf$kt)




APCfit_FRANCEm <- fit(apc(link = "logit"),
                   data = FRANCEmStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit=y.fit.min:y.fit.max,
                   wxt = wxt)

APCfit_FRANCEf <- fit(apc(link = "logit"),
                   data = FRANCEfStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit=y.fit.min:y.fit.max,
                   wxt = wxt)



M7fitFRANCEm <- fit(m7(link = "logit"),
                  data = FRANCEmStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt)

M7fitFRANCEf <- fit(m7(link = "logit"),
                  data = FRANCEfStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt)



PLATfit_FRANCEm <- fit(PLAT,
                    data = FRANCEmStMoMoI,
                    ages.fit = a.min:a.max,
                    years.fit=y.fit.min:y.fit.max,
                    wxt = wxt)

PLATfit_FRANCEf <- fit(PLAT,
                    data = FRANCEfStMoMoI,
                    ages.fit = a.min:a.max,
                    years.fit=y.fit.min:y.fit.max,
                    wxt = wxt)

# Analisi e confronto tra modelli -----------------------------------------

######### Confronto tra modelli

### Log-like

## Popolazione maschile

max((c(LCfit_FRANCEm$loglik, ## -37602.23
       RHfit_FRANCEm$loglik, ##  -28701.93
       APCfit_FRANCEm$loglik, ## -38436.26
       M7fitFRANCEm$loglik,  ## -379891.7
       PLATfit_FRANCEm$loglik))) ## -35276.75

## Popolazione femminile

max(c(LCfit_FRANCEf$loglik, ## -30753.76
      APCfit_FRANCEf$loglik, ## -33640.73
      M7fitFRANCEf$loglik, ##  -270848.5
      PLATfit_FRANCEf$loglik)) ## -31030.02

###### BIC

## Popolazione maschile 

min(c(BIC(LCfit_FRANCEm), ## 77348.76
      BIC(RHfit_FRANCEm), ##  60761.28
      BIC(APCfit_FRANCEm), ## 79367.08
      BIC(M7fitFRANCEm), ## 762286.5
      BIC(PLATfit_FRANCEm))) ## 73466.68

## Popolazione femminile

min(c(BIC(LCfit_FRANCEf), ## 63651.83
      BIC(APCfit_FRANCEf), ## 69776.04
      BIC(M7fitFRANCEf),  ## 544200.2
      BIC(PLATfit_FRANCEf))) ## 64973.22


###### AIC

## Popolazione maschile 

min(c(AIC(LCfit_FRANCEm), ##  75706.45
      AIC(RHfit_FRANCEm), ##  58189.86
      AIC(APCfit_FRANCEm), ##  77456.51
      AIC(M7fitFRANCEm), ##   760369.4
      AIC(PLATfit_FRANCEm))) ## 71235.5

## Popolazione femminile

min(c(AIC(LCfit_FRANCEf), ## 62009.52
      AIC(APCfit_FRANCEf), ## 67865.46
      AIC(M7fitFRANCEf),  ## 542283.1
      AIC(PLATfit_FRANCEf))) ## 62742.04



################### Rappresentazione dei parametri

plot(LCfit_FRANCEm, nCol = 3) ## Lee-Carter per la popolazione maschile

plot(LCfit_FRANCEf, nCol = 3) ## Lee-Carter per la popolazione femminile

plot(RHfit_FRANCEm, nCol = 3) ## Renshaw-Haberman per la popolazione maschile

plot(PLATfit_FRANCEf, parametricbx = F) ## Plat per la popolazione femminile


###### Residui

LCres_FRANCEm <- residuals(LCfit_FRANCEm) ## Residui Lee-Carter per la popolazione maschile

LCres_FRANCEf <- residuals(LCfit_FRANCEf) ## Residui Lee-Carter per la popolazione femminile

RHres_FRANCEm <- residuals(RHfit_FRANCEm) ## Residui Renshaw-Haberman per la popolazione maschile

PLATres_FRANCEf <- residuals(PLATfit_FRANCEf) ## Residui Plat per la popolazione femminile


### rappresentazione scatter dei residui

TYPE = "scatter"
par(mfrow=c(1, 1))

plot(LCres_FRANCEm, type=TYPE) ## scatter dei residui Lee-Carter per la popolazione maschile

plot(LCres_FRANCEf, type=TYPE) ## scatter dei residui Lee-Carter per la popolazione femminile


plot(RHres_FRANCEm, type=TYPE) ## scatter dei residui Renshaw-Haberman per la popolazione maschile

plot(PLATres_FRANCEf, type=TYPE) ## scatter dei residui Plat per la popolazione femminile

### colourmap dei residui


plot(LCres_FRANCEm, type = "colourmap", reslim = c(-3.5, 3.5)) ## colourmap dei residui Lee-Carter

## per la popolazione maschile


plot(LCres_FRANCEf, type = "colourmap", reslim = c(-3.5, 3.5)) ## colourmap dei residui Lee-Carter 

## per la popolazione femminile


plot(RHres_FRANCEm, type = "colourmap", reslim = c(-3.5, 3.5)) ## colourmap dei residui Renshaw-Haberman

## per la popolazione maschile


plot(PLATres_FRANCEf, type = "colourmap", reslim = c(-3.5, 3.5)) ## colourmap dei residui Plat

## per la popolazione femminile


# Proiezione del parametro k_t considerando il modello Lee-Carter

# Orizzonte di previsione: 50 anni

# k_t model: Random Walk with drift

LCfor_FRANCEm <- forecast(LCfit_FRANCEm, h=y.pred) ## Popolazione maschile

plot(LCfor_FRANCEm, only.kt = TRUE) ## Rappresentazione

LCfor_FRANCEf <- forecast(LCfit_FRANCEf, h=y.pred) ## Popolazione femminile

plot(LCfor_FRANCEf, only.kt = TRUE) ## Rappresentazione



# Proiezioni dei parametri k_t e gc considerando il modello Renshaw-Haberman 

# Orizzonte di previsione: 50 anni

# k_t model: Random Walk with drift

# gc model: ARIMA(1,1,0) with drift

RHfor_FRANCEm <- forecast(RHfit_FRANCEm, h=y.pred) ## Popolazione maschile

plot(RHfor_FRANCEm, only.kt = TRUE) ## Rappresentazione del parametro k_t

plot(RHfor_FRANCEm, only.gc = TRUE) ## Rappresentazione del parametro gc


# Proiezioni dei parametri k_t(1), k_t(2) e gc considerando il modello Plat 

# Orizzonte di previsione: 50 anni

# k_t(1) model: Random Walk with drift

# k_t(2) model: Random Walk with drift

# gc model: ARIMA(1,1,0) with drift

PLATfor_FRANCEf <- forecast(PLATfit_FRANCEf, h=y.pred) ## Popolazione femminile

plot(PLATfor_FRANCEf, only.kt=TRUE) ## Rappresentazione dei parametri k_t(1) e k_t(2)

plot(PLATfor_FRANCEf, only.gc=TRUE) ## Rappresentazione del parametro gc

## Miglior modello ARIMA per k_t(1) considerando il modello Lee-Carter: ARIMA(1,1,0) with drift

LCfor_FRANCEmBArima <- forecast(LCfit_FRANCEm, h=y.pred, kt.method = "iarima",  kt.order = NULL) ## Popolazione maschile

LCfor_FRANCEfBArima <- forecast(LCfit_FRANCEf, h=y.pred, kt.method = "iarima",  kt.order = NULL) ## Popolazione femminile


## Miglior modello ARIMA per k_t(1) considerando il modello Renshaw-Haberman: ARIMA(1,2,1)


RHfor_FRANCEmBArima <- forecast(RHfit_FRANCEm, h=y.pred, kt.method = "iarima",  kt.order = NULL) ## Popolazione maschile

## Miglior modello ARIMA per k_t(1) e k_t(2) considerando il modello Plat

# k_t(1): ARIMA(2,1,2) with drift

# k_t(2): ARIMA(0,2,2)

PLATfor_FRANCEfBArima <- forecast(PLATfit_FRANCEf, h=y.pred, kt.method = "iarima",  kt.order = NULL)  # Popolazione femminile


#####  generazione di probabilità dai tassi (osservati e proiettati)

## Lee-Carter, popolazione maschile

rates_LCfit_FRANCEm <- fitted(LCfit_FRANCEm, type = "rates")

rates_LCfor_FRANCEm <- LCfor_FRANCEmBArima$rates


rates_LC_FRANCEm <- cbind(FRANCEmRates,rates_LCfor_FRANCEm)

q_LC_FRANCEm <- 1- exp(-rates_LC_FRANCEm)

#### estensione delle probabilità qxt ad età estreme

q_LC_FRANCEm.ext  <- extrapolation.fit(q_LC_FRANCEm)

## Lee-Carter, popolazione femminile

rates_LCfit_FRANCEf <- fitted(LCfit_FRANCEf, type = "rates")

rates_LCfor_FRANCEf <- LCfor_FRANCEfBArima$rates

rates_LC_FRANCEf <- cbind(FRANCEfRates,rates_LCfor_FRANCEf)

q_LC_FRANCEf <- 1- exp(-rates_LC_FRANCEf)

#### estensione delle probabilità qxt ad età estreme

q_LC_FRANCEf.ext  <- extrapolation.fit(q_LC_FRANCEf)

q_LC.FRANCEm<-as.data.frame(q_LC_FRANCEm.ext)

q_LC.FRANCEf<-as.data.frame(q_LC_FRANCEf.ext)


#####  generazione di probabilità dai tassi (osservati e proiettati)

## Renshaw-Haberman, popolazione maschile

rates_RHfit_FRANCEm <- fitted(RHfit_FRANCEm, type = "rates")

rates_RHfor_FRANCEm <- RHfor_FRANCEmBArima$rates

rates_RH_FRANCEm <- cbind(FRANCEmRates,rates_RHfor_FRANCEm)

q_RH_FRANCEm <- 1- exp(-rates_RH_FRANCEm)

#### estensione delle probabilità qxt ad età estreme

q_RH_FRANCEm.ext  <- extrapolation.fit(q_RH_FRANCEm)

q_RH.FRANCEm<-as.data.frame(q_RH_FRANCEm.ext)

#####  generazione di probabilità dai tassi (osservati e proiettati)

## Plat, popolazione femminile

rates_PLATfit_FRANCEf <- fitted(PLATfit_FRANCEf, type = "rates")

rates_PLATfor_FRANCEf <- PLATfor_FRANCEfBArima$rates

rates_PLAT_FRANCEf <- cbind(FRANCEfRates,rates_PLATfor_FRANCEf)

q_PLAT_FRANCEf <- 1- exp(-rates_PLAT_FRANCEf)

#### estensione delle probabilità qxt ad età estreme

q_PLAT_FRANCEf.ext  <- extrapolation.fit(q_PLAT_FRANCEf)

q_PLAT.FRANCEf<-as.data.frame(q_PLAT_FRANCEf.ext)



# Simulazioni (considerando i migliori modelli ARIMA visti in precedenza)

## Lee-Carter, popolazione maschile
 
LCsim_FRANCEm.Aarima <- simulate(LCfit_FRANCEm, nsim = n.sim, h=y.pred, kt.method = "iarima", kt.order = NULL)

rates_LC_FRANCEm.st <- LCsim_FRANCEm.Aarima$rates

q_LC_FRANCEm.st <- 1- exp(-rates_LC_FRANCEm.st)

q_LC_FRANCEm.st.ext <-  extrapolation.sim(q_LC_FRANCEm.st)


## Lee-Carter, popolazione femminile

LCsim_FRANCEf.Aarima <- simulate(LCfit_FRANCEf, nsim = n.sim, h=y.pred, kt.method = "iarima", kt.order = NULL)


rates_LC_FRANCEf.st <- LCsim_FRANCEf.Aarima$rates

q_LC_FRANCEf.st <- 1- exp(-rates_LC_FRANCEf.st)

q_LC_FRANCEf.st.ext <-  extrapolation.sim(q_LC_FRANCEf.st)


## Renshaw-Haberman, popolazione maschile

RHsim_FRANCEm.Aarima <- simulate(RHfit_FRANCEm, nsim = n.sim, h=y.pred, kt.method = "iarima", kt.order = NULL)

rates_RH_FRANCEm.st <- RHsim_FRANCEm.Aarima$rates

q_RH_FRANCEm.st <- 1- exp(-rates_RH_FRANCEm.st)

q_RH_FRANCEm.st.ext <-  extrapolation.sim(q_RH_FRANCEm.st)


## Plat, popolazione femminile

PLATsim_FRANCEf.Aarima <- simulate(PLATfit_FRANCEf, nsim = n.sim, h=y.pred, kt.method = "iarima", kt.order = NULL)

rates_PLAT_FRANCEf.st <- PLATsim_FRANCEf.Aarima$rates

q_PLAT_FRANCEf.st <- 1- exp(-rates_PLAT_FRANCEf.st)

q_PLAT_FRANCEf.st.ext <-  extrapolation.sim(q_PLAT_FRANCEf.st)


## Fan-charts considerando i modelli Lee-Carter per la popolazione maschile e femminile

library("fanplot")
dev.new(width = 550, height = 330, unit = "px")
par(mfrow=c(1, 1))
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(LCfit_FRANCEm$years, t(q_LC_FRANCEm[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.max+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "male mortality rate (log scale)")
fan(t(LCsim_FRANCEm.Aarima$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(LCsim_FRANCEm.Aarima$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(LCsim_FRANCEm.Aarima$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_LC_FRANCEm[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
matplot(LCfit_FRANCEf$years, t(q_LC_FRANCEf[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.max+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(LCsim_FRANCEf.Aarima$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(LCsim_FRANCEf.Aarima$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(LCsim_FRANCEf.Aarima$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_LC_FRANCEf[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))


## Fan-charts considerando il modello Renshaw-Haberman per la popolazione maschile

library("fanplot")
dev.new(width = 550, height = 330, unit = "px")
par(mfrow=c(1, 1))
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(RHfit_FRANCEm$years, t(q_RH_FRANCEm[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.max+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "male mortality rate (log scale)")
fan(t(RHsim_FRANCEm.Aarima$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(RHsim_FRANCEm.Aarima$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(RHsim_FRANCEm.Aarima$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_LC_FRANCEm[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))



## Fan-charts considerando il modello Plat per la popolazione femminile

library("fanplot")
dev.new(width = 550, height = 330, unit = "px")
par(mfrow=c(1, 1))
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(PLATfit_FRANCEf$years, t(q_PLAT_FRANCEf[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.max+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(PLATsim_FRANCEf.Aarima$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(PLATsim_FRANCEf.Aarima$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(PLATsim_FRANCEf.Aarima$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_LC_FRANCEm[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))



# Bootstrap 

## Lee-Carter, popolazione maschile

LCboot_FRANCEm <- bootstrap(LCfit_FRANCEm, nBoot = n.boot, type = "semiparametric")

par(mfrow=c(1, 1))
plot(LCboot_FRANCEm, nCol = 3) ## Rappresentazione


LCsim_FRANCEm.boot <- simulate(LCboot_FRANCEm, nsim = n.sim/n.boot, h = y.pred, kt.method = "iarima", kt.order = NULL)
rates_LC_FRANCEm.boot.st <- LCsim_FRANCEm.boot$rates
q_LC_FRANCEm.boot.st <- 1- exp(-rates_LC_FRANCEm.boot.st)
q_LC_FRANCEm.boot.st.ext <-  extrapolation.sim(q_LC_FRANCEm.boot.st)


## Lee-Carter, popolazione femminile

LCboot_FRANCEf <- bootstrap(LCfit_FRANCEf, nBoot = n.boot, type = "semiparametric")

par(mfrow=c(1, 1))
plot(LCboot_FRANCEf, nCol = 3) ## Rappresentazione

LCsim_FRANCEf.boot <- simulate(LCboot_FRANCEf, nsim = n.sim/n.boot, h = y.pred, kt.method = "iarima", kt.order = NULL)
rates_LC_FRANCEf.boot.st <- LCsim_FRANCEf.boot$rates
q_LC_FRANCEf.boot.st <- 1- exp(-rates_LC_FRANCEf.boot.st)
q_LC_FRANCEf.boot.st.ext <-  extrapolation.sim(q_LC_FRANCEf.boot.st)


## Renshaw-Haberman, popolazione maschile

RHboot_FRANCEm <- bootstrap(RHfit_FRANCEm, nBoot = n.boot, type = "semiparametric")

par(mfrow=c(1, 1))
plot(RHboot_FRANCEm, nCol = 3) ## Rappresentazione


RHsim_FRANCEm.boot <- simulate(RHboot_FRANCEm, nsim = n.sim/n.boot, h = y.pred, kt.method = "iarima", kt.order = NULL)
rates_RH_FRANCEm.boot.st <- RHsim_FRANCEm.boot$rates
q_RH_FRANCEm.boot.st <- 1- exp(-rates_RH_FRANCEm.boot.st)
q_RH_FRANCEm.boot.st.ext <-  extrapolation.sim(q_RH_FRANCEm.boot.st)


## Plat, popolazione femminile

PLATboot_FRANCEf<- bootstrap(PLATfit_FRANCEf, nBoot = n.boot, type = "semiparametric")

par(mfrow=c(1, 1))
plot(PLATboot_FRANCEf, nCol = 3) ## Rappresentazione

PLATsim_FRANCEf.boot <- simulate(PLATboot_FRANCEf, nsim = n.sim/n.boot, h = y.pred, kt.method = "iarima", kt.order = NULL)
rates_PLAT_FRANCEf.boot.st <- PLATsim_FRANCEf.boot$rates
q_PLAT_FRANCEf.boot.st <- 1- exp(-rates_PLAT_FRANCEf.boot.st)
q_PLAT_FRANCEf.boot.st.ext <-  extrapolation.sim(q_PLAT_FRANCEf.boot.st)


## Lee-Carter, popolazione maschile e femminile

# Intervalli di confidenza 

##### Probabilità di morte entro 1 anno

conf.lev <- 0.9
q_LC_FRANCEm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_FRANCEm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_FRANCEm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_FRANCEm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))

for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_LC_FRANCEm.q95[j,k] <- quantile(q_LC_FRANCEm.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_FRANCEm.q05[j,k] <- quantile(q_LC_FRANCEm.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_LC_FRANCEm.boot.q95[j,k] <- quantile(q_LC_FRANCEm.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_FRANCEm.boot.q05[j,k] <- quantile(q_LC_FRANCEm.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_LC_FRANCEf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_FRANCEf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_FRANCEf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_FRANCEf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))

for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_LC_FRANCEf.q95[j,k] <- quantile(q_LC_FRANCEf.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_FRANCEf.q05[j,k] <- quantile(q_LC_FRANCEf.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_LC_FRANCEf.boot.q95[j,k] <- quantile(q_LC_FRANCEf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_FRANCEf.boot.q05[j,k] <- quantile(q_LC_FRANCEf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_LC.FRANCEm.q95<-as.data.frame(q_LC_FRANCEm.q95)
q_LC.FRANCEm.q05<-as.data.frame(q_LC_FRANCEm.q05)
q_LC.FRANCEm.boot.q95<-as.data.frame(q_LC_FRANCEm.boot.q95)
q_LC.FRANCEm.boot.q05<-as.data.frame(q_LC_FRANCEm.boot.q05)

q_LC.FRANCEf.q95<-as.data.frame(q_LC_FRANCEf.q95)
q_LC.FRANCEf.q05<-as.data.frame(q_LC_FRANCEf.q05)
q_LC.FRANCEf.boot.q95<-as.data.frame(q_LC_FRANCEf.boot.q95)
q_LC.FRANCEf.boot.q05<-as.data.frame(q_LC_FRANCEf.boot.q05)

plot(q_LC_FRANCEm.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1,xlab = "years", ylab = "Mortality rates") # evolution with t (same age)
lines(q_LC_FRANCEm.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_FRANCEm.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_FRANCEm.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_FRANCEm.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)



plot(q_LC_FRANCEf.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1,xlab = "years", ylab = "Mortality rates")
lines(q_LC_FRANCEf.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_FRANCEf.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_FRANCEf.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_FRANCEf.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)



## Renshaw-Haberman, popolazione maschile

conf.lev <- 0.9
q_RH_FRANCEm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_RH_FRANCEm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_RH_FRANCEm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_RH_FRANCEm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))

for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_RH_FRANCEm.q95[j,k] <- quantile(q_RH_FRANCEm.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_RH_FRANCEm.q05[j,k] <- quantile(q_RH_FRANCEm.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_RH_FRANCEm.boot.q95[j,k] <- quantile(q_RH_FRANCEm.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_RH_FRANCEm.boot.q05[j,k] <- quantile(q_RH_FRANCEm.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}


q_RH.FRANCEm.q95<-as.data.frame(q_RH_FRANCEm.q95)
q_RH.FRANCEm.q05<-as.data.frame(q_RH_FRANCEm.q05)
q_RH.FRANCEm.boot.q95<-as.data.frame(q_RH_FRANCEm.boot.q95)
q_RH.FRANCEm.boot.q05<-as.data.frame(q_RH_FRANCEm.boot.q05)


plot(q_RH_FRANCEm.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1, xlab = "years", ylab = "Mortality rates")
lines(q_RH_FRANCEm.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_RH_FRANCEm.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_RH_FRANCEm.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_RH_FRANCEm.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)


## Plat, popolazione femminile

conf.lev <- 0.9
q_PLAT_FRANCEf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_FRANCEf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_FRANCEf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_PLAT_FRANCEf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))

for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_PLAT_FRANCEf.q95[j,k] <- quantile(q_PLAT_FRANCEf.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_PLAT_FRANCEf.q05[j,k] <- quantile(q_PLAT_FRANCEf.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_PLAT_FRANCEf.boot.q95[j,k] <- quantile(q_PLAT_FRANCEf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_PLAT_FRANCEf.boot.q05[j,k] <- quantile(q_PLAT_FRANCEf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}


q_PLAT.FRANCEf.q95<-as.data.frame(q_PLAT_FRANCEf.q95)
q_PLAT.FRANCEf.q05<-as.data.frame(q_PLAT_FRANCEf.q05)
q_PLAT.FRANCEf.boot.q95<-as.data.frame(q_PLAT_FRANCEf.boot.q95)
q_PLAAT.FRANCEf.boot.q05<-as.data.frame(q_PLAT_FRANCEf.boot.q05)


plot(q_PLAT_FRANCEf.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1, xlab = "years", ylab = "Mortality rates")
lines(q_PLAT_FRANCEf.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_FRANCEf.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_PLAT_FRANCEf.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_PLAT_FRANCEf.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

# Valore attuale di rendite vitalizie

## Lee-Carter, popolazione maschile e femminile

p_LC_FRANCEm.st.ext <- 1-q_LC_FRANCEm.st.ext
ann_LC_FRANCEm.st <- annuity.st(q_LC_FRANCEm.st.ext,v_vect,0.9)
ann_LC_FRANCEm.mean <- ann_LC_FRANCEm.st[[1]]
ann_LC_FRANCEm.q95 <- ann_LC_FRANCEm.st[[2]]
ann_LC_FRANCEm.q05 <- ann_LC_FRANCEm.st[[3]]

dev.new(width = 550, height = 330, unit = "px")
par(mfrow=c(1, 2))
plot(ann_LC_FRANCEm.mean["65",], type="l", lty=1) # evoluzione rispetto a t (stessa età x)
lines(ann_LC_FRANCEm.q95["65",], type="l", lty=2)
lines(ann_LC_FRANCEm.q05["65",], type="l", lty=2)
plot(ann_LC_FRANCEm.mean[,"2050"], type="l", lty=1) # evoluzione rispetto a x (stesso anno t)
lines(ann_LC_FRANCEm.q95[,"2050"], type="l", lty=2)
lines(ann_LC_FRANCEm.q05[,"2050"], type="l", lty=2)



p_LC_FRANCEf.st.ext <- 1-q_LC_FRANCEf.st.ext
ann_LC_FRANCEf.st <- annuity.st(q_LC_FRANCEf.st.ext,v_vect,0.9)
ann_LC_FRANCEf.mean <- ann_LC_FRANCEf.st[[1]]
ann_LC_FRANCEf.q95 <- ann_LC_FRANCEf.st[[2]]
ann_LC_FRANCEf.q05 <- ann_LC_FRANCEf.st[[3]]

dev.new(width = 550, height = 330, unit = "px")
par(mfrow=c(1, 2))
plot(ann_LC_FRANCEf.mean["65",], type="l", lty=1) # evoluzione rispetto a t (stessa età x)
lines(ann_LC_FRANCEf.q95["65",], type="l", lty=2)
lines(ann_LC_FRANCEf.q05["65",], type="l", lty=2)
plot(ann_LC_FRANCEf.mean[,"2050"], type="l", lty=1) # evoluzione rispetto a x (stesso anno t)
lines(ann_LC_FRANCEf.q95[,"2050"], type="l", lty=2)
lines(ann_LC_FRANCEf.q05[,"2050"], type="l", lty=2)

## Renshaw-Haberman, popolazione maschile

p_RH_FRANCEm.st.ext <- 1-q_RH_FRANCEm.st.ext
ann_RH_FRANCEm.st <- annuity.st(q_RH_FRANCEm.st.ext,v_vect,0.9)
ann_RH_FRANCEm.mean <- ann_RH_FRANCEm.st[[1]]
ann_RH_FRANCEm.q95 <- ann_RH_FRANCEm.st[[2]]
ann_RH_FRANCEm.q05 <- ann_RH_FRANCEm.st[[3]]

dev.new(width = 550, height = 330, unit = "px")
par(mfrow=c(1, 2))
plot(ann_RH_FRANCEm.mean["65",], type="l", lty=1) # evoluzione rispetto a t (stessa età x)
lines(ann_RH_FRANCEm.q95["65",], type="l", lty=2)
lines(ann_RH_FRANCEm.q05["65",], type="l", lty=2)
plot(ann_RH_FRANCEm.mean[,"2050"], type="l", lty=1) # evoluzione rispetto a x (stesso anno t)
lines(ann_RH_FRANCEm.q95[,"2050"], type="l", lty=2)
lines(ann_RH_FRANCEm.q05[,"2050"], type="l", lty=2)

## Plat, popolazione femminile

p_PLAT_FRANCEf.st.ext <- 1-q_PLAT_FRANCEf.st.ext
ann_PLAT_FRANCEf.st <- annuity.st(q_PLAT_FRANCEf.st.ext,v_vect,0.9)
ann_PLAT_FRANCEf.mean <- ann_PLAT_FRANCEf.st[[1]]
ann_PLAT_FRANCEf.q95 <- ann_PLAT_FRANCEf.st[[2]]
ann_PLAT_FRANCEf.q05 <- ann_PLAT_FRANCEf.st[[3]]

dev.new(width = 550, height = 330, unit = "px")
par(mfrow=c(1, 2))
plot(ann_PLAT_FRANCEf.mean["65",], type="l", lty=1) # evoluzione rispetto a t (stessa età x)
lines(ann_PLAT_FRANCEf.q95["65",], type="l", lty=2)
lines(ann_PLAT_FRANCEf.q05["65",], type="l", lty=2)
plot(ann_PLAT_FRANCEf.mean[,"2050"], type="l", lty=1) # evoluzione rispetto a x (stesso anno t)
lines(ann_PLAT_FRANCEf.q95[,"2050"], type="l", lty=2)
lines(ann_PLAT_FRANCEf.q05[,"2050"], type="l", lty=2)




