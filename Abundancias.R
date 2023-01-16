library(tidyverse)
library(lme4)
library(ggplot2)
library(performance)
library(readxl)
library(lubridate)
library(XNomial)
library(ggpubr)
library(rstatix)
# abundancias <- read_excel("Doctorado/2022/Abundance/Abundancias_2021.xlsx", 
#                                sheet = "primer_filtro", col_types = c("skip", 
#                                                                       "text", "text", "text", "text", "skip", 
#                                                                       "text", "numeric", "numeric", "numeric", 
#                                                                       "text", "text", "skip", "skip", "text", 
#                                                                       "date", "skip", "skip", "date", "skip", 
#                                                                       "skip", "skip", "skip", "numeric", 
#                                                                       "numeric", "text", "numeric", "text", 
#                                                                       "numeric", "numeric", "text"))
setwd("C:\\Users\\lblan\\OneDrive\\Escritorio\\CEAB\\2022\\Abundance_study")
abundancias <- read_excel("C:\\Users\\lblan\\OneDrive\\Escritorio\\CEAB\\2022\\Abundance_study\\Abundancias_2021.xlsx",
                               sheet = "primer_filtro", col_types = c("skip",
                                                                      "text", "text", "text", "text", "skip",
                                                                      "text", "numeric", "numeric", "numeric",
                                                                      "text", "text", "skip", "skip", "text",
                                                                      "date", "skip", "skip", "date", "skip",
                                                                      "skip", "skip", "skip", "numeric",
                                                                      "numeric", "text", "numeric", "text",
                                                                      "numeric", "numeric", "text"))

abundancias$species
abundancias_aedes <- filter(abundancias, species == "Aedes albopictus")
abundancias_aedes <- abundancias_aedes %>%
  rename("start_date" = `Start date`)
abundancias_aedes <- abundancias_aedes %>%
  rename("trap_name" = `trap name`)
abundancias_aedes_bl <- filter(abundancias_aedes, city == "Blanes")

ggplot(data = abundancias_aedes) +
  geom_line(aes(x= start_date, y= nrperspecies, color = trap_name)) #plot todas trampas

ggplot(data = abundancias_aedes_bl) +
  geom_line(aes(x= start_date, y= nrperspecies, color = trap_name)) #plot trampas blanes

#eliminamos la que no estaba en el jardín y la de la entrada (disaster)
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_1")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_2")


abundancias_aedes_bl <- abundancias_aedes_bl %>% #creamos columna mes
  add_column(month = NA)
abundancias_aedes_bl$month <- lubridate::month(ymd(abundancias_aedes_bl$start_date))

abundancias_aedes_bl <- abundancias_aedes_bl %>% #creamos columna semana
  add_column(week = NA)
abundancias_aedes_bl$week <- lubridate::week(ymd(abundancias_aedes_bl$start_date))

ggplot() +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_9")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_10")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_11")) 

abundancias_aedes_bl$trap_name = factor(abundancias_aedes_bl$trap_name) # convert to nominal factor
abs_bl = xtabs( nrperspecies ~ trap_name, data=abundancias_aedes_bl)
abs_bl

#Analysis of variance
checkin <- lm(nrperspecies ~ trap_name, data = abundancias_aedes_bl)
summary(checkin)
plot(checkin)

#normality
ggqqplot(abundancias_aedes_bl$nrperspecies)
ggdensity(abundancias_aedes_bl$nrperspecies, fill = "lightgray") #we can't assume normality

#Shapiro-test for grouped data to check normality
abundancias_aedes_bl %>%
  group_by(trap_name) %>%
  shapiro_test(nrperspecies) #not normal

#GENERALISED LINEAR MODELS IN RSTUDIO

#WHAT KIND OF DISTRIBUTION FOLLOWS OUR DATA?
library(fitdistrplus)
#Gamma
fg<-fitdist(abundancias_aedes_bl$nrperspecies, "gamma")
plot(fg)
summary(fg)
#lognormal
fln<-fitdist(abundancias_aedes_bl$nrperspecies, "lnorm")
plot(fln)
summary(fln)
#negative binomial
fnb<-fitdist(abundancias_aedes_bl$nrperspecies, "nbinom")
plot(fnb)
summary(fnb)
#weibull
fwei <- fitdist(abundancias_aedes_bl$nrperspecies, "weibull")
plot(fwei)
summary(fwei)
#poisson
fpois <- fitdist(abundancias_aedes_bl$nrperspecies, "pois")
plot(fpois)
summary(fpois)

#goodness-of-fit
fit <- list(fg, fln, fnb, fwei, fpois)
gofstat(fit, fitnames = NULL) #exponential gamma with lower AIC

#GENERALISED LINEAR MODELS
library(DHARMa)
library(performance)

#simple models
m0 <- glm(nrperspecies ~ 1, family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m0)            
m1 <- glmer(nrperspecies ~ week + (1|trap_name), family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m1)
plot(m1)
testResiduals(m1)
testZeroInflation(m1)
m2 <- glmer(nrperspecies ~week +`trapping effort` + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m2)
plot(m2)
testResiduals(m2)
testZeroInflation(m2)

m3 <- glmer(nrperspecies ~week + log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m3)
plot(m3)
testResiduals(m2)
testZeroInflation(m2)

AIC(m2,m3) #model 3 is better

#adding clima
clima <- read_excel("C:/Users/lblan/OneDrive/Escritorio/PrimerCap/Chapter-Survival/PrimerCap/Clima_aemet.xlsx", 
                    sheet = "Hoja2", 
                    col_types = c("date",
                                  "numeric", "numeric", "numeric", 
                                  "numeric",
                                  "numeric", "numeric", "numeric", 
                                  "numeric"))

clima <- clima %>%
  rename(start_date = "Date")
clima <- clima %>%
  rename(wind = "wind(km/h)")
library(pollen)
clima <- clima %>% 
  mutate(gdd = gdd(tmax = tmax, tmin = tmin, tbase = 10,
                   tbase_max = 30)) %>% 
  mutate(daily_acc_gdd = c(NA, diff(gdd)))

clima2 <-clima[,c("start_date","gdd","daily_acc_gdd")]

abundancias_aedes_bl <- inner_join(abundancias_aedes_bl, clima2, 
                         by = c("start_date"), all= TRUE) 

for(i in 1:nrow(abundancias_aedes_bl)){
  int <- seq(as.POSIXct(abundancias_aedes_bl$start_date[i], format = '%Y-%m-%d'),
             as.POSIXct(abundancias_aedes_bl$`End date`[i], format = '%Y-%m-%d'),
             by = 'day')

  tmean <- mean(clima$tmean[clima$start_date %in% int])
  tmax <- mean(clima$tmax[clima$start_date %in% int])
  tmin <- mean(clima$tmin[clima$start_date %in% int])
  Rhmean <- mean(clima$Rhmean[clima$start_date %in% int])
  Rhmin <- mean(clima$Rhmin[clima$start_date %in% int])
  Rhmax <- mean(clima$Rhmax[clima$start_date %in% int])
  rainfall <- mean(clima$rainfall[clima$start_date %in% int])
  wind <- mean(clima$wind[clima$start_date %in% int])
  gdd2 <- mean(clima$gdd[clima$start_date %in% int])
  
  abundancias_aedes_bl$tmean[i] <- tmean
  abundancias_aedes_bl$tmax[i] <- tmax
  abundancias_aedes_bl$tmin[i] <- tmin
  abundancias_aedes_bl$Rhmean[i] <- Rhmean
  abundancias_aedes_bl$Rhmin[i] <- Rhmin
  abundancias_aedes_bl$Rhmax[i] <- Rhmax
  abundancias_aedes_bl$rainfall[i] <- rainfall
  abundancias_aedes_bl$wind[i] <- wind
  abundancias_aedes_bl$gdd2[i] <- gdd2
    
}

m2 <- glmer(nrperspecies ~tmean +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m2)
testResiduals(m2)

m3 <- glmer(nrperspecies ~tmin +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m3)
testResiduals(m3)

m4 <- glmer(nrperspecies ~tmax +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m4)
testResiduals(m4)

m5 <- glmer(nrperspecies ~  Rhmax +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m5)
testResiduals(m5)
plot(m5)

m6 <- glmer(nrperspecies ~ Rhmean +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m6)
testResiduals(m6)
plot(m6)

m7 <- glmer(nrperspecies ~ Rhmin +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m7)
testResiduals(m7)
plot(m7)

m8 <- glmer(nrperspecies ~ rainfall +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m8)
testResiduals(m8)
plot(m8)

m9 <- glmer(nrperspecies ~ wind +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m9)
testResiduals(m9)
plot(m9)

m10 <- glmer(nrperspecies ~ scale(gdd) +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m10)
testResiduals(m10)
plot(m10)

m11 <- glmer(nrperspecies ~ scale(daily_acc_gdd) +log(`trapping effort`) + (1|trap_name),family=Gamma(link="log"), data = abundancias_aedes_bl)
summary(m11)
testResiduals(m11)
plot(m11)


abundancias_aedes_bl$doy <- yday(abundancias_aedes_bl$start_date)
modelo_doy<- glmer(nrperspecies ~ scale(doy) + (1|trap_name),
                   family=Gamma(link="log"),
                   control = glmerControl(optimizer = "Nelder_Mead",
                                          optCtrl = list(maxfun = 2000000)), data = abundancias_aedes_bl)
summary(modelo_doy)

#COLINEARITY
modelo<- glmer(nrperspecies ~ tmin + tmean + tmax + Rhmin +
                 Rhmean + Rhmax + `trapping effort` +(1|trap_name),
                   family=Gamma(link="log"), 
               control = glmerControl(optimizer = "Nelder_Mead",
                                      optCtrl = list(maxfun = 2e5)),
               data = abundancias_aedes_bl)
summary(modelo)
car::vif(modelo)

modelo<- glmer(nrperspecies ~ tmin + tmean + tmax + Rhmin +
                 Rhmax  + `trapping effort` +(1|trap_name),
               family=Gamma(link="log"), 
               control = glmerControl(optimizer = "Nelder_Mead",
                                      optCtrl = list(maxfun = 2e5)),
               data = abundancias_aedes_bl)
summary(modelo)
car::vif(modelo)

modelo<- glmer(nrperspecies ~ tmin  + tmax + Rhmin +
                 Rhmax  + `trapping effort` +(1|trap_name),
               family=Gamma(link="log"), 
               control = glmerControl(optimizer = "Nelder_Mead",
                                      optCtrl = list(maxfun = 2e5)),
               data = abundancias_aedes_bl)
summary(modelo)
car::vif(modelo)

modelo<- glmer(nrperspecies ~  tmax + Rhmin +
                 Rhmax  + `trapping effort` +(1|trap_name),
               family=Gamma(link="log"), 
               control = glmerControl(optimizer = "Nelder_Mead",
                                      optCtrl = list(maxfun = 2e5)),
               data = abundancias_aedes_bl)
summary(modelo)
car::vif(modelo)

modelo<- glmer(nrperspecies ~  tmax + Rhmin + scale(doy) +
                 Rhmax  + `trapping effort` +(1|trap_name),
               family=Gamma(link="log"), 
               control = glmerControl(optimizer = "Nelder_Mead",
                                      optCtrl = list(maxfun = 2e5)),
               data = abundancias_aedes_bl)
summary(modelo)
car::vif(modelo)

modelo<- glmer(nrperspecies ~  tmax + Rhmin + scale(doy) + wind +
                 Rhmax  + `trapping effort` +(1|trap_name),
               family=Gamma(link="log"), 
               control = glmerControl(optimizer = "Nelder_Mead",
                                      optCtrl = list(maxfun = 2e5)),
               data = abundancias_aedes_bl)
summary(modelo)
car::vif(modelo)

# modelo<- glmer(nrperspecies ~  tmax + Rhmin + scale(doy) + wind + rainfall+
#                  Rhmax  + `trapping effort` +(1|trap_name),
#                family=Gamma(link="log"), 
#                control = glmerControl(optimizer = "Nelder_Mead",
#                                       optCtrl = list(maxfun = 2e5)),
#                data = abundancias_aedes_bl)
# summary(modelo)
# car::vif(modelo) 

#EN ESTE MODELO HEMOS AÑADIDO RAINFALL Y PARECE QUE PETA... TENDRÍAMOS QUE VER
#SI EN VEZ DE LA MEDIA POR CADA HORQUILLA DE TRAPPING EFFORT CON RAINFALL SERÍA
#MÁS ADECUADO HACER UNA ACUMULACIÓN DE RAIN ???

modelo<- glmer(nrperspecies ~  tmax + Rhmin + scale(doy) + rainfall +
                 Rhmax  + `trapping effort` +(1|trap_name),
               family=Gamma(link="log"), 
               control = glmerControl(optimizer = "Nelder_Mead",
                                      optCtrl = list(maxfun = 2e5)),
               data = abundancias_aedes_bl)
summary(modelo)
car::vif(modelo)
#SI SUSTITUIMOS WIND POR RAINFALL NO PETA!! PERO NO SE SI CREERME EL RESULTADO 
#POR LO QUE HE DICHO ANTES... MEJOR ACUMULACIÓN QUE MEDIA?

modelo<- glmer(nrperspecies ~  scale(gdd) + Rhmin + scale(doy) + wind +
                 Rhmax  + `trapping effort` +(1|trap_name),
               family=Gamma(link="log"), 
               control = glmerControl(optimizer = "Nelder_Mead",
                                      optCtrl = list(maxfun = 2e5)),
               data = abundancias_aedes_bl)
summary(modelo) #con gdd nos peta un poco y encima lo que sale no tiene ningún sentido
car::vif(modelo) #aqui la respuesta: colineal con doy!!!!!

# modelo<- glmer(nrperspecies ~  scale(gdd) + Rhmin  + wind +
#                  Rhmax  + `trapping effort` +(1|trap_name),
#                family=Gamma(link="log"), 
#                control = glmerControl(optimizer = "Nelder_Mead",
#                                       optCtrl = list(maxfun = 2e5)),
#                data = abundancias_aedes_bl)
# summary(modelo)
# car::vif(modelo)
# 
# modelo<- glmer(nrperspecies ~  scale(gdd2)   + wind +
#                  Rhmax  + `trapping effort` +(1|trap_name),
#                family=Gamma(link="log"), 
#                control = glmerControl(optimizer = "Nelder_Mead",
#                                       optCtrl = list(maxfun = 2e5)),
#                data = abundancias_aedes_bl)
# summary(modelo) #no me convence lo de que el gdd esté así relacionado con la abundancia
# car::vif(modelo)

#preguntar fede
#range(abundancias_aedes_bl$gdd)
#range(clima$gdd)

