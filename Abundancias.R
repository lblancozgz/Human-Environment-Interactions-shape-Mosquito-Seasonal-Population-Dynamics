library(tidyverse)
library(lme4)
library(ggplot2)
library(performance)
library(readxl)
library(lubridate)
library(XNomial)
library(ggpubr)
library(rstatix)
library(DHARMa)
library(performance)
library(pollen)
library(meteor)
library(zoo)
library(parallelly)
library(parallel)
library(MuMIn)
library(glmmTMB)
library(GLMMadaptive)
library(car)
library(patchwork)
# library(curl)
# library(jsonlite)
# library(sf)
# library(leaflet)
# library(mgcv)  # Load the mgcv library for splines
library(httr)
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

#Limpieza de datos
abundancias$species
abundancias_aedes <- filter(abundancias, species == "Aedes albopictus")

abundancias_aedes <- abundancias_aedes %>%
  rename("start_date" = `Start date`)
abundancias_aedes <- abundancias_aedes %>%
  rename("end_date" = `End date`)
abundancias_aedes <- abundancias_aedes %>%
  rename("trap_name" = `trap name`)
abundancias_aedes_bl <- filter(abundancias_aedes, city == "Blanes")

abundancias_aedes_bl <- abundancias_aedes_bl %>% 
  rename(trapef = 'trapping effort')
abundancias_aedes_bl <- abundancias_aedes_bl %>% 
  rename(land_use = 'Land use')


#eliminamos la que no estaba en el jardín y la de la entrada (disaster)
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_1")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_2")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_10")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_11")

abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_12")

#añadimos nuevas variables temporales
abundancias_aedes_bl <- abundancias_aedes_bl %>% #creamos columna mes
  add_column(month = NA)
abundancias_aedes_bl$month <- lubridate::month(ymd(abundancias_aedes_bl$start_date))
abundancias_aedes_bl$month <- as.factor(abundancias_aedes_bl$month)
abundancias_aedes_bl <- abundancias_aedes_bl %>% #creamos columna semana
  add_column(week = NA)
abundancias_aedes_bl$week <- lubridate::week(ymd(abundancias_aedes_bl$start_date))
abundancias_aedes_bl$week <- as.factor(abundancias_aedes_bl$week)

abundancias_aedes_bl <- abundancias_aedes_bl %>% #creamos columna semana
  add_column(YEAR = NA)
abundancias_aedes_bl$YEAR <- lubridate::year(ymd(abundancias_aedes_bl$start_date))
abundancias_aedes_bl$doy <- yday(abundancias_aedes_bl$start_date)
abundancias_aedes_bl$YEAR <- as.factor(abundancias_aedes_bl$YEAR)


# Visualización de los datos
ggplot() +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_9")) 

abundancias_aedes_bl$trap_name = factor(abundancias_aedes_bl$trap_name) # convert to nominal factor
str(abundancias_aedes_bl$trap_name)

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

#datos clima estación malgrat
clima_malgrat <- read_csv("Malgrat-de-mar-2020-2022.csv", 
                                     col_types = cols(service = col_skip(), 
                                                      station_id = col_skip(), station_name = col_skip(), 
                                                      station_province = col_skip(), altitude = col_skip(), 
                                                      mean_wind_direction = col_skip(), 
                                                      FW = col_skip(), FT = col_skip(), 
                                                      mwi = col_skip()))
clima_malgrat <- clima_malgrat %>% #creamos columna mes
  add_column(month = NA)
clima_malgrat$month <- lubridate::month(ymd(clima_malgrat$timestamp))

clima_malgrat <- clima_malgrat %>% #creamos columna semana
  add_column(week = NA)
clima_malgrat$week <- lubridate::week(ymd(clima_malgrat$timestamp))

clima_malgrat <- clima_malgrat %>% #creamos columna semana
  add_column(YEAR = NA)
clima_malgrat$YEAR <- lubridate::year(ymd(clima_malgrat$timestamp))
clima_malgrat <- filter(clima_malgrat, !YEAR==2022)
clima_malgrat <- filter(clima_malgrat, !YEAR==2020)

clima_malgrat$doy <- yday(clima_malgrat$timestamp) #calculamos doy
clima_malgrat$month <- as.factor(clima_malgrat$month)
clima_malgrat$week <- as.factor(clima_malgrat$week)
clima_malgrat$YEAR <- as.factor(clima_malgrat$YEAR)

#Visualización clima
climatemp <- clima_malgrat %>%
  ggplot() +
  geom_line(aes(x = doy, y = max_temperature), col= "red")+
  geom_line(aes(x = doy, y = min_temperature), col= "darkblue")+
  geom_line(aes(x = doy, y = mean_temperature), col="darkgreen")+
  facet_wrap(~YEAR) +
  ylab("Temperatures") +
  theme_bw()

climarh <- clima_malgrat %>%
  ggplot() +
  geom_line(aes(x = doy, y = max_relative_humidity), col= "red")+
  geom_line(aes(x = doy, y = min_relative_humidity), col= "darkblue")+
  geom_line(aes(x = doy, y = mean_relative_humidity), col="darkgreen")+
  facet_wrap(~YEAR) +
  ylab("RH") +
  theme_bw()

climalluvia <- clima_malgrat %>%
  group_by(week, month, YEAR) %>% 
  summarise(rainfall = sum(precipitation)) %>% 
  ggplot() +
  geom_col(aes(x = month, y = rainfall))+
  facet_wrap(~YEAR) +
  ylab("Rain") +
  theme_bw()

clima_malgrat <- clima_malgrat %>%
  rename(start_date = "timestamp")

#mwi
clima_malgrat= clima_malgrat %>% mutate(
  FHme = case_when(mean_relative_humidity < 40~0, mean_relative_humidity >95~0, 
                   (mean_relative_humidity >=40 & mean_relative_humidity <= 95)~
                     ((mean_relative_humidity/55)-(40/55)) ),
  FTme = case_when(mean_temperature<=15~0, mean_temperature>30~0, 
                   (mean_temperature>15 & mean_temperature <=20)~
                     (.2*mean_temperature)-3,
                   (mean_temperature>20 & mean_temperature<=25)~1,
                   (mean_temperature>25 & mean_temperature <= 30)~
                     (-.2*mean_temperature)+6),
  mwime = FHme*FTme)
#photoperiod
clima_malgrat$photoperiod <- photoperiod(1:365, 41.676854)
#medias y acumuladas de ventanas temporales previas 7/14/21 días 
clima_malgrat <- clima_malgrat %>% 
  mutate(mintemp7 = rollmean(min_temperature, k = 7, align = "right", fill = NA),
         mintemp14 = rollmean(min_temperature, k = 14, align = "right", fill = NA),
         mintemp21 = rollmean(min_temperature, k = 21, align = "right", fill = NA),
         meantemp7 = rollmean(mean_temperature, k = 7, align = "right", fill = NA),
         meantemp14 = rollmean(mean_temperature, k = 14, align = "right", fill = NA),
         meantemp21 = rollmean(mean_temperature, k = 21, align = "right", fill = NA),
         pho7 = rollmean(photoperiod, k = 7, align = "right", fill = NA),
         pho14 = rollmean(photoperiod, k = 14, align = "right", fill = NA),
         pho21 = rollmean(photoperiod, k = 21, align = "right", fill = NA),
         maxtemp7 = rollmean(max_temperature, k = 7, align = "right", fill = NA),
         maxtemp14 = rollmean(max_temperature, k = 14, align = "right", fill = NA),
         maxtemp21 = rollmean(max_temperature, k=21, align = "right", fill = NA),
         mwi7 = rollmean(mwime, k = 7, align = "right", fill = NA),
         mwi14 = rollmean(mwime, k = 14, align = "right", fill = NA),
         mwi21 = rollmean(mwime, k=21, align = "right", fill = NA),
         minrh7 = rollmean(min_relative_humidity, k = 7, align = "right", fill = NA),
         minrh14 = rollmean(min_relative_humidity, k = 14, align = "right", fill = NA),
         minrh21 = rollmean(min_relative_humidity, k = 21, align = "right", fill = NA),
         meanrh7 = rollmean(mean_relative_humidity, k = 7, align = "right", fill = NA),
         meanrh14 = rollmean(mean_relative_humidity, k = 14, align = "right", fill = NA),
         meanrh21 = rollmean(mean_relative_humidity, k = 21, align = "right", fill = NA),
         maxrh7 = rollmean(max_relative_humidity, k = 7, align = "right", fill = NA),
         maxrh14 = rollmean(max_relative_humidity, k = 14, align = "right", fill = NA),
         maxrh21 = rollmean(max_relative_humidity, k = 21, align = "right", fill = NA),
         rain7 = rollapplyr(precipitation, 7, FUN=sum, align = "right", fill = NA),
         rain14 = rollapplyr(precipitation, 14, FUN=sum, align = "right", fill = NA),
         rain21 = rollapplyr(precipitation, 21, FUN=sum, align = "right", fill = NA))

str(clima_malgrat$YEAR)
range(abundancias_aedes_bl$nrperspecies)

#chequeamos nuevas variables
range(clima_malgrat$mwime, na.rm = T) #debe estar comprendido entre 0 y 1
ggplot(clima_malgrat, aes(x=doy, y=mwi7, color=YEAR)) + geom_line()+ facet_wrap(~YEAR)
ggplot(clima_malgrat, aes(x=doy, y=photoperiod, color=YEAR)) + geom_line() + facet_wrap(~YEAR)
ggplot(clima_malgrat, aes(x=doy, y=pho21, color=YEAR)) + geom_line() + facet_wrap(~YEAR)

#nos aseguramos que las fechas estén en el mismo formato
clima_malgrat$start_date <- as.POSIXct(clima_malgrat$start_date,format = '%Y-%m-%d')
abundancias_aedes_bl$start_date <- as.POSIXct(abundancias_aedes_bl$start_date,format = '%Y-%m-%d')

#Fusionamos clima + datos abundancias
abundancias_jardin_tot <- full_join(abundancias_aedes_bl, clima_malgrat, 
                                   by = c("start_date", "YEAR", "week", "month", "doy"), all= TRUE) 


#importante trapping effort y calcular las medias de cada variable climática
#por el trapping effort

for(i in 1:nrow(abundancias_jardin_tot)){
  x0 <- as.POSIXct(abundancias_jardin_tot$start_date[i], format = '%Y-%m-%d')
  x1 <- as.POSIXct(abundancias_jardin_tot$end_date[i], format = '%Y-%m-%d')
  if(is.na(x0) || is.na(x1)){
    next
  }
  int <- seq(x0, x1, by= "day")


  mean_temperature <- mean(clima_malgrat$mean_temperature[clima_malgrat$start_date %in% int])
  max_temperature <- mean(clima_malgrat$max_temperature[clima_malgrat$start_date %in% int])
  min_temperature <- mean(clima_malgrat$min_temperature[clima_malgrat$start_date %in% int])
  mean_relative_humidity <- mean(clima_malgrat$mean_relative_humidity[clima_malgrat$start_date %in% int])
  min_relative_humidity <- mean(clima_malgrat$min_relative_humidity[clima_malgrat$start_date %in% int])
  max_relative_humidity <- mean(clima_malgrat$max_relative_humidity[clima_malgrat$start_date %in% int])
  precipitation <- sum(clima_malgrat$precipitation[clima_malgrat$start_date %in% int])
  photoperiod <- mean(clima_malgrat$photoperiod[clima_malgrat$start_date %in% int])
  mwime <- mean(clima_malgrat$mwime[clima_malgrat$start_date %in% int])
  rain7 <- clima_malgrat$rain7[clima_malgrat$start_date %in% x0]
  rain14 <- clima_malgrat$rain14[clima_malgrat$start_date %in% x0]
  rain21 <-clima_malgrat$rain21[clima_malgrat$start_date %in% x0]
  mintemp7 <- clima_malgrat$mintemp7[clima_malgrat$start_date %in% x0]
  mintemp14 <- clima_malgrat$mintemp14[clima_malgrat$start_date %in% x0]
  mintemp21 <- clima_malgrat$mintemp21[clima_malgrat$start_date %in% x0]
  minrh7 <- clima_malgrat$minrh7[clima_malgrat$start_date %in% x0]
  minrh14 <- clima_malgrat$minrh14[clima_malgrat$start_date %in% x0]
  minrh21 <- clima_malgrat$minrh21[clima_malgrat$start_date %in% x0]
  meantemp7 <- clima_malgrat$meantemp7[clima_malgrat$start_date %in% x0]
  meantemp14 <- clima_malgrat$meantemp14[clima_malgrat$start_date %in% x0]
  meantemp21 <- clima_malgrat$meantemp21[clima_malgrat$start_date %in% x0]
  meanrh7 <- clima_malgrat$meanrh7[clima_malgrat$start_date %in% x0]
  meanrh14 <- clima_malgrat$meanrh14[clima_malgrat$start_date %in% x0]
  meanrh21 <- clima_malgrat$meanrh21[clima_malgrat$start_date %in% x0]
  maxtemp7 <- clima_malgrat$maxtemp7[clima_malgrat$start_date %in% x0]
  maxtemp14 <- clima_malgrat$maxtemp14[clima_malgrat$start_date %in% x0]
  maxtemp21 <- clima_malgrat$maxtemp21[clima_malgrat$start_date %in% x0]
  maxrh7 <- clima_malgrat$maxrh7[clima_malgrat$start_date %in% x0]
  maxrh14 <- clima_malgrat$maxrh14[clima_malgrat$start_date %in% x0]
  maxrh21 <- clima_malgrat$maxrh21[clima_malgrat$start_date %in% x0]
  pho7 <- clima_malgrat$pho7[clima_malgrat$start_date %in% x0]
  pho14 <- clima_malgrat$pho14[clima_malgrat$start_date %in% x0]
  pho21 <- clima_malgrat$pho21[clima_malgrat$start_date %in% x0]
  mwi7 <- clima_malgrat$mwi7[clima_malgrat$start_date %in% x0]
  mwi14 <- clima_malgrat$mwi14[clima_malgrat$start_date %in% x0]
  mwi21 <- clima_malgrat$mwi21[clima_malgrat$start_date %in% x0]
  
  abundancias_jardin_tot$mean_temperature[i] <- mean_temperature
  abundancias_jardin_tot$max_temperature[i] <- max_temperature
  abundancias_jardin_tot$min_temperature[i] <- min_temperature
  abundancias_jardin_tot$mean_relative_humidity[i] <- mean_relative_humidity
  abundancias_jardin_tot$min_relative_humidity[i] <- min_relative_humidity
  abundancias_jardin_tot$max_relative_humidity[i] <- max_relative_humidity
  abundancias_jardin_tot$precipitation[i] <- precipitation
  abundancias_jardin_tot$mwime[i] <- mwime
  abundancias_jardin_tot$photoperiod[i] <- photoperiod
  abundancias_jardin_tot$pho7[i] <- pho7
  abundancias_jardin_tot$pho14[i] <- pho14
  abundancias_jardin_tot$pho21[i] <- pho21
  abundancias_jardin_tot$rain7[i] <- rain7
  abundancias_jardin_tot$rain14[i] <- rain14
  abundancias_jardin_tot$rain21[i] <- rain21
  abundancias_jardin_tot$mintemp7[i] <- mintemp7
  abundancias_jardin_tot$mintemp14[i] <- mintemp14
  abundancias_jardin_tot$mintemp21[i] <- mintemp21
  abundancias_jardin_tot$minrh7[i] <- minrh7
  abundancias_jardin_tot$minrh14[i] <- minrh14
  abundancias_jardin_tot$minrh21[i] <- minrh21
  abundancias_jardin_tot$meantemp7[i] <- meantemp7
  abundancias_jardin_tot$meantemp14[i] <- meantemp14
  abundancias_jardin_tot$meantemp21[i] <- meantemp21
  abundancias_jardin_tot$meanrh7[i] <- meanrh7
  abundancias_jardin_tot$meanrh14[i] <- meanrh14
  abundancias_jardin_tot$meanrh21[i] <- meanrh21
  abundancias_jardin_tot$maxtemp7[i] <- maxtemp7
  abundancias_jardin_tot$maxtemp14[i] <- maxtemp14
  abundancias_jardin_tot$maxtemp21[i] <- maxtemp21
  abundancias_jardin_tot$maxrh7[i] <- maxrh7
  abundancias_jardin_tot$maxrh14[i] <- maxrh14
  abundancias_jardin_tot$maxrh21[i] <- maxrh21
  abundancias_jardin_tot$mwi7[i] <- mwi7
  abundancias_jardin_tot$mwi14[i] <- mwi14
  abundancias_jardin_tot$mwi21[i] <- mwi21
  
}

abundancias_jardin_tot <- filter(abundancias_jardin_tot, !is.na(trapef))

#GDD
# row by row
# Requires paraell or parallelly package, depends on the R version
this_result <- data.frame() # Mandatory object instead not using
abundances_with_gdd <-  bind_rows(mclapply(1:nrow(abundancias_jardin_tot), function(i){
  # Verbose
  if (i %% 50 == 0) {print(i)}
  
  x = abundancias_jardin_tot[i,]
  this_start_date = x$start_date
  this_end_date = x$end_date
  
  # Exposure week
  wth_v <- clima_malgrat %>%
    filter(start_date >= this_start_date & start_date <= this_end_date)
  
  l0gdd = gdd(tmax = wth_v$max_temperature,
            tmin = wth_v$min_temperature,
            tbase = 10,
            tbase_max = 30)
  
  weather_info <- x %>%
    mutate(
      l0gdd = tail(l0gdd, 1) 
    )
  
  # lag 7 days
  wth_v <- clima_malgrat %>%
    filter(start_date >= this_start_date - days(7) & start_date <= this_end_date) 

  l7gdd = gdd(tmax = wth_v$max_temperature,
              tmin = wth_v$min_temperature,
              tbase = 10,
              tbase_max = 30)
  
  weather_info <- weather_info %>%
    mutate(
      l7gdd = tail(l7gdd, 1)
    )
  
  # lag 14 days
  wth_v <- clima_malgrat %>%
    filter(start_date >= this_start_date - days(14) & start_date <= this_end_date) 
 
  l14gdd = gdd(tmax = wth_v$max_temperature,
               tmin = wth_v$min_temperature,
               tbase = 10,
               tbase_max = 30)
  weather_info <- weather_info %>%
    mutate(
      l14gdd = tail(l14gdd, 1)
    )
  
  # lag 21 days
  wth_v <- clima_malgrat %>%
    filter(start_date >= this_start_date - days(21) & start_date <= this_end_date) 
  
  l21gdd = gdd(tmax = wth_v$max_temperature,
               tmin = wth_v$min_temperature,
               tbase = 10,
               tbase_max = 30)
  weather_info <- weather_info %>%
    mutate(
      l21gdd = tail(l21gdd, 1)
      )
}))




options(na.action = "na.fail") #para hacer el dredge() es necesario
which(is.na(abundances_with_gdd$nrperspecies)) #remove all NAs that can affect negatively to the df

#MODELOS
#function collinearity
max.r <- function(x){
  corm <- cov2cor(vcov(x))
  corm <- as.matrix(corm)
  if (length(corm)==1){
    corm <- 0
    max(abs(corm))
  } else if (length(corm)==4){
    cormf <- corm[2:nrow(corm),2:ncol(corm)]
    cormf <- 0
    max(abs(cormf))
  } else {
    cormf <- corm[2:nrow(corm),2:ncol(corm)]
    diag(cormf) <- 0
    max(abs(cormf))
  }
}


#PRIMERA SERIE DE MODELOS - WEEK SAMPLING - TODAS LAS VARIABLES
m0_poly <-glmer.nb(nrperspecies ~ scale(max_temperature) +
                 scale(min_temperature) +
                scale(mean_temperature) +
                scale(max_relative_humidity) + scale(mean_relative_humidity) + scale(min_relative_humidity) +
                scale(precipitation) + scale(mwime) + log(trapef) + scale(l0gdd) +
                (1|trap_name),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
              data = abundances_with_gdd)
summary(m0_poly)
vif(m0_poly)

#UTILIZAMOS DREDGE PARA ELECCIÓN DE VARIABLES
ms3 <- dredge(m0_poly, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))

NCM3 <- get.models(ms3, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)

NCMDF <- model.sel(NCM3)

#with poly, we got polymintemp, logtrapef, max_rl_hmd, min_relt_humd, precipitation
#without poly, we mintemp, logtrapef, max_rel_hm, min_rel_humd, precipit


#Modelo resultante con las variables climáticas de la semana de muestreo
m0_dredge <- glmer.nb(nrperspecies ~ scale(max_relative_humidity) + scale(min_relative_humidity) + 
                              scale(min_temperature) + scale(precipitation) + log(trapef) + (1|trap_name),
                             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                             data = abundances_with_gdd)
summary(m0_dredge)
testZeroInflation(m0_dredge)
testResiduals(m0_dredge)
testDispersion(m0_dredge)


#SET DE VARIABLES - 7 DÍAS PREVIOS A LA SEMANA DE MUESTREO
m1 <-glmer.nb(nrperspecies ~ scale(maxtemp7) + scale(mintemp7) + scale(meantemp7) +
                scale(maxrh7) + scale(meanrh7) + scale(minrh7) +
                scale(rain7) + scale(mwi7) + log(trapef) + scale(l7gdd) +
                (1|trap_name), 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
              data = abundances_with_gdd)

summary(m1)


ms4 <- dredge(m1, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))
NCM4 <- get.models(ms4, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDF4 <- model.sel(NCM4)
#with poly, we got log(trapef), poly(maxtemp7,2), maxrh7, minrh7, rain7
#without poly, we got log(trapef), maxrh7, maxtemp7, minrh7, rain7

#Modelo con las variables obtenidas después de dredge para 7 días previos
m1_dredge <-glmer.nb(nrperspecies ~  log(trapef)+ scale(maxrh7) + scale(minrh7) +
                              scale(rain7) + scale(maxtemp7)+
                              (1|trap_name), 
                            control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                            data = abundances_with_gdd) 
summary(m1_dredge)
car::vif(m1_dredge)
testZeroInflation(m1_dredge)
testResiduals(m1_dredge)
testDispersion(m1_dredge)

#SET DE VARIABLES - 14 DÍAS PREVIOS A LA SEMANA DE MUESTREO

m2 <-glmer.nb(nrperspecies ~ scale(maxtemp14) + scale(mintemp14) + scale(meantemp14) +
                scale(maxrh14) + scale(meanrh14) + scale(minrh14) +
                scale(rain14) + scale(mwi14) + log(trapef) + scale(l14gdd) +
                (1|trap_name), 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
              data = abundances_with_gdd)

ms5 <- dredge(m2, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))
NCM5 <- get.models(ms5, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDF5 <- model.sel(NCM5)

#Modelo con las variables obtenidas después de dredge para 14 días previos

m2_dredge <-glmer.nb(nrperspecies ~  log(trapef) + scale(maxrh14) + scale(minrh14) + scale(mintemp14) + scale(rain14) +
                       (1|trap_name), 
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                     data = abundances_with_gdd)
summary(m2_dredge)
car::vif(m2_dredge)
testZeroInflation(m2_dredge)
testResiduals(m2_dredge)
testDispersion(m2_dredge)

#SET DE VARIABLES - 21 DÍAS PREVIOS A LA SEMANA DE MUESTREO

m3 <-glmer.nb(nrperspecies ~ scale(maxtemp21) + scale(mintemp21) + scale(meantemp21) +
                scale(maxrh21) + scale(meanrh21) + scale(minrh21) +
                scale(rain21) + scale(mwi21) + log(trapef) + scale(l21gdd) +
                (1|trap_name), 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
              data = abundances_with_gdd)

summary(m3)
ms6 <- dredge(m3, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))
NCM6 <- get.models(ms6, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDF6 <- model.sel(NCM6)

#Modelo con las variables obtenidas después de dredge para 21 días previos

m3_dredge <-glmer.nb(nrperspecies ~ log(trapef) + scale(maxrh21) + scale(minrh21) + scale(mintemp21) + scale(rain21) +
                       (1|trap_name), 
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                     data = abundances_with_gdd)

summary(m3_dredge)
car::vif(m3_dredge)
testZeroInflation(m3_dredge)
testResiduals(m3_dredge)
testDispersion(m3_dredge)

#Full model con todas las variables climáticas (no añadiendo turistas, embornales...)
fullmodel <- glmer.nb(nrperspecies ~ log(trapef) + scale(maxrh21) + scale(minrh21)
                      +scale(mintemp21) + scale(rain21) + scale(maxrh14) + scale(minrh14) +
                        scale(mintemp14) + scale(rain14) + scale(maxrh7) + scale(minrh7) +
                        scale(maxtemp7) + scale(rain7) + scale(max_relative_humidity) + 
                        scale(min_relative_humidity) +
                        poly(scale(min_temperature),2) + scale(precipitation) + (1|trap_name), 
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                      data = abundances_with_gdd) 
mfull <- dredge(fullmodel, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))
NCMfull <- get.models(mfull, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDFfull <- model.sel(NCMfull)
selection_full <- as.data.frame(NCMDFfull)
full_dredge <- glmer.nb(nrperspecies ~ log(trapef) + scale(maxrh14) + scale(mintemp21) +
                          scale(minrh7) + scale(precipitation) +(1|trap_name), 
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                        data = abundances_with_gdd) 
summary(full_dredge)
car::vif(full_dredge)
testZeroInflation(full_dredge)
testResiduals(full_dredge)
testDispersion(full_dredge)
testUniformity(full_dredge)

##winner model
NEWd <- data.frame(minrh7 = abundances_with_gdd$minrh7,
                   mintemp21 = abundances_with_gdd$mintemp21,
                   maxrh14 = abundances_with_gdd$maxrh14,
                   precipitation = abundances_with_gdd$precipitation,
                   trapef = 8,
                   trap_name=  "A_SP_BL_15")

plot(predict(full_dredge, NEWd, type="response", allow.new.levels = TRUE))
predict(full_dredge, NEWd, type="response", allow.new.levels = TRUE)
predictiondata <- as.data.frame(predict(full_dredge, NEWd, type="response", allow.new.levels = TRUE))
predictiondata$predi <- predictiondata$`predict(full_dredge, NEWd, type = "response", allow.new.levels = TRUE)`
abundances <- cbind(abundances_with_gdd, predictiondata)
abundances$week <- as.character(abundances$week)
abundances$week <- as.numeric(abundances$week)

ggplot() +
  geom_line(aes(x= week, y= nrperspecies), data = filter(abundances, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= week, y= nrperspecies), data = filter(abundances, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= week, y= nrperspecies), data = filter(abundances, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= week, y= nrperspecies), data = filter(abundances, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= week, y= nrperspecies), data = filter(abundances, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= week, y= nrperspecies), data = filter(abundances, trap_name == "A_SP_BL_9")) +
  scale_x_continuous(breaks = seq(from = 18, to = 46, by = 3))+
  geom_line(aes(x= week, y= predi),col="red", size = 1.1, data = abundances) +
  xlab("Year 2021") + ylab("Abundance") +
  theme_bw()


# Añadimos datos de turistas

datatourist <- read_excel("Origen i Pernocta 2021_XL.xlsx", 
                                        col_types = c("date", "numeric"))

datatourist$start_date<- strftime(datatourist$start_date, format="%Y-%m-%d")
datatourist$start_date <- as_date(as.POSIXct(datatourist$start_date,format = '%Y-%m-%d'))
abundances_with_gdd$start_date <- as_date(as.POSIXct(abundances_with_gdd$start_date,format = '%Y-%m-%d'))


#Tratamiento de datos de turistas

datatourist <- datatourist %>% 
  group_by(start_date) %>% 
  summarise(nrturist = sum(Tourist))
str(datatourist$nrturist)
sum(duplicated(datatourist$start_date))
datatourist <- datatourist %>% #creamos columna mes
  add_column(month = NA)
datatourist$month <- lubridate::month(ymd(datatourist$start_date))
datatourist$month <- as.factor(datatourist$month)
datatourist <- datatourist %>% #creamos columna semana
  add_column(week = NA)
datatourist$week <- lubridate::week(ymd(datatourist$start_date))
datatourist$week <- as.factor(datatourist$week)
#calculamos con rollmean variables como el número de turistas de los últimos 7,14,21 días
datatourist <- datatourist %>% 
  mutate(nturist7 = rollapplyr(nrturist, 7, FUN=sum, align = "right", fill = NA),
         nturist14 = rollapplyr(nrturist, 14, FUN=sum, align = "right", fill = NA),
         nturist21 = rollapplyr(nrturist, 21, FUN=sum, align = "right", fill = NA))

datatourist$start_date <- as.POSIXct(datatourist$start_date,format = '%Y-%m-%d')
abundances_with_gdd$start_date <- as.POSIXct(abundances_with_gdd$start_date,format = '%Y-%m-%d')

#juntamos clima con datos de turistas 
abundances_with_tourist<- merge(abundances_with_gdd, datatourist, by = c("start_date", "week", "month"), all.x = T) 
for(i in 1:nrow(abundances_with_tourist)){
  x0 <- as.POSIXct(abundances_with_tourist$start_date[i], format = '%Y-%m-%d')
  x1 <- as.POSIXct(abundances_with_tourist$end_date[i], format = '%Y-%m-%d')
  if(is.na(x0) || is.na(x1)){
    next
  }
  int <- seq(x0, x1, by= "day")
  nturist7 <- datatourist$nturist7[datatourist$start_date %in% x0]
  nturist14 <- datatourist$nturist14[datatourist$start_date %in% x0]
  nturist21 <- datatourist$nturist21[datatourist$start_date %in% x0]
  
  abundances_with_tourist$nturist7[i] <- nturist7
  abundances_with_tourist$nturist14[i] <- nturist14
  abundances_with_tourist$nturist21[i] <- nturist21
}
  

### Breeding sites - imbornales
matrizdistancias <- read_excel("matrizdistancias.xlsx", 
                               sheet = "matrizmodelos", col_types = c("text", 
                                                                      "numeric", "numeric", "numeric", 
                                                                      "numeric", "numeric", "numeric", 
                                                                      "numeric"))
abundancias_em <- left_join(abundances_with_tourist, matrizdistancias)
abundancias_em$trap_name <- as.factor(abundancias_em$trap_name)
abundancias_em$total_imb = rowSums(abundancias_em[ , c(77,78,79)], na.rm = T)

modelo_tot1 <- glmer.nb(nrperspecies ~ log(trapef) + scale(n_imbornales_50)  + 
                            scale(n_imbornales_100)+ scale(n_imbornales_150)+ scale(n_imbornales_agua)+
                            scale(n_imbornales_larva)+ scale(total_imb) + (1|trap_name), 
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                          data = abundancias_em) 
mfull_em <- dredge(modelo_tot1, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 2), extra= c(max.r))
NCMfull <- get.models(mfull_em, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDFfull <- model.sel(NCMfull)

#control ----------
TAULA_ATRIBUTS_GENERAL_MARIMURTRA <- read_excel("TAULA ATRIBUTS GENERAL MARIMURTRA.xlsx", 
                                                sheet = "datos_clean", col_types = c("numeric", 
                                                                                     "numeric", "date", "numeric", "text", 
                                                                                     "numeric", "numeric", "text", "text"))

TAULA_ATRIBUTS_GENERAL_MARIMURTRA <- TAULA_ATRIBUTS_GENERAL_MARIMURTRA %>% #creamos columna mes
  add_column(YEAR = NA)
TAULA_ATRIBUTS_GENERAL_MARIMURTRA$YEAR <- lubridate::year(ymd(TAULA_ATRIBUTS_GENERAL_MARIMURTRA$fecha))
datacontrol <- filter(TAULA_ATRIBUTS_GENERAL_MARIMURTRA, YEAR == 2021)  #separamos 2021 para el análisis
datacontrol <- datacontrol %>% 
  group_by(fecha) %>% 
  drop_na() %>% 
  summarise(volumen_Control = sum(volumen_tto))

#Tratamiento de datos
datacontrol <- datacontrol %>%
  rename("start_date" = `fecha`)
datacontrol <- datacontrol %>% #creamos columna mes
  add_column(week = NA)

datacontrol$week <- lubridate::week(ymd(datacontrol$start_date))
datacontrol$week <- as.factor(datacontrol$week)

abundancias_em_con <- left_join(abundancias_em, datacontrol, by = "week")
abundancias_em_con$volumen_Control[is.na(abundancias_em_con$volumen_Control)] <- 0

abundancias_aedes_bl$week <- as.character(abundancias_aedes_bl$week)


abundancias_aedes_bl$week <- as.numeric(abundancias_aedes_bl$week)
datacontrol$week <- as.character(datacontrol$week)
datacontrol$week <- as.numeric(datacontrol$week)
abundances_with_tourist$week <- as.character(abundances_with_tourist$week)
abundances_with_tourist$week <- as.numeric(abundances_with_tourist$week)

##calculate cumulative date control

all_dates<- as.data.frame(seq.Date(from = as_date("2021-04-12"), to = as_date("2021-11-22"), by = "day"))
all_dates <- all_dates %>%
  rename("date" = `seq.Date(from = as_date("2021-04-12"), to = as_date("2021-11-22"), by = "day")`)
abundancias_em_con$start_date.x <- as_date(abundancias_em_con$start_date.x)
abundancias_em_con$end_date <- as_date(abundancias_em_con$end_date)

c1 <- all_dates %>% 
  filter(between(date, as_date('2021-04-12'), as_date('2021-06-16'))) %>% 
  mutate(week_c =rep(1:100, each = 7, length.out=66))
c2 <- all_dates %>% 
  filter(between(date, as_date('2021-06-17'), as_date('2021-08-05'))) %>% 
  mutate(week_c =rep(1:100, each = 7, length.out=50))
c3 <- all_dates %>% 
  filter(between(date, as_date('2021-08-06'), as_date('2021-09-26'))) %>% 
  mutate(week_c =rep(1:100, each = 7, length.out=52))
c4 <- all_dates %>% 
  filter(between(date, as_date('2021-09-27'), as_date('2021-11-22'))) %>% 
  mutate(week_c =rep(1:100, each = 7, length.out=57))
new_table <- rbind(c1,c2,c3,c4)

abundancias_em_con <-  mclapply(1:nrow(abundancias_em_con), function(i){
  # Verbose
  if (i %% 50 == 0) {print(i)}
  print(i)
  
  x = abundancias_em_con[i,]
  this_start_date = x$start_date.x
  this_end_date = x$end_date
  
  # Exposure week
  y <- new_table %>%
    filter(date > this_start_date & date < this_end_date) 
  
  # if(length(levels(as.factor(y$week_c))) == 2) {
  #   x$weeksince = min(y$week_c)
  # }
  x$weeksince = max(y$week_c)
  
  return(x)
}
)

abundancias_em_con <- do.call(rbind, abundancias_em_con)
# p <- abundancias_em_con %>% 
#   select(start_date.x, end_date, start_date.y, weeksince)
# firstc <- p[1:22, ]
# firstc$weeksince <- firstc$weeksince -1
# second <- p[23:60,]
# second$weeksince <- second$weeksince -1
# third <- p[61:101,]
# third$weeksince <- third$weeksince -1
# forth <- p[102:142,]
# p<- rbind(firstc,second,third, forth)
# abundancias_em_con$weeksince <- p$weeksince
# abundancias_step_point <- abundancias_em_con %>% 
#   select(start_date.x, end_date, start_date.y, weeksince)
# test<- abundancias_step_point%>% 
#   mutate(control_f = case_when((weeksince >= 4) ~ 0,
#                             (weeksince < 4) ~ 1)) 
# abundancias_em_con$control_f <- test$control_f
# abundancias_em_con$control_f <- as.factor(abundancias_em_con$control_f)




#Modelo oficial, con todas las variables climáticas previas tras el dredge por fases
#pero añadiendo control, turistas y embornales...
modelo_ofi <- glmer.nb(nrperspecies ~ log(trapef) + scale(maxrh21) + scale(minrh21)
                       +scale(mintemp21) + scale(rain21) + scale(maxrh14) + scale(minrh14) +
                         scale(mintemp14) + scale(rain14) + scale(maxrh7) + scale(minrh7) +
                         scale(maxtemp7) + scale(rain7) + scale(max_relative_humidity) + 
                         scale(min_relative_humidity) +  scale(n_imbornales_100) + scale(total_imb) +
                         scale(n_imbornales_larva) +  scale(nrturist) + scale(n_imbornales_agua) + 
                         scale(poly(weeksince,2))
                       + scale(precipitation) + (1|trap_name), 
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                       data = abundancias_em_con) 


summary(modelo_ofi)
mfull_em <- dredge(modelo_ofi, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 7), extra= c(max.r))
NCMfull <- get.models(mfull_em, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDFfull <- model.sel(NCMfull)
selection_full <- as.data.frame(NCMDFfull)
sw(NCMDFfull)
writexl::write_xlsx(selection_full, path = "C:\\Users\\lblan\\OneDrive\\Escritorio\\CEAB\\2022\\Abundance_study\\selection_full.xlsx")
modelo_ofi <- glmer.nb(nrperspecies ~ log(trapef) +  scale(max_relative_humidity) + scale(mintemp21) +
                         scale(minrh7) + scale(rain21) +
                         scale(nrturist) + scale(n_imbornales_agua) + 
                         scale(poly(weeksince,2))+
                         (1|trap_name),
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                       data = abundancias_em_con)

summary(modelo_ofi)
check_model(modelo_ofi)
newdata <- data.frame(minrh7 = abundancias_em_con$minrh7,
                   mintemp21 = abundancias_em_con$mintemp21,
                   max_relative_humidity = abundancias_em_con$max_relative_humidity,
                   rain21= abundancias_em_con$rain21,
                   nrturist = abundancias_em_con$nrturist,
                   trapef = 7,
                   n_imbornales_agua= 18,
                   weeksince= abundancias_em_con$weeksince,
                   trap_name= abundancias_em_con$trap_name)

predict(modelo_ofi, newdata, type="response")
predictiondata <- as.data.frame(predict(modelo_ofi, newdata, type="response"))
predictiondata$predi <- predictiondata$`predict(modelo_ofi, newdata, type = "response")`
abundances <- cbind(abundancias_em_con, predictiondata)
abundances$week <- as.character(abundances$week)
abundances$week <- as.numeric(abundances$week)
plot_abu <- abundances %>% 
  group_by(week) %>% 
  summarise(sumpredi = sum(predi)) %>% 
  ggplot() +
  geom_line(aes(x=week , y= sumpredi),col="red", size = 1.1) +
  theme_bw()

preDIC <- abundances %>% 
  group_by(start_date.x, month, end_date, trapef) %>% 
  summarise(summos = sum(predi))
preDIC$cumsumos = cumsum(preDIC$summos)

writexl::write_xlsx(abundancias_filter, path = "C:\\Users\\lblan\\OneDrive\\Escritorio\\CEAB\\2022\\Abundance_study\\abundanciasdf.xlsx")


#prediction con todo 2021 para Jesús
clima_malgrat <- read_csv("Malgrat-de-mar-2020-2022.csv", 
                          col_types = cols(service = col_skip(), 
                                           station_id = col_skip(), station_name = col_skip(), 
                                           station_province = col_skip(), altitude = col_skip(), 
                                           mean_wind_direction = col_skip(), 
                                           FW = col_skip(), FT = col_skip(), 
                                           mwi = col_skip()))
clima_malgrat <- clima_malgrat %>% #creamos columna mes
  add_column(month = NA)
clima_malgrat$month <- lubridate::month(ymd(clima_malgrat$timestamp))

clima_malgrat <- clima_malgrat %>% #creamos columna semana
  add_column(week = NA)
clima_malgrat$week <- lubridate::week(ymd(clima_malgrat$timestamp))

clima_malgrat <- clima_malgrat %>% #creamos columna semana
  add_column(YEAR = NA)
clima_malgrat$YEAR <- lubridate::year(ymd(clima_malgrat$timestamp))
clima_malgrat <- filter(clima_malgrat, !YEAR==2022)
clima_malgrat <- filter(clima_malgrat, !YEAR==2020)

clima_malgrat$doy <- yday(clima_malgrat$timestamp) #calculamos doy
clima_malgrat$month <- as.factor(clima_malgrat$month)
clima_malgrat$week <- as.factor(clima_malgrat$week)
clima_malgrat$YEAR <- as.factor(clima_malgrat$YEAR)
#photoperiod
clima_malgrat$photoperiod <- photoperiod(1:365, 41.676854)
#medias y acumuladas de ventanas temporales previas 7/14/21 días 
clima_malgrat <- clima_malgrat %>% 
  mutate(mintemp7 = rollmean(min_temperature, k = 7, align = "right", fill = NA),
         mintemp14 = rollmean(min_temperature, k = 14, align = "right", fill = NA),
         mintemp21 = rollmean(min_temperature, k = 21, align = "right", fill = NA),
         minrh7 = rollmean(min_relative_humidity, k = 7, align = "right", fill = NA),
         minrh14 = rollmean(min_relative_humidity, k = 14, align = "right", fill = NA),
         minrh21 = rollmean(min_relative_humidity, k = 21, align = "right", fill = NA),
         maxrh7 = rollmean(max_relative_humidity, k = 7, align = "right", fill = NA),
         maxrh14 = rollmean(max_relative_humidity, k = 14, align = "right", fill = NA),
         maxrh21 = rollmean(max_relative_humidity, k = 21, align = "right", fill = NA),
         rain7 = rollapplyr(precipitation, 7, FUN=sum, align = "right", fill = NA),
         rain14 = rollapplyr(precipitation, 14, FUN=sum, align = "right", fill = NA),
         rain21 = rollapplyr(precipitation, 21, FUN=sum, align = "right", fill = NA))

clima_malgrat <- clima_malgrat %>% 
  rename('start_date' = "timestamp")
datatourist$start_date <- as_date(datatourist$start_date)

pred_data <- left_join(clima_malgrat, datatourist, by="start_date")
pred_data$trap_name <- "A_SP_BL_12"



sequence <- rep(1:4, each = 7, length.out = 365)


newdata <- data.frame(minrh7 = pred_data$minrh7,
                      mintemp21 = pred_data$mintemp21,
                      max_relative_humidity = pred_data$max_relative_humidity,
                      rain21= pred_data$rain21,
                      nrturist = pred_data$nrturist,
                      trapef = 7,
                      n_imbornales_agua= 18,
                      weeksince= sequence,
                      trap_name= pred_data$trap_name)

predict(modelo_ofi, newdata, type="response", allow.new.levels = T)
predictiondata <- as.data.frame(predict(modelo_ofi, newdata, type="response", allow.new.levels = T))
predictiondata$predi <- predictiondata$`predict(modelo_ofi, newdata, type = "response", allow.new.levels = T)`
abundances <- cbind(pred_data, predictiondata)

plot_abu <- abundances %>% 
  ggplot() +
  geom_line(aes(x=start_date , y= predi),col="red", size = 1.1) +
  theme_bw()

preDIC <- abundances %>% 
  group_by(start_date) %>% 
  select(start_date, predi)
preDIC$cumsumos = cumsum(preDIC$predi)

##PREDICCIÓN JARDÍN REAL TIME -----------------------------------

# DF vacío
marimurtra <-data.frame()

# Para 2022
url <- "https://webapi.www.ecowitt.net/uploads/60907/EasyWeather-WIFIBC89%28202201010000-202212312359%29.xlsx"

mm <- openxlsx::read.xlsx(url, startRow = 2)

marimurtra <- rbind(marimurtra, mm)
marimurtra$X1 <- as.Date(marimurtra$X1)

marimurtra <- marimurtra[-1:-140,]

# Para 2023

url <- "https://webapi.www.ecowitt.net/uploads/60907/EasyWeather-WIFIBC89%28202301010000-202305311121%29.xlsx"
mm <- openxlsx::read.xlsx(url, startRow = 2)
marimurtra <- rbind(marimurtra, mm)
marimurtra <- marimurtra[,-27:-29]

marimurtra <- marimurtra %>% 
  rename("max_temperature" = `Temperature.High(℃)`)
marimurtra <- marimurtra %>% 
  rename("min_temperature" =  `Temperature.Low(℃)`)
marimurtra <- marimurtra %>% 
  rename("mean_temperature" = `Temperature(℃)` )
marimurtra <- marimurtra %>% 
  rename("max_relative_humidity" = `Humidity.High(%)`)
marimurtra <- marimurtra %>% 
  rename("min_relative_humidity" = `Humidity.Low(%)`  )
marimurtra <- marimurtra %>% 
  rename("mean_relative_humidity" = `Humidity(%)`)
marimurtra <- marimurtra %>% 
  rename("precipitation" = 'Daily(mm)')
marimurtra <- marimurtra %>% 
  rename(DATA = "X1")
class(marimurtra$DATA)


marimurtra <- marimurtra %>% 
  mutate(mintemp7 = rollmean(min_temperature, k = 7, align = "right", fill = NA),
         mintemp14 = rollmean(min_temperature, k = 14, align = "right", fill = NA),
         mintemp21 = rollmean(min_temperature, k = 21, align = "right", fill = NA),
         meantemp7 = rollmean(mean_temperature, k = 7, align = "right", fill = NA),
         meantemp14 = rollmean(mean_temperature, k = 14, align = "right", fill = NA),
         meantemp21 = rollmean(mean_temperature, k = 21, align = "right", fill = NA),
         maxtemp7 = rollmean(max_temperature, k = 7, align = "right", fill = NA),
         maxtemp14 = rollmean(max_temperature, k = 14, align = "right", fill = NA),
         maxtemp21 = rollmean(max_temperature, k=21, align = "right", fill = NA),
         minrh7 = rollmean(min_relative_humidity, k = 7, align = "right", fill = NA),
         minrh14 = rollmean(min_relative_humidity, k = 14, align = "right", fill = NA),
         minrh21 = rollmean(min_relative_humidity, k = 21, align = "right", fill = NA),
         meanrh7 = rollmean(mean_relative_humidity, k = 7, align = "right", fill = NA),
         meanrh14 = rollmean(mean_relative_humidity, k = 14, align = "right", fill = NA),
         meanrh21 = rollmean(mean_relative_humidity, k = 21, align = "right", fill = NA),
         maxrh7 = rollmean(max_relative_humidity, k = 7, align = "right", fill = NA),
         maxrh14 = rollmean(max_relative_humidity, k = 14, align = "right", fill = NA),
         maxrh21 = rollmean(max_relative_humidity, k = 21, align = "right", fill = NA),
         rain7 = rollapplyr(precipitation, 7, FUN=sum, align = "right", fill = NA),
         rain14 = rollapplyr(precipitation, 14, FUN=sum, align = "right", fill = NA),
         rain21 = rollapplyr(precipitation, 21, FUN=sum, align = "right", fill = NA)
  )

new_wth <- function(all_dates){
  
  wth_data <-data.frame()
  
  for(date in all_dates){ 
    
    x0 <- as_date(date) - 7
    x1 <- as_date(date)

    int <- seq(x0, x1, by= "day")
    
    max_temperature <- mean(marimurtra$max_temperature[marimurtra$DATA %in% int])
    min_temperature <- mean(marimurtra$min_temperature[marimurtra$DATA %in% int])
    mean_relative_humidity <- mean(marimurtra$mean_relative_humidity[marimurtra$DATA %in% int])
    min_relative_humidity <- mean(marimurtra$min_relative_humidity[marimurtra$DATA %in% int])
    max_relative_humidity <- mean(marimurtra$max_relative_humidity[marimurtra$DATA %in% int])
    precipitation <- sum(marimurtra$precipitation[marimurtra$DATA %in% int])
    rain7 <- marimurtra$rain7[marimurtra$DATA %in% x0]
    rain14 <- marimurtra$rain14[marimurtra$DATA %in% x0]
    rain21 <-marimurtra$rain21[marimurtra$DATA %in% x0]
    mintemp7 <- marimurtra$mintemp7[marimurtra$DATA %in% x0]
    mintemp14 <- marimurtra$mintemp14[marimurtra$DATA %in% x0]
    mintemp21 <- marimurtra$mintemp21[marimurtra$DATA %in% x0]
    maxrh7 <- marimurtra$maxrh7[marimurtra$DATA %in% x0]
    maxrh14 <- marimurtra$maxrh14[marimurtra$DATA %in% x0]
    maxrh21 <- marimurtra$maxrh21[marimurtra$DATA %in% x0]
    minrh7 <- marimurtra$minrh7[marimurtra$DATA %in% x0]
    minrh14 <- marimurtra$minrh14[marimurtra$DATA %in% x0]
    minrh21 <- marimurtra$minrh21[marimurtra$DATA %in% x0]
    
    wth <- data.frame(
      DATA = x1,
      max_temperature = max_temperature,
      min_temperature = min_temperature,
      mean_relative_humidity = mean_relative_humidity,
      min_relative_humidity = min_relative_humidity,
      max_relative_humidity = max_relative_humidity,
      precipitation = precipitation,
      rain7 = rain7,
      rain14 = rain14,
      rain21 = rain21,
      mintemp7 = mintemp7,
      mintemp14 = mintemp14,
      mintemp21 = mintemp21,
      maxrh7 = maxrh7,
      maxrh14 = maxrh14,
      maxrh21 = maxrh21,
      minrh7 = minrh7,
      minrh14 = minrh14,
      minrh21 = minrh21
    )
    
    wth_data <- rbind(wth_data, wth)} 
  
  return(wth_data)
}



all_dates <- seq.Date(from = as_date("2022-08-01"), to = as_date("2023-05-18"), by = "day")

wth_pred <- new_wth(all_dates)
TOURIST23 <- read_excel("TOURIST23.xlsx", 
                        col_types = c("date", "numeric"))
all_dates <- data.frame(Data = seq.Date(from = as_date("2023-01-01"), to = as_date("2023-05-18"), by = "day"))
TOURIST23$Data <- as_date(TOURIST23$Data)
Datatourist23 <-  merge(TOURIST23, all_dates,all = TRUE)
Datatourist23 <- Datatourist23 %>% 
  rename(start_date = "Data")
Datatourist23 <- Datatourist23 %>% 
  rename(nrturist = "Visitants")
Tourist22 <- Tourist22[-1:-212,]
tourist_2223 <- rbind(Tourist22, Datatourist23)

NEWd <- data.frame(minrh7 = wth_pred$minrh7,
                   mintemp21 = wth_pred$mintemp21,
                   maxrh14 = wth_pred$maxrh14,
                   precipitation = wth_pred$precipitation,
                   nrturist = tourist_2223$nrturist,
                   n_imbornales_agua = 18,
                   trapef = 7,
                   weeksince1_2=5,
                   trap_name=  "A_SP_BL_4")


