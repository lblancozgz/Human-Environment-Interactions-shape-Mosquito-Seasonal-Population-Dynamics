#Uploading packages
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
library(effects)
library(corrplot)
library(RColorBrewer)

#set directory and upload dataframe
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


#Filtering and renaming
abundancias_aedes <- filter(abundancias, species == "Aedes albopictus") #we want only Ae. albopictus

#Renaming  variables
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


#removing BG-TRAPS not located in the Jardí Botánic
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_1")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_2")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_10")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_11")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_12")

#Temporal variables:

#Month
abundancias_aedes_bl <- abundancias_aedes_bl %>%
  add_column(month = NA)
abundancias_aedes_bl$month <- lubridate::month(ymd(abundancias_aedes_bl$start_date))
abundancias_aedes_bl$month <- as.factor(abundancias_aedes_bl$month) #is a factor


#week
abundancias_aedes_bl <- abundancias_aedes_bl %>% 
  add_column(week = NA)
abundancias_aedes_bl$week <- lubridate::week(ymd(abundancias_aedes_bl$start_date))
abundancias_aedes_bl$week <- as.factor(abundancias_aedes_bl$week) #is a factor

#Year
abundancias_aedes_bl <- abundancias_aedes_bl %>% 
  add_column(YEAR = NA)
abundancias_aedes_bl$YEAR <- lubridate::year(ymd(abundancias_aedes_bl$start_date))
abundancias_aedes_bl$doy <- yday(abundancias_aedes_bl$start_date) #days of the year also
abundancias_aedes_bl$YEAR <- as.factor(abundancias_aedes_bl$YEAR) #is a factor

# Data visualization: abundance per trap
ggplot() +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date, y= nrperspecies), data = filter(abundancias_aedes_bl, trap_name == "A_SP_BL_9")) 

abundancias_aedes_bl$trap_name = factor(abundancias_aedes_bl$trap_name) # convert trap to factor
str(abundancias_aedes_bl$trap_name)

#normality
ggqqplot(abundancias_aedes_bl$nrperspecies)
ggdensity(abundancias_aedes_bl$nrperspecies, fill = "lightgray") #we can't assume normality
#Shapiro-test for grouped data to check normality
abundancias_aedes_bl %>%
  group_by(trap_name) %>%
  shapiro_test(nrperspecies) #not normal

# Uploading climate data Malgrat de Mar station
clima_malgrat <- read_csv("Malgrat-de-mar-2020-2022.csv", 
                          col_types = cols(service = col_skip(), 
                                           station_id = col_skip(), station_name = col_skip(), 
                                           station_province = col_skip(), altitude = col_skip(), 
                                           mean_wind_direction = col_skip(), 
                                           FW = col_skip(), FT = col_skip(), 
                                           mwi = col_skip()))

#Creating temporal columns as before

#Month
clima_malgrat <- clima_malgrat %>% 
  add_column(month = NA)
clima_malgrat$month <- lubridate::month(ymd(clima_malgrat$timestamp))
clima_malgrat$month <- as.factor(clima_malgrat$month)

#Week
clima_malgrat <- clima_malgrat %>% 
  add_column(week = NA)
clima_malgrat$week <- lubridate::week(ymd(clima_malgrat$timestamp))
clima_malgrat$week <- as.factor(clima_malgrat$week)

#Year
clima_malgrat <- clima_malgrat %>% 
  add_column(YEAR = NA)
clima_malgrat$YEAR <- lubridate::year(ymd(clima_malgrat$timestamp))
clima_malgrat <- filter(clima_malgrat, !YEAR==2022) #we want only 2021
clima_malgrat <- filter(clima_malgrat, !YEAR==2020)
clima_malgrat$doy <- yday(clima_malgrat$timestamp) #days of the year also

clima_malgrat$YEAR <- as.factor(clima_malgrat$YEAR)
clima_malgrat <- clima_malgrat %>% #rename variable
  rename(start_date = "timestamp")
#Visualize climate variables
climatemp <- clima_malgrat %>%
  ggplot() +
  geom_line(aes(x = start_date, y = max_temperature), col= "red")+
  geom_line(aes(x = start_date, y = min_temperature), col= "darkblue")+
  geom_line(aes(x = start_date, y = mean_temperature), col="darkgreen")+
  ylab("Temperatures") +
  theme_bw()

climarh <- clima_malgrat %>%
  ggplot() +
  geom_line(aes(x = start_date, y = max_relative_humidity), col= "red")+
  geom_line(aes(x = start_date, y = min_relative_humidity), col= "darkblue")+
  geom_line(aes(x = start_date, y = mean_relative_humidity), col="darkgreen")+
  ylab("RH") +
  theme_bw()

climalluvia <- clima_malgrat %>%
  group_by(week, month, YEAR) %>% 
  summarise(rainfall = sum(precipitation)) %>% 
  ggplot() +
  geom_col(aes(x = month, y = rainfall))+
  ylab("Rain") +
  theme_bw()



# Calculate MWI
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
#Calculate  photoperiod
clima_malgrat$photoperiod <- photoperiod(1:365, 41.676854)


#Calculate average and accumulated variables in the 1-2-3 weeks prior to sampling 
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
         rain21 = rollapplyr(precipitation, 21, FUN=sum, align = "right", fill = NA),
         rain28 = rollapplyr(precipitation, 21, FUN=sum, align = "right", fill = NA))



#Check some calculations
range(clima_malgrat$mwime, na.rm = T) #it mus be between 0 y 1
ggplot(clima_malgrat, aes(x=start_date, y=mwi7)) + geom_line()
ggplot(clima_malgrat, aes(x=start_date, y=photoperiod)) + geom_line()
ggplot(clima_malgrat, aes(x=start_date, y=pho21)) + geom_line()

#Dates from the two main dataframes in the same format (clima + abundances)
clima_malgrat$start_date <- as.POSIXct(clima_malgrat$start_date,format = '%Y-%m-%d')
abundancias_aedes_bl$start_date <- as.POSIXct(abundancias_aedes_bl$start_date,format = '%Y-%m-%d')

#Join two dataframes into one
abundancias_jardin_tot <- full_join(abundancias_aedes_bl, clima_malgrat, 
                                    by = c("start_date", "YEAR", "week", "month", "doy"), all= TRUE) 

#Important to calculate the averages of each climatic variable in the intervals of
#the trapping effort (sampling period)

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
  rain28 <-clima_malgrat$rain28[clima_malgrat$start_date %in% x0]
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

#keep the dataframe with the data containing abundances
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




options(na.action = "na.fail") #necessary to run dredge() function
which(is.na(abundances_with_gdd$nrperspecies)) #remove all NAs that can affect negatively to the df

#### MODELS ###
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


#First full model, with all the climatic variables from the sampling week window
m0 <-glmer.nb(nrperspecies ~ scale(max_temperature) +
                     scale(min_temperature) +
                     scale(mean_temperature) +
                     scale(max_relative_humidity) + scale(mean_relative_humidity) + scale(min_relative_humidity) +
                     scale(precipitation) + scale(mwime) + log(trapef) + scale(l0gdd) +
                     (1|trap_name),
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                   data = abundances_with_gdd)
summary(m0)
vif(m0)

#DREDGE FOR MODEL SELECTION
ms3 <- dredge(m0, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))

NCM3 <- get.models(ms3, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)

NCMDF <- model.sel(NCM3)

#In the sampling week window, we got log(trapef), max_RH, min_RH, minTemp, precipitation


#Model with the sampling week window after dredge
m0_dredge <- glmer.nb(nrperspecies ~ scale(max_relative_humidity) + scale(min_relative_humidity) + 
                        scale(min_temperature) + scale(precipitation) + log(trapef) + (1|trap_name),
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                      data = abundances_with_gdd)
summary(m0_dredge)
testZeroInflation(m0_dredge)
testResiduals(m0_dredge)
testDispersion(m0_dredge)


# 1-week prior sampling week, we do the same
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
#We got log(trapef), maxrh7, maxtemp7, minrh7, rain7

# 1-week prior sampling week, after dredge
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

#2-week prior to sampling week, we do the same

m2 <-glmer.nb(nrperspecies ~ scale(maxtemp14) + scale(mintemp14) + scale(meantemp14) +
                scale(maxrh14) + scale(meanrh14) + scale(minrh14) +
                scale(rain14) + scale(mwi14) + log(trapef) + scale(l14gdd) +
                (1|trap_name), 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
              data = abundances_with_gdd)

ms5 <- dredge(m2, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))
NCM5 <- get.models(ms5, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDF5 <- model.sel(NCM5)
#MaxRH 14, MinRH 14, Mintemp14, rain14

#2-week prior to sampling week, after dredge

m2_dredge <-glmer.nb(nrperspecies ~  log(trapef) + scale(maxrh14) + scale(minrh14) + scale(mintemp14) + scale(rain14) +
                       (1|trap_name), 
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                     data = abundances_with_gdd)
summary(m2_dredge)
car::vif(m2_dredge)
testZeroInflation(m2_dredge)
testResiduals(m2_dredge)
testDispersion(m2_dredge)

#3-week prior to sampling week
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
# log(trapef, maxrh21, minrh21, mintemp21, rain21)

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

#Full model with all previous chosen variables (excluding water drains, tourist... by the moment)
fullmodel <- glmer.nb(nrperspecies ~ scale(maxrh21) + scale(minrh21)
                      +scale(mintemp21) + scale(rain21) +
                        scale(mintemp14) + scale(rain14) + scale(minrh14) + scale(maxrh14) +
                        scale(maxrh7) + scale(minrh7) +
                        scale(maxtemp7) + scale(rain7) + scale(max_relative_humidity) + scale(min_relative_humidity) +
                        scale(min_temperature) + 
                        (1|trap_name), 
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                      data = abundances_with_gdd) 
mfull <- dredge(fullmodel, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 3), extra= c(max.r))
NCMfull <- get.models(mfull, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDFfull <- model.sel(NCMfull)
sw(NCMfull) #variables that most explain are mintemp-21, and maxrh21. The others are <0.6 in weight
#We will remove now the other RH variables and we will see again to be sure that not other variables are important here

fullmodel <- glmer.nb(nrperspecies ~ scale(maxrh21) + 
                      +scale(mintemp21) + scale(rain21) +
                        scale(rain14) + 
                         scale(rain7) + 
                        (1|trap_name), 
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                      data = abundances_with_gdd) 
mfull <- dredge(fullmodel, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 3), extra= c(max.r))
NCMfull <- get.models(mfull, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDFfull <- model.sel(NCMfull)
sw(NCMfull)

#it makes more sense that rain21 is affecting more the number of posterior adults if
#we take into account the biology and life cycle of the mosquito
full_dredge <- glmer.nb(nrperspecies ~ scale(maxrh21) + scale(mintemp21) +
                          scale(rain21) + (1|trap_name), 
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                        data = abundances_with_gdd) 
summary(full_dredge)
car::vif(full_dredge)
testZeroInflation(full_dredge)
testResiduals(full_dredge)
testDispersion(full_dredge)
testUniformity(full_dredge)

# We add tourist data

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

#Calculating the sum of tourist in the 1,2,3 weeks prior to sampling
datatourist <- datatourist %>% 
  mutate(nturist7 = rollapplyr(nrturist, 7, FUN=sum, align = "right", fill = NA),
         nturist14 = rollapplyr(nrturist, 14, FUN=sum, align = "right", fill = NA),
         nturist21 = rollapplyr(nrturist, 21, FUN=sum, align = "right", fill = NA))

datatourist$start_date <- as.POSIXct(datatourist$start_date,format = '%Y-%m-%d')
abundances_with_gdd$start_date <- as.POSIXct(abundances_with_gdd$start_date,format = '%Y-%m-%d')

ggplot(datatourist, aes(x=start_date, y=nrturist)) + geom_line() #plot tourist
#We join data with tourist, be careful to be linked with the sampling weeks
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


### Water drains
matrizdistancias <- read_excel("matrizdistancias.xlsx", 
                               sheet = "matrizmodelos", col_types = c("text", 
                                                                      "numeric", "numeric", "numeric", 
                                                                      "numeric", "numeric", "numeric", 
                                                                      "numeric"))
abundancias_em <- left_join(abundances_with_tourist, matrizdistancias)
abundancias_em$trap_name <- as.factor(abundancias_em$trap_name)
abundancias_em$total_imb = rowSums(abundancias_em[ , c(78,79,80)], na.rm = T)


#control ----------
TAULA_ATRIBUTS_GENERAL_MARIMURTRA <- read_excel("TAULA ATRIBUTS GENERAL MARIMURTRA.xlsx", 
                                                sheet = "datos_clean", col_types = c("numeric", 
                                                                                     "numeric", "date", "numeric", "text", 
                                                                                     "numeric", "numeric", "text", "text"))

#Year column
TAULA_ATRIBUTS_GENERAL_MARIMURTRA <- TAULA_ATRIBUTS_GENERAL_MARIMURTRA %>% 
  add_column(YEAR = NA)
TAULA_ATRIBUTS_GENERAL_MARIMURTRA$YEAR <- lubridate::year(ymd(TAULA_ATRIBUTS_GENERAL_MARIMURTRA$fecha))

#We filter for 2021 
datacontrol <- filter(TAULA_ATRIBUTS_GENERAL_MARIMURTRA, YEAR == 2021)  
datacontrol <- datacontrol %>% 
  group_by(fecha) %>% 
  drop_na() %>% 
  summarise(volumen_Control = sum(volumen_tto))

#Tratamiento de datos
datacontrol <- datacontrol %>%  #renaming
  rename("start_date" = `fecha`)
datacontrol <- datacontrol %>% #week column
  add_column(week = NA)
datacontrol$week <- lubridate::week(ymd(datacontrol$start_date))
datacontrol$week <- as.factor(datacontrol$week) #week as a factor

abundancias_em_con <- left_join(abundancias_em, datacontrol, by = "week") #joining with our dataset
abundancias_em_con$volumen_Control[is.na(abundancias_em_con$volumen_Control)] <- 0

##calculate WEEK SINCE CONTROL
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
p <- abundancias_em_con %>%
  select(start_date.x, end_date, start_date.y, weeksince) #Check calculations


#We take again the PREVIOUS FULL MODEL WITH ALL VARIABLES CHOSEN in the ESCALATED dredges
#and we add water drains, weeks since control treatments...
modelo_ofi <- glmer.nb(nrperspecies ~ log(trapef) + scale(maxrh21) + scale(minrh21)
                       +scale(mintemp21) + scale(rain21) + scale(maxrh14) + scale(minrh14) +
                         scale(mintemp14) + scale(rain14) + scale(maxrh7) + scale(minrh7) +
                         scale(maxtemp7) + scale(rain7) + scale(max_relative_humidity) + 
                         scale(min_relative_humidity) + 
                         scale(n_imbornales_agua) + 
                         scale(weeksince)
                       + scale(precipitation) + (1|trap_name), 
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                       data = abundancias_em_con) 


summary(modelo_ofi)
mfull_em <- dredge(modelo_ofi, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0,8), extra= c(max.r))
NCMfull <- get.models(mfull_em, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDFfull <- model.sel(NCMfull)

selection_full <- as.data.frame(NCMDFfull)
sw(NCMDFfull) # maximum RH is the humidity variable which appears the most, so we will keep it and we will repeat again the process
#removing humidities, avoiding like this models with same variables explaining also the same

#mintemp21, maxrh, precipitation

modelo_ofi2 <- glmer.nb(nrperspecies ~ log(trapef) +
                       +scale(mintemp21) + scale(rain21) + 
                         scale(mintemp14) + scale(rain14) + 
                         scale(maxtemp7) + scale(rain7) + scale(max_relative_humidity) + 
                         scale(n_imbornales_agua) + 
                         scale(weeksince)
                       + scale(precipitation) + (1|trap_name), 
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                       data = abundancias_em_con) 


summary(modelo_ofi2)

mfull_em2 <- dredge(modelo_ofi2, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0,8), extra= c(max.r))
NCMfull <- get.models(mfull_em2, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDFfull2 <- model.sel(NCMfull)
sw(NCMDFfull2) #we obtain in order of weight: MaxRH > weeksince > mintemp21 >n_imbornales_agua > rain21
selection_full <- as.data.frame(NCMDFfull2)

writexl::write_xlsx(selection_full, path = "C:\\Users\\lblan\\OneDrive\\Escritorio\\CEAB\\2022\\Abundance_study\\selection_full.xlsx")

#Model after dredge
modelo_ofi <- glmer.nb(nrperspecies ~ log(trapef) + scale(mintemp21) + scale(rain21) +
                         scale(max_relative_humidity) +  scale(n_imbornales_agua) + 
                         scale(nrturist) +  scale(weeksince)
                       + (1|trap_name), 
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                       data = abundancias_em_con)

summary(modelo_ofi)
car::Anova(modelo_ofi)


#Our database is HUGE, so we will select the important variables and create
#a more little dataframe
abundancias_df <- abundancias_em_con %>% 
  select(trapef, max_relative_humidity, minrh7, mintemp21, weeksince, n_imbornales_agua, 
         nrturist, start_date.x, end_date, week, nrperspecies, trap_name, rain21,
         rain14, rain7, rain28, precipitation,n_imbornales_agua, total_imb, volumen_Control)

#We will also remove these 2 outliers
abundancias_df <- filter(abundancias_df, !weeksince == "9") 
abundancias_df <- filter(abundancias_df, !weeksince == "10") 


#Is there any interaction between rain and the weeks passed since control treatments??
modelo_ofi <- glmer.nb(nrperspecies ~  scale(mintemp21) + scale(weeksince)*scale(rain21) +
                         scale(max_relative_humidity) +  scale(n_imbornales_agua) + scale(weeksince) +
                         scale(rain21) +
                         scale(nrturist) +
                         + (1|trap_name), 
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                       data = abundancias_df)

summary(modelo_ofi)
car::Anova(modelo_ofi, type = 3)
testResiduals(modelo_ofi)
testQuantiles(modelo_ofi)
testDispersion(modelo_ofi)
testZeroInflation(modelo_ofi)
car::vif(modelo_ofi)
testQuantiles(modelo_ofi)
performance::r2_nakagawa(modelo_ofi)

#PREDICTION
newdata <- data.frame(mintemp21 = abundancias_df$mintemp21,
                      max_relative_humidity = abundancias_df$max_relative_humidity,
                      rain21 = abundancias_df$rain21, 
                      n_imbornales_agua =43,
                      nrturist = abundancias_df$nrturist,
                      weeksince= abundancias_df$weeksince,
                      trap_name= "A_SP_BL_12")

predict(modelo_ofi, newdata, type="response", allow.new.levels = T)
predictiondata <- as.data.frame(predict(modelo_ofi, newdata, type="response",allow.new.levels = T))
predictiondata$predi <- predictiondata$`predict(modelo_ofi, newdata, type = "response", allow.new.levels = T)`
abundances <- cbind(abundancias_df, predictiondata)
abundances$week <- as.character(abundances$week)
abundances$week <- as.numeric(abundances$week)
abundances$start_date.x <- as_date(abundances$start_date.x)

Sys.setlocale(category = "LC_ALL", locale = "EN")
plot_abu_tur <- abundances %>% 
  ggplot() +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_9")) +
  geom_line(aes(x=start_date.x , y= predi),col="darkblue", size = 1.9) +
  scale_x_date("Sampling season", date_breaks = "1 month", date_labels = "%b") +
  ylab("Abundance (average)") + ylim(0,100) +
  # scale_y_continuous(breaks = seq(from = 0, to = 700, by = 20))+
  theme_classic(base_size = 21)
plotb <- plot_abu_tur #plot b in panel Figure 3


#Plot a in panel Figure 3
plota <- ggplot(abundancias_em_con, aes(x=start_date.x)) +
  geom_col( aes(y=nrperspecies,),fill="#DE9826", alpha= 0.7,width =8) + 
  geom_line(aes(y= rain21), color="#3C4C8F", linewidth=1.3) +  labs() + 
  geom_line(aes(y = min_temperature*10),  col="#670012",size = 1.5) +  # Scale y2 by 10 for better visualization
  scale_y_continuous(
    name = "Abundance (total counts) \n Acc, Rainfall 3 weeks (mm)",
    limits = c(0, 300),
    sec.axis = sec_axis(~./10, name = "Min. temperature (Celsius)")
  ) +
  scale_x_date("Sampling season", date_breaks = "1 month", date_labels = "%b") +
  theme_classic(base_size = 21) 


#Plot C in panel Figure 3
plotc<- ggplot(data=abundancias_df, aes(x=factor(weeksince), y=nrperspecies)) + 
  geom_boxplot(col="#652770", linewidth = 0.8) +
  geom_jitter(col="#DE9826", alpha=0.3, size=4.8) +
  geom_smooth(aes(x=factor(weeksince), y=nrperspecies))+
  xlab("Weeks since treatments") + ylab("Average abundance") + 
  theme_classic(base_size = 21)



#effect plot

interactions_weeksince <- plot(allEffects(modelo_ofi, residuals=F)) # With this argument I can plot the partial residuals
plot(effect("scale(weeksince)*scale(rain21)", modelo_ofi))
interacdata <- effect("scale(weeksince)*scale(rain21)",modelo_ofi)
effect_data <- as.data.frame(interacdata)
effect_data <- effect_data %>%
  rename("Number_mosquitoes" = `fit`)
#Heat plot - Figure 4
heat<- ggplot(effect_data, aes(factor(rain21), factor(weeksince))) +
  geom_tile(aes(fill = Number_mosquitoes)) +
  ylab("Weeks since treatments") + xlab("Acc. rainfall 3 weeks")+
  geom_text(aes(label = round(Number_mosquitoes, 1)), size = 9) +
  scale_fill_gradient(low = "white", high = "#652770") +
  theme_classic(base_size = 25)
heat <- heat + theme(legend.title = element_blank())  
heat


#Plot D figure 3
effect_data <- as.data.frame(interacdata)
effect_data$rain21 <- as.factor(effect_data$rain21)

rain_names <- list(
  '0.2'="Rain 0.2mm",
  '54'="Rain 54 mm",
  '110'="Rain 110 mm",
  '160'="Rain 160 mm",
  '210'="Rain 210 mm"
  
)
rain_labeller <- function(variable,value){
  return(rain_names[value])
}

interac<- ggplot(effect_data, aes(x = weeksince, y = fit, color = rain21)) +
  geom_line(linewidth=2, color="#652770") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#652770", color="gray87") +
  facet_wrap(~rain21, nrow = 1, labeller = rain_labeller) +
  labs(x = "Weeks since treatments", y = "Predicted abundance \n of mosquitoes") +
  theme_bw(base_size = 21)
interac <- interac + theme(legend.position = "none",
                           panel.spacing = unit(0, "cm"),
                           strip.text.x = element_text(size = 14),
                           strip.background=element_rect(fill="white"))
interac
plotd <- interac 

#SUPP MATERIAL

#Figure Panel plot

# newdata <- data.frame(mintemp21 = abundancias_df$mintemp21, 
#                       max_relative_humidity = abundancias_df$max_relative_humidity,
#                       rain21 = abundancias_df$rain21, WE FIX acc.RAIN to 0.2,54,110 and 210
#                       n_imbornales_agua =43,
#                       nrturist = abundancias_df$nrturist,
#                       weeksince= abundancias_df$weeksince, WE FIX WEEKSINCE in 1 and 6 to compare
#                       trap_name= "A_SP_BL_12")
# dataplot <- as.data.frame(predictiondata$predi02_1)
# dataplot$predi02_6 <- predictiondata$predi02_6
# dataplot$predi54_1 <- predictiondata$predi54_1
# dataplot$predi54_6 <- predictiondata$predi54_6
# dataplot$predi110_1 <- predictiondata$predi110_1
# dataplot$predi110_6 <- predictiondata$predi110_6
# dataplot$predi210_1 <- predictiondata$predi210_1
# dataplot$predi210_6 <- predictiondata$predi210_6
# dataplot$start_date.x <- abundances$start_date.x

writexl::write_xlsx(dataplot, path = "C:\\Users\\lblan\\OneDrive\\Escritorio\\CEAB\\2022\\Abundance_study\\dataplotsm.xlsx")
dataplot$start_date.x <- as_date(dataplot$start_date.x)
Sys.setlocale(category = "LC_ALL", locale = "EN")
plot02_16 <- abundances %>% 
  ggplot() +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_9")) +
  geom_line(data= dataplot, aes(x=start_date.x , y= `predictiondata$predi02_1`),col="darkblue", size = 1.9) +
  geom_line(data= dataplot, aes(x=start_date.x , y= predi02_6),col="red", size = 1.9) +
  scale_x_date("Sampling season", date_breaks = "1 month", date_labels = "%b") +
  ylab("Abundance (average)") + ylim(0,100) + ggtitle("0.2mm Acc. Rainfall") +
  # scale_y_continuous(breaks = seq(from = 0, to = 700, by = 20))+
  theme_classic(base_size = 18)

plot54_16 <- abundances %>% 
  ggplot() +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_9")) +
  geom_line(data= dataplot, aes(x=start_date.x , y= predi54_1),col="darkblue", size = 1.9) +
  geom_line(data= dataplot, aes(x=start_date.x , y= predi54_6),col="red", size = 1.9) +
  scale_x_date("Sampling season", date_breaks = "1 month", date_labels = "%b") +
  ylab("Abundance (average)") + ylim(0,100) + ggtitle("54mm Acc. Rainfall") +
  # scale_y_continuous(breaks = seq(from = 0, to = 700, by = 20))+
  theme_classic(base_size = 18)

plot110_16 <- abundances %>% 
  ggplot() +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_9")) +
  geom_line(data= dataplot, aes(x=start_date.x , y= predi110_1),col="darkblue", size = 1.9) +
  geom_line(data= dataplot, aes(x=start_date.x , y= predi110_6),col="red", size = 1.9) +
  scale_x_date("Sampling season", date_breaks = "1 month", date_labels = "%b") +
  ylab("Abundance (average)") + ylim(0,100) + ggtitle("110mm Acc. Rainfall") +
  # scale_y_continuous(breaks = seq(from = 0, to = 700, by = 20))+
  theme_classic(base_size = 18)

plot210_16 <- abundances %>% 
  ggplot() +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_9")) +
  geom_line(data= dataplot, aes(x=start_date.x , y= predi210_1),col="darkblue", size = 1.9) +
  geom_line(data= dataplot, aes(x=start_date.x , y= predi210_6),col="red", size = 1.9) +
  scale_x_date("Sampling season", date_breaks = "1 month", date_labels = "%b") +
  ylab("Abundance (average)") + ylim(0,150) + ggtitle("210mm Acc. Rainfall") +
  # scale_y_continuous(breaks = seq(from = 0, to = 700, by = 20))+
  theme_classic(base_size = 18)


library(patchwork)

(plot02_16 + plot54_16) / (plot110_16 + plot210_16)

#FIGURE S4 Supp.material

#We fix Temperature
mean(abundancias_df$mintemp21) #we take the average of the mintemp21
newdata <- data.frame(mintemp21 = 15.84,
                      max_relative_humidity = abundancias_df$max_relative_humidity,
                      rain21 = abundancias_df$rain21, 
                      n_imbornales_agua =43,
                      nrturist = abundancias_df$nrturist,
                      weeksince= abundancias_df$weeksince,
                      trap_name= "A_SP_BL_12")

predict(modelo_ofi, newdata, type="response", allow.new.levels = T)
predictiondata <- as.data.frame(predict(modelo_ofi, newdata, type="response",allow.new.levels = T))
predictiondata$predi <- predictiondata$`predict(modelo_ofi, newdata, type = "response", allow.new.levels = T)`
abundances <- cbind(abundancias_df, predictiondata)
abundances$week <- as.character(abundances$week)
abundances$week <- as.numeric(abundances$week)
abundances$start_date.x <- as_date(abundances$start_date.x)

Sys.setlocale(category = "LC_ALL", locale = "EN")
plot_temp_fix <- abundances %>% 
  ggplot() +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_9")) +
  geom_line(aes(x=start_date.x , y= predi),col="darkblue", size = 1.9) +
  scale_x_date("Sampling season", date_breaks = "1 month", date_labels = "%b") +
  ylab("Abundance (average)") + ylim(0,100) +
  # scale_y_continuous(breaks = seq(from = 0, to = 700, by = 20))+
  theme_classic(base_size = 21)
plot_temp_fix #plot b in panel Figure 3


#Now we fix and limit rain
newdata <- data.frame(mintemp21 = abundancias_df$mintemp21,
                      max_relative_humidity = abundancias_df$max_relative_humidity,
                      rain21 =54, 
                      n_imbornales_agua =43,
                      nrturist = abundancias_df$nrturist,
                      weeksince= abundancias_df$weeksince,
                      trap_name= "A_SP_BL_12")

predict(modelo_ofi, newdata, type="response", allow.new.levels = T)
predictiondata <- as.data.frame(predict(modelo_ofi, newdata, type="response",allow.new.levels = T))
predictiondata$predi <- predictiondata$`predict(modelo_ofi, newdata, type = "response", allow.new.levels = T)`
abundances <- cbind(abundancias_df, predictiondata)
abundances$week <- as.character(abundances$week)
abundances$week <- as.numeric(abundances$week)
abundances$start_date.x <- as_date(abundances$start_date.x)

Sys.setlocale(category = "LC_ALL", locale = "EN")
plot_rain_fix <- abundances %>% 
  ggplot() +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_4")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6, data = filter(abundances, trap_name == "A_SP_BL_5")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col= "#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_6")) +
  geom_line(aes(x= start_date.x, y= nrperspecies), col= "#DE9826", linewidth = 0.65, alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_7")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_8")) +
  geom_line(aes(x= start_date.x, y= nrperspecies),  col="#DE9826", linewidth = 0.65,alpha=0.6,data = filter(abundances, trap_name == "A_SP_BL_9")) +
  geom_line(aes(x=start_date.x , y= predi),col="darkblue", size = 1.9) +
  scale_x_date("Sampling season", date_breaks = "1 month", date_labels = "%b") +
  ylab("Abundance (average)") + ylim(0,150) +
  # scale_y_continuous(breaks = seq(from = 0, to = 700, by = 20))+
  theme_classic(base_size = 21)
plot_rain_fix #plot b in panel Figure 3

(plot_temp_fix + plot_rain_fix) + plot_annotation(tag_levels = "a")
