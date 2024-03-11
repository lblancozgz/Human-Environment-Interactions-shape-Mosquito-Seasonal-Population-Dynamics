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
  rename("end_date" = `End date`)
abundancias_aedes <- abundancias_aedes %>%
  rename("trap_name" = `trap name`)

abundancias_aedes_bl <- filter(abundancias_aedes, city == "Blanes")
abundancias_aedes_bl <- abundancias_aedes_bl %>% 
  rename(trapef = 'trapping effort')
abundancias_aedes_bl <- abundancias_aedes_bl %>% 
  rename(land_use = 'Land use')
ggplot(data = abundancias_aedes) +
  geom_line(aes(x= start_date, y= nrperspecies, color = trap_name)) #plot todas trampas

ggplot(data = abundancias_aedes_bl) +
  geom_line(aes(x= start_date, y= nrperspecies, color = trap_name)) #plot trampas blanes

#eliminamos la que no estaba en el jardín y la de la entrada (disaster)
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_1")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_2")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_10")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_11")
abundancias_aedes_bl <- filter(abundancias_aedes_bl, !trap_name == "A_SP_BL_12")


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


#GENERALISED LINEAR MODELS
library(DHARMa)
library(performance)


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

library(pollen)
library(meteor)

#tenemos que hacer el acúmulo para cada año!!!!!!!!!!!!!! y photoperiod
library(zoo)
clima_malgrat$photoperiod <- photoperiod(1:365, 41.676854)
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
         rain21 = rollapplyr(precipitation, 21, FUN=sum, align = "right", fill = NA)
  )



str(clima_malgrat$YEAR)
range(abundancias_aedes_bl$nrperspecies)
#chequeamos nuevas variables

range(clima_malgrat$mwime, na.rm = T)
range(clima_malgrat$photoperiod)
ggplot(clima_malgrat, aes(x=doy, y=mwi7, color=YEAR)) + geom_line()+ facet_wrap(~YEAR)
ggplot(clima_malgrat, aes(x=doy, y=photoperiod, color=YEAR)) + geom_line() + facet_wrap(~YEAR)
ggplot(clima_malgrat, aes(x=doy, y=pho21, color=YEAR)) + geom_line() + facet_wrap(~YEAR)

clima_malgrat$start_date <- as.POSIXct(clima_malgrat$start_date,format = '%Y-%m-%d')
abundancias_aedes_bl$start_date <- as.POSIXct(abundancias_aedes_bl$start_date,format = '%Y-%m-%d')

#juntamos clima con datos
abundancias_jardin_tot <- full_join(abundancias_aedes_bl, clima_malgrat, 
                                   by = c("start_date", "YEAR", "week", "month", "doy"), all= TRUE) 

str(clima_malgrat$start_date)

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
  rain7 <- clima_malgrat$rain7[clima_malgrat$start_date %in% x1]
  rain14 <- clima_malgrat$rain14[clima_malgrat$start_date %in% x1]
  rain21 <-clima_malgrat$rain21[clima_malgrat$start_date %in% x1]
  mintemp7 <- clima_malgrat$mintemp7[clima_malgrat$start_date %in% x1]
  mintemp14 <- clima_malgrat$mintemp14[clima_malgrat$start_date %in% x1]
  mintemp21 <- clima_malgrat$mintemp21[clima_malgrat$start_date %in% x1]
  minrh7 <- clima_malgrat$minrh7[clima_malgrat$start_date %in% x1]
  minrh14 <- clima_malgrat$minrh14[clima_malgrat$start_date %in% x1]
  minrh21 <- clima_malgrat$minrh21[clima_malgrat$start_date %in% x1]
  meantemp7 <- clima_malgrat$meantemp7[clima_malgrat$start_date %in% x1]
  meantemp14 <- clima_malgrat$meantemp14[clima_malgrat$start_date %in% x1]
  meantemp21 <- clima_malgrat$meantemp21[clima_malgrat$start_date %in% x1]
  meanrh7 <- clima_malgrat$meanrh7[clima_malgrat$start_date %in% x1]
  meanrh14 <- clima_malgrat$meanrh14[clima_malgrat$start_date %in% x1]
  meanrh21 <- clima_malgrat$meanrh21[clima_malgrat$start_date %in% x1]
  maxtemp7 <- clima_malgrat$maxtemp7[clima_malgrat$start_date %in% x1]
  maxtemp14 <- clima_malgrat$maxtemp14[clima_malgrat$start_date %in% x1]
  maxtemp21 <- clima_malgrat$maxtemp21[clima_malgrat$start_date %in% x1]
  maxrh7 <- clima_malgrat$maxrh7[clima_malgrat$start_date %in% x1]
  maxrh14 <- clima_malgrat$maxrh14[clima_malgrat$start_date %in% x1]
  maxrh21 <- clima_malgrat$maxrh21[clima_malgrat$start_date %in% x1]
  pho7 <- clima_malgrat$pho7[clima_malgrat$start_date %in% x1]
  pho14 <- clima_malgrat$pho14[clima_malgrat$start_date %in% x1]
  pho21 <- clima_malgrat$pho21[clima_malgrat$start_date %in% x1]
  mwi7 <- clima_malgrat$mwi7[clima_malgrat$start_date %in% x1]
  mwi14 <- clima_malgrat$mwi14[clima_malgrat$start_date %in% x1]
  mwi21 <- clima_malgrat$mwi21[clima_malgrat$start_date %in% x1]
  
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
  abundancias_jardin_tot$mwi14[i] <- mwi21
  
}

abundancias_jardin_tot <- filter(abundancias_jardin_tot, !is.na(trapef))
#GDD
# row by row
# Requires paraell or parallelly package, depends on the R version
library(parallelly)

library(parallel)

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

str(abundances_with_gdd$mwi14)
range(abundances_with_gdd$mwi21)
abundances_with_gdd$week <- as.character(abundances_with_gdd$week)
abundances_with_gdd$week <- as.numeric(abundances_with_gdd$week)
ggplot(abundances_with_gdd, aes(x=week, y=rain21, color=YEAR)) + geom_line() + facet_wrap(~YEAR)

options(na.action = "na.fail")
lapply(abundances_with_gdd, function(x) any(is.na(x)))
which(is.na(abundances_with_gdd$nrperspecies)) #remove all NAs that can affect negatively to the df


#more plots
library(ggthemes)
abundancias21 <- filter(abundances_with_gdd, YEAR == 2021)
range(abundancias21$maxtemp21)
abundancias21$week <- as.character(abundancias21$week)
abundancias21$week <- as.numeric(abundancias21$week)
str(abundancias21$week)

  ylim.prim <- c(0, 321)   
  ylim.sec <- c(0, 140)
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1]
ggplot(abundances_with_gdd, aes(x=week)) +
    geom_col( aes(y=nrperspecies), fill="#527614", alpha= 0.6) + 
    geom_line(aes(y = a + l0gdd*b), color = "red", linewidth = 1.5) +
    geom_line(aes(y=rain21), color= "blue", linewidth = 1.5) +
    scale_x_continuous(breaks = seq(from = 18, to = 46, by = 3))+
    scale_y_continuous("Abundance", sec.axis = sec_axis(~ (. - a)/b, name = "Maximum temperatures")) +
    theme_few(base_size = 22)

ylim.prim <- c(0, 321)   
ylim.sec <- c(15,30)
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
ggplot(abundancias21, aes(x=week)) +
  geom_col( aes(y=nrperspecies), fill="#527614", alpha= 0.6) + 
  geom_line(aes(y = a + max_temperature*b), color = "red", linewidth = 1.5) +
  geom_line(aes(y=precipitation), color= "blue", linewidth = 1.5) +
  scale_x_continuous(breaks = seq(from = 18, to = 46, by = 3))+
  scale_y_continuous("Abundance", sec.axis = sec_axis(~ (. - a)/b, name = "Maximum temperatures")) +
  theme_few(base_size = 22)

ylim.prim <- c(0, 321)   
ylim.sec <- c(30,100)
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
ggplot(abundancias21, aes(x=week)) +
  geom_col( aes(y=nrperspecies), fill="gray", alpha= 0.6) + 
  geom_line(aes(y = a + max_relative_humidity*b), color = "red", linewidth = 1.5) +
  geom_line(aes(y = a + mean_relative_humidity*b), color = "lightblue", linewidth = 1.5) +
  geom_line(aes(y = a + min_relative_humidity*b), color = "darkgreen", linewidth = 1.5) +
  geom_line(aes(y=precipitation), color= "blue", linewidth = 1.5) +
  scale_x_continuous(breaks = seq(from = 18, to = 46, by = 3))+
  scale_y_continuous("Abundance", sec.axis = sec_axis(~ (. - a)/b, name = "Maximum temperatures")) +
  theme_few(base_size = 22)



library(MuMIn)
library(glmmTMB)
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

# install.packages("GLMMadaptive")
library(GLMMadaptive)
library(car)
m0 <-glmer.nb(nrperspecies ~ scale(max_temperature) + scale(min_temperature) + scale(mean_temperature) +
                  scale(max_relative_humidity) + scale(mean_relative_humidity) + scale(min_relative_humidity) +
                  scale(precipitation) + scale(mwime) + log(trapef) + scale(l0gdd) +
                   (1|trap_name),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                 data = abundances_with_gdd)

summary(m0)
vif(m0)

library(MuMIn)
ms3 <- dredge(m0, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))

NCM3 <- get.models(ms3, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)

NCMDF <- model.sel(NCM3)
model.avg(NCM3)
sw(NCM3)

m0_dredge <- glmer.nb(nrperspecies ~ scale(max_relative_humidity) + scale(min_relative_humidity) + 
                          scale(min_temperature) + scale(precipitation) + log(trapef) + (1|trap_name),
                        control = glmerControl(optimizer = "Nelder_Mead"),
                        data = abundances_with_gdd)
summary(m0_dredge)
car::vif(m0_dredge) # maxtemperature out


testZeroInflation(m0_dredge)
testResiduals(m0_dredge)
testDispersion(m0_dredge)


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

m1_dredge <-glmer.nb(nrperspecies ~  log(trapef)+ scale(maxrh7) + scale(minrh7) +  scale(l7gdd)  +
                (1|trap_name), 
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                data = abundances_with_gdd) #in dredge, also mintemp7 and maxtemp7 were added, but very collinear
summary(m1_dredge)
car::vif(m1_dredge)
testZeroInflation(m1_dredge)
testResiduals(m1_dredge)
testDispersion(m1_dredge)

m2 <-glmer.nb(nrperspecies ~ scale(maxtemp14) + scale(mintemp14) + scale(meantemp14) +
                scale(maxrh14) + scale(meanrh14) + scale(minrh14) +
                scale(rain14) + scale(mwi14) + log(trapef) + scale(l14gdd) +
                (1|trap_name), 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
              data = abundances_with_gdd)

summary(m2)
ms5 <- dredge(m2, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))
NCM5 <- get.models(ms5, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDF5 <- model.sel(NCM5)

m2_dredge <-glmer.nb(nrperspecies ~  log(trapef) + scale(maxrh14) + scale(minrh14) +  scale(l14gdd) + scale(rain14) +
                       (1|trap_name), 
                     control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                     data = abundances_with_gdd)
summary(m2_dredge)
car::vif(m2_dredge)
testZeroInflation(m2_dredge)
testResiduals(m2_dredge)
testDispersion(m2_dredge)


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

m3_dredge <-glmer.nb(nrperspecies ~ log(trapef) + scale(maxrh21) + scale(minrh21) + scale(mwi21) + scale(rain21) +
                       (1|trap_name), 
                     control = glmerControl(optimizer = "Nelder_Mead"),
                     data = abundances_with_gdd)

summary(m3_dredge)
car::vif(m3_dredge)
testZeroInflation(m3_dredge)
testResiduals(m3_dredge)
testDispersion(m3_dredge)

abundances_with_gdd$l7gdd

fullmodel <- glmer.nb(nrperspecies ~ log(trapef) + scale(maxrh21) + scale(minrh21)
                      + scale(mwi21) + scale(rain21) + scale(maxrh14) + scale(minrh14) +
                        scale(l14gdd) + scale(rain14) + scale(maxrh7) + scale(minrh7) +
                        scale(l7gdd) + scale(max_relative_humidity) + scale(min_relative_humidity) +
                        scale(min_temperature) + scale(precipitation) + (1|trap_name), 
                      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                      data = abundances_with_gdd) 
summary(fullmodel)
mfull <- dredge(fullmodel, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 5), extra= c(max.r))
NCMfull <- get.models(mfull, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDFfull <- model.sel(NCMfull)

full_dredge <- glmer.nb(nrperspecies ~ log(trapef) + scale(maxrh21) + scale(minrh21) +
                          scale(mwi21) + scale(rain14) +(1|trap_name), 
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                        data = abundances_with_gdd) 
summary(full_dredge)
car::vif(full_dredge)
testZeroInflation(full_dredge)
testResiduals(full_dredge)
testDispersion(full_dredge)
testUniformity(full_dredge)

##winner model
NEWd <- data.frame(mwi21 = abundances_with_gdd$mwi21,
                   minrh21 = abundances_with_gdd$minrh21,
                   maxrh21 = abundances_with_gdd$maxrh21,
                   rain14 = abundances_with_gdd$rain14,
                   trapef = 8,
                   trap_name=  "A_SP_BL_15")

plot(predict(full_dredge, NEWd, type="response", allow.new.levels = TRUE))
predict(full_dredge, NEWd, type="response", allow.new.levels = TRUE)
predictiondata <- as.data.frame(predict(full_dredge, NEWd, type="response", allow.new.levels = TRUE))
predictiondata$predi <- predictiondata$`predict(full_dredge, NEWd, type = "response", allow.new.levels = TRUE)`
abundances <- cbind(abundances_with_gdd, predictiondata)
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


#TOURIST

datatourist <- read_excel("Origen i Pernocta 2021_XL.xlsx", 
                                        col_types = c("date", "numeric"))

datatourist$start_date<- strftime(datatourist$start_date, format="%Y-%m-%d")
datatourist$start_date <- as_date(as.POSIXct(datatourist$start_date,format = '%Y-%m-%d'))
abundances_with_gdd$start_date <- as_date(as.POSIXct(abundances_with_gdd$start_date,format = '%Y-%m-%d'))

str(datatourist$start_date)
str(abundances_with_gdd$start_date)


#grouping by date

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
library(ggthemes)
# plot_turist<- datatourist %>% 
#   group_by(month) %>% 
#   summarise(count = sum(nrturist)) %>% 
#   drop_na() %>% 
#   ggplot(aes(x = month, y = count, fill = month)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   scale_fill_brewer(name = "Month", palette = "Paired",
#                     labels = c("January", "February", "March", "April", "May", "June",
#                                "July", "August", "September", "October", "November", "December")) +
#   geom_text(aes(label=count), position=position_dodge(width=0.9), vjust=-0.25, size = 4) +
#   scale_x_discrete(labels=c("1" = "January", "2" = "February", "3" = "March", "4" = "April",
#                             "5" = "May", "6" = "June", "7"= "July","8" = "August", "9" = "September",
#                             "10" = "October", "11"= "November", "12" = "December")) +
#   xlab("Month (2021)") + ylab("Number of tourists") +
#   theme_few(base_size = 22) + scale_color_few()
# plot_turist + theme(legend.position = "none")

#juntamos clima con datos
abundances_with_tourist<- merge(abundances_with_gdd, datatourist, by = c("start_date", "week", "month"), all.x = T) 
abundances_with_tourist_plot <- abundances_with_tourist %>% 
  group_by(month) %>% 
  summarise(nrturistmonth = sum(nrturist))

range(abundances_with_tourist_plot$nrturistmonth)
a <- ylim.prim[1] - b*ylim.sec[1]
datatourist$month <- as.character(datatourist$month)
datatourist$month <- as.numeric(datatourist$month)

modelo_turis <- glmer.nb(nrperspecies ~ log(trapef) + scale(nrturist) + scale(maxrh21) + scale(minrh21) +
                           scale(mwi21) + scale(rain14) +(1|trap_name), 
                         control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                         data = abundances_with_tourist) 
summary(modelo_turis)
car::vif(modelo_turis)

mtur <- dredge(modelo_turis, trace = 2, rank = "AICc", REML = FALSE, m.lim=c(0, 6), extra= c(max.r))
NCMtur <- get.models(mtur, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
NCMDFtur <- model.sel(NCMtur)
modelo_turis_dredge <- glmer.nb(nrperspecies ~ log(trapef) + scale(nrturist) + scale(maxrh21) + scale(minrh21) +
                           scale(mwi21) + scale(rain14) +(1|trap_name), 
                         control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                         data = abundances_with_tourist) 
summary(modelo_turis_dredge)
sw(NCMDFtur)
NEWd <- data.frame(mwi21 = abundances_with_tourist$mwi21,
                   minrh21 = abundances_with_tourist$minrh21,
                   maxrh21 = abundances_with_tourist$maxrh21,
                   rain14 = abundances_with_tourist$rain14,
                   nrturist = abundances_with_tourist$nrturist,
                   trapef = 8,
                   trap_name=  "A_SP_BL_15")

plot(predict(modelo_turis, NEWd, type="response", allow.new.levels = TRUE))
predict(modelo_turis, NEWd, type="response", allow.new.levels = TRUE)
predictiondata <- as.data.frame(predict(modelo_turis, NEWd, type="response", allow.new.levels = TRUE))
predictiondata$predi <- predictiondata$`predict(modelo_turis, NEWd, type = "response", allow.new.levels = TRUE)`
abundances <- cbind(abundances_with_tourist, predictiondata)
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







# plot_abun<- abundances_with_tourist %>% 
#   group_by(month, trap_name) %>% 
#   summarise(sumcoun = sum(nrperspecies)) %>% 
#   drop_na() %>% 
#   ggplot(aes(x = month, y = sumcoun, fill = trap_name)) +
#   geom_bar(stat = "identity", position = position_dodge()) + 
#   scale_fill_manual(name = "Trap name", 
#                       values= c("#2C8437", "#71266E", "#00989F", "#AA3F39", "#2C417A", "#FF7100"), 
#                       labels = c("Trap 4", "Trap 5", "Trap 6", "Trap 7", "Trap 8", "Trap 9")) +
#   geom_text(aes(label=sumcoun), position=position_dodge(width=0.9), vjust=-0.25, size = 4) +
#   scale_x_discrete(labels=c("5" = "May", "6" = "June", "7" = "July", "8" = "August",
#                             "9" = "September", "10" = "October", "11"= "November")) +
#   xlab("Month (2021)") + ylab("Abundance") +
#   theme_classic(base_size = 25)
# 
# figturist <- visreg::visreg(modelo_turis, "nrturist",  overlay = TRUE, trans = exp,
#                          points=list(size=4, pch = 16, alpha = 0.35),
#                          line.par=list(col="#71266E"),
#                          whitespace=0.1, gg = TRUE, rug = FALSE, partial = TRUE, legend = F) +
#     theme_few() +
#     scale_colour_manual(values = c("#71266E"))+
#     theme(axis.title.x = element_text(size=18), 
#           axis.title.y = element_text(size = 18), 
#           axis.text.x = element_text(size = 16, face = "italic"),
#           axis.text.y = element_text(size = 16),
#           legend.position = "bottom",
#           legend.title=element_text(size = 18),
#           legend.text = element_text(size = 17.75, face = "italic")) +
#     labs(y = "Abundances",
#          x = "Number of turist") 
# 
# ggplot(abundances_with_gdd, aes(x= week, y= photoperiod)) + geom_line() + 
#   scale_x_continuous(breaks = seq(from = 18, to = 46, by = 3))
