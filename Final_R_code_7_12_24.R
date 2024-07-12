## Set working directory & library ####
library(fishboot)
library(TropFishR)
library(tidyverse)

wd = "C:/Users/06061016/OneDrive - Nord universitet/Desktop/Mulimbwa/Final_Mulimbwa_7_10_24"
setwd(wd)
set.seed(1)

## LFQ ESTIMATION 1987 - 1989 ########################## ######
## 1. Load data ####
load("Stolothrissa-8788.Rdata")
load("Limnothrissa-8788.Rdata")


## As an artifact of digitizing, some data appeared which is wrong. 
## Limits according to the publication from which the data was taken (Mulimbwa and Shirakihara 1994), 
## Limnothrissa: 2.25 - 11.25
## Stolothrissa: 2.25 - 9.25
## turned into 0
limnothrissa9[4,] # This needs to be corrected
limnothrissa9[4,2:19] = 0
limnothrissa9[24,] # This needs to be corrected
limnothrissa9[24,2:19] = 0

stolothrissa9[4,] # That is already fine
stolothrissa9[20,] # This needs to be corrected
stolothrissa9[20,2:19] = 0
stolothrissa9[21,2:19] = 0

# Transform the data from SL to TL. This is done with a and b values from 
# Growth parameter values in the "Growth and mortialities parameters.xls" 
data_std = stolothrissa9 %>% mutate(lengthClass =  0.08729 + 1.1562 * lengthClass)
data_lmd = limnothrissa9 %>% mutate(lengthClass =  round(0.16658 + 1.1873 * lengthClass,2))

data_dig = rbind(mutate(data_std, Species = "Stolo"),
                 mutate(data_lmd, Species = "Limno"))

# Save the data transformed into TL and without the digitalisation mistakes

write.csv(data_dig, "data_final_8788_both_species.csv")
his = data.frame(length = data_dig$lengthClass, Species = data_dig$Species,
                 freq = rowSums(data_dig[,2:(ncol(data_dig)-1)],na.rm = T))

# Big individuals not present at all
ggplot(his, aes(x=length, y=freq)) +
  geom_bar(stat = "identity", width = 2) + facet_wrap(~ Species) +
  ggtitle("Stolothrissa raw data distribution across lengths") +
  theme(plot.title = element_text(hjust = 0.5))

colnames(data_std)[2:ncol(data_std)] <- as.character(as.Date(colnames(data_std)[2:ncol(data_std)]))
colnames(data_lmd)[2:ncol(data_lmd)] <- as.character(as.Date(colnames(data_lmd)[2:ncol(data_lmd)]))

dates_std <- paste(colnames(data_std)[2:ncol(data_std)],sep="")
dates_std <- as.Date(dates_std, "%Y-%m-%d")
dates_lmd <- paste(colnames(data_lmd)[2:ncol(data_lmd)],sep="")
dates_lmd <- as.Date(dates_lmd, "%Y-%m-%d")


data_new_std <- list(dates = dates_std,
                     midLengths = data_std$lengthClass,
                     catch = as.matrix(data_std[,2:ncol(data_std)]))
data_new_lmd <- list(dates = dates_lmd,
                     midLengths = data_lmd$lengthClass,
                     catch = as.matrix(data_lmd[,2:ncol(data_lmd)]))

## 2. Create LFQ objects ####
class(data_new_std) <- "lfq"
class(data_new_lmd) <- "lfq"
lfq0.std1 <- data_new_std
lfq1.std1 <- lfqModify(lfq0.std1, bin_size = 0.8) # Stolo
lfq0.lmd1 <- data_new_lmd
lfq1.lmd1 <- lfqModify(lfq0.lmd1, bin_size = 0.8) # Limno

par(mfrow = c(2,2))
plot(lfq0.std1, Fname = "catch", main = "Stolothrissa 88")
plot(lfq1.std1, Fname = "catch", main = "Stolothrissa 88")
plot(lfq0.lmd1, Fname = "catch", main = "Limnothrissa 88")
plot(lfq1.lmd1, Fname = "catch", main = "Limnothrissa 88")

# Assign a moving average.
ma <- 5 # See youngest cohort
lfq2a.std1 <- lfqRestructure(lfq0.std1, MA = ma, addl.sqrt = FALSE)
lfq2b.std1 <- lfqRestructure(lfq1.std1, MA = ma, addl.sqrt = FALSE)
ma <- 5 # See youngest cohort
lfq2a.lmd1 <- lfqRestructure(lfq0.lmd1, MA = ma, addl.sqrt = FALSE)
lfq2b.lmd1 <- lfqRestructure(lfq1.lmd1, MA = ma, addl.sqrt = FALSE)
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))

plot(lfq2a.std1, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 88")
plot(lfq2b.std1, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 88")
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2a.lmd1, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 88")
plot(lfq2b.lmd1, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 88")

lfq2.std <- lfq2a.std1
lfq2.lmd <- lfq2a.lmd1

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2.std, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 88")
plot(lfq2.lmd, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 88")

## 3. Estimate parameters ####
## Stolothrissa
set.seed(1)
Sys.time()
lfq_GA.std <- fishboot::ELEFAN_GA_boot(
  lfq2.std,  MA = lfq2.std$MA, seasonalised = F, popSize=100, pmutation = 0.5,
  low_par = list(Linf = 6, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 12, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 18/12)
print("Fi del primer")

Sys.time()
lfq_SA.std <- fishboot::ELEFAN_SA_boot(
  lfq2.std,  MA = lfq2.std$MA,seasonalised = F, SA_temp = 10e+05,maxit = 500,
  low_par = list(Linf = 6, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 12, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 18/12)
Sys.time()
print("Fi del segon")
Sys.time()
saveRDS(list(lfq_GA.std, lfq_SA.std), "Stolothrissa_digitalized_87_88_bootstrapped.Rdata")

## Limnothrissa
set.seed(1)
Sys.time()
lfq_GA.lmd <- fishboot::ELEFAN_GA_boot(
  lfq2.lmd,  MA = lfq2.lmd$MA,seasonalised = F, popSize=100, pmutation = 0.5,
  low_par = list(Linf = 8, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 22, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 30/12)
print("Fi del primer")
Sys.time()
lfq_SA.lmd <- fishboot::ELEFAN_SA_boot(
  lfq2.lmd,  MA = lfq2.lmd$MA,seasonalised = F, SA_temp = 10e+05,maxit = 500,
  low_par = list(Linf = 8, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 22, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 30/12)
print("Fi del segon")
Sys.time()
saveRDS(list(lfq_GA.lmd, lfq_SA.lmd), "Limnothrissa_digitalized_87_88_bootstrapped.Rdata")

## LFQ ESTIMATION 1999 - 2001 ########################## ####
## Stolothrissa ####
## 1. Load data ####
data_st2 = readxl::read_excel("Stolothrissa_data_LFQ.xlsx", sheet = "Stolo_99_01")
data_st2[is.na(data_st2)] = 0

## 2. Create LFQ object ####
colnames(data_st2)[1] <- "lengthClass"
colnames(data_st2)[2:ncol(data_st2)] <- as.character(as.Date(as.numeric(colnames(data_st2)[2:ncol(data_st2)]), origin = "1900-01-01"))

which(colSums(data_st2[,2:ncol(data_st2)])==0)
data_st2 = data_st2[,-22]

dates_st2 <- paste(colnames(data_st2)[2:ncol(data_st2)],sep="")
dates_st2 <- as.Date(dates_st2, "%Y-%m-%d")

data.new.st2 <- list(dates = dates_st2,
                     midLengths = data_st2$lengthClass/10,
                     catch = as.matrix(data_st2[,2:ncol(data_st2)]))

data.new.st2
class(data.new.st2) <- "lfq"

plot(data.new.st2, Fname = "catch", main = "Stolothrissa 99-01")

lfq0.st2 <- data.new.st2
lfq1.st2 <- lfqModify(lfq0.st2, bin_size = 0.8) # Stolo

par(mfrow = c(2,1))
plot(lfq0.st2, Fname = "catch", main = "Stolothrissa 99-01")
plot(lfq1.st2, Fname = "catch", main = "Stolothrissa 99-01")

## 3. Assign a moving average. 
ma <- 7 # See youngest cohort August 1999 has 6 bins
lfq2a.st2 <- lfqRestructure(lfq0.st2, MA = ma, addl.sqrt = FALSE)
lfq2b.st2 <- lfqRestructure(lfq1.st2, MA = ma, addl.sqrt = FALSE)

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2a.st2, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 99-01")
plot(lfq2b.st2, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 99-01")

lfq2.st2 <- lfq2b.st2

## 3. Estimate Parameters Stolothrissa ####
Sys.time()
lfq_fin_st_GA <- fishboot::ELEFAN_GA_boot(
  lfq2.st2,  MA = lfq2.st2$MA,seasonalised = F, popSize=100, pmutation = 0.5,
  low_par = list(Linf = 6, K = 1, t_anchor = 0, C = 0, ts = 0),  #t_anchor = 7/12
  up_par  = list(Linf = 12, K = 10, t_anchor = 1, C = 1, ts = 1),#t_anchor = 9/12
  nresamp = 400, agemax = 18/12)

print("Fi del primer")
Sys.time()

lfq_fin_st_SA <- fishboot::ELEFAN_SA_boot(
  lfq2.st2,  MA = lfq2.st2$MA,seasonalised = F, SA_temp = 10e+05,maxit = 500,
  low_par = list(Linf = 6, K = 1, t_anchor = 0, C = 0, ts = 0),#t_anchor = 7/12
  up_par  = list(Linf = 12, K = 10, t_anchor = 1, C = 1, ts = 1),  #t_anchor = 9/12
  nresamp = 400, agemax = 18/12)
print("Fi del segon")
Sys.time()

saveRDS(list(lfq_fin_st_GA,lfq_fin_st_SA), "Stolothrissa_99_01_bootstrapped.Rdata")

## Limnothrissa ####
## 1. Load data #### 
data_lm_99 = readxl::read_excel("Limnothrissa_data_LFQ.xlsx", sheet = "99-01")

## 2. Create LFQ objects ####

colnames(data_lm_99)[1] <- "lengthClass"
colnames(data_lm_99)[2:ncol(data_lm_99)] <- as.character(as.Date(as.numeric(colnames(data_lm_99)[2:ncol(data_lm_99)]), origin = "1900-01-01"))

which(colSums(data_lm_99[,2:ncol(data_lm_99)])==0)

dates_lm_99 <- paste(colnames(data_lm_99)[2:ncol(data_lm_99)],sep="")
dates_lm_99 <- as.Date(dates_lm_99, "%Y-%m-%d")

data_new_lm_99 <- list(dates = dates_lm_99,
                       midLengths = data_lm_99$lengthClass/10,
                       catch = as.matrix(data_lm_99[,2:ncol(data_lm_99)]))

class(data_new_lm_99) <- "lfq"
which(colSums(data_lm_99[,2:ncol(data_lm_99)])==0)

dates_lm_99 <- paste(colnames(data_lm_99)[2:ncol(data_lm_99)],sep="")
dates_lm_99 <- as.Date(dates_lm_99, "%Y-%m-%d")

data_new_lm_99 <- list(dates = dates_lm_99,
                     midLengths = data_lm_99$lengthClass/10,
                     catch = as.matrix(data_lm_99[,2:ncol(data_lm_99)]))

class(data_new_lm_99) <- "lfq"

plot(data_new_lm_99, Fname = "catch", main = "Limnothrissa 99-01")

lfq0.lm2 <- data_new_lm_99
lfq1.lm2 <- lfqModify(lfq0.lm2, bin_size = 0.8) # Limnothrissa

par(mfrow = c(2,1))
plot(lfq0.lm2, Fname = "catch", main = "Limnothrissa 99-01")
plot(lfq1.lm2, Fname = "catch", main = "Limnothrissa 99-01")

#Assign a moving average. 
ma <- 7 # Youngest cohort December 2001
lfq2a.lm2 <- lfqRestructure(lfq0.lm2, MA = ma, addl.sqrt = FALSE)
lfq2b.lm2 <- lfqRestructure(lfq1.lm2, MA = ma, addl.sqrt = FALSE)

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2a.lm2, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 99-01")
plot(lfq2b.lm2, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 99-01")

lfq2.lm_99 <- lfq2b.lm2

## 3. Estimate parameters Limnothrissa ####
## Limnothrissa 99-2001 

Sys.time()
lfq_GA.lm_99 <- fishboot::ELEFAN_GA_boot(
  lfq2.lm_99,  MA = lfq2.lm_99$MA,seasonalised = F, popSize=100, pmutation = 0.5,
  low_par = list(Linf = 8, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 22, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 30/12)
print("Fi del primer")
Sys.time()

lfq_SA.lm_99 <- fishboot::ELEFAN_SA_boot(
  lfq2.lm_99,  MA = lfq2.lm_99$MA,seasonalised = F, SA_temp = 10e+05,maxit = 500,
  low_par = list(Linf = 8, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 22, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 30/12)
print("Fi del segon")
Sys.time()

saveRDS(list(lfq_GA.lm_99, lfq_SA.lm_99), "Limnothrissa_99_01_bootstrapped.Rdata")

## LFQ ESTIMATION 2007 - 2008 ########################## ####
## Stolothrissa ####
## 1. Load data ####
set.seed(1)
data_st_07 = readxl::read_excel("LFQ_stolo_2007_2008.xls", sheet = "Feuil3")
data_st_07[is.na(data_st_07)] = 0

## 2. Create LFQ object ####
colnames(data_st_07)[1] <- "lengthClass"
colnames(data_st_07)[2:ncol(data_st_07)] <- as.character(as.Date(as.numeric(colnames(data_st_07)[2:ncol(data_st_07)]), origin = "1900-01-01"))

which(colSums(data_st_07[,2:ncol(data_st_07)])==0)

dates_st_07 <- paste(colnames(data_st_07)[2:ncol(data_st_07)],sep="")
dates_st_07 <- as.Date(dates_st_07, "%Y-%m-%d")

data.new.st_07 <- list(dates = dates_st_07,
                     midLengths = data_st_07$lengthClass/10,
                     catch = as.matrix(data_st_07[,2:ncol(data_st_07)]))

class(data.new.st_07) <- "lfq"

plot(data.new.st_07, Fname = "catch", main = "Stolothrissa 07-08")

lfq0.st_07 <- data.new.st_07
lfq1.st_07 <- lfqModify(lfq0.st_07, bin_size = 0.8) # Stolo


par(mfrow = c(2,1))
plot(lfq0.st_07, Fname = "catch", main = "Stolothrissa 07-08")
plot(lfq1.st_07, Fname = "catch", main = "Stolothrissa 07-08")


## 3. Assign a moving average. 
ma <- 5 # See youngest cohort August 1999
lfq2a.st_07 <- lfqRestructure(lfq0.st_07, MA = ma, addl.sqrt = FALSE)
lfq2b.st_07 <- lfqRestructure(lfq1.st_07, MA = ma, addl.sqrt = FALSE)

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2a.st_07, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 99-01")
plot(lfq2b.st_07, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 99-01")

lfq2.st_07 <- lfq2b.st_07

## 3. Estimate Parameters Stolothrissa ####
library(fishboot)
Sys.time()
lfq_fin_st_GA_07 <- fishboot::ELEFAN_GA_boot(
  lfq2.st_07,  MA = lfq2.st_07$MA,seasonalised = F, popSize=100, pmutation = 0.5,
  low_par = list(Linf = 6, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 12, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 18/12)

print("Fi del primer")
Sys.time()

lfq_fin_st_SA_07 <- fishboot::ELEFAN_SA_boot(
  lfq2.st_07,  MA = lfq2.st_07$MA,seasonalised = F, SA_temp = 10e+05,maxit = 500,
  low_par = list(Linf = 6, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 12, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 18/12)
print("Fi del segon")
Sys.time()

saveRDS(list(lfq_fin_st_GA_07,lfq_fin_st_SA_07), "Stolothrissa_07_bootstrapped.Rdata")

## Limnothrissa ####
## 1. Load data ####
set.seed(1)
data_lm_07 = readxl::read_excel("LFQ_limno_2007_2008.xls", sheet = "Feuil3")

data_lm_07[is.na(data_lm_07)] = 0

## 2. Create LFQ object ####
colnames(data_lm_07)[1] <- "lengthClass"
colnames(data_lm_07)[2:ncol(data_lm_07)] <- as.character(as.Date(as.numeric(colnames(data_lm_07)[2:ncol(data_lm_07)]), origin = "1900-01-01"))

which(colSums(data_lm_07[,2:ncol(data_lm_07)])==0)

dates_lm_07 <- paste(colnames(data_lm_07)[2:ncol(data_lm_07)],sep="")
dates_lm_07 <- as.Date(dates_lm_07, "%Y-%m-%d")

data.new.lm_07 <- list(dates = dates_lm_07,
                     midLengths = data_lm_07$lengthClass/10,
                     catch = as.matrix(data_lm_07[,2:ncol(data_lm_07)]))

class(data.new.lm_07) <- "lfq"

plot(data.new.lm_07, Fname = "catch", main = "Limnothrissa 99-01")

lfq0.lm_07 <- data.new.lm_07
lfq1.lm_07 <- lfqModify(lfq0.lm_07, bin_size = 0.8) # Stolo


par(mfrow = c(2,1))
plot(lfq0.lm_07, Fname = "catch", main = "Limnothrissa 99-01")
plot(lfq1.lm_07, Fname = "catch", main = "Limnothrissa 99-01")


## 3. Assign a moving average.
ma <- 7 # See youngest cohort August 2007
lfq2a.lm_07 <- lfqRestructure(lfq0.lm_07, MA = ma, addl.sqrt = FALSE)
lfq2b.lm_07 <- lfqRestructure(lfq1.lm_07, MA = ma, addl.sqrt = FALSE)

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2a.lm_07, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 99-01")
plot(lfq2b.lm_07, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 99-01")

lfq2.lm_07 <- lfq2b.lm_07

## 3. Estimate Parameters Limnothrissa ####
library(fishboot)
Sys.time()
lfq_fin_lm_GA_07 <- fishboot::ELEFAN_GA_boot(
  lfq2.lm_07,  MA = lfq2.lm_07$MA,seasonalised = F, popSize=100, pmutation = 0.5,
  low_par = list(Linf = 8, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 22, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 30/12)

print("Fi del primer")
Sys.time()

lfq_fin_lm_SA_07 <- fishboot::ELEFAN_SA_boot(
  lfq2.lm_07,  MA = lfq2.lm_07$MA,seasonalised = F, SA_temp = 10e+05,maxit = 500,
  low_par = list(Linf = 8, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 22, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 30/12)
print("Fi del segon")
Sys.time()

saveRDS(list(lfq_fin_lm_GA_07,lfq_fin_lm_SA_07), "Limnothrissa_07_bootstrapped.Rdata")

## OTOLITH DERIVED ESTIMATIONS ########################## ####
## 1. Load data ####
data_st_o = readxl::read_excel("Allen-method-Stolothrissa L.inf and K.xls", sheet = "R")
data_lm_o = readxl::read_excel("Allen-method-Limnothrissa L.inf and K.xls", sheet = "R")

data_st = list(age=data_st_o$Age/365, length= data_st_o$Length/10)
data_lm = list(age=data_lm_o$Age/365, length= data_lm_o$Length/10)

par(mfrow = c(2,1))
plot(data_st$length ~ data_st$age, main = "Stolothrissa")
plot(data_lm$length ~ data_lm$age, main = "Limnothrissa")
abline(10,0)

otolith_data_stolo = data.frame(Age = data_st$age, Length = data_st$length)
otolith_data_limno = data.frame(Age = data_lm$age, Length = data_lm$length)

write.csv(otolith_data_stolo, "Final_dataset_otolith_2004_Stolothrissa.csv")
write.csv(otolith_data_limno, "Final_dataset_otolith_2004_Limnothrissa.csv")

## 2. Fit data to a model that estimates Lmax and K from Otolith data using Least sum of squares method (LSM) #####

LSM_est_lm = growth_length_age(param = data_lm, method = "LSM", Linf_init = 500, CI = TRUE)
LSM_est_st = growth_length_age(param = data_st, method = "LSM", 
                               Linf_init = 50, CI = TRUE)
## SUMMARY OF ALL RESULTS AND FINAL TABLE ########################## ####

## The "get_values" function takes part of the code from the function
## univariate_density() from fishboot. This code collects the 
## most likely estimation of the value, and the CI for each value
## y = 1 is the Linf, y = 2, is K, and y = 3 is t0. 
## the first element of the list in val is GA, and the second is SA, both
## methods of ELEFAN boost

get_values = function(val,y){
  res = val$bootRaw
  x = ks::kde(res[,y])
  CItxt <- paste0(round(100-95), "%")
  inCI <- rle( x$estimate > x$cont[CItxt] )
  start.idx <- c(1, cumsum(inCI$lengths[-length(inCI$lengths)])+1)
  end.idx <- cumsum(inCI$lengths)
  limCI <- range(x$eval.points[start.idx[min(which(inCI$values))]:end.idx[max(which(inCI$values))]])
  in1 <- which(x$estimate > x$cont["99%"])
  mean1 <- mean(x$eval.points[in1])
  
  limCI = unlist(lapply(limCI, round,digits=2))
  mean1 = round(mean1,2)
  return(c(mean1, limCI[1], limCI[2]))
}

Stolo_88 = readRDS("Stolothrissa_digitalized_87_88_bootstrapped.Rdata")
Limno_88 = readRDS("Limnothrissa_digitalized_87_88_bootstrapped.Rdata")
Stolo_99 = readRDS("Stolothrissa_99_01_bootstrapped.Rdata")
Limno_99 = readRDS("Limnothrissa_99_01_bootstrapped.Rdata")
Stolo_07 = readRDS("Stolothrissa_07_bootstrapped.Rdata")
Limno_07 = readRDS("Limnothrissa_07_bootstrapped.Rdata")


stolo_table_87_GA = c(get_values(val = Stolo_88[[1]],y = 1),
                      get_values(val = Stolo_88[[1]],y = 2),
                      get_values(val = Stolo_88[[1]],y = 3))
stolo_table_87_SA = c(get_values(val = Stolo_88[[2]],y = 1),
                      get_values(val = Stolo_88[[2]],y = 2),
                      get_values(val = Stolo_88[[2]],y = 3))
limno_table_87_GA = c(get_values(val = Limno_88[[1]],y = 1),
                      get_values(val = Limno_88[[1]],y = 2),
                      get_values(val = Limno_88[[1]],y = 3))
limno_table_87_SA = c(get_values(val = Limno_88[[2]],y = 1),
                      get_values(val = Limno_88[[2]],y = 2),
                      get_values(val = Limno_88[[2]],y = 3))

stolo_table_99_GA = c(get_values(val = Stolo_99[[1]],y = 1),
                      get_values(val = Stolo_99[[1]],y = 2),
                      get_values(val = Stolo_99[[1]],y = 3))
stolo_table_99_SA = c(get_values(val = Stolo_99[[2]],y = 1),
                      get_values(val = Stolo_99[[2]],y = 2),
                      get_values(val = Stolo_99[[2]],y = 3))
limno_table_99_GA = c(get_values(val = Limno_99[[1]],y = 1),
                      get_values(val = Limno_99[[1]],y = 2),
                      get_values(val = Limno_99[[1]],y = 3))
limno_table_99_SA = c(get_values(val = Limno_99[[2]],y = 1),
                      get_values(val = Limno_99[[2]],y = 2),
                      get_values(val = Limno_99[[2]],y = 3))

stolo_table_07_GA = c(get_values(val = Stolo_07[[1]],y = 1),
                      get_values(val = Stolo_07[[1]],y = 2),
                      get_values(val = Stolo_07[[1]],y = 3))
stolo_table_07_SA = c(get_values(val = Stolo_07[[2]],y = 1),
                      get_values(val = Stolo_07[[2]],y = 2),
                      get_values(val = Stolo_07[[2]],y = 3))
limno_table_07_GA = c(get_values(val = Limno_07[[1]],y = 1),
                      get_values(val = Limno_07[[1]],y = 2),
                      get_values(val = Limno_07[[1]],y = 3))
limno_table_07_SA = c(get_values(val = Limno_07[[2]],y = 1),
                      get_values(val = Limno_07[[2]],y = 2),
                      get_values(val = Limno_07[[2]],y = 3))

results_together = rbind(stolo_table_87_GA,
                         stolo_table_87_SA,
                         limno_table_87_GA,
                         limno_table_87_SA,
                         stolo_table_99_GA,
                         stolo_table_99_SA,
                         limno_table_99_GA,
                         limno_table_99_SA,
                         stolo_table_07_GA,
                         stolo_table_07_SA,
                         limno_table_07_GA,
                         limno_table_07_SA)

# Now the rows of the Otolith analysis
otolith_st = c(LSM_est_st$estimates[1,2],LSM_est_st$estimates[1,3],LSM_est_st$estimates[1,4],
  LSM_est_st$estimates[2,2],LSM_est_st$estimates[2,3],LSM_est_st$estimates[2,4],
  LSM_est_st$estimates[3,2],LSM_est_st$estimates[3,3],LSM_est_st$estimates[3,4])

otolith_lm = c(LSM_est_lm$estimates[1,2],LSM_est_lm$estimates[1,3],LSM_est_lm$estimates[1,4],
  LSM_est_lm$estimates[2,2],LSM_est_lm$estimates[2,3],LSM_est_lm$estimates[2,4],
  LSM_est_lm$estimates[3,2],LSM_est_lm$estimates[3,3],LSM_est_lm$estimates[3,4])

#Bind everything together & save
results_together = rbind(results_together,otolith_st, otolith_lm)
colnames(results_together) = c("Lmax",	"CI_low",	"CI_high",
                               "K", "CI_low_K",	"CI_high_K",	
                               "t_anchor",	"CI_low_t",	"CI_high_t")

write.csv(results_together, "results_together.csv", row.names = F)
results = results_together
####################################### ####
############# FIGURES ################# ####
### -- PAPER PLOTS -- ### ####
### -- Figure 2: Histogram plots -- ### #####

# Stolothrissa
load("Stolothrissa-8788.Rdata")
stolothrissa9[4,] # That is already fine
stolothrissa9[20,] # This needs to be corrected
stolothrissa9[20,2:19] = 0
stolothrissa9[21,2:19] = 0
data_std88 = stolothrissa9 %>% mutate(lengthClass =(0.08729 + 1.1562 * lengthClass)*10)
data_st99 = readxl::read_excel("Stolothrissa_data_LFQ.xlsx", sheet = "Stolo_99_01")
data_st07 = readxl::read_excel("LFQ_stolo_2007_2008.xls", sheet = "Feuil3")
data_st_o = readxl::read_excel("Allen-method-Stolothrissa L.inf and K.xls", sheet = "R")
colnames(data_st07)[1] = "lengthClass"
ot_st = data_st_o %>% mutate(Length = ceiling(Length) - ceiling(Length) %% 5) %>% group_by(Length) %>% summarise(Freq = n())
his_o_st = data.frame(time = "Otolith data (04-05)",length = ot_st$Length, freq = ot_st$Freq)
his_st_88 = data.frame(time = "87-89",length = data_std88$lengthClass, freq = rowSums(data_std88[,2:ncol(data_std88)],na.rm = T))
his_st_99 = data.frame(time = "99-01",length = data_st99$lengthClass, freq = rowSums(data_st99[,2:ncol(data_st99)],na.rm = T))
his_st_07 = data.frame(time = "07-08",length = data_st07$lengthClass, freq = rowSums(data_st07[,2:ncol(data_st07)],na.rm = T))

hist_st = rbind(his_o_st,his_st_88, his_st_99, his_st_07)
hist_st$length = hist_st$length/10

# Limnothrissa 
load("Limnothrissa-8788.Rdata")
limnothrissa9[4,] # This needs to be corrected
limnothrissa9[4,2:19] = 0
limnothrissa9[24,] # This needs to be corrected
limnothrissa9[24,2:19] = 0

data_lmd88 = limnothrissa9 %>% mutate(lengthClass =  (0.16658 + 1.1873 * lengthClass)*10)

data_lm99 = readxl::read_excel("Limnothrissa_data_LFQ.xlsx", sheet = "99-01")
data_lm07 = readxl::read_excel("LFQ_limno_2007_2008.xls", sheet = "Feuil3")
data_lmot = readxl::read_excel("Allen-method-Limnothrissa L.inf and K.xls", sheet = "R")
colnames(data_lm07)[1] = "lengthClass"

ot_lm = data_lmot %>% mutate(Length = ceiling(Length) - ceiling(Length) %% 5) %>% group_by(Length) %>% summarise(Freq = n())
his_o_lm = data.frame(time = "Otolith data (04-05)",length = ot_lm$Length, freq = ot_lm$Freq)
his_lm_88 = data.frame(time = "87-89",length = data_lmd88$lengthClass, freq = rowSums(data_lmd88[,2:ncol(data_lmd88)],na.rm = T))
his_lm_99 = data.frame(time = "99-01",length = data_lm99$lengthClass, freq = rowSums(data_lm99[,2:ncol(data_lm99)],na.rm = T))
his_lm_07 = data.frame(time = "07-08",length = data_lm07$lengthClass, freq = rowSums(data_lm07[,2:ncol(data_lm07)],na.rm = T))

hist_lm = rbind(his_o_lm,his_lm_88,his_lm_99, his_lm_07)
hist_lm$length = hist_lm$length/10 

hist_lm$Sp = "L.miodon"
hist_st$Sp = "S.tanganicae"
hist = rbind(hist_lm,hist_st)

hist$time = factor(hist$time, levels = c("87-89", "99-01", "Otolith data (04-05)", "07-08"))
hist$Sp = factor(hist$Sp, levels = c("S.tanganicae", "L.miodon"))

levels(hist$time) <- c("1987-1989 (LFA)","1999-2001 (LFA)","2004-2005 (LAA)","2007-2008 (LFA)")
levels(hist$Sp) <- c("S. tanganicae", "L. miodon")

histy = hist %>% group_by(time, Sp) %>%
  mutate(total_freq = sum(freq),
         proportion = round(freq / total_freq,2),
         min_length = min(length[proportion > 0]),
         max_length = max(length[proportion > 0])) 

## Plot together
pdf("FIGURE_2.pdf", width = 12)
ggplot(histy, aes(x=length, y=proportion, fill = Sp)) +
  geom_bar(stat = "identity", width = 0.15) +
  facet_wrap(.~ time + Sp, scale = "free", nrow = 1) +
  xlab("Total length (cm)") +
  ylab("Number of individuals sampled") +
  scale_x_continuous(limits = c(0.1, 18)) +
  scale_y_reverse() +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") + 
  coord_flip() + 
  geom_vline(aes(xintercept=min_length, col=Sp), linetype="dashed") +
  geom_vline(aes(xintercept=max_length, col=Sp), linetype="dashed") +
  theme(strip.text = element_text(face = "italic"))
dev.off()

## Figure 3. In the paper is a plot with the VBGF curves of all de data
### -- Figure 3: VBGF curves -- ### #####

Stolo_88 = readRDS("Stolothrissa_digitalized_87_88_bootstrapped.Rdata")
Limno_88 = readRDS("Limnothrissa_digitalized_87_88_bootstrapped.Rdata")
Stolo_99 = readRDS("Stolothrissa_99_01_bootstrapped.Rdata")
Limno_99 = readRDS("Limnothrissa_99_01_bootstrapped.Rdata")
Stolo_07 = readRDS("Stolothrissa_07_bootstrapped.Rdata")
Limno_07 = readRDS("Limnothrissa_07_bootstrapped.Rdata")

days_s = c(1:548) # Days fpr Stolo
days = c(1:720) # Days for Limno

# S.tanganicae GA
a = data.frame(days = days_s,label = "A",Year = "1987-1989 (LFA)",algorithm = "S. tanganicae GA",
               length=VBGF(list(Linf=results[1,"Lmax"], K=results[1,"K"]/365, t0 = results[1,"t_anchor"]), t = days_s, L = NA, na.rm = FALSE),
               ymin = VBGF(list(Linf=results[1,"CI_low"], K=results[1,"CI_low_K"]/365, t0 = results[1,"CI_low_t"]), t = days_s, L = NA, na.rm = FALSE),
               ymax = VBGF(list(Linf=results[1,"CI_high"], K=results[1,"CI_high_K"]/365,t0 = results[1,"CI_high_t"]), t = days_s, L = NA, na.rm = FALSE))

a = rbind(a,data.frame(days = days_s,label = "A",Year = "1999-2001 (LFA)",algorithm = "S. tanganicae GA",
                       length=VBGF(list(Linf=results[5,"Lmax"], K=results[5,"K"]/365, t0=results[5,"t_anchor"]), t = days_s, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[5,"CI_low"], K=results[5,"CI_low_K"]/365, t0=results[5,"CI_low_t"]), t = days_s, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[5,"CI_high"], K=results[5,"CI_high_K"]/365,t0=results[5,"CI_high_t"]), t = days_s, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days_s,label = "A",Year = "2007-2008 (LFA)",algorithm = "S. tanganicae GA",
                       length=VBGF(list(Linf=results[9,"Lmax"], K=results[9,"K"]/365, t0=results[9,"t_anchor"]), t = days_s, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[9,"CI_low"], K=results[9,"CI_low_K"]/365, t0=results[9,"CI_low_t"]), t = days_s, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[9,"CI_high"], K=results[9,"CI_high_K"]/365,t0=results[9,"CI_high_t"]), t = days_s, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days_s,label = "A",Year = "2004-2005 (LAA)",algorithm = "S. tanganicae GA",
                       length=VBGF(list(Linf=results[13,"Lmax"], K=results[13,"K"]/365, t0 = results[10,"t_anchor"]), t = days_s, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[13,"CI_low"], K=results[13,"CI_low_K"]/365, t0 = results[10,"CI_low_t"]), t = days_s, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[13,"CI_high"], K=results[13,"CI_high_K"]/365, t0 = results[10,"CI_high_t"]), t = days_s, L = NA, na.rm = FALSE)))

# S. tanganicae SA
a = rbind(a,data.frame(days = days_s,label = "B",Year = "1987-1989 (LFA)",algorithm = "S. tanganicae SA",
                       length=VBGF(list(Linf=results[2,"Lmax"], K=results[2,"K"]/365, t0=results[2,"t_anchor"]), t = days_s, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[2,"CI_low"], K=results[2,"CI_low_K"]/365, t0=results[2,"CI_low_t"]), t = days_s, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[2,"CI_high"], K=results[2,"CI_high_K"]/365, t0=results[2,"CI_high_t"]), t = days_s, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days_s,label = "B",Year = "1999-2001 (LFA)",algorithm = "S. tanganicae SA",
                       length=VBGF(list(Linf=results[6,"Lmax"], K=results[6,"K"]/365, t0=results[6,"t_anchor"]), t = days_s, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[6,"CI_low"], K=results[6,"CI_low_K"]/365, t0=results[6,"CI_low_t"]), t = days_s, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[6,"CI_high"], K=results[6,"CI_high_K"]/365, t0=results[6,"CI_high_t"]), t = days_s, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days_s,label = "B",Year = "2007-2008 (LFA)",algorithm = "S. tanganicae SA",
                       length=VBGF(list(Linf=results[10,"Lmax"], K=results[10,"K"]/365, t0=results[10,"t_anchor"]), t = days_s, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[10,"CI_low"], K=results[10,"CI_low_K"]/365, t0=results[10,"CI_low_t"]), t = days_s, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[10,"CI_high"], K=results[10,"CI_high_K"]/365, t0=results[10,"CI_high_t"]), t = days_s, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days_s,label = "B",Year = "2004-2005 (LAA)",algorithm = "S. tanganicae SA",
                       length=VBGF(list(Linf=results[13,"Lmax"], K=results[13,"K"]/365, t0 = results[13,"t_anchor"]), t = days_s, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[13,"CI_low"], K=results[13,"CI_low_K"]/365, t0 = results[13,"CI_low_t"]), t = days_s, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[13,"CI_high"], K=results[13,"CI_high_K"]/365, t0 = results[13,"CI_high_t"]), t = days_s, L = NA, na.rm = FALSE)))

# L. miodon GA
a = rbind(a,data.frame(days = days,label = "C",Year = "1987-1989 (LFA)",algorithm = "L. miodon GA",
                       length=VBGF(list(Linf=results[3,"Lmax"], K=results[3,"K"]/365, t0=results[3,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[3,"CI_low"], K=results[3,"CI_low_K"]/365, t0=results[3,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[3,"CI_high"], K=results[3,"CI_high_K"]/365, t0=results[3,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days,label = "C",Year = "1999-2001 (LFA)",algorithm = "L. miodon GA",
                       length=VBGF(list(Linf=results[7,"Lmax"], K=results[7,"K"]/365, t0=results[7,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[7,"CI_low"], K=results[7,"CI_low_K"]/365, t0=results[7,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[7,"CI_high"], K=results[7,"CI_high_K"]/365, t0=results[7,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days,label = "C",Year = "2007-2008 (LFA)",algorithm = "L. miodon GA",
                       length=VBGF(list(Linf=results[11,"Lmax"], K=results[11,"K"]/365, t0=results[11,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[11,"CI_low"], K=results[11,"CI_low_K"]/365, t0=results[11,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[11,"CI_high"], K=results[11,"CI_high_K"]/365, t0=results[11,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days,label = "C",Year = "2004-2005 (LAA)",algorithm = "L. miodon GA",
                       length=VBGF(list(Linf=results[14,"Lmax"], K=results[14,"K"]/365, t0 = results[10,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[14,"CI_low"], K=results[14,"CI_low_K"]/365, t0 = results[10,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[14,"CI_high"], K=results[14,"CI_high_K"]/365, t0 = results[10,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

# L. miodon SA
a = rbind(a,data.frame(days = days,label = "D",Year = "1987-1989 (LFA)",algorithm = "L. miodon SA",
                       length=VBGF(list(Linf=results[4,"Lmax"], K=results[4,5]/365, t0=results[4,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[4,"CI_low"], K=results[4,"CI_low_K"]/365, t0=results[4,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[4,"CI_high"], K=results[4,"CI_high_K"]/365,t0=results[4,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days,label = "D",Year = "1999-2001 (LFA)",algorithm = "L. miodon SA",
                       length=VBGF(list(Linf=results[8,"Lmax"], K=results[8,"K"]/365, t0=results[8,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[8,"CI_low"], K=results[8,"CI_low_K"]/365, t0=results[8,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[8,"CI_high"], K=results[8,"CI_high_K"]/365,t0=results[8,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days,label = "D",Year = "2007-2008 (LFA)",algorithm = "L. miodon SA",
                       length=VBGF(list(Linf=results[12,"Lmax"], K=results[12,"K"]/365, t0=results[8,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[12,"CI_low"], K=results[12,"CI_low_K"]/365, t0=results[9,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[12,"CI_high"], K=results[12,"CI_high_K"]/365,t0=results[10,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days,label = "D",Year = "2004-2005 (LAA)",algorithm = "L. miodon SA",
                       length=VBGF(list(Linf=results[14,"Lmax"], K=results[14,"K"]/365, t0 = results[14,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[14,"CI_low"], K=results[14,"CI_low_K"]/365, t0 = results[14,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[14,"CI_high"], K=results[14,"CI_high_K"]/365, t0 = results[14,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

a$algorithm = factor(a$algorithm)
a$algorithm = factor(a$algorithm, levels = c("S. tanganicae GA", 
                                             "S. tanganicae SA",
                                             "L. miodon GA",
                                             "L. miodon SA"))
# Plot with all lines
pdf("FIGURE_3.pdf", width = 7, height = 5)
ggplot(a) + 
  geom_line(aes(x = days, y = length, col = Year)) + 
  geom_ribbon(aes(x = days, ymax = ymax, ymin = ymin, fill = Year), alpha = 0.3) + 
  facet_wrap(.~algorithm, scale = "free") +
  ylab("Total length (cm)") + xlab("Time (days)") +
  theme_bw() +
  geom_text(data = a,
            mapping = aes(x = 30,
                          y = 12,
                          label = label)) +
  theme(strip.text = element_text(face = "italic"))

dev.off()

### -- Figure 4: Plots of bins and lines -- ### #####

# With raw data
pdf("FIGURE_4.pdf", height = 8, width = 9)
par(mfrow = c(3,2), mar = c(4,4, 5, 3))

# 1988 S. tanganicae
# raw data
plot(lfq2.std, "catch", main = substitute(paste(italic('S. tanganicae'))), ylab="Total length (cm)", ylim=c(0,11))
abline(lfqFitCurves(lfq2.std, par = list(Linf = results[1,"Lmax"], K = results[1,"K"], t_anchor = results[1,"t_anchor"]), col = "red", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.std, par = list(Linf = results[2,"Lmax"], K = results[2,"K"], t_anchor = results[2,"t_anchor"]), col = "blue", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.std, par = list(Linf = results[13,"Lmax"], K = results[13,"K"], t_anchor = 0.003), col = "black", draw = T,lwd = 1.2))

plot(lfq2.lmd, "catch", main = substitute(paste(italic('L. miodon'))), ylab="Total length (cm)", ylim=c(0,14.5))
abline(lfqFitCurves(lfq2.lmd, par = list(Linf = results[3,"Lmax"], K = results[3,"K"], t_anchor = results[3,"t_anchor"]), col = "red", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.lmd, par = list(Linf = results[4,"Lmax"], K = results[4,"K"], t_anchor = results[4,"t_anchor"]), col = "blue", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.lmd, par = list(Linf = results[14,"Lmax"], K = results[14,"K"], t_anchor = 0.003), col = "black", draw = T, lwd = 1.2))

### 1999
# S. tanganicae
# raw data
plot(lfq2.st2, "catch", ylab="Total length (cm)", ylim=c(0,11))
abline(lfqFitCurves(lfq2.st2, par = list(Linf = results[5,"Lmax"], K = results[5,"K"], t_anchor = results[5,"t_anchor"]), col = "red", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.st2, par = list(Linf = results[6,"Lmax"], K = results[6,"K"], t_anchor = results[6,"t_anchor"]), col = "blue", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.st2, par = list(Linf = results[13,"Lmax"], K = results[13,"K"], t_anchor = 0.003), col = "black", draw = T, lwd = 1.2))

# L. miodon
# With raw data
plot(lfq2.lm_99, "catch", ylab="Total length (cm)", ylim=c(0,14.5))
abline(lfqFitCurves(lfq2.lm_99, par = list(Linf = results[7,"Lmax"], K = results[7,"K"], t_anchor = results[7,"t_anchor"]), col = "red", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.lm_99, par = list(Linf = results[8,"Lmax"], K = results[8,"K"], t_anchor = results[8,"t_anchor"]), col = "blue", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.lm_99, par = list(Linf = results[14,"Lmax"], K = results[14,"K"], t_anchor = 0.003), col = "black", draw = T, lwd = 1.2))

## 2007 
# S. tanganicae
# raw data
plot(lfq2.st_07, "catch", ylab="Total length (cm)", ylim=c(0,11))
abline(lfqFitCurves(lfq2.st_07, par = list(Linf = results[9,"Lmax"], K = results[9,"K"], t_anchor = results[9,"t_anchor"]), col = "red", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.st_07, par = list(Linf = results[10,"Lmax"], K = results[10,"K"], t_anchor = results[10,"t_anchor"]), col = "blue", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.st_07, par = list(Linf = results[13,"Lmax"], K = results[13,"K"], t_anchor = 0.003), col = "black", draw = T,lwd = 1.2))

# L. miodon
# With raw data
plot(lfq2.lm_07, "catch", ylab="Total length (cm)", ylim=c(0,14.5))
abline(lfqFitCurves(lfq2.lm_07, par = list(Linf = results[11,"Lmax"], K = results[11,"K"], t_anchor = results[11,"t_anchor"]), col = "red", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.lm_07, par = list(Linf = results[12,"Lmax"], K = results[12,"K"], t_anchor = results[12,"t_anchor"]), col = "blue", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.lm_07, par = list(Linf = results[14,"Lmax"], K = results[14,"K"], t_anchor = 0.003), col = "black", draw = T, lwd = 1.2))

dev.off()

### -- Figure S1: Density plots -- ### ####

Stolo_88 = readRDS("Stolothrissa_digitalized_87_88_bootstrapped.Rdata")
Limno_88 = readRDS("Limnothrissa_digitalized_87_88_bootstrapped.Rdata")
Stolo_99 = readRDS("Stolothrissa_99_01_bootstrapped.Rdata")
Limno_99 = readRDS("Limnothrissa_99_01_bootstrapped.Rdata")
Stolo_07 = readRDS("Stolothrissa_07_bootstrapped.Rdata")
Limno_07 = readRDS("Limnothrissa_07_bootstrapped.Rdata")

par(mfrow = c(5,5))
pdf("FIGURE_S1.pdf", width = 3, height = 2)
univariate_density(Stolo_88[[1]], cex = 0.5)
univariate_density(Stolo_88[[2]], cex = 0.5)
univariate_density(Limno_88[[1]], cex = 0.5)
univariate_density(Limno_88[[2]], cex = 0.5)

univariate_density(Stolo_99[[1]], cex = 0.5)
univariate_density(Stolo_99[[2]], cex = 0.5)
univariate_density(Limno_99[[1]], cex = 0.5)
univariate_density(Limno_99[[2]], cex = 0.5)

univariate_density(Stolo_07[[1]], cex = 0.5)
univariate_density(Stolo_07[[2]], cex = 0.5)
univariate_density(Limno_07[[1]], cex = 0.5)
univariate_density(Limno_07[[2]], cex = 0.5)
dev.off()

### -- Figure S2: Curves -- ### ######

ot_par = list(Linf = 10.64, K = 3.18, t0 = 0.03) # This I make up, does not change the shape
st_oto_low.conf = list(Linf = 10.32, K = 2.76, t0 = 0.01)
st_oto_up.conf = list(Linf = 11.05, K = 3.63, t0 = 0.05)

lm_par = list(Linf = 12.82, K = 1.73, t0 = 0.02) # This I make up, does not change the shape
lm_oto_low.conf = list(Linf = 12.18, K = 1.47, t0 = 0.01)
lm_oto_up.conf = list(Linf = 13.64, K = 2.01, t0 = 0.05)



data_st_pred_ot = data.frame(age = LSM_est_st$age,
                             length = VBGF(param = ot_par, t = LSM_est_st$age),
                             ymin = VBGF(param = st_oto_low.conf, t = LSM_est_st$age),
                             ymax = VBGF(param = st_oto_up.conf, t = LSM_est_st$age))

data_lm_pred_ot = data.frame(age = LSM_est_lm$age,
                             length = VBGF(param = lm_par, t = LSM_est_lm$age),
                             ymin = VBGF(param = lm_oto_low.conf, t = LSM_est_lm$age),
                             ymax = VBGF(param = lm_oto_up.conf, t = LSM_est_lm$age))


figure_s2 = ggplot() + 
  geom_point(aes(x = LSM_est_st$age, y = LSM_est_st$length), col = "#8B0000", alpha = 0.4)+
  
  geom_line(data = data_st_pred_ot, aes(x = age, y = length), col = "#8B0000") + 
  geom_ribbon(data = data_st_pred_ot, aes(x = age, y = length, 
                                          ymax = ymax, ymin = ymin), fill = "#8B0000", alpha = 0.2)+
  geom_point(aes(x = LSM_est_lm$age, y = LSM_est_lm$length), col = "darkblue", alpha = 0.4)+
  geom_line(data = data_lm_pred_ot, aes(x = age, y = length), col = "darkblue") + 
  geom_ribbon(data = data_lm_pred_ot, aes(x = age, y = length, 
                                          ymax = ymax, ymin = ymin), fill = "darkblue", alpha = 0.1) +
  
  ylab("Total length (cm)") + xlab("Age (years)") + theme(legend.position="top")

pdf("FIGURE_S2.pdf", height = 5, width = 6)
figure_s2
dev.off()
