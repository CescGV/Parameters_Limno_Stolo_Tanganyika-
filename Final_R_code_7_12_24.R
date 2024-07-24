## Set working directory & library ####
library(fishboot)
library(TropFishR)
library(tidyverse)
library(gridExtra)
library(pdftools)
library(magick)
library(grid)

wd = "C:/Users/06061016/OneDrive - Nord universitet/Documents"
set.seed(1)
setwd(wd)

## LFQ ESTIMATION 1987 - 1989 ########################### ######
## 1. Load data ####
stolo_87_raw = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-1987-1989-LFA", col_names = T)
names(stolo_87_raw)[str_detect(names(stolo_87_raw), "\\d{5}")] <- format(as.Date(as.numeric(names(stolo_87_raw)[str_detect(names(stolo_87_raw), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")

limno_87_raw = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Limno-1987-1989-LFA", col_names = T)
names(limno_87_raw)[str_detect(names(limno_87_raw), "\\d{5}")] <- format(as.Date(as.numeric(names(limno_87_raw)[str_detect(names(limno_87_raw), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")

## 2. Create LFQ objects ####

data_std_87 <- list(dates = as.Date(names(stolo_87_raw)[str_detect(names(stolo_87_raw), "-")]),
                     midLengths = stolo_87_raw$Length/10,
                     catch = as.matrix(stolo_87_raw[,2:ncol(stolo_87_raw)]))

data_lmd_87 <- list(dates = as.Date(names(limno_87_raw)[str_detect(names(limno_87_raw), "-")]),
                     midLengths = limno_87_raw$Length/10,
                     catch = as.matrix(limno_87_raw[,2:ncol(limno_87_raw)]))

class(data_std_87) <- "lfq"
class(data_lmd_87) <- "lfq"
lfq0.std1 <- data_std_87
lfq1.std1 <- lfqModify(lfq0.std1, bin_size = 0.8) # Stolo
lfq0.lmd1 <- data_lmd_87
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

lfq2.std <- lfq2b.std1
lfq2.lmd <- lfq2b.lmd1

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

## LFQ ESTIMATION 1999 - 2001 ########################### ####
## 1. Load data ####

stolo_99_raw = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-1999-2001-LFA", col_names = T)
names(stolo_99_raw)[str_detect(names(stolo_99_raw), "\\d{5}")] <- format(as.Date(as.numeric(names(stolo_99_raw)[str_detect(names(stolo_99_raw), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")
stolo_99_raw = stolo_99_raw[,!colSums(is.na(stolo_99_raw)) == nrow(stolo_99_raw)]

limno_99_raw = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Limno-1999-2001-LFA", col_names = T)
names(limno_99_raw)[str_detect(names(limno_99_raw), "\\d{5}")] <- format(as.Date(as.numeric(names(limno_99_raw)[str_detect(names(limno_99_raw), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")
limno_99_raw = limno_99_raw[,!colSums(is.na(limno_99_raw)) == nrow(limno_99_raw)]

## 2. Create LFQ object ####

data_st_99 <- list(dates = as.Date(names(stolo_99_raw)[str_detect(names(stolo_99_raw), "-")]),
                   midLengths = stolo_99_raw$Length/10,
                   catch = as.matrix(stolo_99_raw[,2:ncol(stolo_99_raw)]))

data_lm_99 <- list(dates = as.Date(names(limno_99_raw)[str_detect(names(limno_99_raw), "-")]),
                   midLengths = limno_99_raw$Length/10,
                   catch = as.matrix(limno_99_raw[,2:ncol(limno_99_raw)]))

class(data_st_99) <- "lfq"
class(data_lm_99) <- "lfq"
lfq0.st_99 <- data_st_99
lfq1.st_99 <- lfqModify(lfq0.st_99, bin_size = 0.8) # Stolo
lfq0.lm_99 <- data_lm_99
lfq1.lm_99 <- lfqModify(lfq0.lm_99, bin_size = 0.8) # Limno

par(mfrow = c(2,2))
plot(lfq0.st_99, Fname = "catch", main = "Stolothrissa 99-01")
plot(lfq1.st_99, Fname = "catch", main = "Stolothrissa 99-01")
plot(lfq0.lm_99, Fname = "catch", main = "Limnothrissa 99-01")
plot(lfq1.lm_99, Fname = "catch", main = "Limnothrissa 99-01")

# Assign a moving average.
ma <- 7 # See youngest cohort
lfq2a.st_99 <- lfqRestructure(lfq0.st_99, MA = ma, addl.sqrt = FALSE)
lfq2b.st_99 <- lfqRestructure(lfq1.st_99, MA = ma, addl.sqrt = FALSE)
ma <- 7 # See youngest cohort
lfq2a.lm_99 <- lfqRestructure(lfq0.lm_99, MA = ma, addl.sqrt = FALSE)
lfq2b.lm_99 <- lfqRestructure(lfq1.lm_99, MA = ma, addl.sqrt = FALSE)
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))

plot(lfq2a.st_99, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 99-01")
plot(lfq2b.st_99, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 99-01")
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2a.lm_99, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 99-01")
plot(lfq2b.lm_99, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 99-01")

lfq2.st_99 <- lfq2b.st_99
lfq2.lm_99 <- lfq2b.lm_99

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2.st_99, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 88")
plot(lfq2.lm_99, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 88")

## 3. Estimate parameters ####
## Stolothrissa
set.seed(1)
Sys.time()
lfq_GA.st_99 <- fishboot::ELEFAN_GA_boot(
  lfq2.st_99,  MA = lfq2.st_99$MA, seasonalised = F, popSize=100, pmutation = 0.5,
  low_par = list(Linf = 6, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 12, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 18/12)
print("Fi del primer")

Sys.time()
lfq_SA.st_99 <- fishboot::ELEFAN_SA_boot(
  lfq2.st_99,  MA = lfq2.st_99$MA,seasonalised = F, SA_temp = 10e+05,maxit = 500,
  low_par = list(Linf = 6, K = 1, t_anchor = 0, C = 0, ts = 0),
  up_par  = list(Linf = 12, K = 10, t_anchor = 1, C = 1, ts = 1),
  nresamp = 400, agemax = 18/12)
Sys.time()
print("Fi del segon")
Sys.time()
saveRDS(list(lfq_GA.std_99, lfq_SA.std_99), "Stolothrissa_digitalized_99_01_bootstrapped.Rdata")

## Limnothrissa
set.seed(1)
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
saveRDS(list(lfq_GA.lm_99, lfq_SA.lm_99), "Limnothrissa_digitalized_99_01_bootstrapped.Rdata")

## LFQ ESTIMATION 2007 - 2008 ########################### ####
## 1. Load data ####
stolo_07_raw = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-2007-2008-LFA", col_names = T)
names(stolo_07_raw)[str_detect(names(stolo_07_raw), "\\d{5}")] <- format(as.Date(as.numeric(names(stolo_07_raw)[str_detect(names(stolo_07_raw), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")
stolo_07_raw = stolo_07_raw[,!colSums(is.na(stolo_07_raw)) == nrow(stolo_07_raw)]

limno_07_raw = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Limno-2007-2008-LFA", col_names = T)
names(limno_07_raw)[str_detect(names(limno_07_raw), "\\d{5}")] <- format(as.Date(as.numeric(names(limno_07_raw)[str_detect(names(limno_07_raw), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")
limno_07_raw = limno_07_raw[,!colSums(is.na(limno_07_raw)) == nrow(limno_07_raw)]

## 2. Create LFQ object ####

data_st_07 <- list(dates = as.Date(names(stolo_07_raw)[str_detect(names(stolo_07_raw), "-")]),
                   midLengths = stolo_07_raw$Length/10,
                   catch = as.matrix(stolo_07_raw[,2:ncol(stolo_07_raw)]))

data_lm_07 <- list(dates = as.Date(names(limno_07_raw)[str_detect(names(limno_07_raw), "-")]),
                   midLengths = limno_07_raw$Length/10,
                   catch = as.matrix(limno_07_raw[,2:ncol(limno_07_raw)]))

class(data_st_07) <- "lfq"
class(data_lm_07) <- "lfq"
lfq0.st_07 <- data_st_07
lfq1.st_07 <- lfqModify(data_st_07, bin_size = 0.8) # Stolo
lfq0.lm_07 <- data_lm_07
lfq1.lm_07 <- lfqModify(data_lm_07, bin_size = 0.8) # Limno

par(mfrow = c(2,2))
plot(lfq0.st_07, Fname = "catch", main = "Stolothrissa 07-08")
plot(lfq1.st_07, Fname = "catch", main = "Stolothrissa 07-08")
plot(lfq0.lm_07, Fname = "catch", main = "Limnothrissa 07-08")
plot(lfq1.lm_07, Fname = "catch", main = "Limnothrissa 07-08")

# Assign a moving average.
ma <- 5 # See youngest cohort
lfq2a.st_07 <- lfqRestructure(lfq0.st_07, MA = ma, addl.sqrt = FALSE)
lfq2b.st_07 <- lfqRestructure(lfq1.st_07, MA = ma, addl.sqrt = FALSE)
ma <- 7 # See youngest cohort
lfq2a.lm_07 <- lfqRestructure(lfq0.lm_07, MA = ma, addl.sqrt = FALSE)
lfq2b.lm_07 <- lfqRestructure(lfq1.lm_07, MA = ma, addl.sqrt = FALSE)
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))

plot(lfq2a.st_07, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 99-01")
plot(lfq2b.st_07, Fname = "rcounts", date.axis = "modern", main = "Stolothrissa 99-01")
opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2a.lm_07, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 99-01")
plot(lfq2b.lm_07, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 99-01")

lfq2.st_07 <- lfq2b.st_07
lfq2.lm_07 <- lfq2b.lm_07

opar <- par(mfrow = c(2,1), mar = c(2,5,2,3), oma = c(2,0,0,0))
plot(lfq2.st_07, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 88")
plot(lfq2.lm_07, Fname = "rcounts", date.axis = "modern", main = "Limnothrissa 88")

## 3. Estimate parameters ####

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

## Limnothrissa
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
data_st_o = readxl::read_excel("Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-2004-2005-LAA")
data_lm_o = readxl::read_excel("Mulimbwa et al 2024 Data.xlsx", sheet = "Limno-2004-2005-LAA")
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

## 2. Fit data  #####

LSM_est_lm = growth_length_age(param = data_lm, method = "LSM", Linf_init = 500, CI = TRUE)
LSM_est_st = growth_length_age(param = data_st, method = "LSM", 
                               Linf_init = 50, CI = TRUE)
## RESULTS SUMMARY AND FINAL TABLE ###################### ####

## The "get_values" function takes part of the code from the function
## univariate_density() from fishboot. This code collects the 
## most likely estimation of the value, and the CI for each value
## y = 1 is the Linf, y = 2, is K, and y = 3 is t0. 
## the first element of the list in val is GA, and the second is SA, both
## methods of ELEFAN boost
setwd("C:/Users/06061016/OneDrive - Nord universitet/Desktop/Mulimbwa/Final_Mulimbwa_7_10_24")
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

write.csv(results_together, "results_together2.csv", row.names = T)
results = results_together
####################################### ####
############# FIGURES ################# ####
### -- PAPER PLOTS -- ### ####
### -- Figure 2: Histogram plots -- ### #####

# Stolothrissa

data_std88 = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-1987-1989-LFA", col_names = T)
names(data_std88)[str_detect(names(data_std88), "\\d{5}")] <- format(as.Date(as.numeric(names(data_std88)[str_detect(names(data_std88), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")

data_st99 = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-1999-2001-LFA", col_names = T)
names(data_st99)[str_detect(names(data_st99), "\\d{5}")] <- format(as.Date(as.numeric(names(data_st99)[str_detect(names(data_st99), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")
data_st99 = data_st99[,!colSums(is.na(data_st99)) == nrow(data_st99)]

data_st07 = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-2007-2008-LFA", col_names = T)
names(data_st07)[str_detect(names(data_st07), "\\d{5}")] <- format(as.Date(as.numeric(names(data_st07)[str_detect(names(data_st07), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")
data_st07 = data_st07[,!colSums(is.na(data_st07)) == nrow(data_st07)]

data_st_o = readxl::read_excel("Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-2004-2005-LAA")

ot_st = data_st_o %>% mutate(Length = ceiling(Length) - ceiling(Length) %% 5) %>% group_by(Length) %>% summarise(Freq = n())
his_o_st = data.frame(time = "Otolith data (04-05)",length = ot_st$Length, freq = ot_st$Freq)
his_st_88 = data.frame(time = "87-89",length = data_std88$Length, freq = rowSums(data_std88[,2:ncol(data_std88)],na.rm = T))
his_st_99 = data.frame(time = "99-01",length = data_st99$Length, freq = rowSums(data_st99[,2:ncol(data_st99)],na.rm = T))
his_st_07 = data.frame(time = "07-08",length = data_st07$Length, freq = rowSums(data_st07[,2:ncol(data_st07)],na.rm = T))

hist_st = rbind(his_o_st,his_st_88, his_st_99, his_st_07)
hist_st$length = hist_st$length/10

# Limnothrissa 


data_lmd88 = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Limno-1987-1989-LFA", col_names = T)
names(data_lmd88)[str_detect(names(data_lmd88), "\\d{5}")] <- format(as.Date(as.numeric(names(data_lmd88)[str_detect(names(data_lmd88), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")

data_lm99 = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-1999-2001-LFA", col_names = T)
names(data_lm99)[str_detect(names(data_lm99), "\\d{5}")] <- format(as.Date(as.numeric(names(data_lm99)[str_detect(names(data_lm99), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")
data_lm99 = data_lm99[,!colSums(is.na(data_lm99)) == nrow(data_lm99)]

data_lm07 = read_excel(path = "Mulimbwa et al 2024 Data.xlsx", sheet = "Stolo-2007-2008-LFA", col_names = T)
names(data_lm07)[str_detect(names(data_lm07), "\\d{5}")] <- format(as.Date(as.numeric(names(data_lm07)[str_detect(names(data_lm07), "\\d{5}")]), origin = "1899-12-30"), "%Y-%m-%d")
data_lm07 = data_lm07[,!colSums(is.na(data_lm07)) == nrow(data_lm07)]

data_lm_o = readxl::read_excel("Mulimbwa et al 2024 Data.xlsx", sheet = "Limno-2004-2005-LAA")

ot_lm = data_lm_o %>% mutate(Length = ceiling(Length) - ceiling(Length) %% 5) %>% group_by(Length) %>% summarise(Freq = n())
his_o_lm = data.frame(time = "Otolith data (04-05)",length = ot_lm$Length, freq = ot_lm$Freq)
his_lm_88 = data.frame(time = "87-89",length = data_lmd88$Length, freq = rowSums(data_lmd88[,2:ncol(data_lmd88)],na.rm = T))
his_lm_99 = data.frame(time = "99-01",length = data_lm99$Length, freq = rowSums(data_lm99[,2:ncol(data_lm99)],na.rm = T))
his_lm_07 = data.frame(time = "07-08",length = data_lm07$Length, freq = rowSums(data_lm07[,2:ncol(data_lm07)],na.rm = T))

hist_lm = rbind(his_o_lm,his_lm_88, his_lm_99, his_lm_07)
hist_lm$length = hist_lm$length/10

# Join both species with a new column to distinguish them

hist_lm$Sp = "L.miodon"
hist_st$Sp = "S.tanganicae"
hist = rbind(hist_lm,hist_st)

hist$time = factor(hist$time, levels = c("87-89", "99-01", "Otolith data (04-05)", "07-08"))
hist$Sp = factor(hist$Sp, levels = c("S.tanganicae", "L.miodon"))

levels(hist$time) <- c("1987-1989 (LFA)","1999-2001 (LFA)","2004-2005 (LAA)","2007-2008 (LFA)")
levels(hist$Sp) <- c("S. tanganicae", "L. miodon")

histy1 = hist %>% group_by(time, Sp) %>%
  mutate(total_freq = sum(freq),
         proportion = round(freq / total_freq,2),
         min_length = min(length[proportion > 0]),
         max_length = max(length[proportion > 0])) 

histy2 = hist %>% group_by(time, Sp) %>%
  mutate(total_freq = sum(freq),
         proportion = freq / total_freq,
         min_length = min(length[proportion > 0]),
         max_length = max(length[proportion > 0])) 

histy <- histy2

## Plot together
pdf("Mulimbwa et al 2024 Figure 2_ULT2.pdf", width = 14)
ggplot(histy, aes(x=length, y=proportion, fill = Sp)) +
  geom_bar(stat = "identity", width = 0.15) +
  facet_wrap(.~ time + Sp, scale = "free", nrow = 1) +
  xlab("Total length (cm)") +
  ylab("Proportion of individuals sampled") +
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
                       length=VBGF(list(Linf=results[13,"Lmax"], K=results[13,"K"]/365, t0 = results[13,"t_anchor"]), t = days_s, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[13,"CI_low"], K=results[13,"CI_low_K"]/365, t0 = results[13,"CI_low_t"]), t = days_s, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[13,"CI_high"], K=results[13,"CI_high_K"]/365, t0 = results[13,"CI_high_t"]), t = days_s, L = NA, na.rm = FALSE)))

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
                       length=VBGF(list(Linf=results[14,"Lmax"], K=results[14,"K"]/365, t0 = results[14,"t_anchor"]), t = days, L = NA, na.rm = FALSE), # corrected t0 from row 10 to row 14
                       ymin = VBGF(list(Linf=results[14,"CI_low"], K=results[14,"CI_low_K"]/365, t0 = results[14,"CI_low_t"]), t = days, L = NA, na.rm = FALSE), # corrected t0 from row 10 to row 14
                       ymax = VBGF(list(Linf=results[14,"CI_high"], K=results[14,"CI_high_K"]/365, t0 = results[14,"CI_high_t"]), t = days, L = NA, na.rm = FALSE))) # corrected t0 from row 10 to row 14

# L. miodon SA
a = rbind(a,data.frame(days = days,label = "D",Year = "1987-1989 (LFA)",algorithm = "L. miodon SA",
                       length=VBGF(list(Linf=results[4,"Lmax"], K=results[4,"K"]/365, t0=results[4,"t_anchor"]), t = days, L = NA, na.rm = FALSE), # corrected K
                       ymin = VBGF(list(Linf=results[4,"CI_low"], K=results[4,"CI_low_K"]/365, t0=results[4,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[4,"CI_high"], K=results[4,"CI_high_K"]/365,t0=results[4,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days,label = "D",Year = "1999-2001 (LFA)",algorithm = "L. miodon SA",
                       length=VBGF(list(Linf=results[8,"Lmax"], K=results[8,"K"]/365, t0=results[8,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[8,"CI_low"], K=results[8,"CI_low_K"]/365, t0=results[8,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[8,"CI_high"], K=results[8,"CI_high_K"]/365,t0=results[8,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

a = rbind(a,data.frame(days = days,label = "D",Year = "2007-2008 (LFA)",algorithm = "L. miodon SA",
                       length=VBGF(list(Linf=results[12,"Lmax"], K=results[12,"K"]/365, t0=results[12,"t_anchor"]), t = days, L = NA, na.rm = FALSE),
                       ymin = VBGF(list(Linf=results[12,"CI_low"], K=results[12,"CI_low_K"]/365, t0=results[12,"CI_low_t"]), t = days, L = NA, na.rm = FALSE),
                       ymax = VBGF(list(Linf=results[12,"CI_high"], K=results[12,"CI_high_K"]/365,t0=results[12,"CI_high_t"]), t = days, L = NA, na.rm = FALSE)))

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
pdf("Mulimbwa et al 2024 Figure 3_ULT.pdf", width = 7, height = 5)
ggplot(a) + 
  geom_line(aes(x = days, y = length, col = Year)) + 
  geom_ribbon(aes(x = days, ymax = ymax, ymin = ymin, fill = Year), alpha = 0.3) + 
  facet_wrap(.~algorithm, scale = "free") +
  ylab("Total length (cm)") + xlab("Time (days)") +
  theme_bw() +
  geom_text(data = a,
            mapping = aes(x = ifelse(algorithm %in% c("S. tanganicae GA","S. tanganicae SA"),24,30),
                          y = ifelse(algorithm %in% c("S. tanganicae GA","S. tanganicae SA"),10.5,12.5),
                          label = label)) +
  theme(strip.text = element_text(face = "italic"))

dev.off()

### -- Figure 4: Plots of bins and lines -- ### #####

# With raw data
pdf("Mulimbwa et al 2024 Figure 4_ULT2.pdf", height = 7.5, width = 9)
par(mfrow = c(3,2), mar = c(4, 4, 3, 3))

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
plot(lfq2.st_99, "catch", ylab="Total length (cm)", ylim=c(0,11))
abline(lfqFitCurves(lfq2.st_99, par = list(Linf = results[5,"Lmax"], K = results[5,"K"], t_anchor = results[5,"t_anchor"]), col = "red", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.st_99, par = list(Linf = results[6,"Lmax"], K = results[6,"K"], t_anchor = results[6,"t_anchor"]), col = "blue", draw = T, lwd = 1.2))
abline(lfqFitCurves(lfq2.st_99, par = list(Linf = results[13,"Lmax"], K = results[13,"K"], t_anchor = 0.003), col = "black", draw = T, lwd = 1.2))

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

pdf("Mulimbwa et al 2024 Figure S2_ULT.pdf", width = 6, height = 4)
univariate_density(Stolo_88[[1]], cex = 1)
univariate_density(Stolo_88[[2]], cex = 1)
univariate_density(Limno_88[[1]], cex = 1)
univariate_density(Limno_88[[2]], cex = 1)

univariate_density(Stolo_99[[1]], cex = 1)
univariate_density(Stolo_99[[2]], cex = 1)
univariate_density(Limno_99[[1]], cex = 1)
univariate_density(Limno_99[[2]], cex = 1)

univariate_density(Stolo_07[[1]], cex = 1)
univariate_density(Stolo_07[[2]], cex = 1)
univariate_density(Limno_07[[1]], cex = 1)
univariate_density(Limno_07[[2]], cex = 1)
dev.off()

# Path to your PDF file
pdf_path <- "C:/Users/06061016/OneDrive - Nord universitet/Desktop/Mulimbwa/Final_Mulimbwa_7_10_24/Mulimbwa et al 2024 Figure S2_ULT.pdf"

# Read the PDF file
pdf <- image_read_pdf(pdf_path)

# Convert each PDF page to an image
images <- lapply(seq_along(pdf), function(i) {
  image <- image_convert(pdf[i], format = "png")
  image
})

# Convert magick images to grobs for grid.arrange
image_grobs <- lapply(images, function(img) {
  rasterGrob(img, interpolate = TRUE)
})

pdf("Mulimbwa et al 2024 Figure S2_ARRANGED.pdf", width = 8.27, height = 11.69)
grid.arrange(grobs = image_grobs, ncol = 2, nrow = 6)
dev.off()

### -- Figure S2: Curves -- ### ######

ot_par = list(Linf = 10.64, K = 3.18, t0 = 0.03) # This I make up, does not change the shape
st_oto_low.conf = list(Linf = 10.32, K = 2.76, t0 = 0.01)
st_oto_up.conf = list(Linf = 11.05, K = 3.63, t0 = 0.05)

lm_par = list(Linf = 12.82, K = 1.73, t0 = 0.02) # This I make up, does not change the shape
lm_oto_low.conf = list(Linf = 12.18, K = 1.47, t0 = -0.01)
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

pdf("Mulimbwa et al 2024 Figure S3_ULT.pdf", height = 5, width = 6)
figure_s2
dev.off()
