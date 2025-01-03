
# Load packages
library(tidyverse)
library(data.table)
library(cmdstanr)
library(posterior)
set_cmdstan_path("~/cmdstan/")
library(lubridate)

# Create the output directory
var_name <- c("H3K4me1")
out_99 <- paste0("01_STM_forloop/99/", var_name, "/")
for(i in 1:length(var_name)){
  if(file.exists(out_99[i])==F){
    dir.create(out_99[i], recursive=T)
  }
}
out_mcmc <- paste0("01_STM_forloop/mcmc/", var_name, "/")
for(i in 1:length(var_name)){
  if(file.exists(out_mcmc[i])==F){
    dir.create(out_mcmc[i], recursive=T)
  }
}

# Load functions
source("functions/stm.R")

# Load data
H3K4me1 <- as.data.frame(fread(paste0("data/", "MZs_H3K4me1_GENE_log2rpkm.csv")))

# Join dummy data before and after the data frames for H3K4me1
for(i in 1){
  
  eval(parse(text = paste0(
    "df <- ", var_name[i]
  )))
  
  data_names <- names(df)[11:ncol(df)]
  date <- as.character(unique(str_sub(data_names, 3, 8)))
  date_2006 <- ymd(paste0("06", str_sub(date, 3, 6)))
  
  # Second date of dummy1
  diff <- difftime(date_2006[1], date_2006, units = "days")
  #diff[diff > 182] <- - (365 - diff[diff > 182])
  idx2_dummy1 <- tail(which(diff > 10), 1)
  name2_dummy1 <- date[idx2_dummy1]
  # First date of dummy1
  idx1_dummy1 <- idx2_dummy1 - 11
  name1_dummy1 <- date[idx1_dummy1]
  # Start and end of dummy1
  st_dummy1 <- str_which(names(df), name1_dummy1)[1]
  en_dummy1 <- tail(str_which(names(df), name2_dummy1), 1)
  
  # First date of dummy2
  diff <- difftime(date_2006[length(date_2006)], date_2006, units = "days") 
  diff[diff > 182] <- - (365 - diff[diff > 182])
  idx1_dummy2 <- head(which(diff < -10), 1)
  name1_dummy2 <- date[idx1_dummy2]
  # Second date of dummy2
  idx2_dummy2 <- idx1_dummy2 + 11
  name2_dummy2 <- date[idx2_dummy2]
  # Start and end of dummy2
  st_dummy2 <- str_which(names(df), name1_dummy2)[1]
  en_dummy2 <- tail(str_which(names(df), name2_dummy2), 1)
  
  # Rename columns of dummys
  dummy1 <- df[,st_dummy1:en_dummy1]
  names(dummy1) <- str_replace(names(dummy1), "mz18", "mz16")
  dummy2 <- df[,st_dummy2:en_dummy2]
  names(dummy2) <- str_replace(names(dummy2), "mz16", "mz18")
  names(dummy2) <- str_replace(names(dummy2), "mz17", "mz19")
  
  # Join df and dummys
  df <- cbind(df[,1:10], dummy1, df[,11:ncol(df)], dummy2)
  print(names(df))
  
  eval(parse(text = paste0(
    var_name[i], " <- df"
  )))
}


# Prepare data list
df_list <- list()
for(i in 1:length(var_name)){
  eval(parse(text = paste0(
    "df_list[[i]] <- ", var_name[i]
  )))
}


# Load data with full time points
test <- H3K4me1


## Execution of smooth trend model
# for loop of genes for a variable in the STM_function was very slow
# thus, for loop of genes for all variables was adopted
for(i in 1:100){  # for loop of genes
  for(j in 1:length(var_name)){  # for loop of variables
    stm(data = df_list[[j]], gene_idx = i, out_99 = out_99[j], out_mcmc = out_mcmc[j], var_name = var_name[j])
  }
}

