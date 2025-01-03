
## Quantile function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.5, 0.975, 0.995), names = TRUE)
}


## Smooth trend model function
stm <- function(data, gene_idx, out_99, out_mcmc, var_name){
  
  # Load stan model
  model <- cmdstan_model("stan_model/STM.stan")
  #model <- cmdstan_model("/Volume3/hnishio/R/CREST_seasonal/stan_model/STM_mRNA_prior_211207.stan")
  
  # Name of time points
  data_names <- names(data)[11:ncol(data)]
  mzdate <- str_sub(data_names, 1, 8)
  #date <- as.POSIXct(unique(str_sub(data_names, 8, 13)), format="%y%m%d")
  date <- as.character(unique(str_sub(data_names, 3, 8)))
  
  # Get index of observed values
  date_full <- as.character(unique(str_sub(names(test)[11:ncol(test)], 3, 8)))
  mzdate_full <- paste0("mz", date_full)
  na_date <- setdiff(date_full, date)
  na_pos <- NULL
  for(j in 1:length(na_date)){
    na_pos <- c(na_pos, which(date_full == na_date[j]))
  }
  if(length(na_pos) > 0){
    obsidx <- (1:length(date_full))[-na_pos]
  }else{
    obsidx <- 1:length(date_full)
  }
  
  # Prepare data_list
  data_list <- list(
    Ndate_full = length(date_full),
    Ndate = length(unique(mzdate)),
    Nrep = as.numeric(table(mzdate)),
    sumNrep = sum(as.numeric(table(mzdate))),
    Y = as.numeric(data[gene_idx, 11:ncol(data)]),
    obsidx = obsidx
  )
  
  # Observed data
  d <- data.frame(
    date = mzdate,
    obs = as.numeric(data[gene_idx, 11:ncol(data)])
  )
  for(j in 1:4){ 
    eval(parse(text = 
                 paste0(  "d", j, "<- d %>% group_by(date) %>% slice(", j, ")"  )
    ))
  }
  d1d2 <- full_join(d1, d2, by = "date")
  names(d1d2) <- c("date", "rep1", "rep2")
  d1d2d3 <- full_join(d1d2, d3, by = "date")
  d1d2d3d4 <- full_join(d1d2d3, d4, by = "date")
  names(d1d2d3d4) <- c("date", "rep1", "rep2", "rep3", "rep4")
  
  # Run MCMC
  fit <- model$sample(
    data = data_list,
    #init = function() { list(mu = as.numeric(tapply(d$obs, d$date, mean, na.rm=T))) },
    seed = 1,
    iter_warmup = 1000,
    iter_sampling = 1000,
    chains = 4,
    parallel_chains = 4,
    refresh = 0,
    show_messages = F,
    sig_figs = 4,
    output_dir = "/tmp",
    output_basename = paste0("gene_", gene_idx)
  )
  
  # Summarization of the results
  
  # Obtain MCMC samples
  outcsv_name <- list.files("/tmp")
  outcsv_name <- outcsv_name[grep(paste0("gene_", gene_idx, "-"), outcsv_name)]
  tmp_csv_mu <- NULL
  for(j in 1:length(outcsv_name)){
    tmp_csv <- as.data.frame(fread(cmd = paste0("grep -v '^#' ", "/tmp/", outcsv_name[j])))
    tmp_csv_mu <- rbind(tmp_csv_mu, tmp_csv[,str_starts(names(tmp_csv), "mu.")])
  }
  
  # Remove tmp files
  file.remove(paste0("/tmp/", outcsv_name))
  
  data_st <- which(date_full == "161122")
  data_en <- which(date_full == "181127")
  
  # Save MCMC samples of mu
  fwrite(tmp_csv_mu[,data_st:data_en],
         file = paste0(out_mcmc, "MZs_", var_name, "_log2rpkm_mcmc_", data$Ahal_ID[gene_idx], ".csv"))
  
  # 99% Bayesian credible intervals
  df <- as.data.frame(t(apply(tmp_csv_mu, 2, quantile99)))
  #df <- as.data.frame(round(df, digits = 2))
  df <- df %>%
    mutate(date = mzdate_full) %>%
    relocate(date)
  
  # Add observed data
  df_rev <- full_join(df, d1d2d3d4, by = "date")
  
  fwrite(df_rev[data_st:data_en,],
         file = paste0(out_99, "MZs_", var_name, "_log2rpkm_99_", data$Ahal_ID[gene_idx], ".csv"))
  print(paste0(var_name, "_", data$Ahal_ID[gene_idx], " has finished!"))
  
  
  # Resetting memory
  rm(tmp_csv, tmp_csv_mu, df, d, df_rev, d1, d1d2, d1d2d3, d1d2d3d4, d2, d3, d4, 
     data_list, fit)
  #gc(reset = T); gc(reset = T)
}
