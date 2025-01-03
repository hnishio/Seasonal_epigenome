
# Load packages
library(tidyverse)
library(patchwork)
library(data.table)
library(cmdstanr)
library(posterior)
set_cmdstan_path("~/cmdstan/")
library(network)
library(sna)
library(GGally)

# Create the output directory
out <- paste0("02_MM_allbins_150bp/")
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

# Quantile function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.05, 0.5, 0.95, 0.975, 0.995), names = TRUE)
}
norm_c <- function(x){(x-min(x))/(max(x)-min(x))}



### Markov model

# Load class
df_class <- as.data.frame(fread(file = paste0("data/", "MZs_150bp_log2rpkm_splined_kmean_kNN_tp1.csv")))
head(df_class)
df_class_rev <- df_class

df_class_rev[,-(1:6)][df_class[,-(1:6)]=="E1"] <- 1
df_class_rev[,-(1:6)][df_class[,-(1:6)]=="E2"] <- 2
df_class_rev[,-(1:6)][df_class[,-(1:6)]=="I1"] <- 3
df_class_rev[,-(1:6)][df_class[,-(1:6)]=="I2"] <- 4
df_class_rev[,-(1:6)][df_class[,-(1:6)]=="F1"] <- 5
df_class_rev[,-(1:6)][df_class[,-(1:6)]=="F2"] <- 6
df_class_rev[,-(1:6)][df_class[,-(1:6)]=="H1"] <- 7
df_class_rev[,-(1:6)][df_class[,-(1:6)]=="H2"] <- 8
df_class_rev[,-(1:6)][df_class[,-(1:6)]=="N"] <- 9

df_class_rev2 <- df_class_rev

# Set data_list
data_list <- list(
  K = 9, # No. of categories
  T = ncol(df_class_rev2)-6, # No. of time
  N = nrow(df_class_rev2), # No. of bins
  z = df_class_rev2 [,-(1:6)], # Category matrix (row: no. of bins, column: no. of time)
  alpha = rep(1, 9) # Prior of transition probabilities
)

# Load stan model
model <- cmdstan_model("stan_model/MM.stan")

# Execute MCMC
fit <- model$sample(
  data = data_list,
  seed = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  parallel_chains = 4,
  refresh = 200,
  #sig_figs = 4,
  output_dir = "/tmp",
  output_basename = paste0("MallM")
)

# Obtain MCMC samples
outcsv_name <- list.files("/tmp")
outcsv_name <- outcsv_name[grep(paste0("MallM"), outcsv_name)]
tmp_csv_theta <- NULL
for(j in 1:length(outcsv_name)){
  tmp_csv <- as.data.frame(fread(cmd = paste0("grep -v '^#' ", "/tmp/", outcsv_name[j])))
  tmp_csv_theta <- rbind(tmp_csv_theta, tmp_csv[,str_starts(names(tmp_csv), "theta.")])
}

# Remove tmp files
file.remove(paste0("/tmp/", outcsv_name))

# Save MCMC samples of mu
fwrite(tmp_csv_theta,
       file = paste0(out, "MZs_150bp_log2rpkm_splined_kmean_kNN_MM_mcmc.csv"))

# 99% Bayesian credible intervals
df_theta <- as.data.frame(t(apply(tmp_csv_theta, 2, quantile99)))
df_theta_rev <- df_theta %>%
  mutate(from = str_sub(row.names(df_theta),7,7),
         to = str_sub(row.names(df_theta),9,9)) %>%
  relocate(from, to)

fwrite(df_theta_rev,
       file = paste0(out, "MZs_150bp_log2rpkm_splined_kmean_kNN_MM_99.csv"))



### Visualization of network

# Edit data
df_theta <- as.data.frame(fread(file = paste0(out, "MZs_150bp_log2rpkm_splined_kmean_kNN_MM_99.csv")))
df_net <- df_theta[,c("from","to","50%")]
df_net_rev <- df_net
df_net_rev[,1:2][df_net[,1:2]==1] <- "E1"
df_net_rev[,1:2][df_net[,1:2]==2] <- "E2"
df_net_rev[,1:2][df_net[,1:2]==3] <- "I1"
df_net_rev[,1:2][df_net[,1:2]==4] <- "I2"
df_net_rev[,1:2][df_net[,1:2]==5] <- "F1"
df_net_rev[,1:2][df_net[,1:2]==6] <- "F2"
df_net_rev[,1:2][df_net[,1:2]==7] <- "H1"
df_net_rev[,1:2][df_net[,1:2]==8] <- "H2"
df_net_rev[,1:2][df_net[,1:2]==9] <- "N"
df_net_rev2 <- df_net_rev %>%
  slice(-c(1,11,21,31,41,51,61,71,81)) %>%
  rename(prob = `50%`) %>%
  mutate(prob = norm_c(prob)*2+0.1)

# Plotting (- prob < thre[i])
thre <- seq(0, 1, by = 0.2)
for(i in 1:length(thre)){
  df_net_rev3 <- df_net_rev2[which(df_net_rev2$prob > thre[i]),]
  n <- as.network(df_net_rev3, directed = TRUE)
  set.seed(1)
  g <- ggnet2(n, edge.size = "prob",
              size = 6, label.size = 9/ggplot2::.pt, 
              label = TRUE, arrow.size = 5, arrow.gap = 0.035)
  ggsave(paste0(out, "MZs_150bp_log2rpkm_splined_kmean_kNN_MM_50_prob>", thre[i], ".pdf"),
         g, height = 70, width = 70, units = "mm")
}

