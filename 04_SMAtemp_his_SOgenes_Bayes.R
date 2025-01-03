
# Load packages
library(tidyverse)
library(data.table)
library(patchwork)
library(scales)
library(cmdstanr)
library(posterior)
set_cmdstan_path("~/cmdstan/")

# Create the output directory
out <- "04_SMAtemp_his_SOgenes_Bayes/"
if(file.exists(out)==F){
  dir.create(out, recursive=T)
}

var_name <- c("mRNA", "H3K4me1", "H3K4me2", "H3K4me3", "H4ac", "H3K36me3", 
              "H3K27me3", "H2AZ", "H3K9me2", "H3")


# Quantile function
quantile99 <- function(x){
  quantile(x, probs = c(0.005, 0.025, 0.5, 0.975, 0.995), names = TRUE)
}

# Load stan model
model <- cmdstan_model("stan_model/ordering.stan")

# Load data
thre_r2_all <- c(0.7)

# for loop of threshold R2 values
for(i in 1:length(thre_r2_all)){
  
  # Load data
  thre_r2 <- thre_r2_all[i]
  df_day_max <- as.data.frame(fread(paste0("data/SMAtemp_his_2016-2018_dayR2max_R2thre", thre_r2, ".csv")))
  
  # extract genes with more than two marks
  mat_TorF <- apply(df_day_max[,-1], 1, is.na)
  num_mark <- as.numeric(apply(!(mat_TorF), 2, sum))
  df_day_max2 <- df_day_max[num_mark >= 2,]
  
  # Initialize an empty vector to store the constraints
  constraints_all <- character()
  
  # for loop of genes
  for(j in 1:nrow(df_day_max2)){
    
    df_gene <- df_day_max2[j,-1]
    # Create the dataframe
    df <- data.frame(var = var_name[!is.na(df_gene)], 
                     level = df_gene[!is.na(df_gene)])
    
    # Sort the dataframe by 'level'
    # df_sorted <- df[order(df$level), ]
    
    # Initialize an empty vector to store the constraints
    constraints <- character()
    
    # Iterate over all unique pairs of variables
    for (i in 1:(nrow(df) - 1)) {
      for (j in (i + 1):nrow(df)) {
        if (df$level[i] == df$level[j]) {
          # Equal levels: use '='
          constraints <- c(constraints, paste0(df$var[i], " = ", df$var[j]))
        } else if (df$level[i] < df$level[j]) {
          # var[i] level is less than var[j] level: use '<'
          constraints <- c(constraints, paste0(df$var[i], " < ", df$var[j]))
        } else {
          # var[i] level is greater than var[j] level: use '<'
          constraints <- c(constraints, paste0(df$var[j], " < ", df$var[i]))
        }
      }
    }
    
    # Integrate all constraints
    constraints_all <- c(constraints_all, constraints)
  }# for loop of genes
  
  
  ## Remove "=" from constraints
  idx_equal <- grep(" = ", constraints_all)
  constraints_all2 <- constraints_all[-idx_equal]
  df_constraints <- data.frame(constraint = constraints_all2)
  
  fwrite(df_constraints, paste0(out, "constraints_", thre_r2, ".csv"))
  
  
  # Read the CSV file
  df_constraints <- as.data.frame(fread(paste0(out, "constraints_", thre_r2, ".csv"), sep=","))
  
  # Process the constraints
  # Split each constraint into two parts: var1 and var2
  constraints <- str_split_fixed(df_constraints$constraint, " < ", 2)
  colnames(constraints) <- c("var1", "var2")
  constraints <- as.data.frame(constraints)
  
  # Create a list for Stan data
  data_list <- list(N_constraints = nrow(constraints), 
                    var1 = as.integer(factor(constraints$var1, levels = var_name)),
                    var2 = as.integer(factor(constraints$var2, levels = var_name)),
                    sigma = 1) # Example standard deviation
  
  # Run MCMC
  fit <- model$sample(
    data = data_list,
    seed = 1,
    iter_warmup = 1000,
    iter_sampling = 1000,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    show_messages = F,
    sig_figs = 4,
    output_dir = "/tmp",
    output_basename = "ordering"
  )
  
  # Obtain MCMC samples
  outcsv_name <- list.files("/tmp")
  outcsv_name <- outcsv_name[grep("ordering", outcsv_name)]
  tmp_csv_means <- NULL
  tmp_csv_diff <- NULL
  for(j in 1:length(outcsv_name)){
    tmp_csv <- as.data.frame(fread(cmd = paste0("grep -v '^#' ", "/tmp/", outcsv_name[j])))
    tmp_csv_means <- rbind(tmp_csv_means, tmp_csv[,str_starts(names(tmp_csv), "means.")])
    tmp_csv_diff <- rbind(tmp_csv_diff, tmp_csv[,str_starts(names(tmp_csv), "diff.")])
  }
  fwrite(tmp_csv_means, file = paste0(out, "mcmc_means_", thre_r2, ".csv"))
  fwrite(tmp_csv_diff, file = paste0(out, "mcmc_diff_", thre_r2, ".csv"))
  
  
  # 99% Bayesian credible intervals
  df_means <- as.data.frame(t(apply(tmp_csv_means, 2, quantile99)))
  df_diff <- as.data.frame(t(apply(tmp_csv_diff, 2, quantile99)))
  df_q <- cbind(data.frame(var = var_name), 
                df_means, 
                rbind(NA, df_diff))
  names(df_q) <- paste0(c("", rep("means_",5), rep("diff_",5)), names(df_q))
  fwrite(df_q, file = paste0(out, "quantile_", thre_r2, ".csv"))
  
  # Order of var
  df_q$var[order(df_q$`means_50%`)]
  df_q$var[order(df_q$`diff_50%`)[-nrow(df_q)]]
  
  # Remove tmp files
  file.remove(paste0("/tmp/", outcsv_name))
  
}# for loop of threshold R2 values





##### Box plot of dayR2max

# mRNA, H3K4me1, H3K4me2, H3K4me3, H4ac, H3K36me3, H3K27me3, H2AZ, H3K9me2, H3
cols <- c("#000000", "#FFA600", "#FF0000", "#8C0000", "#2d8c0a", 
          "#FF00FF", "#0000FF", "#66CCFF", "#000080", "#808080")

# Load data
thre_r2_all <- c(0.7)

# for loop of threshold R2 values
for(i in 1:length(thre_r2_all)){
  
  thre_r2 <- thre_r2_all[i]
  
  df_q <- as.data.frame(fread(paste0(out, "quantile_", thre_r2, ".csv")))
  df_q_rev <- df_q %>%
    na.omit() %>%
    select(var, `diff_2.5%`, `diff_50%`, `diff_97.5%`)

  df_q_rev$var <- factor(df_q_rev$var, levels = var_name)
  
  
  # significant difference
  loc_sig <- sort(c(which(apply(df_q_rev[,2:4]>0, 1, sum) == 3),
                    which(apply(df_q_rev[,2:4]<0, 1, sum) == 3)))
  
  ymin <- min(df_q_rev[,2:4])
  ymax <- max(df_q_rev[,2:4])
  yrange <- ymax - ymin
  
  # Plotting
  g <- ggplot(data = df_q_rev, aes(x = var, color = var)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(aes(ymin = `diff_2.5%`, ymax = `diff_97.5%`),
                  linewidth = 0.3, width = 0.3) +
    geom_point(aes(y = `diff_50%`), size = 0.6) +
    annotate("text", x=loc_sig, y=ymax+yrange*0.05, label="*", size=7/ggplot2::.pt) +
    coord_cartesian(ylim = c(ymin, ymax+yrange*0.1)) +
    scale_color_manual(values = cols[-1]) +
    theme_bw() +
    theme(#plot.title=element_text(size=6),
      legend.position = "none",
      axis.title=element_text(size=6), 
      axis.title.x=element_blank(), 
      axis.text=element_text(size=6),
      axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(#title=var_name[j], 
      #x="Distance between\nadjacent genes (kb)", 
      y="Relative time delay\nfrom mRNA")
  
  ggsave(paste0(out, "diff_", thre_r2, ".pdf"),
         g, height = (6/7)*45, width = (6/7)*60, units = "mm")
}  

