
import numpy as np
import pandas as pd
import os
# import matplotlib.pyplot as plt
import time


# Set output and input directories
out = "03_Distance_chrom_pattern_150bp_weighted/"
out_dist = out + "dist/"
out_minidx = out + "minidx/"
os.makedirs(out_dist, exist_ok=True)
os.makedirs(out_minidx, exist_ok=True)


##### Calculate the distance between genes #####

# Prepare lag series
num = range(0, 52, 1) #7:58

for i in range(len(num)):
  if i == 0:
    ref = np.arange(num[0], 52).reshape(52, 1)
  else:
    temp = np.append(np.arange(num[i], 52), np.arange(0, num[i])).reshape(52, 1)
    ref = np.append(ref, temp, axis = 1)


# Prerequisite
tp = 1
df_class_rev = pd.read_csv("data/MZs_150bp_log2rpkm_splined_kmean_kNN_tp" + str(tp) + ".csv", usecols=['Ahal_ID'])
uniq_gene = set(df_class_rev["Ahal_ID"])
test = pd.read_csv("data/" + "MZs_H3K4me1_GENE_log2rpkm.csv")
test = test.iloc[:10,:]

# Transition probability of Markov model
df_theta = pd.read_csv("data/" + "transition_dist.csv")



### for loop iteration
all_comb = test.shape[0] * test.shape[0] / 2 + test.shape[0] / 2
each_comb = all_comb / 1 # 1 cores

# No. of colums in each row
no_col = np.sort(np.arange(0, test.shape[0])+1)[::-1]

# Calculate the iteration of for loop
sum_no_col = 0
st_loop = np.array([0])
en_loop = np.array([0])
for i in np.arange(test.shape[0]):
  sum_no_col = sum_no_col + no_col[i]
  if sum_no_col >= each_comb:
    st_loop = np.append(st_loop, en_loop[en_loop.shape[0]-1]+1)
    en_loop = np.append(en_loop, i)
    sum_no_col = 0

# Adjust the result
st_loop[1] = 0
st_loop = np.append(st_loop, en_loop[en_loop.shape[0]-1]+1)
st_loop = np.delete(st_loop, 0)
en_loop = np.append(en_loop, test.shape[0]-1)
en_loop = np.delete(en_loop, 0)

# Confirmation
no_row = en_loop - st_loop + 1
sum(no_row)

  



### Distance between gene i and j

# np.arange(st_loop[29], en_loop[29]+1)
for i in np.arange(st_loop[0], en_loop[0]+1):
  
  gene_id1 = test["Ahal_ID"][i]
  
  if gene_id1 in finished_genes:
    continue # when gene has been finished in saizo, skip
  
  if gene_id1 not in uniq_gene:
    continue # when gene locations are overlapped, second genes are not included in "30_kmeans_150bp_230912/"MZs_150bp_log2rpkm_splined_kmean_clus_AhalID_TEID_tp", tp, ".csv""
  df_gene1 = np.loadtxt("data/" + "MZs_150bp_log2rpkm_splined_kmean_kNN_" + gene_id1 + "_tp1_strand_15bin.csv", delimiter=',', skiprows=1, usecols=range(6,58), dtype = "object")
  
  mat_dist = np.empty(no_col[st_loop[0]])
  mat_dist[:] = np.nan
  mat_minidx = np.empty(no_col[st_loop[0]])
  mat_minidx[:] = np.nan
  
  # np.arange(i, test.shape[0])
  for j in np.arange(i, test.shape[0]):
    
    gene_id2 = test["Ahal_ID"][j]
    if gene_id2 not in uniq_gene:
      continue # when gene locations are overlapped, second genes are not included in "30_kmeans_150bp_230912/"MZs_150bp_log2rpkm_splined_kmean_clus_AhalID_TEID_tp", tp, ".csv""
    df_gene2 = np.loadtxt("data/" + "MZs_150bp_log2rpkm_splined_kmean_kNN_" + gene_id2 + "_tp1_strand_15bin.csv", delimiter=',', skiprows=1, usecols=range(6,58), dtype = "object")
    
    dists = np.empty(52)
    vec_gene1 = df_gene1.flatten()
    for k in range(52):
      vec_gene2 = df_gene2[:,ref[:,k]].flatten()
      df_transition = pd.DataFrame({"transition": vec_gene1 + "_" + vec_gene2})
      df_merge = pd.merge(df_transition, df_theta, on='transition')
      dists[k] = sum(df_merge["dist"])
    
    mat_dist[j-st_loop[0]] = min(dists)
    if np.var(dists) != 0:
      mat_minidx[j-st_loop[0]] = np.argmin(dists)
  
  np.savetxt(out_dist + 'dist_' + gene_id1 + '.txt', mat_dist, fmt="%0.2f")
  np.savetxt(out_minidx + 'minidx_' + gene_id1 + '.txt', mat_minidx, fmt="%0.2f")

