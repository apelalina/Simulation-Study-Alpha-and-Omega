n_sim <- 6

cond <- 6
results <- t(sapply(1:n_sim, FUN = function(i){
  simulate_reliability_trial(condition_grid[cond, "loading_type"], condition_grid[cond, "n"], condition_grid[cond, "n_items"])
}))

for(cond in 1:nrow(condition_grid)){
  results <- t(sapply(1:n_sim, FUN = function(i){
    simulate_reliability_trial(condition_grid[cond, "loading_type"], condition_grid[cond, "n"], condition_grid[cond, "n_items"])
  }))
}

library(pbapply)
library(parallel)
cl <- makeCluster(6)
clusterExport(cl = cl, varlist = c("simulate_reliability_trial", "n_sim", "cond", "condition_grid", "N_BOOT", "alpha_boot", "omega_boot", "lavaan_omega"))
clusterEvalQ(cl = cl, expr = {library(lavaan)
  library(psych)
  library(boot)
  })
result <- t(pbsapply(1:n_sim, FUN = function(i){
  simulate_reliability_trial(condition_grid[cond, "loading_type"], condition_grid[cond, "n"], condition_grid[cond, "n_items"])
}, cl = cl))
stopCluster(cl)
