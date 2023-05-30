#' Run stacking model
#'
#' This function matches summary stats with loaded plink files using bigsnpr::snp_matchs. Matches alleles and flips if needed.
#'
#' @param bigSNP.obj - bigsnp object containing genotypes, window_size - size of blocks/windows
#' @return dataframe containing local ancestry inference in windows 
#' @export
run_slaPRS <- function(bigSNP.obj, matched_betas_pvals, window_anc_df, full_model, ind.train, ind.test){
  
  ### Iterate over p value threshold (5e-4, 5e-6, 5e-8)
  rsq_results <- as.data.frame(matrix(nrow = 0, ncol=4))
  best_rsq <- 0
  
  for(thres in c(5e-2, 5e-4, 5e-6, 5e-8)){
    print(thres)
    #1) Compute local PRS in each window
    local_prs_df <- make_local_PRS_df(bigSNP.obj, matched_betas_pvals, window_size_num = 5000000, r_sq = 0.10, thres = thres, out_file_path = 'asdf')
    
    #2) Make df for local population specific PRS in each window, local ancestry in each window
    y <- bigSNP.obj$fam$affection
    
    #If full model specific, need to compute interaction terms in each window
    if(full_model == TRUE){
      windows_vector <- sapply(strsplit(colnames(window_anc_df[,-1]), "_"), "[", 3) #Get list of window names
      temp_all_data <- cbind(local_prs_df[, -1], window_anc_df) #df with just local prs and ancestry
      local_data <- as.data.frame(matrix(nrow=nrow(local_prs_df), ncol=0))
      for(window in windows_vector){
        temp_df <-  temp_all_data[ ,grepl(window, names(temp_all_data)) ]
        temp_df$eur_PRS_ancestry <- temp_df[,1]*temp_df[,3]
        temp_df$afr_PRS_ancestry <- temp_df[,2]*temp_df[,3]
        colnames(temp_df) <- c(paste0("eur_PRS_train_20:", window), paste0("afr_PRS_train_20:", window),
                               paste0("anc_window_:", window), 
                               paste0("eur_PRS_ancestry_train_20:", window),
                               paste0("afr_PRS_ancestry_train_20:", window))
        local_data <- cbind(local_data, temp_df)
      }
      all_local_data <- cbind(used.bigSNP$fam$affection, local_data)
    } else {
      all_local_data <- cbind(used.bigSNP$fam$affection, local_prs_df, window_anc_df)
    }
    
    #3) Split into training and testing. Then Run stacking model. 
    #First define phenotype vector from bigSNP object then use corresponding all_local_data df
    colnames(all_local_data)[1] <- 'y'
    
    all_train <- all_local_data[ind.train, ]; all_test <- all_local_data[ind.test, ]
    
    lambda <- 10^seq(-5, 3, length = 20); alpha <- seq(0, 10, length=20)*(0.1) #hyper parameters for elastic net
    
    en_model <- train(
      y ~ ., data = all_train, method = "glmnet", trControl = trainControl("cv", number = 4),
      tuneGrid = expand.grid(alpha = alpha, lambda = lambda)
    )
  
    max_idx <- which.max(en_model$results$Rsquared)
    rsq_results <- rbind(rsq_results, c(thres, en_model$results$Rsquared[max_idx], en_model$results$alpha[max_idx], en_model$results$lambda[max_idx]))
    
    top_rsq <- en_model$results$Rsquared[max_idx]
    
    #Save weights for best model fit across p, alpha, lambda
    if(top_rsq >= best_rsq){
      coefs <- as.data.frame(as.matrix(coef(en_model$finalModel, en_model$bestTune$lambda)))
      coefs$var <- row.names(coefs)
      colnames(coefs) <- c("Weight", "Var")
    }
    
  }
  
  #Get best p value, alpha, and lambda; Not really needed but curious to look at
  colnames(rsq_results) <- c("pval", "rsquared", "alpha", "lambda") 
  best_idx <- which.max(rsq_results$rsquared) 
  best_p <- rsq_results$pval[best_idx]; best_alpha = rsq_results$alpha[best_idx]; best_lambda=rsq_results$lambda[best_idx]
  
  print(paste0("best p-val thres: ", best_p)); print(paste0("best alpha: ", best_alpha))
  print(paste0("best lambda: ", best_lambda))
  
  #Use weights from best parameter combination on testing data
  intercept <- coefs$Weight[1]; coefs_matrix <- as.matrix(coefs$Weight[-1])
  prs_matrix <- as.matrix(all_test[, -1]) #remove phenotype column
  stacked_PRS <- prs_matrix %*% coefs_matrix
  
  return(stacked_PRS)
}



#' Run C+T with single population GWAS
#'
#' This function matches summary stats with loaded plink files using bigsnpr::snp_matchs. Matches alleles and flips if needed.
#'
#' @param bigSNP.obj - bigsnp object containing genotypes, window_size - size of blocks/windows
#' @return dataframe containing local ancestry inference in windows 
#' @export
run_CT <- function(bigSNP.obj, matched_betas_pvals, rsquared, ind.train, ind.test){
  G <- bigSNP.obj$genotypes
  
  pop1_idx <- snp_clumping(G, infos.chr = bigSNP.obj$map$chromosome, ind.row = ind.train, S = -log10(matched_betas_pvals$pop1_pval),
                          thr.r2 = rsquared, size = 250, infos.pos = bigSNP.obj$map$physical.pos,ncores = 1)
  
  #pop1_stats_keep <- sumstats[pop1_idx, ]
  
  pop2_idx <- snp_clumping(G, infos.chr = bigSNP.obj$map$chromosome, ind.row = ind.train, S = -log10(matched_betas_pvals$pop2_pval),
                          thr.r2 = rsquared, size = 250, infos.pos = bigSNP.obj$map$physical.pos,ncores = 1)
  
  #pop2_stats_keep <- sumstats[afr_idx, ]
  
  THR <- seq_log(1, 8, length.out = 20)
  pop1_CTs <- snp_PRS(G, matched_betas_pvals$pop1_beta[pop1_idx], ind.test = ind.train, ind.keep = pop1_idx, same.keep = rep(TRUE, length(pop1_idx)), 
                     lpS.keep = matched_betas_pvals$pop1_pval[pop1_idx], thr.list = THR )
  
  pop2_CTs <- snp_PRS(G, matched_betas_pvals$pop2_beta[pop2_idx], ind.test = ind.train, ind.keep = pop2_idx, same.keep = rep(TRUE, length(pop2_idx)), 
                     lpS.keep = matched_betas_pvals$pop1_pval[pop2_idx], thr.list = THR)
  
  #Get best p value to use in testing data 
  train_y <- as.numeric(bigSNP.obj$fam$affection)[ind.train]
  pop1_cors <- apply(pop1_CTs, 2, cor, train_y)
  plot(THR, pop1_cors, xlab = "-log10(p-value)", ylab = "cors", pch = 20)
  pop1_best_thres <- THR[which.max(pop1_cors)]
  
  pop2_cors <- apply(pop2_CTs, 2, cor, train_y)
  
  plot(THR, pop2_cors, xlab = "-log10(p-value)", ylab = "cors", pch = 20)
  pop2_best_thres <- THR[which.max(pop2_cors)]
  
  #Get C+T in testing data now 
  pop1_CT <- as.vector(snp_PRS(G, matched_betas_pvals$pop1_beta[pop1_idx], ind.test = ind.test, ind.keep = pop1_idx, 
                                         lpS.keep = matched_betas_pvals$pop1_pval[pop1_idx], thr.list = pop1_best_thres )
  )
  
  pop2_CT <- as.vector(snp_PRS(G, matched_betas_pvals$pop2_beta[pop2_idx], ind.test = ind.test, ind.keep = pop2_idx, 
                                         lpS.keep = matched_betas_pvals$pop2_pval[pop2_idx], thr.list = pop2_best_thres ))
  
  return(list("pop1_CT" = pop1_CT, "pop2_CT" = pop2_CT))
}
