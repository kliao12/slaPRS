#' Load summary statistics 
#'
#' This function matches summary stats with loaded plink files using bigsnpr::snp_matchs. Matches alleles and flips if needed.
#'
#' @param bigSNP.obj - imported rds file of admixed plink files, pop1_stats_df - loaded  summary statistics for pop 1, pop2_stat_df - loaded summary statistics for pop 2
#' @return List of pop 1 matched betas, pop 1 matched p vals, pop 2 matched betas, pop 2 matched p vals
#' @export
match_stats <- function(bigSNP.obj, pop1_stats_df, pop2_stats_df){
  map <- bigSNP.obj$map[,-(2:3)]
  names(map) <- c("chr", "pos", "a0", "a1")
  
  #Check alleles and get vector of Population 1 betas and pvals
  pop1_info_snp <- bigsnpr::snp_match(pop1_stats_df, map)
  pop1_beta <- rep(0, ncol(G2)); pop1_beta[pop1_info_snp$`_NUM_ID_`] <- pop1_info_snp$beta
  pop1_pval <- rep(0, ncol(G2)); pop1_pval[pop1_info_snp$`_NUM_ID_`] <- pop1_info_snp$pvalue_neg_log10_GC
  
  #Repeat for Population 2 GWAS data
  pop2_info_snp <- bigsnpr::snp_match(pop2_stats_df, map)
  pop2_beta <- rep(0, ncol(G2)); pop2_beta[pop2_info_snp$`_NUM_ID_`] <- pop2_info_snp$beta
  pop2_pval <- rep(0, ncol(G2)); pop2_pval[pop2_info_snp$`_NUM_ID_`] <- pop2_info_snp$pvalue_neg_log10_GC
  return(list("pop1_beta" = pop1_beta, "pop1_pval" = pop1_pval, "pop2_beta" = pop2_beta, "pop2_pval" = pop2_pval)) 
}
