#' Get local population-specific PRS in a window of determined start and end
#'
#' This function matches summary stats with loaded plink files using bigsnpr::snp_matchs. Matches alleles and flips if needed.
#'
#' @param bigSNP.obj - bigsnp object containing genotypes, window_size - size of blocks/windows
#' @return dataframe containing local population-specific PRS in regions of window_size
#' @export
get_local_PRS <- function(chr_bigSNP, ind.idx, start, end, r2, thres, matched_betas_pvals){
  chr_POS <- chr_bigSNP$map$physical.pos 
  window_bigSNP <- snp_attach(snp_subset(chr_bigSNP, ind.row = ind.idx, ind.col = which(start <= chr_POS & chr_POS < window_end)))
  window_G <- window_bigSNP$genotypes; window_G <- window_G$copy(code = c(0, 1, 2, rep(0, 253)))
  window_map <- window_bigSNP$map
  
  sumstat_idx <- which((map$chr == chrom & map$pos >= start) & (map$chr == chrom & map$pos < end) )
  window_eur_betas <- eur_beta[sumstat_idx]; window_eur_pvals <- eur_pval[sumstat_idx] #these are already -log10 transformed
  window_eur_betas[is.na(window_eur_betas)] <- 0; window_eur_pvals[is.na(window_eur_pvals)] <- 0 #snp_PRS doesnt like NA values 
  
  window_afr_betas <- afr_beta[sumstat_idx]; window_afr_pvals <- afr_pval[sumstat_idx] #these are already -log10 transformed
  window_afr_betas[is.na(window_afr_betas)] <- 0; window_afr_pvals[is.na(window_afr_pvals)] <- 0 #snp_PRS doesnt like NA values 
  
  #Get local PRS via C+T for both eur and african stats
  #europeans
  eur_clumped_snp_idx <- snp_clumping(window_G, infos.chr = window_bigSNP$map$chromosome, ind.row = rows_along(window_G), 
                                      S = window_eur_pvals,
                                      thr.r2 = r2, size = 250, infos.pos = window_bigSNP$map$physical.pos,ncores = 1)
  
  eur_window_PRS <- snp_PRS(window_G, window_eur_betas[eur_clumped_snp_idx], ind.test = rows_along(window_G), ind.keep = eur_clumped_snp_idx, 
                            lpS.keep = window_eur_pvals[eur_clumped_snp_idx], 
                            #same.keep = ref_check,
                            thr.list = -log10(thres))
  
  #africans
  afr_clumped_snp_idx <- snp_clumping(window_G, infos.chr = window_bigSNP$map$chromosome, ind.row = rows_along(window_G),
                                      S = window_afr_pvals,
                                      thr.r2 = r2, size = 250, infos.pos = window_bigSNP$map$physical.pos,ncores = 1)
  
  afr_window_PRS <- snp_PRS(window_G, window_afr_betas[afr_clumped_snp_idx], ind.test = rows_along(window_G), ind.keep = afr_clumped_snp_idx, 
                            lpS.keep = window_afr_pvals[afr_clumped_snp_idx], thr.list = -log10(thres))
  
  prs_list <- list("afr" = afr_window_PRS, "eur" = eur_window_PRS)
  return(prs_list) 
}


#' Make dataframe containing local population-specific PRS in each window of determined size looping through entire chromosome file 
#'
#' @param bigSNP.obj - bigsnp object containing genotypes, window_size - size of blocks/windows
#' @return dataframe containing local population-specific PRS in regions of window_size
#' @export
make_local_PRS_df <- function(bigSNP.obj, pop1_stats_df, pop2_stats_df, window_size_num, r_sq, thres, out_file_path){
  chr_POS <- bigSNP.obj$map$physical.pos
  start <- chr_POS[1]; end <- chr_POS[length(chr_POS)]
  
  while(start < end){
    print(start)
    window_end <- start + window_size_num
    if(length(which(start <= chr_POS & chr_POS < window_end)) == 0){
      start <- start + window_size_num
      next
    } #Check if any SNPs in window. SKip if not
    
    train_list_prs <- get_local_PRS(bigSNP.obj, ind.train, start, start+window_size_num, r_sq, thres)
    train_all_prs[ , ncol(train_all_prs) + 1] <- train_list_prs$eur; colnames(train_all_prs)[ncol(train_all_prs)] <- paste0("pop1_PRS_", chrom, ":", start, "-", start+window_size_num)
    train_all_prs[ , ncol(train_all_prs) + 1] <- train_list_prs$afr; colnames(train_all_prs)[ncol(train_all_prs)] <- paste0("pop2_PRS_", chrom, ":", start, "-", start+window_size_num)
    
    start <- start + window_size_num
  }
  
  #Output files
  train_all_prs2 <- cbind(train_ids, train_all_prs); colnames(train_all_prs2)[1] <- "eid"
  write.table(train_all_prs2, paste0("./output/real/", phenotype, "/", thres, "/", thres, "_AA_20k_reversed_train_stack_local_prs_", window_size, "_chr", chrom, ".txt"), row.names = FALSE, col.names=TRUE, quote=FALSE)
}  