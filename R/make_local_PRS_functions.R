#' Get local population-specific PRS in a window of determined start and end
#'
#' This function matches summary stats with loaded plink files using bigsnpr::snp_matchs. Matches alleles and flips if needed.
#'
#' @param bigSNP.obj - bigsnp object containing genotypes, window_size - size of blocks/windows
#' @return dataframe containing local population-specific PRS in regions of window_size
#' @export
get_local_PRS <- function(bigSNP.obj, ind.idx, start, end, chrom, r2, thres, matched_betas_pvals){
  chr_POS <- bigSNP.obj$map$physical.pos
  map <- bigSNP.obj$map[,-(2:3)]; names(map) <- c("chr", "pos", "a0", "a1")
  
  window_bigSNP <- snp_attach(snp_subset(bigSNP.obj, ind.row = ind.idx, ind.col = which(start <= chr_POS & chr_POS < end)))
  window_G <- window_bigSNP$genotypes; window_G <- window_G$copy(code = c(0, 1, 2, rep(0, 253)))

  sumstat_idx <- which((map$chr == chrom & map$pos >= start) & (map$chr == chrom & map$pos < end) )
  pop1_beta <- matched_betas_pvals$pop1_beta; pop1_pval <- matched_betas_pvals$pop1_pval
  window_pop1_betas <- pop1_beta[sumstat_idx]; window_pop1_pvals <- pop1_pval[sumstat_idx] #these are already -log10 transformed
  window_pop1_betas[is.na(window_pop1_betas)] <- 0; window_pop1_pvals[is.na(window_pop1_pvals)] <- 0 #snp_PRS doesnt like NA values 
  
  pop2_beta <- matched_betas_pvals$pop2_beta; pop2_pval <- matched_betas_pvals$pop1_pval
  window_pop2_betas <- pop2_beta[sumstat_idx]; window_pop2_pvals <- pop2_pval[sumstat_idx] #these are already -log10 transformed
  window_pop2_betas[is.na(window_pop2_betas)] <- 0; window_pop2_pvals[is.na(window_pop2_pvals)] <- 0 #snp_PRS doesnt like NA values 
  
  #Get local PRS via C+T for both pop1 and pop2ican stats
  #pop1opeans
  pop1_clumped_snp_idx <- snp_clumping(window_G, infos.chr = window_bigSNP$map$chromosome, ind.row = rows_along(window_G), 
                                      S = window_pop1_pvals,
                                      thr.r2 = r2, size = 250, infos.pos = window_bigSNP$map$physical.pos,ncores = 1)
  
  pop1_window_PRS <- snp_PRS(window_G, window_pop1_betas[pop1_clumped_snp_idx], ind.test = rows_along(window_G), ind.keep = pop1_clumped_snp_idx, 
                            lpS.keep = window_pop1_pvals[pop1_clumped_snp_idx], 
                            #same.keep = ref_check,
                            thr.list = -log10(thres))
  
  #pop2icans
  pop2_clumped_snp_idx <- snp_clumping(window_G, infos.chr = window_bigSNP$map$chromosome, ind.row = rows_along(window_G),
                                      S = window_pop2_pvals,
                                      thr.r2 = r2, size = 250, infos.pos = window_bigSNP$map$physical.pos,ncores = 1)
  
  pop2_window_PRS <- snp_PRS(window_G, window_pop2_betas[pop2_clumped_snp_idx], ind.test = rows_along(window_G), ind.keep = pop2_clumped_snp_idx, 
                            lpS.keep = window_pop2_pvals[pop2_clumped_snp_idx], thr.list = -log10(thres))
  
  prs_list <- list("pop2" = pop2_window_PRS, "pop1" = pop1_window_PRS)
  return(prs_list) 
}


#' Make dataframe containing local population-specific PRS in each window of determined size looping through entire chromosome file 
#'
#' @param bigSNP.obj - bigsnp object containing genotypes, window_size - size of blocks/windows
#' @return dataframe containing local population-specific PRS in regions of window_size
#' @export
make_local_PRS_df <- function(bigSNP.obj, matched_betas_pvals, window_size_num, r_sq, thres, out_file_path){
  chr_POS <- bigSNP.obj$map$physical.pos
  start <- chr_POS[1]; end <- chr_POS[length(chr_POS)]
  map <- bigSNP.obj$map[,-(2:3)]; names(map) <- c("chr", "pos", "a0", "a1")
  chrom <- unique(map$chr)
  
  all_prs_df <- data.frame(matrix(nrow = N, ncol = 0))
  
  while(start < end){
    #print(start)
    window_end <- start + window_size_num
    if(length(which(start <= chr_POS & chr_POS < window_end)) == 0){
      start <- start + window_size_num
      next
    } #Check if any SNPs in window. SKip if not
    
    list_prs <- get_local_PRS(bigSNP.obj, 1:N, start, start+window_size_num, chrom=20, r2 = 0.10, thres = 5e-4, matched_betas_pvals)
    all_prs_df[ , ncol(all_prs_df) + 1] <- (list_prs$pop1); colnames(all_prs_df)[ncol(all_prs_df)] <- paste0("pop1_PRS_", chrom, ":", start, "-", start+window_size_num)
    all_prs_df[ , ncol(all_prs_df) + 1] <- (list_prs$pop2); colnames(all_prs_df)[ncol(all_prs_df)] <- paste0("pop2_PRS_", chrom, ":", start, "-", start+window_size_num)
    
    start <- start + window_size_num
  }
  
  #Output files
  all_prs_df2 <- cbind(bigSNP.obj$fam$sample.ID, all_prs_df); colnames(all_prs_df2)[1] <- "eid"
  return(all_prs_df2)
  #write.table(train_all_prs2, paste0("./output/real/", phenotype, "/", thres, "/", thres, "_AA_20k_reversed_train_stack_local_prs_", window_size, "_chr", chrom, ".txt"), row.names = FALSE, col.names=TRUE, quote=FALSE)
}  
