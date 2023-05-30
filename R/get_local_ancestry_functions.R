#' Get local population-specific PRS in a window of determined start and end
#'
#' This function matches summary stats with loaded plink files using bigsnpr::snp_matchs. Matches alleles and flips if needed.
#'
#' @param bigSNP.obj - bigsnp object containing genotypes, window_size - size of blocks/windows
#' @return dataframe containing local ancestry inference in windows 
#' @export
get_la_example <- function(la_df, window_size){
  #Create df of window start and end 
  window_df <- data.frame()
  while(start < end){
    window_df <- rbind(window_df, c(start, start+window_size))
    start <- start + window_size
  }
  colnames(window_df) <- c("start", "end")
  
  #Subtact 1 to be able to make values = {0, 1} and get means for faster computation
  num_haps = dim(la_df)[2]
  la_df[,3:num_haps] <- la_df[, 3:num_haps] - 1
  
  local_anc_haps <- as.data.frame(la_df[,-c(1,2)])
  samp_ids <- unique(substr(colnames(local_anc_haps), 1, 8)) #Get list of sample ids
  
  all_window <- as.data.frame(samp_ids)
  
  for(i in 1:nrow(window_df)){
    print(i)
    start <- window_df[i,1]; end <- window_df[i,2]
    local_window <- as.data.frame(subset(la_df, la_df$pos > start & la_df$pos < end))    
    
    samp_anc <- vector()
    for(samp in samp_ids){
      samp_haps <- local_window[ , grepl(samp, colnames(local_window)) ]
      global_anc <- mean(as.matrix(samp_haps))
      samp_anc <- c(samp_anc, global_anc)
    }
    
    window_anc <- as.data.frame(samp_anc)
    colnames(window_anc) <- c(paste0("anc_window_", start, "-", end))
    all_window <- cbind(all_window, window_anc)
  }
 
  return(all_window)   
}  
