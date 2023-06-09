% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_slaPRS_functions.R
\name{run_CT}
\alias{run_CT}
\title{Run C+T with single population GWAS}
\usage{
run_CT(bigSNP.obj, matched_betas_pvals, rsquared, ind.train, ind.test)
}
\arguments{
\item{bigSNP.obj}{\itemize{
\item bigsnp object containing genotypes, window_size - size of blocks/windows
}}
}
\value{
dataframe containing local ancestry inference in windows
}
\description{
This function matches summary stats with loaded plink files using bigsnpr::snp_matchs. Matches alleles and flips if needed.
}
