% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_stats_functions.R
\name{match_stats}
\alias{match_stats}
\title{Load summary statistics}
\usage{
match_stats(bigSNP.obj, pop1_stats_df, pop2_stats_df)
}
\arguments{
\item{bigSNP.obj}{\itemize{
\item imported rds file of admixed plink files, pop1_stats_df - loaded  summary statistics for pop 1, pop2_stat_df - loaded summary statistics for pop 2
}}
}
\value{
List of pop 1 matched betas, pop 1 matched p vals, pop 2 matched betas, pop 2 matched p vals
}
\description{
This function matches summary stats with loaded plink files using bigsnpr::snp_matchs. Matches alleles and flips if needed.
}
