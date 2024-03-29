% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netscore.R
\name{NetScore}
\alias{NetScore}
\title{Evaluating and scoring microbial subnetworks}
\usage{
NetScore(df1)
}
\arguments{
\item{df1}{the network topological property indicators of subnetworks ranked in the top k percentile by node count (k=20) are selected and evaluated by ten network topological property indicators}
}
\value{
This function returns a list with the following components:
\item{subnetworks_score.csv}{the topological attributes scores and rank of the chosen subnetwork using ten network topological property indicators.}
}
\description{
Evaluating and scoring microbial subnetworks
}
