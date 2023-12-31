% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oncoplot.R
\name{oncoplot}
\alias{oncoplot}
\title{main function of oncoplot}
\usage{
oncoplot(mat, varis_color = NULL, clinical, clinical_color = NULL, na.value)
}
\arguments{
\item{mat}{matrix, each row is a gene,
and each column is a sample}

\item{varis_color}{namesd vector, color of varis}

\item{clinical}{dataframe of clinical data,
the first column is sample ID.}

\item{clinical_color}{namesd vector, color of clinical varis.}

\item{na.value}{color of na value.}
}
\value{
an aplot object
}
\description{
main function of oncoplot
}
