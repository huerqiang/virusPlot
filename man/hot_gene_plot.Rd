% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hot_gene.R
\name{hot_gene_plot}
\alias{hot_gene_plot}
\title{Title}
\usage{
hot_gene_plot(
  hot_gene_result,
  hot_gene_host = 5,
  hot_gene_virus = 8,
  break_y = TRUE,
  observed_color = "grey30",
  expected_color = "grey60"
)
}
\arguments{
\item{hot_gene_result}{result of get_hot_gene()}

\item{hot_gene_host}{hot insert genes of host to plot, can be a number of a vector of genes}

\item{hot_gene_virus}{hot insert genes of virus to plot, can be a number of a vector of genes}

\item{break_y}{If ture (the default), use ggbreak to set an y axis break.}

\item{observed_color}{color of observed value.}

\item{expected_color}{color of expected value.}
}
\value{
gg object
}
\description{
Title
}
\examples{
data(insert_info)
virus_info <- data.frame(
      gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
      start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
      end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))
hot_gene <- get_hot_gene(virus_info, insert_info)
insert_plot <- hot_gene_plot(hot_gene)
}
