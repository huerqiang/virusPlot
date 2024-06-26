% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/virus_circle.R
\name{circle_virus}
\alias{circle_virus}
\title{Circle plot of virus insert}
\usage{
circle_virus(virus_info, insert_num, label = "HPV16", bar_color = "#366CB4")
}
\arguments{
\item{virus_info}{The virus information data frame contains three columns.
The first column is the gene name, the second column is the start site,
and the third column is the stop site.}

\item{insert_num}{The data frame of the number of breakpoints in each segment
contains three columns. The first column is the starting site of the segment,
the second column is the ending site of the segment, and the third column is
the number of insertion sites.}

\item{label}{a character of virus name.}

\item{bar_color}{color of bar plot.}
}
\value{
a gg object
}
\description{
Circle plot of virus insert
}
\examples{
virus_info <- data.frame(
         gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
         start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
         end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))

insert_num <- data.frame(start = seq(1, 7801, 100),
    end = seq(101, 7901, 100),
    num = sample(1:28, 79, replace = TRUE))
p <- circle_virus(virus_info, insert_num)
}
