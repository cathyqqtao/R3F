% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rrf_times.R
\name{rrf_times}
\alias{rrf_times}
\title{Estimate relative node times.}
\usage{
rrf_times(
  tree.name = "",
  type = c("NEWICK", "NEXUS"),
  outgroup = "",
  filename = ""
)
}
\arguments{
\item{tree.name}{a file name of the branch length tree.}

\item{type}{specify the format of branch length tree. The default is NEWICK.}

\item{outgroup}{a character list containing names of tips belonging to the rooting outgroup, which will be removed in the calculation. If outgroup = "" (the default), the input tree must be rooted and no tips will be removed.}

\item{filename}{a file name specified for the output files.}
}
\value{
A table of relative node times (<filename>_RRF_times.csv).

A timetree of relative node times (<filename>_RRF_timets.nwk).
}
\description{
This function estimates relative node times from a branch length tree using the relative rate framework (RRF).
}
\examples{
rrf_times(tree.name = "example.nwk", type= "NEWICK", outgroup = c("Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus"), filename = "example")
}
\references{
K. Tamura et al. Mol. Biol. Evol. (2018) 35:1770-1782. doi:10.1093/molbev/msy044.
}
\author{
Qiqing Tao (qiqing.tao@temple.edu) and Sudhir Kumar
}
