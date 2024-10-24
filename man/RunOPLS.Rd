% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModulonR_functions.R
\name{RunOPLS}
\alias{RunOPLS}
\title{Run OPLS}
\usage{
RunOPLS(
  mat,
  annotation,
  ident.1,
  ident.2 = NULL,
  scale_weights = TRUE,
  scale_vipVn = TRUE,
  ...
)
}
\arguments{
\item{mat}{Matrix with features and samples in rows and columns, respectively.}

\item{annotation}{Character vector with the class of the samples.}

\item{ident.1}{Character, query class.}

\item{ident.2}{Character vector with background classes. If NULL, all the classes in annotation but the query class will be included, Default: NULL.}

\item{scale_weights}{Boolean; Scale weightMN and weightStarMN.}

\item{scale_vipVn}{Boolean; Scale vipVn.}
}
\value{
Data frame with a summary of the \code{ropls::opls()}. The columns are:

\itemize{
\item weightMN: OPLS weights. (see ? ropls::opls)
}

\itemize{
\item weightStarMN: OPLS projections. (see ? ropls::opls)
}
}
\description{
Wrapper function to run the ropls::opls().
}
\seealso{
\code{\link[ropls]{opls}}
}
