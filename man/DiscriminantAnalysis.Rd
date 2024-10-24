% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModulonR_functions.R
\name{DiscriminantAnalysis}
\alias{DiscriminantAnalysis}
\title{Discriminant Analysis}
\usage{
DiscriminantAnalysis(
  mat,
  method = "OPLS",
  annotation,
  BackgroundClasses = NULL,
  QueryClasses = NULL
)
}
\arguments{
\item{mat}{Matrix with features and samples in rows and columns, respectively.}

\item{method}{Character; the discriminant analysis method. Possible values are: "OPLS" (see \code{FindMarkersOPLS()}).}

\item{annotation}{Character vector with the class of the samples.}

\item{BackgroundClasses}{Classes to compare with the query class. If NULL, all the classes in annotation except the query class will be included as background. Default: NULL}

\item{QueryClasses}{Character vector with the classes to compare with respect to the background. If NULL, all the classes in annotation will be analyzed (separately) with respect to the background. Default: NULL}
}
\value{
Data frame with the discriminant analysis results for each analyzed class.
}
\description{
This function computes the discriminant analysis (OPLS-DA from the \code{ropls} library), comparing every class in a categorical variable against the background (i.e., the other classes).
}
