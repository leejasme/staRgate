% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMarkerIndicators.R
\name{getMarkerIndicators}
\alias{getMarkerIndicators}
\title{Get indicators of marker positive/negative in numerator and denominator for cell subpopulations}
\usage{
getMarkerIndicators(dat, col_nm = "marker", markers_nm)
}
\arguments{
\item{dat}{dataframe to add indicators to
expects rows to correspond to a subpopulation}

\item{col_nm}{string indicating column in \code{dat} (i.e., \code{dat$col_nm}) to find the subpopulation names to apply indicators on
Default is 'marker'
Expectation that the naming convention is \if{html}{\out{<num>}}\emph{OF}\if{html}{\out{<denom>}}}

\item{markers_nm}{string vector for names of markers to generate indicators for
No default currently}
}
\value{
dataframe \code{dat} returned with additional columns for indicators of positivity
Each value of \code{markers_nm} will have 2 columns to indicate positivity in numerator or denominator
Naming conventions are \verb{tolower(marker_name)_pos} for positivity indicator in numerator
and \verb{tolower(marker_name)_pos_d} for positivity in denominator
}
\description{
Useful for matching subpopulation names ignoring order and naming of markers
}
\details{
Inputs:
}
