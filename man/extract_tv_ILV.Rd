% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_tv_ILV.R
\name{extract_tv_ILV}
\alias{extract_tv_ILV}
\title{extract_tv_ILV()}
\usage{
extract_tv_ILV(nbda_object, ILV_type)
}
\arguments{
\item{nbda_object}{NBDAdata object}

\item{ILV_type}{character indicating types of ILVs: @asocILVdata = asocial, @intILVdata = social, and @multiILVdata = multiplicative (affects both social and asocial equally)}
}
\value{
List of vectors for each ILV
}
\description{
private function to extract the tv ILV matrices and transform them into ILV_tv dataframe for further processing
This is a really silly way of doing it:
\enumerate{
\item the original TV matrix is impossible to reconstruct from NBDA objs, since NBDA only saves values for when individuals are naive.
This shouldn't matter for basic functionality, but if you have a use-case where you wanted to do counterfactuals etc, one needs that data.
\item this basically converts vector > matrix > vector, which will then be converted again into a matrix by the import_NBDA_STb function..
However it works, will come back later to make it more efficient.
}
}
