
#' Metadata associated with the example pedigree
#'
#' A dataframe.
#'
#' @format ## `pedMeta`
#' A 'data.frame' of 7100 individuals (rows) with 12 variables (cols):
#' \describe{
#'   \item{id}{Integer individual ID}
#'   \item{population}{Population code. A, B or AB}
#'   \item{generation}{Generation of the individual}
#'   \item{mid}{dam ID}
#'   \item{fid}{sire ID}
#'   \item{gv1}{genetic value}
#'   \item{pv1}{phenotypic value}
#'   \item{gv2}{genetic value}
#'   \item{pv2}{phenotypic value}
#'   \item{gv}{genetic value}
#'   \item{pv}{phenotypic value}
#'   \item{generationPlotShift}{for plotting}
#' }
#' @source Simulation
"pedMeta"




#' Example pedigree L inverse matrix
#'
#' An L inverse matrix generated from an AlphaSimR simulation of 50 generations.
#' An original population splits into sub-populations A and B. After a number of
#' generations, crossbreeding starts.
#'
#'
#' @format ## `pedLinv`
#' Matrix object of class 'spam' of dimension 7100x7100,
#'     with 21040 (row-wise) nonzero elements.
#'     Density of the matrix is 0.0417%.
#'     Class 'spam' (32-bit)
#' @source Simulation
"pedLinv"
