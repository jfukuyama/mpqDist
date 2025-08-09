#' Gentry Forest Plot Dataset as a Phyloseq Object
#'
#' A dataset containing species abundance data, phylogenetic tree, sample metadata, and taxonomy
#' information for forest plots from the Gentry dataset. This dataset has been processed to
#' aggregate species counts by unique trinomial names, mapped onto a phylogenetic tree using
#' the S.PhyloMaker tool, and organized into a \code{\link[phyloseq]{phyloseq-class}} object.
#'
#' @format A \code{\link[phyloseq]{phyloseq-class}} object containing:
#' \describe{
#'   \item{OTU table}{Species abundance counts, with rows as samples (sites) and columns as species identified by their trinomial names.}
#'   \item{phylogenetic tree}{A phylo object representing the evolutionary relationships among species in the dataset.}
#'   \item{sample data}{Metadata for each forest plot site, including geographic coordinates and other site attributes.}
#'   \item{taxonomy table}{Taxonomic classifications (family, genus, species) matched using the \code{\link[taxonlookup]{taxonlookup}} package.}
#' }
#'
#' @details
#' This dataset was constructed by combining raw species presence data from the Gentry forest plot dataset
#' with phylogenetic information and cleaned taxonomic names. The species names were standardized by removing
#' uncertain annotations such as "cf." and "aff." and mapped to a phylogenetic tree using S.PhyloMaker
#' with scenario S1. Species not matched to the phylogeny (less than 2% of counts) were removed from the
#' final dataset.
#'
#' The data can be used for analyses involving phylogenetic diversity,
#' community ecology, and evolutionary relationships among species in
#' forest plots.
#'
#' @source
#' Raw data originally from:
#' \url{https://github.com/zdealveindy/anadat-r/tree/master/data}
#'
#' Phylogenetic tree and nodes from:
#' \url{https://github.com/jinyizju/S.PhyloMaker}
#'
#' Taxonomic lookups performed using the \code{\link[taxonlookup]{taxonlookup}} package.
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data(gentry)
#'
#' # View the OTU table (species counts)
#' otu_table(gentry)
#'
#' # Plot the phylogenetic tree
#' plot(phy_tree(gentry))
#'
#' # Summarize sample metadata
#' sample_data(gentry)
#'
#' # Extract taxonomic classifications
#' tax_table(gentry)
#' }
#'
#' @keywords datasets phyloseq community ecology forest
"gentry"

#' Subsampled OTU Table from the Gentry Dataset
#'
#' A small subset of the OTU (species abundance) table from the \code{gentry} phyloseq object,
#' containing counts for 50 randomly selected taxa (out of the original 585). This object
#' is primarily intended for use in examples and demonstrations, and is not recommended
#' for substantive ecological or phylogenetic analyses.
#'
#' @format A numeric matrix with 50 columns (taxa) and rows corresponding to samples/sites.
#'
#' @details
#' The subset was created by randomly sampling 50 taxa indices and subsetting
#' the OTU table from the \code{gentry} dataset using \code{subset_taxa()}.
#' The resulting OTU table was coerced to a matrix for convenience.
#'
#' @seealso \code{\link{gentry}}
#'
#' @examples
#' \dontrun{
#' dim(small_otutab)
#' head(small_otutab)
#' }
#'
#' @keywords datasets example data
"small_otutab"


#' Subsampled Phylogenetic Tree from the Gentry Dataset
#'
#' A phylogenetic tree corresponding to the 50 randomly selected taxa used in \code{small_otutab},
#' extracted from the \code{gentry} phyloseq object. This tree is intended solely for
#' illustrative or testing purposes and should not be used for formal phylogenetic analysis.
#'
#' @format A \code{\link[ape]{phylo}} object representing the phylogenetic relationships of the
#' 50 selected taxa.
#'
#' @details
#' The tree was extracted by subsetting the original phylogenetic tree in \code{gentry}
#' to include only the taxa corresponding to the sampled subset used in \code{small_otutab}.
#'
#' @seealso \code{\link{gentry}}, \code{\link{small_otutab}}
#'
#' @examples
#' \dontrun{
#' plot(small_tree)
#' }
#'
#' @keywords datasets example data phylogenetics
"small_tree"
