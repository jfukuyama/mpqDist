## Functions to compute MPQ distances and to make plots

## We will want one function that takes data and a tree and returns a list of distances

## we will want another function to do the animated plots https://plotly.com/ggplot2/animations/

#' Compute MPQr distances for a range of values of r
#' @param X A \eqn{n \times p} data matrix.
#' @param Q A \eqn{p \times p} similarity matrix on the variables
#'     defining an inner product on the rows of \code{X}, can also be
#'     given as an eigendecomposition (formatted as the output from
#'     \code{eigen}).
#' @return A list containing distances between the rows of \code{X},
#'     and the value of \eqn{r} that was used (\code{r}).
#' @examples
#' data(AntibioticSmall)
#' out.agpca = adaptivegpca(AntibioticSmall$X, AntibioticSmall$Q, k = 2)
#' @export
get_mpq_distances <- function(X, Q, rvec = r_transform(0:100/100)) {
    if(is.matrix(Q)) {
        Qeig = eigen(Q, symmetric = TRUE)
    } else if(is.list(Q) & !is.null(Q$vectors) & !is.null(Q$values)) {
        Qeig = Q
    } else {
        stop("Q is not formatted correctly")
    }
    # normalize so that the trace of Q is the same as the trace of the identity matrix
    Qeig$values = ncol(X) * Qeig$values / sum(Qeig$values)
    evecs = Qeig$vectors
    XV = X %*% evecs
    distanceList = list()
    for(r in rvec) {
        if(r == 0) {
            evals = Qeig$values
        } else if (r == 1) {
            evals = rep(1, ncol(X))
        } else {
            evals = (rep(1/(1-r), ncol(X)) + r^(-1) * Qeig$values^(-1))^(-1)
            evals = ncol(X) * evals / sum(evals)
        }
        XVDsqrt = XV * matrix(evals^(.5), nrow = nrow(XV), ncol = ncol(XV), byrow = TRUE)
        distanceList[[as.character(r)]] = dist(XVDsqrt)
    }
    out = list(distances = distanceList, r = rvec)
    return(out)
}

#' A transformation that tends to work better than a linear grid between 0 and 1
#' @param r A value between 0 and 1
#' @return Another value between 0 and 1.
r_transform <- function(r) {
    .5 - .5 * cos(r * pi)
}

#' Sets up a data frame for plotly to make an animation
#' @param distance_list A list of distance objects, for example, the
#'     \code{distance} element from the output of
#'     \code{get_mpq_distances}.
#' @param fn A function to apply to each of the distance objects that
#'     creates a data frame with the variables in the animation. Its
#'     first argument should be a distance object. For example,
#'     \code{get_avg_distances_to_set}.
#' @param ... Extra arguments that are passed to fn.
#' @return A data frame that has one column called \code{frame} (which
#'     can be used by ggplot2 and plotly to make an animation).
#' @export
make_animation_df <- function(distance_list, fn, ...) {
    animation_df_list = list()
    for(i in 1:length(distance_list)) {
        animation_df_list[[i]] = data.frame(frame = i, fn(distance_list[[i]], ...))
    }
    return(dplyr::bind_rows(animation_df_list))
}

#' Average distance between each individual point and a set of base
#' samples
#' @param d_obj A distance object describing the distances between all
#'     the samples
#' @param base_samples The samples against which all the others are
#'     measured.
#' @param sample_data (optional) Extra data that can be added to the
#'     data frame containing the average distances. Useful when you
#'     want to plot average distances against other covariates.
#' @return A data frame with one column containing the average
#'     distances between each sample and the base
#'     samples. Specifically, for sample i, the average distance is
#'     sum_{bs in base_samples} d(i, bs) / |base_samples|. Remaining
#'     columns in the data frame are sample names and sample_data (if
#'     it was provided.)
#' @export
get_avg_distances_to_set <- function(d_obj, base_samples, sample_data = NULL) {
    d_mat = as.matrix(d_obj)
    d_mat_restricted = d_mat[base_samples,]
    avg_distances = colMeans(d_mat_restricted)
    out = data.frame(avg_dist = avg_distances, site = colnames(d_mat_restricted))
    rownames(out) = out$site
    out = cbind(out, sample_data)
    return(out)
}
