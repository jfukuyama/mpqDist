## Functions to compute MPQ distances and to make plots

## We will want one function that takes data and a tree and returns a list of distances

## we will want another function to do the animated plots https://plotly.com/ggplot2/animations/

#' Compute MPQr distances for a range of values of r
#' @param X A \eqn{n \times p} data matrix.
#' @param tr A tree (of class \code{phylo}) with p leaves.
#' @return A list containing distances between the rows of \code{X},
#'     and the value of \eqn{r} that was used (\code{r}).
#' @examples
#' data(AntibioticSmall)
#' out.agpca = adaptivegpca(AntibioticSmall$X, AntibioticSmall$Q, k = 2)
#' @export
get_mpq_distances <- function(X, tr, rvec = r_transform(0:100/100)) {
    Q = ape::vcv(tr, scale = FALSE)
    Qeig = eigen(Q)
    # normalize so that the trace of Q is the same as the trace of the identity matrix
    Qeig$values = Qeig$values / sum(Qeig$values) * ncol(Q)
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
        }
        evals = evals / sum(evals)
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
make_animation_df <- function(mpq_distances, fn, means_and_vars = NULL,  ...) {
    distance_list = mpq_distances$distances
    r_vec = mpq_distances$r
    animation_df_list = list()
    for(i in 1:length(distance_list)) {
        if(is.null(means_and_vars)) {
            animation_df_list[[i]] = data.frame(frame = r_vec[i], fn(distance_list[[i]], ...))
        } else {
            animation_df_list[[i]] = data.frame(frame = r_vec[i],
                                                null_mean = means_and_vars$mean[i],
                                                null_sd = means_and_vars$sds[i],
                                                fn(distance_list[[i]], ...))
        }
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
    d_mat_restricted = d_mat[base_samples,, drop = FALSE]
    avg_distances = colMeans(d_mat_restricted)
    out = data.frame(avg_dist = avg_distances, site = colnames(d_mat_restricted))
    rownames(out) = out$site
    out = cbind(out, sample_data)
    return(out)
}


#' Makes an animated plot showing the the distances between each
#' sample and a set of reference samples
#' @param physeq A phyloseq object containing the abundances, a
#'     phylogenetic tree, and sample data
#' @param reference_samples Either the names or the indices of the set
#'     of reference samples.
#' @param x_variable A variable for plotting the distances against on
#'     the x axis. Should be a string giving the name of one of the
#'     columns of `sample_data(phy_seq)`.
#' @param return_ggplot If TRUE, the function will return a ggplot
#'     object which can be made into an animated plot by calling
#'     ggplotly on the ggplot object. If FALSE, the function will
#'     create the animated plot. The option of returning the ggplot
#'     object is primarily to allow for more customization by the
#'     user.
#' @param ... Extra objects to be passed to `aes_string` in ggplot.
#' @return Either an animated plot created by `ggplotly` (if
#'     `return_ggplot = FALSE`) or a ggplot object (if `return_ggplot
#'     = TRUE`)
#' @export
plot_dist_to_reference_samples <- function(physeq, reference_samples, x_variable, return_ggplot = FALSE, ...) {
    distances = get_mpq_distances(otu_table(physeq), phy_tree(physeq))
    animation_df = make_animation_df(distances,
                                     get_avg_distances_to_set,
                                     reference_samples,
                                     sample_data(physeq))
    p = ggplot(aes_string(x = x_variable, y = "avg_dist", ...),
               data = animation_df) +
        geom_point(aes(frame = frame))
    if(return_ggplot) {
        return(p)
    }
    ggplotly(p)
}

#' Creates a data frame containing MPQr and sample distances
#' @param mpqr_distances A distance object containing the MPQr
#'     distances between the samples.
#' @param sample_distances A distance object containing the distances
#'     between samples.
mpq_and_sample_distances <- function(d_obj, sample_distances) {
    df = data.frame(mpq_dist = as.vector(d_obj), sample_dist = as.vector(sample_distances))
    return(df)
}

#' Creates an animated plot of MPQ distances against some other
#' measure of distance between sites
#' @param physeq A phyloseq object.
#' @param distance_variables Names of columns in sample_data(physeq)
#'     to use to compute distances
#' @param distance_fn A function that takes a matrix and computes
#'     distances between all pairs of rows. It should return an object
#'     of class 'dist'.
#' @param return_ggplot If TRUE, the function will return a ggplot
#'     object which can be made into an animated plot by calling
#'     ggplotly on the ggplot object. If FALSE, the function will
#'     create the animated plot. The option of returning the ggplot
#'     object is primarily to allow for more customization by the
#'     user.
#' @return Either an animated plot created by `ggplotly` (if
#'     `return_ggplot = FALSE`) or a ggplot object (if `return_ggplot
#'     = TRUE`)
#' @export
plot_all_pairs <- function(physeq, sample_distances,
                           return_ggplot = FALSE, ...) {
    mpq_distances = get_mpq_distances(otu_table(physeq), phy_tree(physeq))
    if(is.matrix(sample_distances) &&
       nrow(sample_distances) == nsamples(physeq) &&
       ncol(sample_distances) == nsamples(physeq)) {
        sample_distances = as.dist(sample_distances)
    } else if(class(sample_distances) != "dist") {
        stop("'sample_distances' must be either a nsamples x nsamples matrix or an object of class 'dist'")
    }
    animation_df = make_animation_df(mpq_distances = mpq_distances,
                                     fn = mpq_and_sample_distances,
                                     sample_distances)
    p = ggplot(aes(x = sample_dist, y = mpq_dist, frame = frame), data = animation_df) +
        geom_point()
    if(return_ggplot) {
        return(p)
    }
    ggplotly(p)
}

#' @export
get_null_mean_and_variance <- function(physeq, rvec = r_transform(0:100/100)) {
    Q = ape::vcv(phy_tree(physeq), scale = FALSE)
    Qeig = eigen(Q)
    p = ntaxa(physeq)
    evecs = Qeig$vectors
    sds = apply(otu_table(physeq), 2, sd)
    means = numeric(length(rvec))
    sd_vec = numeric(length(rvec))
    Qeig$values = Qeig$values / sum(Qeig$values) * ncol(Q)
    for(i in seq_along(rvec)) {
        r = rvec[i]
        if(r == 0) {
            evals = Qeig$values
        } else if (r == 1) {
            evals = rep(1, p)
        } else {
            evals = (rep(1/(1-r), p) + r^(-1) * Qeig$values^(-1))^(-1)
        }
        evals = evals / sum(evals)
        Drsqrt = sqrt(evals)
        sigmaV = evecs * sds
        sigmaVDrsqrt = t(t(sigmaV) * Drsqrt)
        svd_sigmaVDrsqrt = svd(sigmaVDrsqrt)
        weights = svd_sigmaVDrsqrt$d^2
        ## mean of mpq^2
        mu2 = sum(weights) * 2
        means[i] = sqrt(mu2)
        ## var of mpq^2
        var2 = sum(weights^2) * 2
        sd_vec[i] = sqrt(var2) / (2 * sqrt(mu2))
        
    }
    return(list(means = means, sds = sd_vec))
}
