#' Compute MPQr distances for a range of values of r
#' @param X A n x p data matrix.
#' @param tr A tree (of class \code{phylo}) with p leaves.
#' @param rvec A vector of values between 0 and 1. For each value in
#'     this vector, the function will compute a distance matrix
#'     corresponding to the MPQr distance for that value of
#'     r. Defaults to a set of values of r that are spaced more
#'     closely close to 0 and 1 and more widely around .5 by applying
#'     \code{r_transform} to equally spaced values between 0 and 1.
#' @return A list of length 2 whose elements are named \code{distances} and \code{r}. \code{distances} is a list of length \code{length(rvec)}, and each element contains a distance object giving the MPQr distances between the rows of \code{X} for one value of r. The second element of the list (\code{r}) gives the values of r that were used for computing the MPQr distances. These values are given in the same order as the distance objects in the \code{distance} element of the list.
#' @examples
#' data(small_otutab)
#' data(small_tree)
#' mpq_distances <- get_mpq_distances(small_otutab, small_tree, rvec = c(0,.5,1))
#' class(mpq_distances$distance[[1]])
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
        evals = get_dr(Qeig, r)
        XVDsqrt = XV * matrix(evals^(.5), nrow = nrow(XV), ncol = ncol(XV), byrow = TRUE)
        distanceList[[as.character(r)]] = dist(XVDsqrt)
    }
    out = list(distances = distanceList, r = rvec)
    return(out)
}


#' Compute normalized eigenvalues for the MPQr distance
#'
#' Given an object containing the eigendecomposition of the matrix \eqn{\T} used for the PQ distance, , this function computes the vector of normalized weights corresponding to the parameter \eqn{r \in [0,1]}.
#'
#' @param eig An 'eigen' object. Should be obtained by calling 'eigen' on the covariance matrix obtained from a phylogenetic tree T.
#' @param r Numeric scalar in \eqn{[0,1]} specifying the tuning parameter of the \eqn{\MPQr} distance.
#'
#' @return A numeric vector of length \code{length(eig$values)}.
#'
#' @keywords internal
#' @examples
#' eig <- eigen(matrix(c(2,1,1,2), 2, 2))
#' # when r = 0, the eigenvalues are equal
#' mpqDist:::get_dr(eig, r = 0)
#' mpqDist:::get_dr(eig, r = 0.5)
#' # when r = 1, the eigenvalues are proportional to the original eigenvalues
#' mpqDist:::get_dr(eig, r = 1)
get_dr <- function(eig, r) {
    eig$values = eig$values / sum(eig$values) * length(eig$values)
    if(r == 1) {
        evals = eig$values
    } else if (r == 0) {
        evals = rep(1, length(eig$values))
    } else {
        evals = (rep(1/r, length(eig$values)) + (1-r)^(-1) * eig$values^(-1))^(-1)
    }
    evals = evals / sum(evals)
    return(evals)
}


#' A transformation that tends to work better than a linear grid between 0 and 1
#' @param r A value between 0 and 1
#' @return Another value between 0 and 1.
#' @keywords internal
r_transform <- function(r) {
    .5 - .5 * cos(r * pi)
}

#' Set up a data frame to make a plotly animation
#' @param mpq_distances The output from calling
#'     \code{get_mpq_distances}. Should be a list with elements
#'     \code{distances} and \code{r}, giving a list of distance
#'     objects and the corresponding values of r for the MPQr
#'     distance, respectively.
#' @param distance_list A list of distance objects, for example, the
#'     \code{distance} element from the output of
#'     \code{get_mpq_distances}.
#' @param fn A function to apply to each of the distance objects that
#'     creates a data frame with the variables in the animation. Its
#'     first argument should be a distance object. For example,
#'     \code{get_avg_distances_to_set}.
#' @param means_and_vars Optional, defaults to NULL. If non-null,
#'     should be the output from the function
#'     \code{get_null_mean_and_variance}. It will add columns
#'     \code{null_mean}, \code{null_sd}, \code{median}, \code{lower},
#'     and \code{upper}, corresponding to those statistics of the null
#'     distribution for each value of r. Can be used to augment the
#'     plots with expected values under a null model.
#' @param ... Extra arguments that are passed to fn.
#' @return A data frame that has one column called \code{frame} (which
#'     can be used by ggplot2 and plotly to make an animation).
#' @examples
#' data(small_otutab)
#' data(small_tree)
#' mpq_distances <- get_mpq_distances(small_otutab, small_tree)
#' # we give make_animation_df the mpq distances and a function we
#' # would like to apply to every distance matrix. The function
#' # requires an extra argument in addition to the distance matrix, so
#' # we pass that as `base_samples = c(1,2)` to specify that we want
#' # to call `get_avg_distances_to_set` on each distance matrix with
#' # the base samples set to be 1 and 2.
#' animation_df = make_animation_df(mpq_distances = mpq_distances, fn = get_avg_distances_to_set, base_samples = c(1,2))
#' # the resulting data frame has one column (frame) that corresponds
#' # to the different values of r and which can be used for animating
#' # the plot. The remaining columns are the columns created by `fn`
#' # (in this example, `get_avg_distances_to_set`.)
#' head(animation_df)
#' 
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
                                                median = means_and_vars$median[i],
                                                lower = means_and_vars$lower[i],
                                                upper = means_and_vars$upper[i],
                                                fn(distance_list[[i]], ...))
        }
    }
    return(dplyr::bind_rows(animation_df_list))
}

#' Get average distance to set
#'
#' Computes the average distance between each individual point and a
#' set of base samples.
#'
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
#'     \ifelse{html}{
#'        \out{\(\sum_{bs \in \text{base_samples}} d(i, bs) / |\text{base_samples}|\)}
#'     }{
#'        \eqn{sum_{bs \in \text{base_samples}} d(i, bs) / |\text{base_samples}|}
#'     }. Remaining
#'     columns in the data frame are sample names and sample_data (if
#'     it was provided.)
#' @examples
#' data(small_otutab)
#' data(small_tree)
#' mpqr_distances = get_mpq_distances(small_otutab, small_tree, rvec = c(0, .5, 1))
#' # get the MPQr distances with r = .5
#' d_obj = mpqr_distances$distances[["0.5"]]
#' ## for each site, get the average of the distances between that site and the base sites (sites 1 and 2)
#' get_avg_distances_to_set(d_obj, base_samples = c(1,2)) |> head()
#' ## compare to doing this by hand for site 3:
#' (d_3_12 = as.matrix(d_obj)[3, c(1,2)])
#' mean(d_3_12)
#' mean(d_3_12) == get_avg_distances_to_set(d_obj, base_samples = c(1,2))$avg_dist[3]
#' @export
get_avg_distances_to_set <- function(d_obj, base_samples, sample_data = NULL) {
    d_mat = as.matrix(d_obj)
    d_mat_restricted = d_mat[base_samples,, drop = FALSE]
    avg_distances = colMeans(d_mat_restricted)
    out = data.frame(avg_dist = avg_distances, site = colnames(d_mat_restricted))
    rownames(out) = out$site
    if(!is.null(sample_data)) {
        out = cbind(out, sample_data)
    }
    return(out)
}


#' Makes an animated plot showing the the distances between each
#' sample and a set of reference samples
#' @param physeq A phyloseq object containing species abundances, a
#'     phylogenetic tree, and sample data
#' @param reference_samples Either the names or the indices of the set
#'     of reference samples.
#' @param x_variable A string giving the variable for plotting the
#'     distances against on the x axis. Should be the name of one of
#'     the columns of `sample_data(phy_seq)`.
#' @param return_ggplot If TRUE, the function will return a ggplot
#'     object which can be made into an animated plot by calling
#'     ggplotly on it. If FALSE, the function will create the animated
#'     plot. The option of returning the ggplot object is primarily to
#'     allow for more customization by the user.
#' @param ... Extra objects to be passed to `aes_string` in ggplot.
#' @return Either an animated plot created by `ggplotly` (if
#'     `return_ggplot = FALSE`) or a ggplot object (if `return_ggplot
#'     = TRUE`)
#' @examples
#' data(gentry)
#' gentry_log = gentry
#' otu_table(gentry_log) = log(otu_table(gentry) + 1)
#' equatorial_samples = sample_names(gentry_log)[which(sample_data(gentry_log)$Lat >= -10 & sample_data(gentry_log)$Lat <= 10)]
#' plot_dist_to_reference_samples(physeq = gentry_log, reference_samples = equatorial_samples,x_variable = "Lat")
#' @importFrom ggplot2 ggplot aes_string geom_point
#' @importFrom plotly ggplotly
#' @export
plot_dist_to_reference_samples <- function(physeq, reference_samples, x_variable, return_ggplot = FALSE, ...) {
    distances = get_mpq_distances(otu_table(physeq), phy_tree(physeq))
    animation_df = make_animation_df(distances,
                                     get_avg_distances_to_set,
                                     means_and_vars = NULL,
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
#' @param d_obj A distance object containing the MPQr distances
#'     between the samples.
#' @param sample_distances A distance object containing the distances
#'     between samples.
#' @keywords internal
#' @return A data frame
mpq_and_sample_distances <- function(d_obj, sample_distances) {
    df = data.frame(mpq_dist = as.vector(d_obj), sample_dist = as.vector(sample_distances))
    return(df)
}

#' Creates an animated plot of MPQ distances against some other
#' measure of distance between sites
#' @param physeq A phyloseq object containing species abundances, a
#'     phylogenetic tree, and sample data.
#' @param sample_distances Either a distance object or a (n_samples x
#'     n_samples) matrix containing distances between the samples.
#' @param return_ggplot If TRUE, the function will return a ggplot
#'     object which can be made into an animated plot by calling
#'     ggplotly on it. If FALSE, the function will create the animated
#'     plot. The option of returning the ggplot object is primarily to
#'     allow for more customization by the user.
#' @return Either an animated plot created by `ggplotly` (if
#'     `return_ggplot = FALSE`) or a ggplot object (if `return_ggplot
#'     = TRUE`)
#' @examples
#' data(gentry)
#' gentry_small = subset_samples(gentry, Precip > 4000)
#' longitude_distance = dist(sample_data(gentry_small)$Long)
#' plot_all_pairs(gentry_small, longitude_distance)
#' @importFrom phyloseq otu_table phy_tree nsamples
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom plotly ggplotly
#' @export
plot_all_pairs <- function(physeq, sample_distances,
                           return_ggplot = FALSE) {
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
                                     means_and_vars = NULL,
                                     sample_distances)
    p = ggplot(aes(x = sample_dist, y = mpq_dist, frame = frame), data = animation_df) +
        geom_point()
    if(return_ggplot) {
        return(p)
    }
    ggplotly(p)
}


#' Compute null distribution statistics for the MPQr distance
#'
#' For a given dataset, encoded as a phyloseq object, and a set of r
#' values, this function computes approximate means, standard
#' deviations, medians, and selected quantiles of the MPQr distance
#' under a null model in which sample values are independent Gaussian
#' random variables with the same variance as in the observed data.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq}} object containing:
#'   - a phylogenetic tree (accessible via \code{\link[phyloseq]{phy_tree}})
#'   - an OTU/ASV abundance table (accessible via \code{\link[phyloseq]{otu_table}})  
#'   The abundance table is used to compute per-taxon standard deviations.
#'
#' @param rvec A numeric vector of r values in [0,1] for which to
#'   compute null distribution summaries. Defaults to \code{r_transform(0:100/100)}.
#'
#' @return
#' A named list of numeric vectors (each of length \code{length(rvec)}). Elements of the list give some values related to the distribution of the MPQr distance between two samples under a null distribution.
#' Under the null distribution, each sample has independent, normally-distributed elements with standard deviation equal to the empirical standard deviation of the observed variables.
#' The statistics are as follows:
#' - \code{means}: Approximate means of the MPQr distance between two samples in the independent-normal null model.
#' - \code{sds}: Approximate standard deviations of the MPQr distance between two samples in the independent normal null model.
#' - \code{medians}: Median of the null distribution.
#' - \code{lower}: Lower 95% confidence bound (2.5% quantile of the null distribution)
#' - \code{upper}: Upper 95% confidence bound (97.5% quantile of the null distribution)
#'
#' @examples
#' data(gentry)
#' get_null_mean_and_variance(gentry, rvec = c(0, 0.5, 1))
#'
#' @importFrom ape vcv
#' @importFrom phyloseq phy_tree ntaxa otu_table
#' @export
get_null_mean_and_variance <- function(physeq, rvec = r_transform(0:100/100)) {
    Q = vcv(phy_tree(physeq), scale = FALSE)
    Qeig = eigen(Q)
    p = ntaxa(physeq)
    evecs = Qeig$vectors
    sds = apply(otu_table(physeq), 2, sd)
    means = numeric(length(rvec))
    sd_vec = numeric(length(rvec))
    medians = numeric(length(rvec))
    lower = numeric(length(rvec))
    upper = numeric(length(rvec))
    for(i in seq_along(rvec)) {
        r = rvec[i]
        evals = get_dr(Qeig, r)
        Drsqrt = sqrt(evals)
        sigmaV = evecs * sds
        sigmaVDrsqrt = t(t(sigmaV) * Drsqrt)
        svd_sigmaVDrsqrt = svd(sigmaVDrsqrt)
        weights = svd_sigmaVDrsqrt$d^2
        ## mean of mpq^2
        mu2 = sum(weights) * 2
        ## var of mpq^2
        var2 = sum(weights^2) * 4
        means[i] = sqrt(mu2) - .125 * mu2^(-3/2) * var2
        sd_vec[i] = sqrt(var2) / (2 * sqrt(mu2))
        medians[i] = sqrt(gx2inv(.5, weights, h = rep(2, length(weights))))
        lower[i] = sqrt(gx2inv(.025, weights, h = rep(2, length(weights))))
        upper[i] = sqrt(gx2inv(.975, weights, h = rep(2, length(weights))))
    }
    return(list(means = means, sds = sd_vec, medians = medians, lower = lower, upper = upper))
}

#' @importFrom CompQuadForm davies
#' @importFrom stats uniroot
gx2inv <- function(p, lambda, h) {
    mu = sum(lambda * h)
    var = 2 * sum(h * lambda^2)
    root = uniroot(
                      function(x) 1 - davies(x, lambda, h)$Qq - p,
                      interval = c(mu - 5 * sqrt(var), mu + 5 * sqrt(var))
                  )$root
    return(root)
}

