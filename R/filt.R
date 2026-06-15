
#' create a filtered biplot
#' @import ggplot2
#' @param pca prcomp output
#' @param pc 2-vector, which pcs to use
#' @param top_n numeric, number of features to use for directional arrows over projection
#' @param min_norm lower bound on L2 norm of feature loadings to use
#' @export
biplot_filtered <- function(pca, pc = c(1, 2), top_n = NULL, min_norm = NULL) {
    stopifnot(inherits(pca, "prcomp"))

    scores   <- as.data.frame(pca$x[, pc])
    loadings <- as.data.frame(pca$rotation[, pc])
    colnames(scores) <- colnames(loadings) <- c("PC1", "PC2")

    # Scale loadings to scores range (standard biplot scaling)
    scale_factor <- max(abs(scores)) / max(abs(loadings)) * 0.8
    loadings_scaled <- loadings * scale_factor
    loadings_scaled$feature <- rownames(loadings)
    loadings_scaled$norm <- sqrt(loadings_scaled$PC1^2 + loadings_scaled$PC2^2)

    # Filter: top_n by norm, or threshold by min_norm
    if (!is.null(top_n))
        loadings_scaled <- loadings_scaled[order(-loadings_scaled$norm)[seq_len(top_n)], ]
    else if (!is.null(min_norm))
        loadings_scaled <- loadings_scaled[loadings_scaled$norm >= min_norm, ]

    ggplot(scores, aes(PC1, PC2)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_segment(
            data = loadings_scaled,
            aes(x = 0, y = 0, xend = PC1, yend = PC2),
            arrow = arrow(length = unit(0.2, "cm")),
            colour = "firebrick"
        ) +
        geom_text(
            data = loadings_scaled,
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = feature),
            size = 3, colour = "firebrick"
        ) +
        coord_equal() +
        theme_bw()
}

