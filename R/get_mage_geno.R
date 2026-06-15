#' Validate PLINK .bed magic bytes
#'
#' A well-formed PLINK 1.9 .bed file begins with 0x6c 0x1b 0x01.
#'
#' @param path Character. Local path to the .bed file.
#' @return Invisibly returns \code{TRUE} on success; stops with an informative
#'   error on failure.
#' @keywords internal
.validate_bed_magic <- function(path) {
    con <- file(path, open = "rb")
    on.exit(close(con))
    magic <- readBin(con, what = "raw", n = 3L)
    expected <- as.raw(c(0x6c, 0x1b, 0x01))
    if (!identical(magic, expected)) {
        stop(
            "PLINK .bed magic byte check failed for: ", path, "\n",
            "  Expected: ", paste(expected, collapse = " "), "\n",
            "  Got:      ", paste(magic,    collapse = " "), "\n",
            "The file may be corrupt or a PLINK 1 (pre-1.9) binary."
        )
    }
    invisible(TRUE)
}

#' Download a URL to a local path with optional staleness check
#'
#' Used internally to populate the PLINK subdirectory. Checks the remote
#' Last-Modified / ETag headers when the file already exists; re-downloads
#' only when needed.
#'
#' @param url    Remote URL.
#' @param dest   Destination file path (must already be writable).
#' @param verbose Logical.
#' @return Invisibly returns \code{dest}.
#' @keywords internal
.download_if_needed <- function(url, dest, verbose = TRUE) {
    if (file.exists(dest)) {
        # HEAD request to compare Last-Modified with local mtime
        needs <- tryCatch({
            h <- curlGetHeaders(url)          # base R, no extra deps
            lm_line <- grep("last-modified", h, ignore.case = TRUE, value = TRUE)
            if (length(lm_line) == 0L) {
                if (verbose)
                    message("  (no Last-Modified header; using cached copy)")
                FALSE
            } else {
                remote_time <- as.POSIXct(
                    sub(".*last-modified:\\s*", "", lm_line[[1L]], ignore.case = TRUE),
                    format = "%a, %d %b %Y %H:%M:%S", tz = "GMT"
                )
                local_time <- file.mtime(dest)
                remote_time > local_time
            }
        }, error = function(e) {
            warning("Could not check remote freshness for '", basename(dest),
                    "': ", conditionMessage(e), "\nUsing cached copy.")
            FALSE
        })

        if (!needs) {
            if (verbose) message("[ cache hit   ] ", basename(dest))
            return(invisible(dest))
        }
        if (verbose) message("[ cache stale ] Re-downloading: ", basename(dest))
    } else {
        if (verbose) message("[ cache miss  ] Downloading:    ", basename(dest))
    }

    utils::download.file(url, destfile = dest, mode = "wb", quiet = !verbose)
    invisible(dest)
}

#' Cache PLINK files for MAGE chr17 genotype data
#'
#' Ensures the three PLINK files (\code{.fam}, \code{.bim}, \code{.bed}) from
#' the BiocMAGE17geno OSN bucket are present in a \strong{shared subdirectory}
#' of the BiocFileCache root, so that tools like \code{BEDMatrix} can locate
#' sibling files by stripping the \code{.bed} extension.
#'
#' @param bfc          A \code{BiocFileCache} object. Defaults to the
#'   user-level cache via \code{BiocFileCache()}.
#' @param verbose      Logical. Emit progress messages? Default \code{TRUE}.
#' @param validate_bed Logical. Check \code{.bed} magic bytes? Default
#'   \code{TRUE}.
#'
#' @return A named list with elements \code{fam}, \code{bim}, and \code{bed},
#'   each a character scalar giving the local path to the cached file. All
#'   three share the same directory and stem (\code{CCDG_mage_chr17}), so
#'   \code{BEDMatrix::BEDMatrix(result$bed)} works without further arguments.
#'
#' @importFrom BiocFileCache BiocFileCache bfccache
#' @export
#'
#' @examples
#' plink <- cache_mage_chr17_plink()
#' bed   <- BEDMatrix::BEDMatrix(plink$bed)
cache_mage_chr17_plink <- function(
    bfc          = BiocFileCache::BiocFileCache(),
    verbose      = TRUE,
    validate_bed = TRUE
) {
    base_url  <- paste0(
        "https://mghp.osn.xsede.org/bir190004-bucket01/",
        "BiocMAGE17geno/CCDG_mage_chr17"
    )
    stem      <- "CCDG_mage_chr17"
    exts      <- c(fam = ".fam", bim = ".bim", bed = ".bed")

    # Place all three files in a dedicated subdirectory of the BFC root so
    # their basenames are bare (e.g. CCDG_mage_chr17.bed) with no UUID prefix.
    plink_dir <- file.path(BiocFileCache::bfccache(bfc), "CCDG_mage_chr17_plink")
    dir.create(plink_dir, showWarnings = FALSE, recursive = TRUE)

    paths <- lapply(names(exts), function(ext) {
        url  <- paste0(base_url, exts[[ext]])
        dest <- file.path(plink_dir, paste0(stem, exts[[ext]]))
        .download_if_needed(url, dest, verbose = verbose)
        dest
    })
    names(paths) <- names(exts)

    if (validate_bed) {
        if (verbose) message("[ validate    ] Checking .bed magic bytes ...")
        .validate_bed_magic(paths$bed)
        if (verbose) message("[ validate    ] .bed magic bytes OK.")
    }

    paths
}

