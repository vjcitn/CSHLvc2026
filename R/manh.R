make_manh_df = function(sid, catalog) {
  kp = which(catalog$`STUDY ACCESSION` == sid)
  stopifnot(length(kp)>2)
  catdf = as.data.frame(catalog[kp])
  catdf$pos = as.numeric(catdf$CHR_POS)
  catdf$chr = as.numeric(catdf$CHR_ID)
  catdf$pval = catdf$P.VALUE
  catdf
}

#' create a manhattan plot given a study accession and a catalog
#' @param sid string STUDY ACCESSION
#' @param catalog instance of gwasWloc from gwascat package
#' @param \dots passed to ggmanh manhattan_plot
#' @examples
#' data("gwc_110626", package="CSHLvc2026")
#' make_manh("GCST90002322". gwc_110626)
#' @export
make_manh = function(sid, catalog, ...) {
   df = make_manh_df(sid, catalog)
   ggmanh::manhattan_plot(df, ...)
} 

