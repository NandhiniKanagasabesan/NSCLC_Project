
# The source code for runDiff used to perform diffferential expression analysis was altered. The "runDiff" in Cytotree was used to
# perform the differential expression of markers in branches, but it is changed to find the differential expression of markers in clusters

#### rundiff with clusters ####
runDiff <- function(object, cluster.id = NULL, cluster.id.2 = NULL, verbose = FALSE) {

  if (verbose) message(Sys.time(), " Calculating differentially expressed markers.")
  if (missing(object)) stop(Sys.time(), " CYT object is missing.")
  if (!"cluster.id" %in% colnames(object@meta.data)) stop(Sys.time(), " cluster.id is missing, please run buildTree first.")

  all.cluster.ids <- unique(object@meta.data$cluster.id)

  total.deg.list <- NULL
  cluster.contrast <- NULL
  ga <- go <- NULL
  if (length(all.cluster.ids) == 1) {
    stop(Sys.time(), " There is only one cluster in the tree.")
  } else {
    pdata <- object@meta.data[which(object@meta.data$dowsample == 1), c("cell", "cluster.id")]
    edata <- object@log.data[which(object@meta.data$dowsample == 1), object@markers.idx]
    if (is.null(cluster.id) & is.null(cluster.id.2)) {
      if (verbose) message(Sys.time(), " All clusteres will be calculated.")
      for (bid in all.cluster.ids) {
        pdata$contrast <- "go"
        pdata$contrast[which(pdata$cluster.id == bid)] = "ga"
        design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
        colnames(design) <- stringr::str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
        fit <- lmFit(t(edata), design)
        contrast <- makeContrasts(ga_go = ga - go,
                                  levels = design)
        fits <- contrasts.fit(fit, contrast)
        ebFit <- eBayes(fits)

        deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = 'fdr', number = Inf)
        deg_sig_list$cluster.contrast <- paste0(bid, "_vs_other")
        deg_sig_list$Gene <- rownames(deg_sig_list)
        total.deg.list <- rbind(total.deg.list, deg_sig_list)
      }
    } else if (is.null(cluster.id.2)) {
      if (verbose) message(Sys.time(), " Some of clusteres will be calculated.")
      pdata$contrast <- "go"
      pdata$contrast[pdata$cluster.id %in% cluster.id] = "ga"
      design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
      colnames(design) <- str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
      fit <- lmFit(t(edata), design)
      contrast <- makeContrasts(ga_go = ga - go,
                                levels = design)
      fits <- contrasts.fit(fit, contrast)
      ebFit <- eBayes(fits)

      deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = 'fdr', number = Inf)
      deg_sig_list$cluster.contrast <- paste0(paste0(cluster.id, collapse = "-"), "_vs_other")
      deg_sig_list$Gene <- rownames(deg_sig_list)
      total.deg.list <- deg_sig_list
    } else {
      if (verbose) message(Sys.time(), " Some of clusteres will be calculated.")
      pdata$contrast <- "gz"
      pdata$contrast[pdata$cluster.id %in% cluster.id] = "ga"
      pdata$contrast[pdata$cluster.id %in% cluster.id.2] = "go"
      design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
      colnames(design) <- str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
      fit <- lmFit(t(edata), design)
      contrast <- makeContrasts(ga_go = ga - go,
                                levels = design)
      fits <- contrasts.fit(fit, contrast)
      ebFit <- eBayes(fits)

      deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = 'fdr', number = Inf)
      deg_sig_list$cluster.contrast <- paste0(paste0(cluster.id, collapse = "-"), "_vs_", paste0(cluster.id.2, collapse = "-"))
      deg_sig_list$Gene <- rownames(deg_sig_list)
      total.deg.list <- deg_sig_list
    }
  }
  if (verbose) message(Sys.time(), " Calculating differentially expressed markers completed")
  return(total.deg.list)
}


