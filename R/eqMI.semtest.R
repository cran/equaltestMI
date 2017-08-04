#' Measurement invariance tests using Lavaan
#'
#' Conventional multiple-group SEM to test measurement invariance. A sequence of chi-square and chi-square difference tests will be conducted.
#'
#' @param ... The same arguments as for any lavaan model. See \code{lavaan::cfa} for more information.
#'
#' Users must explicitly specify the name of the input elements for this function to catch. For example, specify 'data = HolzingerSwineford' instead just 'HolzingerSwineford'.
#' @param output If the function prints out results of covariance structure, mean structure, or both. The value of \code{output} must be \code{mean}, \code{covariance}, or \code{both}. When the tests involve mean structure (\code{output = 'mean' or 'both'}), both the strong and the strict tests of measurement invariance will be conducted.
#' @param quiet If FALSE (default), a summary is printed out containing an overview of the different models that are fitted, together with some model comparison tests. If TRUE, no summary is printed but results will be stored in the object.
#' @return A list is returned with:
#' \describe{
#' \item{\code{LavaanOut}}{A sublist in convention.sem. Contains Lavann style output and results for each chi-square and chi-square difference.}
#' \item{\code{Mean.part}}{A sublist in convention.sem. Contains test statistics and fit measures on invariance tests of mean structure.}
#' \item{\code{Cov.part}}{A sublist in convention.sem. Contains test statistics and fit measures on invariance tests of covariacne structure.}
#' }
#' @details This is a wrapper function of the function \code{measurementInvariance} in package \code{semTools} with the following default options: \code{std.lv = FALSE, fit.measures = "default"}, and \code{method = "satorra.bentler.2001"}. See \code{semTools::measurementInvariance} for more information.
#' @seealso \code{\link[lavaan]{sem}}, \code{\link[semTools]{measurementInvariance}}
#' @references Yuan, K. H., & Chan, W. (2016). Measurement invariance via multigroup SEM: Issues and solutions with chi-square-difference tests. Psychological methods, 21(3), 405-426.
#' @references Yves Rosseel (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36.
#' @references semTools Contributors. (2016). semTools: Useful tools for structural equation modeling. R package version 0.4-14. Retrieved from https://CRAN.R-project.org/package=semTools
#' @importFrom lavaan sem
#' @importFrom lavaan partable
#' @importFrom lavaan anova
#' @export
#' @examples
#' data(HolzingerSwineford)
#' semmodel<-'
#' L1 =~ V1 + V2 + V3
#' L2 =~ V4 + V5 + V6
#' L3 =~ V7 + V8
#' L4 =~ V9 + V10 + V11
#' '
#' \dontrun{
#' run.sem <- eqMI.semtest(model = semmodel, data = HolzingerSwineford,
#'            group = "school", meanstructure = TRUE)
#'}
eqMI.semtest <- function(..., output = 'both', quiet = FALSE) {

  lavaansem <- function(...) { lavaan::sem(...)}
  # check for a group.equal argument in ...
  dotdotdot <- list(...)
  if(!is.null(dotdotdot$group.equal)){
    stop("lavaan ERROR: group.equal argument should not be used")
  }

  # obtain samplestat by group
  if(is.null(dotdotdot$sample.cov) & is.null(dotdotdot$data))
    stop("sample covariances must be provided")
  group1 <- group2 <- dotdotdot
  group1$group <- group2$group <- NULL
  if (!is.null(dotdotdot$sample.cov)){
    group1$sample.cov <- dotdotdot$sample.cov[[1]]
    group2$sample.cov <- dotdotdot$sample.cov[[2]]
    group1$sample.nobs <- dotdotdot$sample.nobs[[1]]
    group2$sample.nobs <- dotdotdot$sample.nobs[[2]]
  }
  if (!is.null(dotdotdot$sample.mean)){
    group1$sample.mean <- dotdotdot$sample.mean[[1]]
    group2$sample.mean <- dotdotdot$sample.mean[[2]]
    }
  if (!is.null(dotdotdot$data)) {
    group.ind <- match(dotdotdot$group, colnames(dotdotdot$data))
    dat.by.group <- split(dotdotdot$data[,-group.ind], f = dotdotdot$data[,group.ind])
    group1$data <- dat.by.group[[1]]
    group2$data <- dat.by.group[[2]]
    }

  res <- list()

  configural <- dotdotdot
  configural$group.equal <- ""
  template <- do.call(lavaansem, configural)
  pttemplate <- lavaan::partable(template)
  varnames <- unique(pttemplate$rhs[pttemplate$op == "=~"])
  facnames <- unique(pttemplate$lhs[(pttemplate$op == "=~") & (pttemplate$rhs %in% varnames)])
  ngroups <- max(pttemplate$group)
  if(ngroups <= 1) stop("Well, the number of groups is 1. Measurement invariance across 'groups' cannot be done.")

  # base-line model: configural invariance
  res$fit.configural.g1 <- do.call(lavaansem, group1)
  res$fit.configural.g2 <- do.call(lavaansem, group2)
  res$fit.combine.groups <- template

  # fix loadings across groups
  loadings <- dotdotdot
  loadings$group.equal <- c("loadings")
  res$fit.metric <- do.call("cfa", loadings)


  # fix loadings and residuals across groups
  residuals <- dotdotdot
  residuals$group.equal <- c("loadings", "residuals")
  res$fit.residuals <- do.call("cfa", residuals)


  # fix loadings, residuals, and factor covariances across groups
  res.var <- dotdotdot
  res.var$group.equal <- c("loadings", "residuals", "lv.variances", "lv.covariances")
  res$fit.varfactor <- do.call("cfa", res.var)

  # fix loadings + intercepts across groups
  intercepts <- dotdotdot
  intercepts$group.equal <- c("loadings", "intercepts")
  res$fit.scalar <- do.call(lavaansem, intercepts)


  # fix loadings + intercepts + means across groups (strong MI)
  means <- dotdotdot
  means$group.equal <- c("loadings", "intercepts", "means")
  res$fit.strong.means <- do.call(lavaansem, means)

  # fix loadings + intercepts + residuals + means across groups (strict MI)
  residuals <- dotdotdot
  residuals$group.equal <- c("loadings", "intercepts", "residuals")
  res$fit.strict.residuals <- do.call(lavaansem, residuals)

  # fix loadings + residuals + intercepts + means
  means <- dotdotdot
  means$group.equal <- c("loadings", "intercepts", "residuals", "means")
  res$fit.strict.means <- do.call(lavaansem, means)


  Mean.part <- Cov.part <- NULL
  fit.measures <- "default"; method <- "satorra.bentler.2001"
  if(output == 'covariance'){
    if(!quiet) {cat('\n', '------------------ Covariance Structure ------------------ ', '\n')}
    Cov.part <- printInvarianceResult(res[1:6], fit.measures, method, quiet)
  } else if(output == 'mean'){
    if(!quiet) {cat('\n', '------------------ Mean Structure ------------------ ', '\n')}
    Mean.part <- printInvarianceResult(res[-c(5:6)], fit.measures, method, quiet)
  } else if(output == 'both'){
    if(!quiet) {cat('\n', '------------------ Covariance Structure ------------------ ', '\n')}
    Cov.part <- printInvarianceResult(res[1:6], fit.measures, method, quiet)
    if(!quiet) {cat('\n', '------------------ Mean Structure ------------------ ', '\n')}
    #cat('NOTE: the chisquare difference test of <fit.strict.residuals> and <fit.strong.means> is not meaningful because they are not nested.', '\n')
    Mean.part <- printInvarianceResult(res[-c(1:2, 5:6)], fit.measures, method, quiet)
  }

    invisible(res)
    return(list(LavaanOut = res, Mean.part = Mean.part, Cov.part = Cov.part))
}

#' Internal function to format the results of measurement invariance tests
#'
#' @param FIT parameter needed by printInvarianceResult().
#' @param fit.measures parameter needed by printInvarianceResult().
#' @param method parameter needed by printInvarianceResult(). Default at "satorra.bentler.2001".
#' @param quiet If \code{quiet=FALSE} (default), a summary is printed out containing an overview of the different models that are fitted, together with some model comparison tests and fit measures. The results of equivalence testing will also be printed if equivalence testing is used. If \code{quiet=TRUE}, no summary is printed but results will be stored in the object.

printInvarianceResult <- function(FIT, fit.measures, method, quiet = FALSE) {

  if ('fit.strict.residuals'%in%names(FIT)){
    keep1 <- anova(FIT$fit.strict.residuals, FIT$fit.scalar)
  }

  # compare models
  NAMES <- names(FIT); names(FIT) <- NULL
  lavaanLavTestLRT <- function(...) { lavaan::lavTestLRT(...) }
  TABLE <- do.call(lavaanLavTestLRT, c(FIT, list(model.names = NAMES, method = method)))

  #not meaningful comparision
  ind <- match(c('fit.configural.g1', 'fit.configural.g2', 'fit.combine.groups'), rownames(TABLE))
  TABLE[ind[!is.na(ind)], c('Chisq diff', 'Df diff', 'Pr(>Chisq)')] <- NA
  if ('fit.strict.residuals'%in%rownames(TABLE)){
    TABLE[match('fit.strict.residuals', rownames(TABLE)), ] <- keep1[2,]
  }

  if(length(fit.measures) == 1L && fit.measures == "default") {
    if(length(lavaan::lavInspect(FIT[[1]], "test")) > 1L) {
      fit.measures <- c("cfi.scaled", "rmsea.scaled")
    } else {
      fit.measures <- c("cfi", "rmsea")
    }
  }

  # add some fit measures
  if(length(fit.measures)) {

    FM <- lapply(FIT, lavaan::fitMeasures, fit.measures)
    FM.table1 <- sapply(fit.measures, function(x) sapply(FM, "[[", x))
    if(length(FM) == 1L) {
      FM.table1 <- rbind( rep(as.numeric(NA), length(fit.measures)),
                          FM.table1 )
    }
    if(length(FM) > 1L) {
      FM.table2 <- rbind(as.numeric(NA),
                         abs(apply(FM.table1, 2, diff)))
      colnames(FM.table2) <- paste(colnames(FM.table2), ".delta", sep="")
      FM.TABLE <- as.data.frame(cbind(FM.table1, FM.table2))
    } else {
      FM.TABLE <- as.data.frame(FM.table1)
    }
    rownames(FM.TABLE) <- rownames(TABLE)
    class(FM.TABLE) <- c("lavaan.data.frame", "data.frame")
  }

  #not meaningful comparision
  ind <- match(c('fit.configural.g1', 'fit.configural.g2', 'fit.combine.groups'), rownames(FM.TABLE))
  FM.TABLE[ind[!is.na(ind)], c('cfi.delta', 'rmsea.delta')] <- NA
  if ('fit.strict.residuals'%in%rownames(FM.TABLE)){
    ind1 <- match('fit.scalar', rownames(FM.TABLE))
    ind2 <- match('fit.strict.residuals', rownames(FM.TABLE))
    FM.TABLE[ind2, 3:4] <- abs(FM.TABLE[ind1, 1:2] - FM.TABLE[ind2, 1:2])
  }

  if(!quiet){
    cat("\n")
    cat("Measurement invariance models \n\n")
    cat(paste(paste("Model", seq_along(FIT), ":", NAMES), collapse = "\n"))
    cat("\n\n")

    print(TABLE)
  }

  if(length(fit.measures)) {
    if(!quiet){
      cat("\n", "Fit measures \n\n")
      pFM <- FM.TABLE
      pFM[] <- sprintf("%8.3f",unlist(FM.TABLE))
      pFM[is.na(FM.TABLE)] <- ''
      print(pFM)
      cat("\n")
    }
    return(data.frame(TABLE, FM.TABLE))
  } else {
    return(TABLE)
  }
}

