#' Runs the GPCM model for use in the Information Matrix Test.
#' @param data A dataset with J columns and n rows.
#' @param constraint Constraint is either "1PL", "rasch" or "gpcm" depending on which parameter constraints should be run.
#' @param IRT.param logical; if TRUE then the usual IRT parametrization is used.
#' @param start.val If not Null, a list of starting values for the parameter estimates
#' @param na.action the na.action to be used on the data
#' @param control See gpcm function in ltm package for details.
#' @return A GPCM object.
#' @examples
#' data(dataset)
#' model = gpcm_IMT(dataset, constraint = "rasch")
#' # See vignette("IMT-vignette") for more examples
#' @import ltm
#' @importFrom stats as.formula binomial coef dnorm glm.fit model.matrix nlminb optim qlogis rnorm
#' @export
gpcm_IMT <-
function (data, constraint = c("gpcm", "1PL", "rasch"), IRT.param = TRUE,
          start.val = NULL, na.action = NULL, control = list()) {
    cl <- match.call()
    if ((!is.data.frame(data) & !is.matrix(data)) || ncol(data) == 1)
        stop("'data' must be either a numeric matrix or a data.frame, with at least two columns.\n")
    constraint <- match.arg(constraint)
    X <- if (!is.data.frame(data)) as.data.frame(data) else data
    X[] <- lapply(X, factor)
    ncatg <- as.vector(sapply(X, function (x) length(levels(x))))
    X <- sapply(X, unclass)
    if (!is.null(na.action))
       X <- na.action(X)
    colnamsX <- colnames(X)
    dimnames(X) <- NULL
    p <- ncol(X)
    pats <- apply(X, 1, paste, collapse = "/")
    freqs <- table(pats)
    nfreqs <- length(freqs)
    obs <- as.vector(freqs)
    X <- unlist(strsplit(cbind(names(freqs)), "/"))
    X[X == "NA"] <- as.character(NA)
    X <- matrix(as.numeric(X), nfreqs, p, TRUE)
    XX <- lapply(1:p, function (j) outer(X[, j], seq(1, ncatg[j] - 1), ">") * 1)
    con <- list(iter.qN = 150, GHk = 21, optimizer = "nlminb", optimMethod = "BFGS", numrDeriv = "fd",
        epsHes = 1e-06, parscale = NULL, verbose = getOption("verbose"))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    GH <- GHpoints(data ~ z1, con$GHk)
    Z <- GH$x[, 2]
    GHw <- GH$w
    start.val = NULL
    IRT.param = TRUE
    init.thetas <- start.val.gpcm(start.val, X, obs, constraint, ncatg, IRT.param)
    environment(loglikgpcm) <- environment(scoregpcm) <- environment(scoregpcm_IMT) <- environment()
    res.qN <- if (con$optimizer == "optim") {
        if (is.null(con$parscale) || length(con$parscale) != length(init.thetas))
            con$parscale <- rep(0.5, length(init.thetas))
        optim(init.thetas, loglikgpcm, scoregpcm, constraint = constraint, method = con$optimMethod, gr="BFGS", control = list(maxit = con$iter.qN, parscale = con$parscale, trace = as.numeric(con$verbose)))
    } else {
        nlb <- nlminb(init.thetas, objective = loglikgpcm, gradient = scoregpcm, constraint = constraint,
                      control = list(iter.max = con$iter.qN, eval.max = con$iter.qN, trace = as.numeric(con$verbose),
                      rel.tol = sqrt(.Machine$double.eps)))
        names(nlb) <- c("par", "value", "convergence", "message", "iterations", "counts")
        nlb
    }
    Hess <- fd.vec(res.qN$par, scoregpcm, constraint = constraint, eps = con$epsHes)

    eps=1e-12
    x=res.qN$par
    f0 <- scoregpcm_IMT(x, constraint = constraint)
    np = length(res.qN$par)
    Hessres <- array(0, c(dim(data)[[1]],np,np))
    ex <- pmax(abs(x), 1)
    for (i in 1:np) {
      x. <- x
      x.[i] <- x[i] + eps * ex[i]
      diff.f <- scoregpcm_IMT(x., constraint = constraint) - f0
      diff.x <- x.[i] - x[i]
      Hessres[,,i] <- diff.f / diff.x
    }


    res.qN$hessian <- 0.5 * (Hess + t(Hess))

    score = scoregpcm(res.qN$par, constraint = constraint)
    score_IMT = scoregpcm_IMT(res.qN$par, constraint = constraint)

    if (all(!is.na(res.qN$hessian) & is.finite(res.qN$hessian))) {
        ev <- eigen(res.qN$hessian, TRUE, TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1])))
            warning("Hessian matrix at convergence is not positive definite; unstable solution.\n")
    }
    thetas <- betas.gpcm(res.qN$par, p, ncatg, constraint)


    names(thetas) <- if (!is.null(colnamsX)) colnamsX else paste("Item", 1:p)
    thetas <- lapply(thetas, function (x) { names(x) <- c(paste("Catgr.", seq(1, length(x) - 1), sep = ""), "Dscrmn"); x })
    max.sc <- max(abs(scoregpcm(res.qN$par, constraint)), na.rm = TRUE)
    fit <- list(coefficients = thetas, log.Lik = -res.qN$value, convergence = res.qN$conv, hessian = res.qN$hessian,
                counts = res.qN$counts, patterns = list(X = X, obs = obs), GH = list(Z = Z, GHw = GHw),
                max.sc = max.sc, constraint = constraint, IRT.param = IRT.param, X = data, control = con,
                na.action = na.action, score = score, score_IMT = score_IMT, Hess_IMT=Hessres, call = cl)
    class(fit) <- "gpcm"
    fit
}
