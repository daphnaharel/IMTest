#' @import ltm
scoregpcm_IMT <-
function (thetas, constraint) {
    betas <- betas.gpcm(thetas, p, ncatg, constraint)
    eta <- eta1 <- log.crf <- num <- den <- vector("list", p)
    log.p.xz <- matrix(0, nfreqs, length(Z))
    na.ind <- is.na(X)
    for (j in 1:p) {
        bt <- betas[[j]]
        nbt <- length(bt)
        eta1[[j]] <- if (IRT.param) t(outer(Z, bt[-nbt], "-")) else matrix(Z, nbt - 1, length(Z), TRUE)
        eta[[j]] <- if (IRT.param) bt[nbt] * eta1[[j]] else outer(bt[-nbt], bt[nbt] * Z , "+")
        num[[j]] <- exp(apply(eta[[j]], 2, cumsum))
        if (!is.matrix(num[[j]]))
            num[[j]] <- t(num[[j]])
        den[[j]] <- 1 + colSums(num[[j]])
        crf <- rbind(1/den[[j]], num[[j]] / rep(den[[j]], each = ncatg[j] - 1))
        eps <- .Machine$double.eps^(1/2)
        if (any(ind <- crf == 1))
            crf[ind] <- 1 - eps
        if (any(ind <- crf == 0))
            crf[ind] <- eps
        log.pr <- log(crf)
        xj <- X[, j]
        log.pr <- log.pr[xj, ]
        if (any(na.ind[, j]))
            log.pr[na.ind[, j], ] <- 0
        log.p.xz <- log.p.xz + log.pr
    }
    p.xz <- exp(log.p.xz)
    p.x <- c(p.xz %*% GHw)
    p.zx <- p.xz / p.x
    check.alpha <- constraint == "gpcm" || constraint == "1PL"
    if (check.alpha){
        scores.alpha <- numeric(p)
        scores.alpha2 = matrix(NA, nrow = length(xj), ncol = p)
        }
    scores.deltas <- vector("list", p)
    scores.deltas2 <- vector("list", p)
    for (j in 1:p) {
        if (check.alpha) {
            xj <- X[, j]
            etaj <- eta[[j]]
            eta1j <- apply(eta1[[j]], 2, cumsum)
            I1 <- rbind(0, eta1j) - rep(colSums(num[[j]] * eta1j) / den[[j]], each = ncatg[j])
            scores.alpha[j] <- - sum((p.zx * I1[xj, ] * obs) %*% GHw, na.rm = TRUE)
            scores.alpha2[,j] = - (p.zx * I1[xj, ]) %*% GHw
        }
        alpha <- betas[[j]][ncatg[j]]
        ii <- seq(1, ncatg[j] - 1)
        ind1 <- unlist(sapply(ii, ":", b = ncatg[j] - 1))
        ind2 <- rep(ii, rev(ii))
        Part1 <- if (IRT.param) - alpha * XX[[j]] else XX[[j]]
        Part2 <- if (IRT.param) - alpha * rowsum(num[[j]][ind1, , drop = FALSE], ind2) else rowsum(num[[j]][ind1, , drop = FALSE], ind2)
        scores.deltas[[j]] <- lapply(ii, function (i) {
            I2 <- outer(Part1[, i], Part2[i, ] / den[[j]], "-")
            if (any(na.ind[, j]))
                I2[na.ind[, j], ] <- 0
            - sum((p.zx * I2) %*% GHw)
        })
        scores.deltas2[[j]] <- lapply(ii, function (i) {
            I2 <- outer(Part1[, i], Part2[i, ] / den[[j]], "-")
            if (any(na.ind[, j]))
                I2[na.ind[, j], ] <- 0
                - (p.zx*I2) %*% GHw})
    }
	names.alphas = paste("alpha", 1:p)
	inds = rep(1:(ncatg[1]-1), p)
	names.betas = paste(paste("beta", rep(1:p, each = (ncatg[1]-1))), inds)

	scoremat = numeric(0)
	for(r in 1:length(scores.deltas2)){
		scoremat = cbind(scoremat,matrix(unlist(scores.deltas2[[r]]),
                                     nrow = length(xj), ncol = (ncatg[r]-1), byrow=F))
    if(constraint=="gpcm") scoremat = cbind(scoremat,scores.alpha2[,r]);
	}
	scoremat2 = scoremat[rep(1:dim(scoremat)[[1]], obs),]
     scoremat2

}
