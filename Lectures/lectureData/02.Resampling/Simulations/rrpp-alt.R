rrpp.alt <- function(f1, iter = 999, seed = NULL, int.first = FALSE,
                    RRPP = TRUE, SS.type = c("I", "II", "III"),
                    data = NULL, Cov = NULL,
                    print.progress = FALSE, Parallel = FALSE, ...) {
  
  L <- c(as.list(environment()), list(...))
  names(L)[which(names(L) == "f1")] <- "formula"
  
  if(int.first) ko = TRUE else ko = FALSE
  SS.type <- match.arg(SS.type)
  
  dots <- list(...)
  if(length(dots) > 0) {
    w <- dots$weights
    o <- dots$offset
  } else w <- o <- NULL
  
########################
  fit.full <- lm.rrpp(L$formula, data = data, SS.type = SS.type, Cov = Cov,
                print.progress = print.progress, iter= iter)
  if (!is.null(Cov)){
    Pcov <- fit.full$LM$Pcov
    y.rrpp.alt <- lapply(1:(iter+1), function(j) Pcov %*% fit.full$LM$gls.fitted +
                           Pcov %*% fit.full$LM$gls.residuals[fit.full$PermInfo$perm.schedule[[j]]]  )
  }
  else{y.rrpp.alt <- lapply(1:(iter+1), function(j) fit.full$LM$fitted +
                 fit.full$LM$residuals[fit.full$PermInfo$perm.schedule[[j]]]  ) }
    
  Models <-fit.full$Models
  Qr <- lapply(Models$reduced, function(x) x$qr)
  Qf <- lapply(Models$full, function(x) x$qr)
  k <- length(Qf)
  Ur <- lapply(Qr, qr.Q)
  Uf <- lapply(Qf, qr.Q)
  n <- fit.full$LM$n
  p <- fit.full$LM$p.prime
  dfs <- Map(function(ur, uf) ncol(uf) - ncol(ur), Ur, Uf)
  dfe <- n - ncol(Uf[[k]])
  
  fastFit <- RRPP:::fastFit
  fastLM <- RRPP:::fastLM
  
  # fast function
  getFs <- function(y, k, dfs, dfe, Ur, Uf) {
    MS <- unlist(Map(function(ur, uf, d) sum((fastFit(ur, y, n, p) - fastFit(uf, y, n, p))^2)/d, 
                     Ur, Uf, dfs))
    MSE <- sum(fastLM(Uf[[k]], as.matrix(y))$residuals^2)/dfe
    MS / MSE
  }

  res <- sapply(1:(iter+1), function(j) {
    getFs(y.rrpp.alt[[j]], k, dfs, dfe, Ur, Uf)
  })  

res  
}