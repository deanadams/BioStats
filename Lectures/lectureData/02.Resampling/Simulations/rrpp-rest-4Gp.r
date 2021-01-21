####  Specific to 4 groups (and equal membership and specific order...)
rrpp.rest <- function(f1, iter = 999, seed = "random", int.first = FALSE,
                    RRPP = TRUE, SS.type = c("I", "II", "III"),
                    data = NULL, Cov = NULL,
                    print.progress = FALSE, Parallel = FALSE, ...) {
  
  L <- c(as.list(environment()), list(...))
#  L <- as.list(environment())
  names(L)[which(names(L) == "f1")] <- "formula"
  
  if(int.first) ko = TRUE else ko = FALSE
  SS.type <- match.arg(SS.type)
  
  dots <- list(...)
  if(length(dots) > 0) {
    w <- dots$weights
    o <- dots$offset
  } else w <- o <- NULL
  
########################  Run bits for each main effect
  fit.full <- lm.rrpp(L$formula, data = data, SS.type = SS.type, Cov = Cov,
                print.progress = print.progress, iter= 1)
  
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
  #specific to THIS example only (4 gp; perfect)
  n.gp <- n/4
  it1 <- RRPP:::perm.index(n.gp, iter,seed = seed)
  it2 <- RRPP:::perm.index(n.gp, iter,seed = seed)
  it3 <- RRPP:::perm.index(n.gp, iter,seed = seed)
  it4 <- RRPP:::perm.index(n.gp, iter,seed = seed)
  it5 <- RRPP:::perm.index(n.gp, iter,seed = seed)
  it6 <- RRPP:::perm.index(n.gp, iter,seed = seed)
  it7 <- RRPP:::perm.index(n.gp, iter,seed = seed)
  it8 <- RRPP:::perm.index(n.gp, iter,seed = seed)
  
  it.x2 <- lapply(1:(iter+1), function(j) c(it1[[j]],n.gp+it2[[j]],
                  (2*n.gp)+it3[[j]],(3*n.gp)+it4[[j]]) )
  it.x1 <- lapply(1:(iter+1), function(j) c(rbind(it5[[j]],n.gp+it6[[j]],
                  (2*n.gp)+it7[[j]],(3*n.gp)+it8[[j]])))
  it.x1[[1]] <- it.x2[[1]]

  if (!is.null(Cov)){
    Pcov <- fit.full$LM$Pcov
    y.rrpp.alt1 <- lapply(1:(iter+1), function(j) Pcov %*% L$data$y[it.x2[[j]]] )
    y.rrpp.alt2 <- lapply(1:(iter+1), function(j) Pcov %*% L$data$y[it.x1[[j]]] )
  }
  
  #Shuffle by x2 to test x1, and x1 for s2
  else{ 
    y.rrpp.alt1 <- lapply(1:(iter+1), function(j) L$data$y[it.x2[[j]]] )
    y.rrpp.alt2 <- lapply(1:(iter+1), function(j) L$data$y[it.x1[[j]]] )
  }

  res1 <- sapply(1:(iter+1), function(j) {
    getFs(y.rrpp.alt1[[j]], k, dfs, dfe, Ur, Uf)[1]
  })  
  res2 <- sapply(1:(iter+1), function(j) {
    getFs(y.rrpp.alt2[[j]], k, dfs, dfe, Ur, Uf)[2]
  })  
  res <-  rbind(res1,res2)
  res
}