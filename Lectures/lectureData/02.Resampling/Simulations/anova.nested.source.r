# aov.single.model Spitting out Fs
aov.single.model.F <- function(object, effect.type = "F",
                             error = NULL) {
  x <- object$ANOVA
  df <- x$df
  k <- length(df)-2
  SS <- x$SS
  MS <- x$MS 
  RSS <-x$RSS
  TSS <- x$TSS
  perms <- object$PermInfo$perms
  pm <- object$PermInfo$perm.method
  trms <- object$LM$term.labels
  
  if(!is.null(error)) {
    if(!inherits(error, "character")) stop("The error description is illogical.  It should be a string of character values matching ANOVA terms.",
                                           call. = FALSE)
    kk <- length(error)
    if(kk != k) stop("The error description should match in length the number of ANOVA terms (not including Residuals)",
                     call. = FALSE)
    MSEmatch <- match(error, c(trms, "Residuals"))
    if(any(is.na(MSEmatch))) stop("At least one of the error terms is not an ANOVA term",
                                  call. = FALSE)
  } else MSEmatch <- NULL
  if(k >= 1) {
    Fs <- x$Fs
    
    if(!is.null(MSEmatch)){
      Fmatch <- which(MSEmatch <= k)
      Fs[Fmatch,] <- MS[Fmatch,]/MS[MSEmatch[Fmatch],]
      F.effect.adj <- apply(Fs, 1, RRPP:::effect.size)
    }
    
    out <- Fs
  } 
  out
}

# Permuted nested only, not nested + overall error
aov.single.model.F2 <- function(object, effect.type = "F",
                               error = NULL) {
  x <- object$ANOVA
  df <- x$df
  k <- length(df)-2
  SS <- x$SS
  MS <- x$MS 
  RSS <-x$RSS
  TSS <- x$TSS
  perms <- object$PermInfo$perms
  pm <- object$PermInfo$perm.method
  trms <- object$LM$term.labels
  
  if(!is.null(error)) {
    if(!inherits(error, "character")) stop("The error description is illogical.  It should be a string of character values matching ANOVA terms.",
                                           call. = FALSE)
    kk <- length(error)
    if(kk != k) stop("The error description should match in length the number of ANOVA terms (not including Residuals)",
                     call. = FALSE)
    MSEmatch <- match(error, c(trms, "Residuals"))
    if(any(is.na(MSEmatch))) stop("At least one of the error terms is not an ANOVA term",
                                  call. = FALSE)
  } else MSEmatch <- NULL
  if(k >= 1) {
    Fs <- x$Fs
    
    if(!is.null(MSEmatch)){
      Fmatch <- which(MSEmatch <= k)
      Fs[Fmatch,] <- MS[Fmatch,]/MS[MSEmatch[Fmatch],]
      F.effect.adj <- apply(Fs, 1, RRPP:::effect.size)
    }
    
    out <- Fs
  } 
  out
}