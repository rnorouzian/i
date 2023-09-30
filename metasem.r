#==============================================================================

make_nested <- function (x) 
{
  if (length(x) == 1) 
    return(x)
  y <- x
  for (i in 2:length(x)) {
    names(y)[i] <- paste(names(x)[1:i], collapse = "/")
    y[i] <- do.call(paste, c(x[1:i], sep = "/"))
  }
  y
}
#==============================================================================

rma_clusters <- function (obj) 
{
  level_dat <- vector(mode = "integer")
  cluster_dat <- data.frame(row.names = 1:obj$k)
  if (obj$withG) {
    level_dat[[obj$g.names[[2]]]] <- obj$g.nlevels[[2]]
    cluster_dat[[obj$g.names[[2]]]] <- obj$mf.g$outer
  }
  if (obj$withH) {
    level_dat[[obj$h.names[[2]]]] <- obj$h.nlevels[[2]]
    cluster_dat[[obj$h.names[[2]]]] <- obj$mf.h$outer
  }
  if (obj$withS) {
    s_levels <- obj$s.nlevels
    names(s_levels) <- obj$s.names
    level_dat <- c(level_dat, s_levels)
    mf_r <- lapply(obj$mf.r, make_nested)
    mf_all <- do.call(cbind, mf_r)
    mf_s <- mf_all[obj$s.names]
    cluster_dat <- cbind(cluster_dat, mf_s)
    cluster_dat <- droplevels(cluster_dat)
  }
  list(level_dat = level_dat, cluster_dat = cluster_dat)
}

#==============================================================================

coef2mat <- function(coefs, sep="[^[:alnum:]]+"){
  row.col <- do.call(rbind, strsplit(names(coefs), sep))
  nm <- sort(unique(c(row.col)))
  corr <- matrix(NA, length(nm), length(nm), dimnames=list(nm, nm))
  diag(corr) <- 1
  corr[rbind(row.col, row.col[,2:1])] <- coefs
  return(corr)
}

#==============================================================================

vcov_match <- function(r_mat, v_mat){
  
  sep <- unique(gsub("[a-zA-Z0-9]", "", rownames(v_mat)))
  
  col.row <- outer(colnames(r_mat), rownames(r_mat), FUN=paste, sep=sep)
  
  v.names <- strsplit(rownames(v_mat), "[^[:alnum:]]+") |>
    lapply(sort, decreasing=TRUE) |>
    sapply(paste, collapse=sep)              
  
  ord <- match(col.row[lower.tri(r_mat)], v.names)
  if(all(is.na(ord))) ord <- match(col.row[lower.tri(r_mat)], lapply(strsplit(v.names, "\\W"), function(x) paste(rev(x), collapse=sep)))
  return(v_mat[ord,ord])
}

#==============================================================================

metasem <- function(rma_fit, sem_model, n_name, cor_var=NULL, n=NULL, 
                    n_fun=mean, cluster_name=NULL, 
                    nearpd=FALSE, tran=NULL, 
                    sep="[^[:alnum:]]+", ...){
  
  if(!inherits(rma_fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  
  dat <- get_data_(rma_fit)
  
  JziLw._ <- if(is.null(rma_fit$formula.yi)) as.character(rma_fit$call$yi) else 
    .all.vars(rma_fit$formula.yi)[1]
  
  cor_var <- if(is.null(cor_var)) 
    as.formula(paste("~",(.all.vars(rma_fit$formula.mods)[1]), collapse = "~"))
  else cor_var
  
  ok <- as.character(cor_var)[2] %in% names(dat)
  if(!ok) stop("Select an accurate 'cor_var= ~VARIABLE_NAME'.", call.=FALSE)
  
  cluster_name <- if(is.null(cluster_name)){ 
    mod_struct <- rma_clusters(rma_fit)
    names(mod_struct$level_dat)[which.min(mod_struct$level_dat)]
  } else cluster_name
  
  ok <- trimws(cluster_name) %in% names(dat)
  if(!ok) stop("'cluster_name=' is incorrect.", call.=FALSE)
  
  n_name <- trimws(n_name)
  ok <- n_name %in% names(dat)
  if(!ok) stop("'n_name=' is incorrect.", call.=FALSE)  
  
  n <- if(is.null(n)) sum(sapply(group_split(dplyr::filter(dat, !is.na(!!!JziLw._) & !is.na(!!!n_name)), 
                                             !!!cluster_name), function(i) 
                                               n_fun(unique(i[[n_name]])))) else n
  
  post <- post_rma(rma_fit, cor_var, tran=tran, type="response")
  
  Rs <- coef(post)
  
  Cov <- coef2mat(Rs, sep = sep)
  aCov <- vcov(post)
  
  aCov <- vcov_match(Cov, aCov)
  
  RAM <- lavaan2RAM(sem_model, rownames(Cov))
  
  Cov <- if(is.pd(Cov)) Cov else 
    if(nearpd) Matrix::nearPD(Cov, corr=TRUE) else 
      stop("r matrix not positive definite: Don't remove NAs or/and use 'nearpd=TRUE'.")
  
  aCov <- if(is.pd(aCov)) aCov else 
    if(nearpd) Matrix::nearPD(aCov) else 
      stop("Sampling covariance matrix not positive definite: Don't remove NAs or/and use 'nearpd=TRUE'.")
  
  wls(Cov=Cov, aCov=aCov, n=n, RAM=RAM, ...)  
  
}
#==============================================================================

#==============================================================================
                                 
source("https://raw.githubusercontent.com/rnorouzian/i/master/3m.r")
                                 
needzzsf <- c('lavaan', 'semPlot', 'metaSEM', 'Matrix')      

not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)

suppressWarnings(
  suppressMessages({ 
    
    invisible(lapply(needzzsf, base::require, character.only = TRUE))
    
  }))
