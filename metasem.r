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

metasem <- function(rma_fit, sem_model, n_name, cor_var=NULL, n=NULL, 
                    n_fun=mean, obs_names=NULL, cluster_name=NULL, 
                    nearpd=FALSE, clean_data=TRUE, tran=NULL, ...){
  
  if(!inherits(rma_fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  
  JziLw._ <- if(is.null(rma_fit$formula.yi)) as.character(rma_fit$call$yi) else 
    .all.vars(rma_fit$formula.yi)[1]
  
  cor_var <- if(is.null(cor_var)) 
    as.formula(paste("~",(.all.vars(rma_fit$formula.mods)[1]), collapse = "~"))
  else cor_var
  
  cluster_name <- if(is.null(cluster_name)){ 
    mod_struct <- rma_clusters(rma_fit)
    names(mod_struct$level_dat)[which.min(mod_struct$level_dat)]
  } else cluster_name
  
  dat <- get_data_(rma_fit)
  
  dat <- if(clean_data) full_clean(dat) else dat
  
  n <- if(is.null(n)) sum(sapply(group_split(dplyr::filter(dat, !is.na(!!!JziLw._) & !is.na(!!!n_name)), 
                                             !!!cluster_name), function(i) 
                                               n_fun(unique(i[[n_name]])))) else n
  
  RAM <- lavaan2RAM(sem_model, obs_names)
  
  post <- post_rma(rma_fit, cor_var, tran=tran, type="response")
  
  rs <- coef(post)
  
  Cov <- if(is.null(obs_names)) vec2mat(rs) else vec2mat(rs, dimnames=obs_names)
  aCov <- vcov(post)
   
  Cov <- if(is.pd(Cov)) Cov else 
    if(nearpd) Matrix::nearPD(Cov, corr = TRUE) else 
      stop("r matrix not positive definite: Don't remove NAs or/and use 'nearpd=TRUE'.")
  
  aCov <- if(is.pd(aCov)) aCov else 
    if(nearpd) Matrix::nearPD(aCov) else 
      stop("Sampling covariance matrix not positive definite: Don't remove NAs or/and use 'nearpd=TRUE'.")
  
  wls(Cov=Cov, aCov=aCov, n=n, RAM=RAM, ...)  
  
}

#==============================================================================
source("https://raw.githubusercontent.com/rnorouzian/i/master/3m.r")
                                 
needzzsf <- c('lavaan', 'semPlot', 'metaSEM', 'Matrix')      

not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)

suppressWarnings(
  suppressMessages({ 
    
    invisible(lapply(needzzsf, base::require, character.only = TRUE))
    
  }))
