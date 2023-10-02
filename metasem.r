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
                    sep="[^[:alnum:]]+", data=NULL, tol=1e-06, ...){
  
  if(!inherits(rma_fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  
  dat <- if(is.null(data)) get_data_(rma_fit) else data
  
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
  
  post <- post_rma(fit=rma_fit, specs=cor_var, tran=tran, type="response", data=data)
  
  Rs <- coef(post)
  
  Cov <- coef2mat(Rs, sep = sep)
  aCov <- vcov(post)
  
  aCov <- vcov_match(Cov, aCov)
  
  RAM <- lavaan2RAM2(sem_model, rownames(Cov))

 ok_Cov <- !inherits(try(solve(Cov), silent=TRUE), "try-error")
ok_aCov <- !inherits(try(solve(aCov), silent=TRUE), "try-error")                                 
                                  
  Cov <- if(is.pd(Cov, tol=tol) & ok_Cov) Cov else 
    if(nearpd) as.matrix(Matrix::nearPD(Cov, corr=TRUE)$mat) else 
      stop("r matrix not positive definite: 
           1) If no moderator or subsetting's involved, Don't remove NAs or/and use 'nearpd=TRUE'.
           2) If a moderator or subsetting's involved, available data is insufficient for moderator analysis.")
  
  aCov <- if(is.pd(aCov, cor.analysis=FALSE, tol=tol) & ok_aCov) aCov else 
    if(nearpd) as.matrix(Matrix::nearPD(aCov)$mat) else 
      stop("Sampling covariance matrix not positive definite: 
           1) If no moderator or subsetting's involved, Don't remove NAs or/and use 'nearpd=TRUE'.
           2) If a moderator or subsetting's involved, available data is insufficient for moderator analysis.")
  
  wls(Cov=Cov, aCov=aCov, n=n, RAM=RAM, ...)  
  
}
                                  
#==============================================================================
                                  
metasem_3m <- function(rma_fit, sem_model, n_name, cor_var=NULL, n=NULL, 
                       n_fun=mean, cluster_name=NULL, 
                       nearpd=FALSE, tran=NULL, ngroups = 1L,
                       sep="[^[:alnum:]]+", moderator=NULL, data=NULL, tol=1e-06, ...){
  
  no_mod <- is.null(moderator)
  
  out <- if(no_mod) { 
    
    metasem(rma_fit=rma_fit, sem_model=sem_model, 
            n_name=n_name, cor_var=cor_var, n=n, 
            n_fun=n_fun, cluster_name=cluster_name, 
            nearpd=nearpd, tran=tran, sep=sep, data=data, tol=tol, ...)
    
  } else {
    
    
    dat_ <- if(is.null(data)) get_data_(rma_fit) else data
    
    mod <- .all.vars(moderator)[1]
    ok <- mod %in% names(dat_)
    if(!ok) stop("'moderator=' not found in the data.", call.=FALSE)
    
    pp <- lavaanify(sem_model, auto.var=TRUE, std.lv=TRUE, fixed.x=FALSE)
    
    regs <- subset(pp, free!=0 & op %in% "~")
    cors <- subset(pp, free!=0 & op %in% "~~")
    lat <- subset(pp, free!=0 & op %in% "=~")
    
    mod_lvls <- as.vector(na.omit(unique(dat_[[mod]])))
    
    mod_list <- lapply(mod_lvls, function(i) 
      suppressWarnings(update.rma(rma_fit, subset = get(mod) == i, data = dat_)))
    
    
    mod_list <- lapply(1:length(mod_list), function(i) 
    { mod_list[[i]]$data <- filter(dat_, !!sym(mod) == mod_lvls[i]); 
    return(mod_list[[i]]) })
    
    mod_list <- lapply(mod_list, function(x) {x$call$subset <- NULL; return(x)})
    
    
    ll = setNames(lapply(1:length(mod_lvls), function(i) transform(pp, label = 
                                                                     ifelse(free!=0 & op %in% "~", 
                                                                            paste0("b",1:nrow(regs),letters[i]),
                                                                            ifelse(free!=0 & op %in% "~~",
                                                                                   paste0("r", 1:nrow(cors),letters[i]),
                                                                                   ifelse(free!=0 & op %in% "=~",
                                                                                          paste0("v", 1:nrow(lat),letters[i]),label))))),
                  mod_lvls
    )
    
    
    mod_lvls <- str_remove(mod_lvls, "[^[:alnum:]]+")
    
    wls_list <- lapply(1:length(ll), function(i) metasem(rma_fit=mod_list[[i]], sem_model=ll[[i]], 
                                                         n_name=n_name, cor_var=cor_var, n=n, data=data,
                                                         n_fun=n_fun, cluster_name=cluster_name, tol=tol,
                                                         nearpd=nearpd, tran=tran, sep=sep, model=mod_lvls[i], run=FALSE, ...=...))  
    
    wls_model <- mxModel(model="combined", wls_list, mxFitFunctionMultigroup(mod_lvls))
    
    wls_fit <- mxRun(wls_model, intervals=TRUE)
    
    ss <- summary(wls_fit)
    
    res <- cbind(ss$parameters[-c(3:4,7:10)], ss$CI[c("lbound","ubound")])
    
    names(res)[4:6] <- c("SE", "Lower", "Upper") 
    
    res
    
  }
  
  return(out)
  
}
#==============================================================================
                                  
lavaan2RAM2 <- function (model, obs.variables = NULL, A.notation = "ON", S.notation = "WITH", 
                         M.notation = "mean", A.start = 0.1, S.start = 0.5, M.start = 0, 
                         auto.var = TRUE, std.lv = TRUE, ngroups = 1, ...) 
{
  my.model <- if(inherits(model,"data.frame")) model else lavaan::lavaanify(model, fixed.x = FALSE, auto.var = auto.var, 
                                std.lv = std.lv, ngroups = ngroups, ...)
  max.gp <- max(my.model$group)
  out <- list()
  for (gp in seq_len(max.gp)) {
    mod <- my.model[my.model$group == gp, ]
    if (any((mod$op == "=~" | mod$op == "~") & is.na(mod$ustart))) {
      mod[(mod$op == "=~" | mod$op == "~") & is.na(mod$ustart), 
      ]$ustart <- A.start
    }
    if (any(mod$op == "~1" & is.na(mod$ustart))) {
      mod[mod$op == "~1" & is.na(mod$ustart), ]$ustart <- M.start
    }
    if (any(mod$op == "~~" & is.na(mod$ustart) & (mod$lhs == 
                                                  mod$rhs))) {
      mod[mod$op == "~~" & is.na(mod$ustart) & (mod$lhs == 
                                                  mod$rhs), ]$ustart <- S.start
    }
    if (any(mod$op == "~~" & is.na(mod$ustart) & (mod$lhs != 
                                                  mod$rhs))) {
      mod[mod$op == "~~" & is.na(mod$ustart) & (mod$lhs != 
                                                  mod$rhs), ]$ustart <- 0
    }
    all.var <- unique(c(mod$lhs, mod$rhs))
    latent <- unique(mod[mod$op == "=~", ]$lhs)
    observed <- all.var[!(all.var %in% latent)]
    observed <- observed[observed != ""]
    if (!is.null(obs.variables)) {
      if (!identical(sort(observed), sort(obs.variables))) {
        stop("Names in \"obs.variables\" do not agree with those in model.\n")
      }
      else {
        observed <- obs.variables
      }
    }
    if (length(latent) > 0) {
      all.var <- c(observed, latent)
    }
    else {
      all.var <- observed
    }
    no.lat <- length(latent)
    no.obs <- length(observed)
    no.all <- no.lat + no.obs
    Amatrix <- matrix(0, ncol = no.all, nrow = no.all, dimnames = list(all.var, 
                                                                       all.var))
    Smatrix <- matrix(0, ncol = no.all, nrow = no.all, dimnames = list(all.var, 
                                                                       all.var))
    Mmatrix <- matrix(0, nrow = 1, ncol = no.all, dimnames = list(1, 
                                                                  all.var))
    for (i in seq_len(nrow(mod))) {
      if (mod[i, ]$label == "") {
        switch(mod[i, ]$op, `=~` = mod[i, ]$label <- paste0(mod[i, 
        ]$rhs, A.notation, mod[i, ]$lhs), `~` = mod[i, 
        ]$label <- paste0(mod[i, ]$lhs, A.notation, 
                          mod[i, ]$rhs), `~~` = mod[i, ]$label <- paste0(mod[i, 
                          ]$lhs, S.notation, mod[i, ]$rhs), `~1` = mod[i, 
                          ]$label <- paste0(mod[i, ]$lhs, M.notation))
      }
    }
    key <- with(mod, ifelse(free == 0, yes = ustart, no = paste(ustart, 
                                                                label, sep = "*")))
    for (i in seq_len(nrow(mod))) {
      my.line <- mod[i, ]
      switch(my.line$op, `=~` = Amatrix[my.line$rhs, my.line$lhs] <- key[i], 
             `~` = Amatrix[my.line$lhs, my.line$rhs] <- key[i], 
             `~~` = Smatrix[my.line$lhs, my.line$rhs] <- Smatrix[my.line$rhs, 
                                                                 my.line$lhs] <- key[i], `~1` = Mmatrix[1, my.line$lhs] <- key[i])
    }
    Fmatrix <- create.Fmatrix(c(rep(1, no.obs), rep(0, no.lat)), 
                              as.mxMatrix = FALSE)
    dimnames(Fmatrix) <- list(observed, all.var)
    out[[gp]] <- list(A = Amatrix, S = Smatrix, F = Fmatrix, 
                      M = Mmatrix)
  }
  names(out) <- seq_along(out)
  if (length(grep("^\\.", my.model$lhs)) > 0) {
    my.model <- my.model[-grep("^\\.", my.model$lhs), ]
  }
  if (any(my.model$group == 0)) {
    mxalgebra <- list()
    con_index <- 1
    y <- my.model[my.model$group == 0, , drop = FALSE]
    for (i in seq_len(nrow(y))) {
      switch(y[i, "op"], `:=` = {
        eval(parse(text = paste0(y[i, "lhs"], "<- mxAlgebra(", 
                                 y[i, "rhs"], ", name=\"", y[i, "lhs"], "\")")))
        eval(parse(text = paste0("mxalgebra <- c(mxalgebra, ", 
                                 y[i, "lhs"], "=", y[i, "lhs"], ")")))
      }, if (y[i, "op"] %in% c("==", ">", "<")) {
        eval(parse(text = paste0("constraint", con_index, 
                                 " <- mxConstraint(", y[i, "lhs"], y[i, "op"], 
                                 y[i, "rhs"], ", name=\"constraint", con_index, 
                                 "\")")))
        eval(parse(text = paste0("mxalgebra <- c(mxalgebra, constraint", 
                                 con_index, "=constraint", con_index, ")")))
        con_index <- con_index + 1
      })
    }
    out[[1]] <- list(A = out[[1]]$A, S = out[[1]]$S, F = out[[1]]$F, 
                     M = out[[1]]$M, mxalgebras = mxalgebra)
  }
  if (max.gp == 1) {
    out <- out[[1]]
  }
  out
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
