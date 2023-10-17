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

metasem_ <- function(rma_fit, sem_model, n_name, cor_var=NULL, n=NULL, 
                     n_fun=mean, cluster_name=NULL, model.name=NULL,
                     nearpd=FALSE, tran=NULL, run=TRUE,
                     sep="[^[:alnum:]]+", data=NULL, tol=1e-06,
                     std.lv=TRUE, auto.var=TRUE, RAM=NULL, ...){
  
  if(!inherits(rma_fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  
  dat <- if(is.null(data)) get_data_(rma_fit) else data
  
  JziLw._ <- if(is.null(rma_fit$formula.yi)) as.character(rma_fit$call$yi) else 
    .all.vars(rma_fit$formula.yi)[1]
  
  cor_var <- if(is.null(cor_var)) 
    as.formula(paste("~",(.all.vars(rma_fit$formula.mods)[1]), collapse = "~"))
  else cor_var
  
  ok <- as.character(cor_var)[2] %in% names(dat)
  if(!ok) stop("Select an accurate 'cor_var= ~VARIABLE_NAME' from the data.", call.=FALSE)
  
  cr <- is_crossed(rma_fit)
  
  ok <- !any(cr)
  
  if(!ok & is.null(cluster_name)) stop("Specify the 'cluster_name=' (usually the 'study' variable).", call.=FALSE)
  
  cluster_name <- if(is.null(cluster_name)){ 
    mod_struct <- rma_clusters(rma_fit)
    cl_nm <- names(mod_struct$level_dat)[which.min(mod_struct$level_dat)]
    message(paste0("NOTE: ",dQuote(cl_nm), " was selected as 'cluster_name=' (usually the 'study' variable). If incorrect, please change it.\n"))
    cl_nm
  } else cluster_name
  
  ok1 <- trimws(cluster_name) %in% names(dat) 
  ok2 <- trimws(cluster_name) %in% names(cr)
  if(!ok1) stop("'cluster_name=' not found in the data.", call.=FALSE)
  if(!ok2) stop("'cluster_name=' not found in the 'rma_fit'.", call.=FALSE)
  
  n_name <- trimws(n_name)
  ok <- n_name %in% names(dat)
  if(!ok) stop("'n_name=' not found in the data.", call.=FALSE)  
  
  n <- if(is.null(n)) sum(sapply(group_split(dplyr::filter(dat, !is.na(!!sym(JziLw._)) & !is.na(!!sym(n_name))), 
                                             !!sym(cluster_name)), function(i) 
                                               n_fun(unique(i[[n_name]])))) else n
  
  post <- post_rma(fit=rma_fit, specs=cor_var, tran=tran, type="response", data=data)
  
  Rs <- coef(post)
  
  Cov <- coef2mat(Rs, sep = sep)
  aCov <- vcov(post)
  
  aCov <- vcov_match(Cov, aCov)
  
  RAM <- if(is.null(RAM)) lavaan2RAM2(sem_model, colnames(Cov), std.lv=std.lv, auto.var=auto.var) else RAM
  
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
  
  if(is.null(model.name)) model.name <- "TSSEM2 Correlation"
  
  out <- wls(Cov=Cov, aCov=aCov, n=n, RAM=RAM, model.name=model.name, run=run, ...)  
  
  status <- if(run) out$mx.fit@output$status[[1]] else mxRun(out,suppressWarnings=TRUE)@output$status$code
  
  if(!status %in% 0:1) warning("Lack of convergence: Try 'rerun(output_of_this_function, extraTries=15)'.",call.=FALSE)
  
  if(run) out <- append(out, list(rma_fit=rma_fit, post_rma_fit=post, 
                                  n_name=n_name, cor_var=cor_var, RAM=RAM,
                                  sem_model=sem_model, cluster_name=cluster_name, 
                                  sep=sep, model.name=model.name, n_fun=n_fun, 
                                  nearpd=nearpd, tran=tran, run=run, 
                                  status=status, data=data))
  
  if(run) class(out) <- "wls"
  return(out) 
}
                                  
#==============================================================================
                                  
metasem_3m <- function(rma_fit, sem_model, n_name, cor_var=NULL, n=NULL, 
                       n_fun=mean, cluster_name=NULL, model.name=NULL,
                       nearpd=FALSE, tran=NULL, ngroups=1L, run=TRUE,
                       sep="[^[:alnum:]]+", moderator=NULL, data=NULL, 
                       std.lv=TRUE, auto.var=TRUE, RAM=NULL, ...){
  
out <- if(is.null(moderator)) { 
    
    metasem_(rma_fit=rma_fit, sem_model=sem_model, 
            n_name=n_name, cor_var=cor_var, n=n, model.name=model.name,
            n_fun=n_fun, cluster_name=cluster_name, run=run,
            nearpd=nearpd, tran=tran, sep=sep, data=data, 
            std.lv=std.lv, auto.var=auto.var, RAM=RAM, ...)
    
  } else {
    
    dat_ <- if(is.null(data)) get_data_(rma_fit) else data
    
    mod <- .all.vars(moderator)[1]
    ok <- mod %in% names(dat_)
    if(!ok) stop("'moderator=' not found in the data.", call.=FALSE)
    
    mod_lvls <- as.vector(na.omit(unique(dat_[[mod]])))
    
    mod_list <- setNames(lapply(mod_lvls, function(i) 
      try(suppressWarnings(update.rma(rma_fit, subset = get(mod) == i, data = dat_)),
          silent=TRUE)), mod_lvls)
    
    mod_list[sapply(mod_list, inherits, what="try-error")] <- NULL
    
    if(length(mod_list)==0) stop("Likely, insufficient data for moderator analysis.
     Try setting 'nearpd=TRUE'.", call.=FALSE)
    
    lost <- setdiff(mod_lvls, names(mod_list))
    if(length(lost)!=0) message(toString(dQuote(lost))," dropped due to lack of data at 1st stage.\n")                             
    
    mod_lvls <- names(mod_list) 
    
    mod_list <- lapply(1:length(mod_list), function(i) 
    { mod_list[[i]]$data <- filter(dat_, !!sym(mod) == mod_lvls[i]); 
    return(mod_list[[i]]) })
    
    mod_list <- lapply(mod_list, function(x) {x$call$subset <- NULL; return(x)})
    
    mod_lvls <- str_remove(mod_lvls, "[^[:alnum:]]+")
    
    wls_list <- setNames(lapply(1:length(mod_lvls), 
                                function(i) try(metasem_(rma_fit=mod_list[[i]], 
                                                        sem_model=sem_model, 
                                                        n_name=n_name, cor_var=cor_var, n=n, data=data,
                                                        n_fun=n_fun, cluster_name=cluster_name, 
                                                        nearpd=nearpd, tran=tran, sep=sep, run=run,
                                                        model.name=mod_lvls[i], std.lv=std.lv, 
                                                        auto.var=auto.var, RAM=RAM, ...=...), silent=TRUE)), mod_lvls)  
    
    
    wls_list[sapply(wls_list, inherits, what="try-error")] <- NULL
    
    if(length(wls_list)==0) stop("Likely, insufficient data for moderator analysis.
     Try setting 'nearpd=TRUE'.", call.=FALSE)
    
    lost <- setdiff(mod_lvls, names(wls_list))
    if(length(lost)!=0) message(toString(dQuote(lost))," dropped due to lack of data at 2nd stage.
                           Try: rerun(output_of_this_function, extraTries=15)") 
    
    if(run) class(wls_list) <- "wls.cluster"
    
    return(wls_list)
    
  }
  return(out)
}
#==============================================================================

plot_sem3m <- function(x, main=NA, reset=TRUE, 
                       index=NULL, line=NA, 
                       cex.main=1, mfrow=NULL, ...)
{
  
  if (!requireNamespace("semPlot", quietly = TRUE)) 
    stop("\"semPlot\" package is required for this function.", call. = FALSE)
  
  if (!inherits(x, c("wls","wls.cluster","list","character"))) 
    stop("\"x\" must be an object of class \"wls\", \"wls.cluster\", \"list\", or \"character\".", call. = FALSE)
  
  if(inherits(x, c("wls","character"))) x <- list(x)
  
  if(inherits(x, c("wls.cluster","list")) & !is.null(index)) {
    
    LL <- length(x)
    
    index <- index[index <= LL & index >= 1]
    if(length(index)==0) index <- NULL
    
    x <- if(is.null(index)) x else x[index]
    
  }
  
  if(reset){
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
  }
  
  h <- length(x)
  if(h>1) { par(mfrow = if(is.null(mfrow)) n2mfrow(h) else mfrow, mgp = c(1.5,.5,0), mar = c(8,.5,.5,.5)+.1, 
                tck = -.02, xpd = FALSE) }
  
  ff <- function(x, main, line, cex.main, ...) { 
    plot(x=x, ...)
    graphics::title(main=main, line=line, 
                    cex.main=cex.main) 
    }
  
  x_nm <- names(x)
  
  cls <- sapply(x, class)
  
  main <- if(anyNA(main) & length(x_nm)!=0) x_nm else if (
    anyNA(main) & length(x_nm)==0 & !"character" %in% cls) unname(sapply(x, '[[', "model.name")) else main
  
  invisible(Map(ff, x=x, main=main, line=line, cex.main=cex.main, ...))
}
                
#===============================================================================
                                
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
