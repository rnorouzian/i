# H==============================================================================

coef2mat <- function(coefs, sep="[^[:alnum:]]+"){
  row.col <- do.call(rbind, strsplit(names(coefs), sep))
  nm <- sort(unique(c(row.col)))
  corr <- matrix(NA, length(nm), length(nm), dimnames=list(nm, nm))
  corr[rbind(row.col, row.col[,2:1])] <- coefs
  diag(corr) <- 1
  return(corr)
}

# H==============================================================================

vcov_match <- function(r_mat, v_mat){
  
  sep <- unique(gsub("[a-zA-Z0-9]", "", rownames(v_mat)))
  
  col.row <- outer(colnames(r_mat), rownames(r_mat), FUN=paste, sep=sep)
  
  v.names <- strsplit(rownames(v_mat), "[^[:alnum:]]+") %>%
    lapply(sort, decreasing=TRUE) %>%
    sapply(paste, collapse=sep)              
  
  ord <- match(col.row[lower.tri(r_mat)], v.names)
  if(all(is.na(ord))) ord <- match(col.row[lower.tri(r_mat)], lapply(strsplit(v.names, "\\W"), function(x) paste(rev(x), collapse=sep)))
  return(v_mat[ord,ord])
}


# H==============================================================================
# Optional covariance-weighted admissibility constraint for non-PD pooled R
#
# This implements the deterministic multi-start hyperspherical-Cholesky
# projection investigated in the dependency-aware MASEM simulations.
# It is deliberately opt-in. If the ordinary pooled R is already PD, the
# identity branch is used and nothing is changed.
#
# When the constraint is active, uncertainty is estimated with the projected-z
# bootstrap. Simulation work indicates that this is useful but can be
# anti-conservative in active/boundary cases, so Stage-2 intervals should be
# interpreted cautiously.

.symm_ <- function(M) (M + t(M))/2

.min_eig_ <- function(M){
  min(eigen(.symm_(M), symmetric=TRUE, only.values=TRUE)$values)
}

.clip_cor_ <- function(x, eps=1e-10){
  pmin(pmax(x, -1 + eps), 1 - eps)
}

.stable_ginv_ <- function(V, tol=1e-10){
  V <- .symm_(as.matrix(V))
  ee <- eigen(V, symmetric=TRUE)
  sc <- max(1, max(abs(ee$values)))
  if(min(ee$values) < -1e-8*sc)
    stop("Sampling covariance matrix is materially non-PSD; constrained admissibility cannot be applied.", call.=FALSE)
  keep <- ee$values > tol*max(ee$values)
  if(!any(keep))
    stop("Sampling covariance matrix has no estimable positive subspace.", call.=FALSE)
  Q <- ee$vectors[,keep,drop=FALSE]
  .symm_(Q %*% diag(1/ee$values[keep], sum(keep)) %*% t(Q))
}

.rmvn_psd_ <- function(n, mu, Sigma){
  Sigma <- .symm_(as.matrix(Sigma))
  ee <- eigen(Sigma, symmetric=TRUE)
  sc <- max(1, max(abs(ee$values)))
  if(min(ee$values) < -1e-8*sc)
    stop("Bootstrap covariance matrix is materially non-PSD.", call.=FALSE)
  A <- ee$vectors %*% diag(sqrt(pmax(ee$values,0)), length(ee$values))
  Z <- matrix(rnorm(n*length(mu)), nrow=n)
  sweep(Z %*% t(A), 2, mu, "+")
}

.eta_to_cor_ <- function(eta, p, angle_eps=1e-6, pd_floor=1e-8){
  phi <- angle_eps + (pi - 2*angle_eps)*plogis(eta)
  L <- matrix(0,p,p)
  L[1,1] <- 1
  idx <- 1L
  for(i in 2:p){
    ang <- phi[idx:(idx+i-2L)]
    idx <- idx+i-1L
    prod_sin <- 1
    for(j in 1:(i-1L)){
      L[i,j] <- prod_sin*cos(ang[j])
      prod_sin <- prod_sin*sin(ang[j])
    }
    L[i,i] <- prod_sin
  }
  C <- .symm_(L %*% t(L))
  diag(C) <- 1
  R <- (1-pd_floor)*C + pd_floor*diag(p)
  diag(R) <- 1
  .symm_(R)
}

.cor_to_eta_ <- function(R, angle_eps=1e-6){
  L <- t(chol(.symm_(R)))
  p <- nrow(R)
  phis <- numeric()
  for(i in 2:p){
    prod_sin <- 1
    for(j in 1:(i-1L)){
      den <- max(prod_sin,1e-12)
      cc <- min(max(L[i,j]/den,-1),1)
      ph <- min(max(acos(cc),angle_eps*2),pi-angle_eps*2)
      phis <- c(phis,ph)
      prod_sin <- prod_sin*sin(ph)
    }
  }
  u <- (phis-angle_eps)/(pi-2*angle_eps)
  qlogis(pmin(pmax(u,1e-10),1-1e-10))
}

.rvec_from_eta_ <- function(eta, p){
  R <- .eta_to_cor_(eta,p)
  R[lower.tri(R)]
}

.zvec_from_eta_ <- function(eta, p){
  atanh(.clip_cor_(.rvec_from_eta_(eta,p)))
}

.make_pd_start_ <- function(R_raw, target_min=.03){
  p <- nrow(R_raw)
  I <- diag(p)
  dimnames(I) <- dimnames(R_raw)
  if(.min_eig_(R_raw) > target_min) return(R_raw)
  emin <- .min_eig_(R_raw)
  lambda <- (target_min-emin)/(1-emin)
  lambda <- min(max(lambda+1e-6,0),.999999)
  R0 <- (1-lambda)*R_raw + lambda*I
  diag(R0) <- 1
  if(.min_eig_(R0) <= 0)
    R0 <- as.matrix(Matrix::nearPD(R_raw,corr=TRUE,keepDiag=TRUE)$mat)
  dimnames(R0) <- dimnames(R_raw)
  R0
}

.make_multistarts_ <- function(default_start, eta_start=NULL,
                               n_perturb=4L, scale=1.5){
  q <- length(default_start)
  starts <- list()
  if(!is.null(eta_start) && length(eta_start)==q && all(is.finite(eta_start)))
    starts[[length(starts)+1L]] <- eta_start
  starts[[length(starts)+1L]] <- default_start
  starts[[length(starts)+1L]] <- rep(0,q)
  if(n_perturb>0L){
    for(s in seq_len(n_perturb)){
      pat <- ((seq_len(q)+s) %% 3L)-1L
      if(all(pat==0)) pat[seq(1L,q,by=2L)] <- 1
      den <- sqrt(sum(pat^2)); if(den==0) den <- 1
      starts[[length(starts)+1L]] <- default_start + scale*(as.numeric(pat)/den)
    }
  }
  keys <- vapply(starts,function(x) paste(round(x,10),collapse="|"),character(1))
  starts[!duplicated(keys)]
}

.project_cor_ <- function(z_target, Vz, R_template, eta_start=NULL,
                          pd_tol=1e-8, opt_maxit=1500L, opt_reltol=1e-8,
                          n_perturb=4L, perturb_scale=1.5){
  z_target <- as.numeric(z_target)
  Vz <- .symm_(as.matrix(Vz))
  p <- nrow(R_template)
  r_raw <- tanh(z_target)
  R_raw <- R_template
  R_raw[lower.tri(R_raw)] <- r_raw
  R_raw[upper.tri(R_raw)] <- t(R_raw)[upper.tri(R_raw)]
  diag(R_raw) <- 1
  raw_min <- .min_eig_(R_raw)

  if(raw_min > pd_tol){
    return(list(z=z_target,r=r_raw,R=R_raw,raw_R=R_raw,
                raw_min_eigen=raw_min,min_eigen=raw_min,active=FALSE,
                convergence=0L,objective=0,eta=NULL,optimizer="identity",
                n_starts=0L,best_start=NA_integer_))
  }

  W <- .stable_ginv_(Vz)
  default_start <- .cor_to_eta_(.make_pd_start_(R_raw))

  objective <- function(eta){
    dz <- z_target - .zvec_from_eta_(eta,p)
    as.numeric(crossprod(dz,W %*% dz))
  }

  starts <- .make_multistarts_(default_start,eta_start,n_perturb,perturb_scale)
  candidates <- vector("list",length(starts))

  for(ss in seq_along(starts)){
    oo <- try(optim(starts[[ss]],objective,method="BFGS",
                    control=list(maxit=opt_maxit,reltol=opt_reltol)),silent=TRUE)
    if(inherits(oo,"try-error") || !is.finite(oo$value)) next
    eta_hat <- oo$par; conv <- oo$convergence; obj <- oo$value
    optimizer <- "BFGS"

    if(conv!=0L){
      nn <- try(nlminb(eta_hat,objective,
                       control=list(iter.max=opt_maxit,eval.max=opt_maxit*2L,
                                    rel.tol=opt_reltol,x.tol=opt_reltol)),silent=TRUE)
      if(!inherits(nn,"try-error") && is.finite(nn$objective)){
        eta_hat <- nn$par; conv <- nn$convergence; obj <- nn$objective
        optimizer <- "nlminb_fallback"
      }
    }
    if(conv!=0L || !is.finite(obj) || any(!is.finite(eta_hat))) next
    candidates[[ss]] <- list(eta=eta_hat,objective=obj,convergence=as.integer(conv),
                             optimizer=optimizer,start_id=ss)
  }

  good <- which(vapply(candidates,function(x)!is.null(x),logical(1)))
  if(!length(good)) stop("All constrained multi-start projection attempts failed.",call.=FALSE)
  objs <- vapply(candidates[good],function(x)x$objective,numeric(1))
  best <- candidates[[good[which.min(objs)]]]

  polish <- try(nlminb(best$eta,objective,
                       control=list(iter.max=opt_maxit,eval.max=opt_maxit*2L,
                                    rel.tol=opt_reltol,x.tol=opt_reltol)),silent=TRUE)
  if(!inherits(polish,"try-error") && polish$convergence==0L &&
     is.finite(polish$objective) && polish$objective <= best$objective+1e-10){
    best$eta <- polish$par; best$objective <- polish$objective
    best$convergence <- as.integer(polish$convergence)
    best$optimizer <- paste0(best$optimizer,"+nlminb_polish")
  }

  z_star <- .zvec_from_eta_(best$eta,p)
  r_star <- tanh(z_star)
  R_star <- R_template
  R_star[lower.tri(R_star)] <- r_star
  R_star[upper.tri(R_star)] <- t(R_star)[upper.tri(R_star)]
  diag(R_star) <- 1

  list(z=z_star,r=r_star,R=R_star,raw_R=R_raw,
       raw_min_eigen=raw_min,min_eigen=.min_eig_(R_star),active=TRUE,
       convergence=best$convergence,objective=best$objective,eta=best$eta,
       optimizer=best$optimizer,n_starts=length(starts),best_start=best$start_id)
}

.projected_z_boot_ <- function(fit,Vz,R_template,B=499L,seed=20261005L,
                               pd_tol=1e-8,opt_maxit=1500L,opt_reltol=1e-8,
                               n_perturb=4L,perturb_scale=1.5){
  if(!fit$active)
    return(list(Vz=Vz,z_draws=NULL,failures=0L,B_used=0L,method="identity"))

  set.seed(seed)
  raw <- .rmvn_psd_(B,fit$z,Vz)
  zz <- t(vapply(seq_len(B),function(b){
    ff <- try(.project_cor_(raw[b,],Vz,R_template,eta_start=fit$eta,
                            pd_tol=pd_tol,opt_maxit=opt_maxit,opt_reltol=opt_reltol,
                            n_perturb=n_perturb,perturb_scale=perturb_scale),silent=TRUE)
    if(inherits(ff,"try-error") || ff$convergence!=0L || any(!is.finite(ff$z)))
      rep(NA_real_,length(fit$z)) else ff$z
  },numeric(length(fit$z))))

  good <- complete.cases(zz)
  if(sum(good) < max(30L,floor(.80*B)))
    stop("Too many projected-z bootstrap failures.",call.=FALSE)
  zz <- zz[good,,drop=FALSE]
  list(Vz=.symm_(cov(zz)),z_draws=zz,failures=sum(!good),
       B_used=nrow(zz),method="projected_z")
}

.constrain_pooled_cor_ <- function(Cov,aCov,B=499L,seed=20261005L,
                                   pd_tol=1e-8,opt_maxit=1500L,opt_reltol=1e-8,
                                   n_perturb=4L,perturb_scale=1.5){
  r0 <- Cov[lower.tri(Cov)]
  r0 <- .clip_cor_(r0)

  # aCov is on the correlation scale in metasem_(). Transform it to Fisher-z.
  J <- diag(1/(1-r0^2),length(r0))
  Vz0 <- .symm_(J %*% as.matrix(aCov) %*% J)

  fit <- .project_cor_(atanh(r0),Vz0,Cov,pd_tol=pd_tol,
                       opt_maxit=opt_maxit,opt_reltol=opt_reltol,
                       n_perturb=n_perturb,perturb_scale=perturb_scale)

  boot <- .projected_z_boot_(fit,Vz0,Cov,B=B,seed=seed,pd_tol=pd_tol,
                             opt_maxit=opt_maxit,opt_reltol=opt_reltol,
                             n_perturb=n_perturb,perturb_scale=perturb_scale)

  # Return uncertainty to the correlation metric expected by metaSEM::wls().
  D <- diag(1-fit$r^2,length(fit$r))
  Vr_star <- .symm_(D %*% boot$Vz %*% D)

  list(Cov=fit$R,aCov=Vr_star,fit=fit,bootstrap=boot,
       raw_Cov=Cov,raw_aCov=aCov,Vz_raw=Vz0)
}


# H==============================================================================

metasem_ <- function(rma_fit, sem_model, n_name, cor_var=NULL, n=NULL, 
                     n_fun=mean, cluster_name=NULL, model.name=NULL,
                     nearpd=FALSE, admissibility=NULL,
                     constrained_B=499L, constrained_seed=20261005L,
                     tran=NULL, run=TRUE, weights="proportional",
                     sep="[^[:alnum:]]+", data=NULL, tol=1e-06,
                     std.lv=TRUE, auto.var=TRUE, RAM=NULL, ...){
  
  if(!inherits(rma_fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  
  dat <- if(is.null(data)) get_data_(rma_fit) else data
  
  nm_dat <- names(dat) 
  
  JziLw._ <- if(is.null(rma_fit$formula.yi)) as.character(rma_fit$call$yi) else 
    .all.vars(rma_fit$formula.yi)[1]
  
  cor_var <- if(is.null(cor_var)) 
    as.formula(paste("~",(.all.vars(rma_fit$formula.mods)[1]), collapse = "~"))
  else cor_var
  
  if(!is_bare_formula(cor_var, lhs=FALSE) || length(.all.vars(cor_var))>1) {
    stop("Select an accurate formula 'cor_var= ~VARIABLE_NAME' consisting of only 1 variable.", call. = FALSE)
  }
  
  ok <- .all.vars(cor_var) %in% nm_dat
  if(!ok) stop("'cor_var= ~VARIABLE_NAME' not found in the data.", call.=FALSE)
  
  cr <- is_crossed(rma_fit)
  
  ok <- !any(cr)
  
  if(!ok & is.null(cluster_name)) stop("Specify the 'cluster_name=' ('study' variable).", call.=FALSE)
  
  cluster_name <- if(is.null(cluster_name)){ 
    mod_struct <- rma_clusters(rma_fit)
    cl_nm <- names(mod_struct$level_dat)[which.min(mod_struct$level_dat)]
    message(paste0("NOTE: ",dQuote(cl_nm), " was selected as 'cluster_name=' (usually the 'study' variable). If incorrect, please change it.\n"))
    cl_nm
  } else cluster_name
  
  cluster_name <- trimws(cluster_name)
  
  ok1 <- cluster_name %in% nm_dat 
  ok2 <- cluster_name %in% names(cr)
  if(!ok1) stop("'cluster_name=' not found in the data.", call.=FALSE)
  if(!ok2) stop("'cluster_name=' not found in the 'rma_fit'.", call.=FALSE)
  
  n_name <- trimws(n_name)
  ok <- n_name %in% nm_dat
  if(!ok) stop("'n_name=' not found in the data.", call.=FALSE)  
  
  n <- if(is.null(n)) sum(sapply(group_split(dplyr::filter(dat, !is.na(!!sym(JziLw._)) & !is.na(!!sym(n_name))), 
                                             !!sym(cluster_name)), function(i) 
                                               n_fun(unique(i[[n_name]])))) else n
  
  post <- post_rma(fit=rma_fit, specs=cor_var, tran=tran, type="response", 
                   weights=weights, data=data)
  
  Rs <- coef(post)
  
  Cov <- coef2mat(Rs, sep = sep)
  aCov <- vcov(post)
  
  aCov <- vcov_match(Cov, aCov)
  
  RAM <- if(is.null(RAM)) lavaan2RAM2(sem_model, colnames(Cov), std.lv=std.lv, auto.var=auto.var) else RAM
  
  # Backward compatibility: nearpd=TRUE retains the old behavior unless
  # admissibility is supplied explicitly.
  if(is.null(admissibility)) admissibility <- if(nearpd) "nearPD" else "stop"
  admissibility <- match.arg(admissibility, c("stop","constrained","nearPD"))

  ok_Cov <- !inherits(try(solve(Cov), silent=TRUE), "try-error")
  ok_aCov <- !inherits(try(solve(aCov), silent=TRUE), "try-error")

  raw_Cov <- Cov
  raw_aCov <- aCov
  admissibility_fit <- NULL

  if(!(is.pd(Cov,tol=tol) && ok_Cov)){
    if(admissibility=="constrained"){
      admissibility_fit <- .constrain_pooled_cor_(
        Cov=Cov, aCov=aCov, B=constrained_B, seed=constrained_seed,
        pd_tol=min(tol,1e-8)
      )
      Cov <- admissibility_fit$Cov
      aCov <- admissibility_fit$aCov

      warning(
        paste0(
          "The pooled correlation matrix was non-positive-definite and the ",
          "covariance-weighted constrained admissibility estimator was applied. ",
          "Simulation evidence supports stable admissible point estimation, but ",
          "projected-z uncertainty can be anti-conservative when the constraint ",
          "is active. Interpret Stage-2 confidence intervals cautiously."
        ),
        call.=FALSE
      )
    } else if(admissibility=="nearPD"){
      Cov <- as.matrix(Matrix::nearPD(Cov,corr=TRUE)$mat)
    } else {
      stop(
        paste0(
          "Pooled correlation matrix is not positive definite. ",
          "Use admissibility='constrained' for the covariance-weighted ",
          "constrained estimator, or admissibility='nearPD' for the geometric ",
          "nearPD sensitivity repair."
        ),
        call.=FALSE
      )
    }
  }

  # The constrained branch jointly supplies its transformed aCov. Otherwise,
  # retain the original aCov check/repair logic.
  if(is.null(admissibility_fit)){
    ok_aCov <- !inherits(try(solve(aCov),silent=TRUE),"try-error")
    if(!(is.pd(aCov,cor.analysis=FALSE,tol=tol) && ok_aCov)){
      if(admissibility=="nearPD"){
        aCov <- as.matrix(Matrix::nearPD(aCov)$mat)
      } else {
        stop(
          paste0(
            "Sampling covariance matrix is not positive definite. ",
            "The constrained admissibility method repairs a non-PD pooled ",
            "correlation matrix; it is not a generic repair for an indefinite ",
            "sampling covariance matrix. Use admissibility='nearPD' only as a ",
            "sensitivity repair if desired."
          ),
          call.=FALSE
        )
      }
    }
  } else {
    ok_aCov_star <- !inherits(try(solve(aCov),silent=TRUE),"try-error")
    if(!(is.pd(aCov,cor.analysis=FALSE,tol=tol) && ok_aCov_star)){
      stop(
        paste0(
          "The constrained estimator produced a non-invertible Stage-2 aCov. ",
          "Do not silently nearPD-repair this covariance because it would break ",
          "the estimator/uncertainty pairing."
        ),
        call.=FALSE
      )
    }
  }

  if(is.null(model.name)) model.name <- "TSSEM2 Correlation"
  
  out <- wls(Cov=Cov, aCov=aCov, n=n, RAM=RAM, model.name=model.name, run=run, ...)  
  
  status <- if(run) out$mx.fit@output$status[[1]] else mxRun(out,suppressWarnings=TRUE)@output$status$code
  
  if(!status %in% 0:1) warning("Lack of convergence: Try 'rerun(output_of_this_function, extraTries=15)'.",call.=FALSE)
  
  if(run) out <- append(out, list(rma_fit=rma_fit, post_rma_fit=post, 
                                  n_name=n_name, cor_var=cor_var, RAM=RAM,
                                  sem_model=sem_model, cluster_name=cluster_name, 
                                  sep=sep, model.name=model.name, n_fun=n_fun, 
                                  nearpd=nearpd, admissibility=admissibility,
                                  admissibility_fit=admissibility_fit,
                                  raw_Cov=raw_Cov, raw_aCov=raw_aCov,
                                  constrained_B=constrained_B,
                                  constrained_seed=constrained_seed,
                                  tran=tran, run=run, 
                                  status=status, data=data))
  
  if(run) class(out) <- "wls"
  return(out) 
}

# M==============================================================================

metasem_3m <- function(rma_fit, sem_model, n_name, cor_var=NULL, n=NULL, 
                       n_fun=mean, cluster_name=NULL, model.name=NULL,
                       nearpd=FALSE, admissibility=NULL,
                       constrained_B=499L, constrained_seed=20261005L,
                       tran=NULL, run=TRUE, weights="proportional",
                       sep="[^[:alnum:]]+", moderator=NULL, data=NULL, 
                       std.lv=TRUE, auto.var=TRUE, RAM=NULL, ...){
  
  out <- if(is.null(moderator)) { 
    
    metasem_(rma_fit=rma_fit, sem_model=sem_model, 
             n_name=n_name, cor_var=cor_var, n=n, model.name=model.name,
             n_fun=n_fun, cluster_name=cluster_name, run=run,
             nearpd=nearpd, admissibility=admissibility,
             constrained_B=constrained_B, constrained_seed=constrained_seed,
             tran=tran, sep=sep, data=data, weights=weights,
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
                                                         n_fun=n_fun, cluster_name=cluster_name, weights=weights,
                                                         nearpd=nearpd, admissibility=admissibility,
                                                         constrained_B=constrained_B,
                                                         constrained_seed=constrained_seed+i-1L,
                                                         tran=tran, sep=sep, run=run,
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

                               
# M==============================================================================

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
                
# H===============================================================================
                                
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
