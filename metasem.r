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
# Credit to: https://groups.google.com/g/lavaan/c/1ueFiue9qLM/m/e6CMxegiAQAJ?pli=1                                 
# Convert a 'lavaan' Parameter
# Table to a 'lavaan' Model Syntax
                                 
ptable_to_syntax <- function(object) {
    if (inherits(object, "lavaan")) {
        if (!ptable_to_syntax_check_fit(object)) {
            stop("The model in 'fit' is not supported.")
          }
        ptable <- lavaan::parameterTable(object)
        names_lv <- lavaan::lavNames(object, "lv.regular")
        names_eqs_y <- lavaan::lavNames(object, "eqs.y")
      } else {
        ptable <- object
        if (!ptable_to_syntax_check_ptable(ptable)) {
            stop("The parameter table in 'fit' is not supported.")
          }
        names_lv <- pt_lavNames_lv_regular(ptable)
        names_eqs_y <- pt_lavNames_eqs_y(ptable)
      }

    names_def <- unname(ptable[ptable$op == ":=",
                        "lhs", drop = TRUE])

    # lv
    if (length(names_lv) != 0) {
        out_lv <- sapply(names_lv,
                         mod_ind,
                         ptable = ptable)
      } else {
        out_lv <- character(0)
      }

    # y
    if (length(names_eqs_y) != 0) {
        out_y <- sapply(names_eqs_y,
                        mod_y,
                        ptable = ptable)

      } else {
        out_y <- character(0)
      }

    # user
    if (length(names_def) != 0) {
        out_def <- sapply(names_def,
                          mod_def,
                          ptable = ptable)
      } else {
        out_def <- character(0)
      }

    # constraints
    out_eq <- mod_eq(ptable)

    # covariances and variances
    out_cov <- mod_cov(ptable)

    # intercepts
    out_int <- mod_int(ptable)

    out <- c(out_lv,
             out_y,
             out_def,
             out_eq,
             out_cov,
             out_int)
    out <- paste(out,
                 collapse = "\n")
    out
  }

#' @describeIn ptable_to_syntax Compare two lavaan parameter tables.
#' @order 2
#' @export

compare_ptables <- function(object1,
                            object2) {
    if (inherits(object1, "lavaan")) {
        ptable1 <- lavaan::parameterTable(object1)
      } else {
        ptable1 <- object1
      }
    if (inherits(object2, "lavaan")) {
        ptable2 <- lavaan::parameterTable(object2)
      } else {
        ptable2 <- object2
      }
    c0 <- c("lhs", "op", "rhs",
            "group", "free", "ustart",
            "exo", "label", "start",
            "plabel")
    pt1 <- try(ptable1[, c0],
               silent = TRUE)
    if (inherits(pt1, "try-error")) {
        stop("Check whether object1 is a lavaan object or a parameter table?")
      }
    pt2 <- try(ptable2[, c0],
               silent = TRUE)
    if (inherits(pt2, "try-error")) {
        stop("Check whether object2 is a lavaan object or a parameter table?")
      }
    plabel_match <- match_plabels(pt1, pt2)
    pt2 <- update_plabel(pt2, plabel_match)
    # Check lhs, op, rhs
    tmp1 <- pt1[order(pt1$group,
                      pt1$op, pt1$lhs, pt1$rhs), ]
    tmp2 <- pt2[order(pt2$group,
                      pt2$op, pt2$lhs, pt2$rhs), ]
    rownames(tmp1) <- NULL
    rownames(tmp2) <- NULL
    if (!identical(tmp1[, c("group", "lhs", "op", "rhs")],
                   tmp2[, c("group", "lhs", "op", "rhs")])) {
        return(FALSE)
      }
    if (!identical(tmp1$free > 0,
                   tmp2$free > 0)) {
        return(FALSE)
      }
    if (!identical(tmp1$label,
                   tmp2$label)) {
        return(FALSE)
      }
    if (!identical(tmp1$exo,
                   tmp2$exo)) {
        return(FALSE)
      }
    tmp1u <- tmp1[which(tmp1$ustart > 0), ]
    tmp2u <- tmp2[which(tmp2$ustart > 0), ]
    if (!identical(tmp1u$start,
                   tmp2u$start)) {
        return(FALSE)
      }
    return(TRUE)
  }

#' @noRd

ptable_to_syntax_check_fit <- function(fit) {
    if (lavaan::lavTech(fit, "ngroups") != 1) {
        stop("Multigroup models not supported.")
      }
    if (lavaan::lavTech(fit, "nlevels") != 1) {
        stop("Multilevel models not supported.")
      }
    if (!all(lavaan::lavTech(fit, "nclusters") == 1)) {
        stop("Models with clusters not supported.")
      }
    if (length(lavaan::lavNames(fit, "ov.ord")) > 0) {
        stop("Models with ordinal variables not supported.")
      }
    if (any(lavaan::parameterTable(fit)$op == "<~")) {
        stop("Models with operator '<~' not supported.")
      }
    ptable <- lavaan::parameterTable(fit)
    ptmp <- ptable
    ptmp$label <- ""
    coef_names <- lavaan::lav_partable_labels(ptmp)
    user_labels <- setdiff(ptable$label, "")
    if (any(user_labels %in% coef_names)) {
        stop("Does not support constraints imposed by 'equal()'.")
      }
    if (any(sapply(ptable$label, any_space))) {
        stop("Does not support labels with spaces.")
      }
    if (any(sapply(ptable$label, any_op))) {
        stop("Does not support labels with syntax operators.")
      }
    TRUE
  }

#' @noRd

ptable_to_syntax_check_ptable <- function(ptable) {
    if (!is.null(ptable$group)) {
        if (max(ptable$group) != 1) {
            stop("Multigroup models not supported.")
          }
      }
    if (!is.null(ptable$level)) {
        if (max(ptable$level) != 1)
        stop("Multilevel models not supported.")
      }
    # if (!all(lavaan::lavTech(fit, "nclusters") == 1)) {
    #     stop("Models with clusters not supported.")
    #   }
    if ("|" %in% ptable$op) {
        stop("Models with ordinal variables not supported.")
      }
    if ("<~" %in% ptable$op) {
        stop("Models with operator '<~' not supported.")
      }
    ptmp <- ptable
    ptmp$label <- ""
    coef_names <- lavaan::lav_partable_labels(ptmp)
    user_labels <- setdiff(ptable$label, "")
    if (any(user_labels %in% coef_names)) {
        stop("Does not support constraints imposed by 'equal()'.")
      }
    if (any(sapply(ptable$label, any_space))) {
        stop("Does not support labels with spaces.")
      }
    if (any(sapply(ptable$label, any_op))) {
        stop("Does not support labels with syntax operators.")
      }
    TRUE
  }


#' @noRd

match_plabels <- function(pt1, pt2) {
    pt01 <- pt1[pt1$plabel != "",
                c("lhs", "op", "rhs", "group", "plabel")]
    pt02 <- pt2[pt2$plabel != "",
                c("lhs", "op", "rhs", "group", "plabel")]
    tmp <- merge(pt01, pt02,
                 by = c("lhs", "op", "rhs", "group"),
                 suffixes = c(".1", ".2"))
    two2one <- tmp$plabel.1
    names(two2one) <- tmp$plabel.2
    two2one
  }

#' @noRd

update_plabel <- function(pt,
                          target) {
    # TODOs:
    # - Simplify the code
    tmp <- match(pt$lhs, names(target))
    if (any(!is.na(tmp))) {
        pt$lhs[which(!is.na(tmp))] <-
          target[tmp[!is.na(tmp)]]
      }
    tmp <- match(pt$rhs, names(target))
    if (any(!is.na(tmp))) {
        pt$rhs[which(!is.na(tmp))] <-
          target[tmp[!is.na(tmp)]]
      }
    tmp <- match(pt$plabel, names(target))
    if (any(!is.na(tmp))) {
        pt$plabel[which(!is.na(tmp))] <-
          target[tmp[!is.na(tmp)]]
      }
    pt
  }

#' @noRd

mod_int <- function(ptable) {
    pt0 <- ptable[(ptable$op == "~1"), , drop = FALSE]
    k <- nrow(pt0)
    out0 <- character(0)
    if (k == 0) {
        return(character(0))
      }
    for (i in seq_len(k)) {
        pt_i <- pt0[i, ]
        if (pt_i$free == 0) {
            if (is.na(pt_i$ustart)) {
                out0 <- c(out0,
                          character(0))
              } else {
                out0 <- c(out0,
                          paste0(pt_i$lhs,
                                 " ~ ",
                                 pt_i$ustart,
                                 "*1"))
              }
          } else {
            if (pt_i$label == "") {
                outi <- paste0("1")
              } else {
                outi <- paste0(pt_i$label,
                               "*1")
              }
            if (!is.na(pt_i$ustart)) {
                outi <- paste0("start(",
                               pt_i$ustart,
                               ")*",
                               outi)
              }
            out0 <- c(out0,
                      paste0(pt_i$lhs,
                             " ~ ",
                             outi))
          }
      }
    out0
  }

#' @noRd

mod_cov <- function(ptable) {
    pt0 <- ptable[(ptable$op == "~~"), , drop = FALSE]
    k <- nrow(pt0)
    out0 <- character(0)
    if (k == 0) {
        return(character(0))
      }
    for (i in seq_len(k)) {
        pt_i <- pt0[i, ]
        if (pt_i$free == 0) {
            if (is.na(pt_i$ustart)) {
                out0 <- c(out0,
                          character(0))
              } else {
                out0 <- c(out0,
                          paste0(pt_i$lhs,
                                 " ~~ ",
                                 pt_i$ustart,
                                 "*",
                                 pt_i$rhs))
              }
          } else {
            if (pt_i$label != "") {
                outi <- paste0(pt_i$label,
                               "*",
                               pt_i$rhs)
              } else {
                outi <- pt_i$rhs
              }
            if (!is.na(pt_i$ustart)) {
                outi <- paste0("start(",
                               pt_i$ustart,
                               ")*",
                               outi)
              }
            out0 <- c(out0,
                      paste0(pt_i$lhs,
                             " ~~ ",
                             outi))
          }
      }
    out0
  }

#' @noRd

mod_eq <- function(ptable) {
    pt0 <- ptable[(ptable$op == "=="), , drop = FALSE]
    plabels <- setdiff(ptable$plabel, "")
    labels <- setdiff(ptable$label, "")
    k <- nrow(pt0)
    out0 <- character(0)
    if (k == 0) {
        return(character(0))
      }
    i1 <- (pt0$lhs %in% plabels) |
          (pt0$rhs %in% plabels)
    i2 <- (pt0$lhs %in% labels) |
          (pt0$rhs %in% labels)
    pt0 <- pt0[!i1 & i2, ]
    k <- nrow(pt0)
    if (k == 0) {
        return(character(0))
      }
    out0 <- c(out0,
              paste(pt0$lhs, pt0$op, pt0$rhs))
    out0
  }

#' @noRd

mod_def <- function(ptable, def) {
    pt0 <- ptable[(ptable$op == ":="), , drop = FALSE]
    k <- nrow(pt0 > 0)
    out0 <- character(0)
    if (k > 0) {
        for (i in seq_len(k)) {
            out0 <- c(out0,
                      paste(pt0[i, "lhs"],
                            ":=",
                            pt0[i, "rhs"]))
          }
      } else {
        return(out0)
      }
    return(out0)
  }

#' @noRd

mod_y <- function(ptable, eqs_y) {
    pt0 <- ptable[(ptable$lhs == eqs_y) &
                  (ptable$op == "~"), , drop = FALSE]
    k <- nrow(pt0 > 0)
    out0 <- character(0)
    if (k > 0) {
        for (i in seq_len(k)) {
            pt_i <- pt0[i, ]
            if (pt_i$free == 0) {
                out0 <- c(out0,
                          paste0(pt_i$ustart,
                                  "*",
                                  pt_i$rhs))
              } else {
                if (pt_i$label != "") {
                    outi <- paste0(pt_i$label,
                                   "*",
                                   pt_i$rhs)
                  } else {
                    outi <- pt_i$rhs
                  }
                if (!is.na(pt_i$ustart)) {
                    outi <- paste0("start(",
                                   pt_i$ustart,
                                   ")*",
                                   outi)
                  }
                out0 <- c(out0,
                          outi)
              }
          }
      } else {
        return(out0)
      }
    if (length(out0) > 1) {
        out1 <- paste0(out0,
                       collapse = " + ")
      } else {
        out1 <- out0
      }
    out2 <- paste(eqs_y, "~", out1)
    out2
  }


#' @noRd

mod_ind <- function(ptable, lv) {
    pt0 <- ptable[(ptable$lhs == lv) &
                  (ptable$op == "=~"), , drop = FALSE]
    k <- nrow(pt0 > 0)
    out0 <- character(0)
    if (k > 0) {
        for (i in seq_len(k)) {
            pt_i <- pt0[i, ]
            if (pt_i$free == 0) {
                # Must set to 1 because
                # whether it is 1 depends on
                # the call, if omiitted.
                out0 <- c(out0,
                          paste0(pt_i$ustart,
                                  "*",
                                  pt_i$rhs))
              } else {
                if (pt_i$label != "") {
                    outi <- paste0(pt_i$label,
                                   "*",
                                   pt_i$rhs)
                  } else {
                    outi <- pt_i$rhs
                  }
                if (!is.na(pt_i$ustart)) {
                    outi <- paste0("start(",
                                   pt_i$ustart,
                                   ")*",
                                   outi)
                  }
                out0 <- c(out0,
                          outi)
              }
          }
      } else {
        return(out0)
      }
    if (length(out0) > 1) {
        out1 <- paste0(out0,
                       collapse = " + ")
      } else {
        out1 <- out0
      }
    out2 <- paste(lv, "=~", out1)
    out2
  }

#' @noRd

any_space <- function(x) {
    x1 <- gsub(" ", "", x, fixed = TRUE)
    if (nchar(x1) < nchar(x)) {
        return(TRUE)
      }
    return(FALSE)
  }

#' @noRd

any_op <- function(x) {
    ops <- c("~",
             "=~",
             "|",
             "~*~",
             "<~",
             ":=",
             "==")
    # "~~" will be treated as "~"
    # OK for now.
    chk <- sapply(ops,
                  grepl,
                  x = x,
                  fixed = TRUE)
    if (any(chk)) {
        return(TRUE)
      }
    return(FALSE)
  }

#' @noRd

pt_lavNames_lv_regular <- function(ptable) {
    out0 <- unique(ptable[ptable$op == "=~", "lhs",
                          drop = TRUE])
    out0
  }

#' @noRd

pt_lavNames_eqs_y <- function(ptable) {
    out0 <- unique(ptable[ptable$op == "~", "lhs",
                          drop = TRUE])
    out0
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
