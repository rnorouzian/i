
Break = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = " R programs for meta-analyzing complex L2 studies.
 Copyright (C) 2022-present  Reza Norouzian, rnorouzian@gmail.com\n"

message(Break, notice, Break)

Break = "\n***********************************************************************************\n"

cite <- "To cite this package, use:\n\nNorouzian, R., & Bui, G. (2024). Meta-analysis of Second Language Research with\nComplex Research Designs. Studies in Second Languages Acquisition, 46(1), 251-276."

message(Break, cite, Break)

# H===============================================================================================================================

trim_ <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}

# H===============================================================================================================================

abber_case <- function(X, abb_except = NULL, abb_length = 10){
  y <- names(Filter(function(i) is.character(i) | is.factor(i), X[setdiff(names(X), abb_except)]))
  X[y] <- lapply(X[y], abbreviate, minlength=abb_length, named=FALSE)
  return(X)
}                    

# H===============================================================================================================================                    

numerize_case <- function(X, num_except = NULL, num_zero = FALSE){
  y <- names(Filter(function(i) is.character(i) | is.factor(i), X[setdiff(names(X), num_except)]))
  X[y] <- lapply(X[y], function(k) as.integer(as.factor(k))-if(num_zero) 1 else 0)
  return(X)
} 

# H================================================================================================================================

rm.allrowNA <- function(X) { 
  
  X[rowSums(is.na(X) | X == "" | X == "NA" | X == "NULL" | is.null(X) ) != ncol(X), , drop = FALSE]  
  
}

# H===============================================================================================================================

rm.allcolNA <- function(X) { 
  
  X[, colSums(is.na(X) | X == "" | X == "NA" | X == "NULL" | is.null(X) ) != nrow(X), drop = FALSE]
  
}

# H===============================================================================================================================

rm.colrowNA <- function(X){
  
  rm.allcolNA(rm.allrowNA(X))  
  
}

# M===============================================================================================================================
                 
# Functions for intact class assignment in studies (Hedges, 2007)
  
g_cluster <- function(g, n_cluster, 
                      Nt, Nc, 
                      icc = .15) 
{ 
  
  N_tot <- Nt + Nc
  n_bar <- N_tot / n_cluster    
  
  g * sqrt(1 - ((2 * (n_bar - 1) * icc) / 
                  (n_cluster * n_bar - 2)))  
  
}

# M=================================================================================================================================================
  
#g_vi_cluster <- function(g, n_class, 
#                          Nt, Nc, 
#                          icc = .15){

#  N_tot <- Nt + Nc
#  n_bar <- N_tot / n_class 
  
#  eta <-  1 + ((n_bar- 1)*icc)
  
#  DF <- (( (N_tot-2)-((2 * ((N_tot/n_class)-1))*icc)   )^2) /
#    ((N_tot-2)*((1-icc)^2))+((N_tot/n_class)*(N_tot-(2*(N_tot/n_class)) )*icc^2)+
#    (2* ( N_tot-(2*(N_tot/n_class)) )*icc*(1-icc) )
  
#  w <- 1 # w_factor(DF)
  
#  z <- (((N_tot-2)*(1-icc)^2) + (n_bar*(N_tot-2*n_bar)*icc^2) + (2*(N_tot-2*n_bar)*icc*(1-icc))) / 
#   (2*((N_tot-2) - 2*(n_bar-1)*icc)^2)
  
# (w*sqrt( (((Nt+Nc)/(Nt*Nc))*eta) + g^2*z))^2
#}
                 
# M=================================================================================================================================================

g_vi_cluster <- function(g, n_cluster, 
                         Nt, Nc, 
                         icc = .15){
  
  N_tot <- Nt + Nc
  n_bar <- N_tot / n_cluster 
  
  ((Nt+Nc)/(Nt*Nc))*(1 + ((n_bar- 1)*icc)) +
    ( g^2 * ( 
      (((N_tot -2)*(1-icc)^2 ) + (n_bar*(N_tot - 2*n_bar)*icc^2) +
         (2* (N_tot - 2*n_bar) * icc * (1 - icc)) ) /
        ((2* (N_tot-2)) * ( (N_tot-2) - (2* (n_bar-1)*icc) )) 
    )  )
  
}
                 
# H=================================================================================================================================================

w_factor <- function (df) 
{
  w <- ifelse(df <= 1, NA_real_, exp(lgamma(df/2) - log(sqrt(df/2)) - 
                                         lgamma((df - 1)/2)))
  return(w)
}      
                 
# M=================================================================================================================================================
  
gen_vi_cluster <- function(vi, n_cluster, 
                           Nt, Nc, 
                           icc = .15)
{
  
  N_tot <- Nt + Nc
  n_bar <- N_tot / n_cluster 
  
  DEF <- (n_bar - 1) * icc + 1
  DEF*vi
}

# End of functions for intact class assignment in studies    

# H=================================================================================================================================================
                 
add_sig_funnel <- function(funnel, level=.05, col="magenta", sig_cex=1,
                           pch=21, bg="cyan", digits=2, refline = 0){
  
  right <- level/2   
  
  x <- funnel$x
  y <- funnel$y
  crit <- qnorm(right, mean=refline, lower.tail = FALSE)
  
  ci1 <- crit*y
  x1 <- x[x > ci1]
  y1 <- y[x > ci1]
  
  points(x1, y1, col=col, pch=pch, bg=bg, cex=sig_cex)
  
  ci2 <- -crit*y
  x2 <- x[x < ci2]
  y2 <- y[x < ci2]
  
  total <- c(x1,x2)  
  
  points(x2, y2, col=col, pch=pch, bg=bg, cex=sig_cex)
  
  out <- data.frame(total = length(total), total_perc=length(total)/length(x)*100, left = length(x2), perc_left = length(x2) / length(x)*100,
                    right = length(x1), perc_right = length(x1) / length(x)*100, sig = level)
  
  round(setNames(out, c("Total","Total(%)","Left","Left(%)","Right","Right(%)","Sig.")), digits = digits)
}


# M=================================================================================================================================================

contour_funnel <- function(fit = NULL, x, vi, level = c(95),
                           shade = c("white"),
                           xlab = "Effect Size", yaxis = "sei",
                           sig = FALSE, sig_level=.05,
                           sig_col = "magenta", sig_cex=NULL,
                           sig_pch=21, sig_bg="cyan", 
                           sig_digits=2, refline = 0, cex=1, ...){
  
  yaxis <- if(sig) "sei" else yaxis
  x <- if(!is.null(fit)) fit$yi else x
  y <- if(!is.null(fit)) fit$vi else vi
  
  if(refline!=0 & sig) {
    message("'refline' reset to '0'.")
    refline <- 0
  }
  
  
  if(is.null(sig_cex)) sig_cex <- cex
  
  f1 <- metafor::funnel.default(    x = x,
                                    vi = y,
                                    level = level, 
                                    shade = shade,
                                    xlab = xlab, 
                                    yaxis = yaxis,
                                    refline = refline,
                                    cex = cex, ...)
  
  if(sig) add_sig_funnel(f1, refline=refline, level=sig_level, col=sig_col, 
                         pch=sig_pch, bg=sig_bg, sig_cex=sig_cex, digits=sig_digits)
}
                 
# H===============================================================================================================================

get_vars_ <- function(gls_fit, as_fml = TRUE){ 
  
  m <- sub("[^:]+\\(([^)]+).*", "\\1",attr(terms(gls_fit), "term.labels"))
  
  if(as_fml) sapply(paste0("~",m), as.formula) else m
}

# H===============================================================================================================================

is_qdrg <- function(rma_fit){ is.null(rma_fit$call$yi) }

# H===============================================================================================================================

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
                 
# H===============================================================================================================================
                 
rma_clusters <- function (obj) 
{
  level_dat <- vector(mode = "integer")
  cluster_dat <- data.frame(row.names = 1:obj$k)
  
  if (obj$withG) {
    g_nm <- tail(obj$g.names,1)
    level_dat[[g_nm]] <- obj$g.nlevels[[2]]
    cluster_dat[[g_nm]] <- obj$mf.g$outer
  }
  if (obj$withH) {
    h_nm <- tail(obj$h.names,1)
    level_dat[[h_nm]] <- obj$h.nlevels[[2]]
    cluster_dat[[h_nm]] <- obj$mf.h$outer
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

# H===============================================================================================================================

is_nested <- function (cluster, fac) 
{
  if (is.list(fac)) {
    res <- sapply(fac, is_nested, cluster = cluster)
    return(res)
  }
  groupings <- tapply(cluster, fac, function(x) length(unique(x)))
  all(groupings == 1L)
}
                      
# M================================================================================================================================

add_blank_row <- function(data, n_blank_row = 1, by = "study", 
                          file_name = NULL, na = "")
{
  
  data <- full_clean(data)
  dat <- group_split(data, !!rlang::sym(by)) %>% 
    map_dfr(~ .x[1:(nrow(.x) + n_blank_row),])
  
  if(!is.null(file_name)){
    file_name <- paste0(file_name, ".csv")
    write_csv(dat, file_name, na = na)
  }
  return(dat)
} 

# M================================================================================================================================                 

 mask <- function(data, what, full = FALSE){
  
  data[] <- lapply(data, function(x) type.convert(as.character(x), as.is = TRUE))
  
  f1 <- function(x) as.numeric(factor(x, levels = unique(x)))
  f2 <- function(x) {
    temp <- substr(x, 1, 1)
    paste0(temp, ave(x, temp, FUN = function(y) match(y, unique(y))))
  }
  cols <- names(data)[sapply(data, is.numeric)]
  num.cols <- cols[cols %in% what]
  cols <- names(data)[sapply(data, is.character)]
  char.cols <- cols[cols %in% what]
  
  if(!full){
    
    if(length(num.cols))  data[num.cols] <- lapply(data[num.cols], f1)
    if(length(char.cols)) data[char.cols] <- lapply(data[char.cols], f2)
    
  }else{
    
    data[what] <- lapply(data[what], f1)
  }
  return(data)
}                  
                 
# H===============================================================================================================================

tran_detect <- function(rma_fit){
  
  if (is.element(rma_fit$measure, 
                 c("RR","OR","PETO","IRR","ROM",
                   "D2OR","D2ORL","D2ORN","CVR",
                   "VR","PLN","IRLN","SDLN","MNLN",
                   "CVLN","ROMC","CVRC","VRC","REH","HRR"))) {
    "log"
  } else
    
    if (is.element(rma_fit$measure, "PLO")) {
      "logit"
    } else

      if (is.element(rma_fit$measure, "PAS")) {
        emmeans::make.tran("asin.sqrt", 1)
        
      } else
        
        if (is.element(rma_fit$measure, "IRS")) {
          "sqrt"
          
        } else
          
          if (is.element(rma_fit$measure, c("ZPHI","ZTET","ZPB","ZBIS","ZCOR","ZPCOR","ZSPCOR"))) {
            
            r2z_tran
            
          } else 
            
            if (is.element(rma_fit$measure, "AHW")) {
              
              AHW_tran 
              
            } else 
              
              if (is.element(rma_fit$measure, "ABT")) {
                
                ABT_tran
              
              } else
                if (is.element(rma_fit$measure, "ZR2")){
                  
                  ZR2_tran  
                  
                } else FALSE
  
}

# H===============================================================================================================================

sigma_detect <- function(rma_fit){
  
  if (!inherits(rma_fit, c("rma.mv"))) {
    sqrt(rma_fit$tau2)
  } else {
    0
  }
} 

# H===============================================================================================================================

df_detect <- function(rma_fit){
if(is.na(rma_fit$ddf[1])) Inf else min(rma_fit$ddf, na.rm = TRUE)
}                 
                 
# H===============================================================================================================================

get_data_ <- function(fit){
  
  f <- function (object) {
  
  if ("data" %in% names(object)) {
    data <- object$data
  } else {
    dat_call <- object$call$data
    envir_names <- sys.frames()
    data <- simpleError("start")
    i <- 1L
    while (inherits(data, "simpleError") & i <= length(envir_names)) {
      data <- tryCatch(eval(dat_call, envir = envir_names[[i]]), error = function(e) e)
      i <- i + 1L
    }
  }
  
  if (inherits(data, "simpleError")) return(NULL)
  
  naAct <- object[["na.action"]]
  if (!is.null(naAct)) {
    data <- if (inherits(naAct, "omit")) {
      data[-naAct, ]
      
    } else if (inherits(naAct, "exclude")) {
      data
    } else eval(object$call$na.action)(data)
  }
  
  subset <- object$call$subset
  if (!is.null(subset)) {
    subset <- eval(asOneSidedFormula(subset)[[2]], data)
    data <- data[subset, ]
  }
  
  data
}                 
  
  data_ <- try(f(fit), silent = TRUE)
  if(!inherits(data_, "try-error")) data_ else recover_data(fit) 
}

# H=============================================================================================================================== 

odds_. <- function(x) subset(x, x %% 2 != 0)

# H===============================================================================================================================

print.post_rma <- function(post_rma_call){
  print(post_rma_call$table)
}                 

# H=============================================================================================================================== 

print.con_rma <- function(con_rma_call){
  print(con_rma_call$table)
}                  

# H=============================================================================================================================== 
                  
print.contrast_rma <- function(contrast_rma_call){
  print(contrast_rma_call$table)
}                  
                  
# H=============================================================================================================================== 

full_clean <- function(data) rm.colrowNA(trim_(data))           

# H===============================================================================================================================

odiag <- function(x) suppressMessages(x[col(x) != row(x)])                   

# H===============================================================================================================================

shift_rows <- function(data, user_index, up = TRUE){
  
  indx <- seq_len(nrow(data))
  remain <- indx[!indx %in% user_index]
  
  data[if(up) c(user_index, remain) else c(remain, user_index), ]
}                    

# H===============================================================================================================================

get_error_rho <- function(fit){
  
  is_V <- any(odiag(fit$V) != 0)  
  
  if(is_V) {
    
    u <- unique(odiag(round(cov2cor(fit$V), 7)))
    u[u!=0]
    
  } else 0
  
}                    

# H===============================================================================================================================                    

is_crossed <- function(obj){
  
  if(!inherits(obj, "rma.mv")) return(FALSE)
  mod_struct <- rma_clusters(obj)
  highest_cluster <- names(mod_struct$level_dat)[which.min(mod_struct$level_dat)]
  cluster <- mod_struct$cluster_dat[[highest_cluster]]
  out <- !is_nested(cluster, fac = mod_struct$cluster_dat)
  return(out)
}      

# H===============================================================================================================================

pluralify_ <- function (x, keep.original = FALSE, 
                        irregular = lexicon::pos_df_irregular_nouns) {
  
  stopifnot(is.data.frame(irregular))
  
  hits <- match(tolower(x), tolower(irregular[[1]]))
  
  ends <- "(sh?|x|z|ch)$"
  plural_ <- ifelse(grepl(ends, x), "es", "s")
  out <- gsub("ys$", "ies", paste0(x, plural_))
  out[which(!is.na(hits))] <- irregular[[2]][hits[which(!is.na(hits))]]
  
  c(if (keep.original) {
    x
  }, out)
}           

# M================================================================================================================================     

cat_pattern <- function(data, cluster, ..., blank_sign = "*"){
  
  data <- full_clean(data)
  cluster <- rlang::ensym(cluster)
  cat_mod <- rlang::ensyms(...)
  cat_nms <- purrr::map_chr(cat_mod, rlang::as_string)
  c_nm <- as_string(cluster)
  nms <- c(c_nm, cat_nms)
  
  idx <- nms %in% names(data)
  if(!all(idx)) stop(toString(dQuote(nms[!idx]))," not found in the 'data'.", call. = FALSE)
  
  setNames(purrr::map(cat_mod,  ~ {
    
    studies_cats <- 
      data %>%
      dplyr::group_by(!!cluster, !!.x) %>%
      dplyr::summarise(effects = n(), .groups = 'drop')
    nm1 <- rlang::as_string(.x)
    cat_names <- paste0(nm1, c(".x", ".y"))
    
    studies_cats <- 
      studies_cats %>%
      dplyr::inner_join(studies_cats, by = rlang::as_string(cluster)) %>%
      dplyr::group_by(!!!rlang::syms(cat_names)) %>%
      dplyr::summarise(
        studies = n(),
        effects = sum(effects.x), .groups = 'drop') %>% 
      dplyr::mutate(n = paste0(studies, " (", effects, ")") )
    
    out1 <- studies_cats %>%
      dplyr::select(-studies, -effects) %>%        
      tidyr::pivot_wider(names_from = cat_names[2], 
                         values_from = n, names_sort = TRUE) %>%
      dplyr::rename_with(~nm1,  cat_names[1]) %>%
      dplyr::arrange(dplyr::across(tidyselect::all_of(nm1))) %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(nm1), as.character))  
    
    out2 <- out1[-1]
    out2[upper.tri(out2)] <- blank_sign
    dplyr::bind_cols(out1[1],out2)
    
  }), cat_nms)
} 

# M===============================================================================================================================
       
bw <- function(data, cluster_name, var_names){
  
  all_names <- trimws(names(data))
  
  id <- grep("study", all_names, value=TRUE, ignore.case=TRUE)
  
  if(!all(cluster_name %in% all_names)) { 
    
    stop(paste(toString(dQuote(cluster_name)), "not found for 'cluster_name' in the 'data'.", if(length(id)>0) 
      paste("\nPossibly you meant to use", toString(dQuote(id)), "as 'cluster_name'?")), call. = FALSE) 
    
  }
  
  ok <- var_names %in% all_names
  
  if(!all(ok)) message(paste("\n", toString(dQuote(var_names[!ok])), "not found in the 'data' thus ignored."))
  
  var_names <- var_names[ok] 
  
  dum_vars <- all_names[sapply(data, function(i)is.character(i)|is.factor(i))]
  
  dum_names <- var_names[var_names %in% dum_vars]
  
  is_dum <- length(dum_names) > 0
  
  num_names <- var_names[!(var_names %in% dum_vars)]
  
  is_num <- length(num_names) > 0
  
  
  if(is_num){
    data <- data %>%
      group_by(across(all_of(cluster_name))) %>%                          
      mutate(across(all_of(num_names), list(wthn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
      as.data.frame()
  }
  
  
  if(is_dum){
    data <- data %>%
      dummy_cols(select_columns = dum_names) %>% 
      group_by(across(all_of(cluster_name))) %>%                          
      mutate(across(starts_with(paste0(dum_names, "_")), list(wthn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
      as.data.frame()
  }
  
  return(data)
}

# M================================================================================================================================

lo_ave_up <- function(var_names, x, at_vals = NA, digits = 0){
  
  data <- if(inherits(x, c("rma.uni", "rma.mv"))) { get_data_(x) }
  
  else if(inherits(x, c("data.frame","tibble"))){ x } else {
    
    stop("'x' must be either a 'data.frame' or an 'rma.uni'/'rma.mv' object.")
  }
  
  data <- full_clean(data)
  
  var_names <- if(is_bare_formula(var_names, lhs=FALSE)) .all.vars(var_names) else if(is.character(var_names)) var_names 
  else stop("'var_names=' can be either a character vector (ex. 'year') or a one-sided formula (ex. ~year).", call. = FALSE)
  
  idx <- var_names %in% names(data)
  if(!all(idx)) stop(toString(dQuote(var_names[!idx])," not found in the data."), call. = FALSE)
  
  num_var_names <- names(Filter(is.numeric, data[var_names]))
  
  if(all(is.na(at_vals)) & length(at_vals) != length(var_names)) at_vals <- rep(at_vals, length(var_names))
  
  out <- ifelse(is.na(at_vals) & var_names %in% num_var_names,
                
                sapply(num_var_names, function(x)
                  round(setNames(mean(data[[x]], na.rm=TRUE) + c(-1, 0, 1)*sd(data[[x]], na.rm=TRUE),
                                 paste0(x, c('-1SD', '.Mean', '+1SD'))), digits), simplify = FALSE)
                , 
                
                type.convert(lapply(seq_along(var_names), function(i) at_vals[i]), as.is=TRUE)
                
  )
  
  out <- if(length(out)==0) NULL else setNames(out, var_names)
  
  out <- Filter(Negate(anyNA), out)
  
 lapply(out, function(i) if(inherits(i, "list")) unlist(i, recursive=FALSE) else i)
  
}
                                                
# H================================================================================================================================ 

data.tree_ <- function(data, toplab = NULL, cex = 1, rowcount = FALSE, cex_top = 1, ...){
  
  toplab <- if(is.null(toplab)) names(data) else toplab
  
  sizetree2(data, toplab = toplab, stacklabels = FALSE, border = 0, base.cex = cex, showcount = rowcount, cex_top = cex_top, ...)
}           

# M===============================================================================================================================

UN_cor_fixer <- function(data, inner, outer, min_cooccur=4) { 
  
  data <- full_clean(data)
  
  mat <- crossprod(table(data[c(outer,inner)])>0)
  lower <- mat[lower.tri(mat)]
  ifelse(lower < min_cooccur, 0, NA)
  
}
                               
# M===============================================================================================================================                               

meta_tree <- function(data, ..., effect = TRUE, highest_level_name = NULL,  
                      structure = c("simple","typical","complex"), toplab = NULL, 
                      main = NULL, main_extra = NULL, rowcount = FALSE, 
                      abb_names = FALSE, abb_length = 6, abb_except = NULL, 
                      num_names = FALSE, num_except = NULL, num_zero = FALSE, 
                      panel_label = TRUE, cex = 1, cex_main = 1, rev_order = TRUE, 
                      rev_page = FALSE, reset = TRUE, index = NULL, cex_top = 1, 
                      detail=TRUE, subset) 
{
  
  data <- full_clean(data) %>%
    dplyr::mutate(effect = dplyr::row_number())
  
  if(!missing(subset)) {
    
    s <- substitute(subset)
    data <- filter(data, eval(s))
  }
  
  dot_cols <- rlang::ensyms(...)
  if(length(dot_cols)==0) stop("Input at least one variable from the 'data'.", call. = FALSE)
  if(length(dot_cols)==1) effect <- TRUE
  dot_cols <- if(effect) append(dot_cols, rlang::sym("effect")) else dot_cols
  str_cols <- purrr::map_chr(dot_cols, rlang::as_string)
  
  ss <- dot_cols[[1]]
  sss <- str_cols[1]
  
  idx <- str_cols %in% names(data)
  if(!all(idx)) stop(toString(dQuote(str_cols[!idx])," not found in the data."), call. = FALSE)
  
  main_org <- main
  
  if(!is.null(highest_level_name) & abb_names) abb_except <- c(sss, abb_except)
  if(!is.null(highest_level_name) & num_names) num_except <- c(sss, num_except)
  
  if(abb_names) data <- abber_case(data, abb_length = abb_length, abb_except = abb_except) 
  if(num_names) data <- numerize_case(data, num_except = num_except, num_zero = num_zero)
  
  if(reset){
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
  }
  
  data <- data %>%
    dplyr::select(!!!dot_cols)
  
  if(is.null(highest_level_name)){
    
    struc <- match.arg(structure) 
    
    hlist <- data %>%
      dplyr::group_by(!!sym(sss)) %>%
      dplyr::mutate(grp = dplyr::across(tidyselect::all_of(str_cols[-1]), if(detail) ~ {
        tmp <- dplyr::n_distinct(.)
        dplyr::case_when(tmp == 1 ~ 1, tmp == n() ~ 2, tmp > 1 & tmp < n() ~ 3,  TRUE ~ 4)
      } else { ~ n_distinct(.) == 1 }) %>%
        purrr::reduce(stringr::str_c, collapse = "")) %>%
      dplyr::ungroup(.) %>%
      dplyr::group_split(grp, .keep = FALSE)
    
    
    hlist <- if(rev_order) rev(hlist) else hlist
    
    res <- Filter(NROW, hlist)
    
    LL <- length(res)
    
    if(!is.null(index)){ 
      
      index <- index[index <= LL & index >= 1]
      if(length(index)==0) index <- NULL
    } 
    
    res <- if(!is.null(index)) res[index] else res
    
    main_no. <- sapply(res, function(i) length(unique(i[[sss]])))
    
    typic <- function(vec) round(median(vec)) # vec[ceiling(length(vec)/2)]
    
    nms <- lapply(res, function(i){
      nr <- sapply(split(i, i[[sss]]), nrow);
      study_type <- if(struc == "typical") {typic(as.numeric(names(table(nr))))
      } else if(struc == "simple") {min(as.numeric(names(table(nr))))
      } else {max(as.numeric(names(table(nr))))};
      names(nr)[nr == study_type][1]
    })
    
    list2plot <- lapply(seq_along(res),function(i) filter(res[[i]], !!ss == nms[i]))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { 
      
      dev <- if(!rev_page) n2mfrow(LL) else rev(n2mfrow(LL))
      par(mfrow = dev) 
      
    }
    
    main <- if(is.null(main)) stringr::str_to_title(ifelse(main_no. > 1, pluralify_(sss), sss)) else main
    
    main <- if(is.null(main_org)) paste(main_no., main) else main
    
    if(panel_label) {
      
      pan_lab <- make.unique(rep(LETTERS,1e1))[seq_along(list2plot)]
      
      main <- paste0("(",pan_lab,") ", main)
    }
    
    if(!is.null(main_extra)) main <- paste0(main, " [",main_extra,"]")
    
    invisible(lapply(seq_along(list2plot), function(i) data.tree_(list2plot[[i]], main = main[i], toplab, cex, rowcount, cex.main = cex_main, cex_top = cex_top)))
    
    invisible(if(panel_label) setNames(res, pan_lab) else res)
    
  } else {
    
    highest_level_name <- trimws(highest_level_name)
    highest_level_names <- unique(data[[sss]])
    
    idx <- highest_level_name %in% highest_level_names 
    
    if(!all(idx)) stop(toString(dQuote(highest_level_name[!idx]))," not found in variable ", paste0("'",sss,"'", "."), call. = FALSE)
    
    list2plot <- lapply(highest_level_name, function(i) filter(data, !!ss == i))
    
    LL <- length(list2plot)
    
    if(!is.null(index)){ 
      
      index <- index[index <= LL & index >= 1]
      if(length(index)==0) index <- NULL
    } 
    
    list2plot <- if(!is.null(index)) list2plot[index] else list2plot
    
    LL <- length(list2plot)
    
    if(LL > 1L) { 
      
      dev <- if(!rev_page) n2mfrow(LL) else rev(n2mfrow(LL))
      par(mfrow = dev) 
      
    }
    
    invisible(lapply(seq_along(list2plot), function(i) data.tree_(list2plot[[i]], main = main[i], toplab, cex, rowcount, cex.main = cex_main, cex_top = cex_top)))
    
    invisible(list2plot)
  }
}
                        
# M================================================================================================================================================

interactive_outlier <- function(fit, cook = NULL, st_del_res_z = NULL, 
                                cex_add_point = .5, cex_multi_point = 1.2,
                                whisker_coef = 2.5, cex_text_outlier = .6,
                                cex_main = .9, parallel = "no", ncpus = 1, 
                                reestimate = FALSE, 
                                file_name_cook = NULL,
                                file_name_res_z = NULL,
                                view = 1, pos = 2)
{
  
  if(!inherits(fit,c("rma.mv"))) stop("Model is not 'rma.mv'.", call. = FALSE)
  datziola <- get_data_(fit) %>%
    mutate(obsss = fit$slab)
  
  hat <- hatvalues.rma.mv(fit)
  
  if(is.null(cook)){
    
    cook <- cooks.distance.rma.mv(fit, progbar=TRUE,
                                  parallel = parallel, 
                                  ncpus = ncpus, 
                                  reestimate = reestimate)
    
    if(!is.null(file_name_cook)){
      
      filenm <- paste0(file_name_cook,".rds")
      saveRDS(cook, filenm)
      
      message("\nNote: Check folder '", basename(getwd()),"' for the ", dQuote(filenm)," file.\n") 
    }
  }
  
  if(is.null(st_del_res_z)){
    
    st_del_res_z <- rstudent.rma.mv(fit, progbar=TRUE,
                                    parallel = parallel, 
                                    ncpus = ncpus, 
                                    reestimate = reestimate)$z
    
    st_del_res_z <- setNames(st_del_res_z, seq_along(st_del_res_z))
    
    if(!is.null(file_name_res_z)){
      
      filenm <- paste0(file_name_res_z,".rds")
      saveRDS(st_del_res_z, filenm)
      
      message("\nNote: Check folder '", basename(getwd()),"' for the ", dQuote(filenm)," file.\n") 
    }
  }
  
  # Make visual size of effects proportional to their hat/cook's distance (estimate influence)
  cex <- cex_add_point+cex_multi_point*sqrt(if(view == 1)cook else hat)
  
  outlier_limits <- qnorm(c(.025,.5,.975))
  ylim <- range(c(outlier_limits, st_del_res_z))
  
  # Plot Leverage against Studentized residuals proportioned on cook's distances
  plot(if(view == 1) hat else cook, st_del_res_z, cex=cex, las=1, mgp=c(1.5,.3,0),
       xlab = if(view == 1) "Leverage (Hat Value)" else "Effect Influence (Cook's Dis.)", 
       ylab = "Outlier (Standardized Del. Value)",pch=19,cex.axis = .9,tcl = -.3,
       col = adjustcolor(1, .5),
       ylim = ylim)
  
  title(if(view == 1) "Size of points denote \nestimate-influencing effects\n (Cook's distances)" 
        else "Size of points denote \nleverage effects\n (Hat value)", 
        cex.main = cex_main, line = .3)
  
  abline(h=outlier_limits, lty=c(3,1,3), lwd=c(1,2,1))
  
  max_hat <- max(mean(range(hat)), boxplot.stats(hat, coef = whisker_coef)$stats[5])
  
  max_cook <- max(mean(range(cook)), boxplot.stats(cook, coef = whisker_coef)$stats[5])
  
  abline(v = if(view == 1) max_hat else max_cook, col=2)
  
  # To be outlier, an estimate must simultaneously (a) be outlying (per studentized value)
  # (b) have high leverage (per hat value), and (c) high model influence (per cook's distance)
  
  i <- abs(st_del_res_z) > outlier_limits[3]  
  j <- hat > max_hat
  k <- cook > max_cook
  L <- which(i & j & k)
  
  if(length(L)==0) { 
    
    message("Note: No interactive outlier detected.") 
    
    return(NA)
  }
  
  u <- par()$usr
  
  if(any(st_del_res_z[L]>0)) rect(if(view == 1) max_hat else max_cook, outlier_limits[3], u[2], u[4], col = adjustcolor(2, .2), border = NA)
  if(any(st_del_res_z[L]<0)) rect(if(view == 1) max_hat else max_cook, outlier_limits[1], u[2], u[3], col = adjustcolor(2, .2), border = NA)
  
  # Show which effects meet all the three conditions
  text(if(view == 1) hat[L] else cook[L], st_del_res_z[L], labels = names(L), pos = pos, col = "magenta", cex = cex_text_outlier, xpd = NA)
  points(if(view == 1) hat[L] else cook[L], st_del_res_z[L], cex=cex[L])
  
  LL <- names(L)
  
  removed <- filter(datziola, obsss %in% LL)
  new_data <- filter(datziola, !obsss %in% LL)
  
  return(list(removed = removed, 
              new_data = new_data))
}      
  
# H=================================================================================================================================================

term_names_ <- function(post_rma_fit, sep = get_emm_option("sep"), na.rm =TRUE){
  
  if(!inherits(post_rma_fit, c("post_rma","contrast_rma"))) stop("post_rma_fit is not 'post_rma()' or 'contrast_rma()'.", call. = FALSE)  
  
  ems <- if(inherits(post_rma_fit, "post_rma")) post_rma_fit$ems else 
    post_rma_fit$con 
  
  if(na.rm) nas <- is.na(predict(ems))
  
   disp <- if(is.null(ems@misc$display)) seq_len(nrow(ems@linfct)) else ems@misc$display
   if(na.rm) disp <- disp[!nas]
   largs <- as.list(ems@grid[disp, seq_along(ems@levels), drop = FALSE])
   largs$sep <- sep
   do.call(paste, largs)  
}       
# H=================================================================================================================================================

fixed_form_rma <- function(fit){ 
  
  a <- fit$formula.yi
  b <- fit$formula.mods
  y <- fit$call$yi
  
  if(!is.null(a) & !is.null(b)) a 
  else if(!is.null(a) & is.null(b)) a
  else if(is.null(a) & !is.null(b)) as.formula(paste(as.character(y), paste(as.character(b), collapse = "")))
  else as.formula(paste(as.character(y), "~ 1"))
}

# H=================================================================================================================================================  

shorten_ <- function(vec, n = 3) { 
  
  gsub("\\s+", "", 
       sub(sprintf("^(([^,]+, ){%s}).*, ([^,]+)$", n), "\\1...,\\3", toString(vec)))  
} 

          
# H=================================================================================================================================================

random_left <- function(random_fml) {
  
  as.formula(as.character(parse(text = sub("\\|.*", "", random_fml))))
  
}

# H=================================================================================================================================================

rma2gls <- function(fit){
  
  fit <- if(!is_qdrg(fit)){ fit } else {
  
  cl <- fit$call
    
  cl[[1]] <- quote(escalc)
  
  escalc_ <- suppressWarnings(eval(cl))
  
  form_ <- if(fit$int.only) { yi~1 } else { as.formula(paste0("yi", paste(as.character(fit$formula.mods),collapse = "")))}
  
  rma(form_, vi=vi, data=escalc_, test = fit$test, method = fit$method)
  
  }
  
  data_. <- get_data_(fit)
  form_. <- fixed_form_rma(fit)
  
  rownames(fit$b)[rownames(fit$b) %in% "intrcpt"] <- "(Intercept)"
  rownames(fit$beta)[rownames(fit$beta) %in% "intrcpt"] <- "(Intercept)"
  rownames(fit$vb)[rownames(fit$vb) %in% "intrcpt"] <- "(Intercept)"
  colnames(fit$vb)[colnames(fit$vb) %in% "intrcpt"] <- "(Intercept)"
  
  fit2 <- nlme::gls(form_., data = data_., na.action = "na.omit",
                    control = glsControl(singular.ok = TRUE))
  
  fit2$call$model <- form_.
  fit2$call$data <- data_.
  
  fit2$coefficients <- coef.rma(fit) 
  fit2$varBeta <- fit$vb
  
  return(fit2)
}  

# H================================================================================================================================================

clean_reg <- function(fm, nm, uniq = TRUE) {
  
  vars <- vapply(attr(terms(fm), "variables"), deparse, "")[-1L]
  subpat <- paste0(gsub("([()])", "\\\\\\1", vars), collapse = "|")
  l <- rapply(strsplit(nm, ":"), sub, how = "list",
              perl = TRUE,
              pattern = sprintf("^(?!(%1$s)$)(%1$s)(.+)$", subpat),
              replacement = "\\3")
  vec <- vapply(l, paste0, "", collapse = ":")
  if(uniq) vec <- make.unique(vec)
  vec[vec=="intrcpt"] <- "(Intercept)"
  return(vec)
}     

# H=================================================================================================================================================     

any_num_vec <- function(vec){
  any(grepl("^[0-9]{1,}$", gsub(":", "", vec)))
}

# M=================================================================================================================================================    

results_rma <- function(fit, digits = 3, robust = FALSE, blank_sign = "", 
                        cat_shown = 1, shift_up = NULL, shift_down = NULL, 
                        drop_rows = NULL, drop_cols = NULL, QM = TRUE, 
                        QE = FALSE, sig = TRUE, clean_names = NULL, 
                        tidy = FALSE, tol_large = 1e4, random_only = FALSE){
  
  if(!inherits(fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  
  fixed_eff <- is.null(fit$random)
  cr <- if(!fixed_eff) is_crossed(fit) else FALSE
  
  lm_fit <- lm(fixed_form_rma(fit), data = get_data_(fit), na.action = "na.omit")
  
  cl <- clean_reg(lm_fit, names(coef(lm_fit)))
  
  if(is.null(clean_names)){
    
    if(any_num_vec(cl)) {
      
      clean_names <- FALSE
      
    } else { 
      
      clean_names <- TRUE
      
    }
  } 
  
  if(clean_names) names(lm_fit$coefficients) <- cl
  if(clean_names) fit <- clean_reg_names(fit)
  
  lm_coef <- coef(lm_fit)
  
  is_singular <- anyNA(lm_coef)
  
  if(is_singular) message("Note:",dQuote(toString(names(lm_coef)[is.na(lm_coef)])), " dropped due to lack of data.\n")
  
  if(robust & any(cr) || robust & fixed_eff) { 
    
    robust <- FALSE
    message("Note: Robust estimation not available for models with", if(any(cr))" crossed random-" else " only fixed-", "effects.\n")
  }
  
  res <- if(!robust) { 
    
    a <- coef(summary(fit))
    
    nm <- c("Estimate","SE","t","Df","p-value","Lower","Upper")
    
    nm <- if(fit$test == "t") { nm } else { nm[3] <- "z";
    nm[-4] }
    
    setNames(a, nm)
    
  } else {
    
    a <- suppressWarnings(as.data.frame(conf_int(fit, vcov = "CR2"))[-1])
    b <- suppressWarnings(coef_test(fit, vcov = "CR2"))
    a$p <- b$p_Satt
    a$t <- b$tstat
    
    a <- a[c(1:2,7,3,6,4:5)]
    
    setNames(a, c("Estimate","SE","t","Df","p-value","Lower","Upper"))
  }
  
  rn <- rownames(res)[1]
  
  if(rn == "intrcpt") rownames(res)[1] <- "(Intercept)"
  
  res_org <- res
  res <- na.omit(res)
  
  if(robust & nrow(res) != nrow(res_org)) message("Note:",dQuote(toString(setdiff(rownames(res_org),rownames(res)))), " dropped due to inestimablity under Robust estimation.\n")
  
  if(random_only) 
  { QE <- FALSE  
  QM <- FALSE
  sig <- FALSE
  drop_cols <- 1:7
  drop_rows <- 1:(nrow(res)+1)
  }
  
  if(QE){
    qe <- data.frame(Estimate = fit$QE, Df = nobs.rma(fit), 
                     pval = fit$QEp, row.names = "QE") %>%
      dplyr::rename("p-value"="pval") 
    
    res <- bind_rows(res, qe)
  }
  
  if(QM){
    
    mc <- try(clubSandwich::Wald_test(fit, constrain_zero(fit$btt), "CR2"), silent = TRUE)   
    
    bad <- inherits(mc,"try-error")
    
    if(robust){
      
      if(bad || !bad && is.na(mc$p_val)) { 
        robust <- FALSE
        message("Note: Robust QM unavailable,likely: \n1- Some moderators in <2 clusters OR/AND \n2- High # of coefficients vs. # of highest clusters.\n3- Some combination of variables in the interactive model are missing.\nQM results are model-based.\n")
      }
      
      if(!bad && mc$Fstat>tol_large) { message("Note: Robust QM is unreasonably large, likely robust estimation is unfit for the model (use 'robust=FALSE').") }
    }
    
    qm <- if(robust) {
      
      data.frame(Estimate = mc$Fstat, Df = mc$df_num, 
                 pval = mc$p_val, row.names = "QM") %>%
        dplyr::rename("p-value"="pval") 
      
    } else {
      
      data.frame(Estimate = fit$QM, Df = fit$QMdf[1], 
                 pval = fit$QMp, row.names = "QM") %>%
        dplyr::rename("p-value"="pval") 
      
    }
    res <- bind_rows(res, qm)
  }
  
  u <- get_error_rho(fit)
  cte <- length(u) == 1
  
  d6 <- data.frame(r = if(cte) u else mean(u, na.rm = TRUE), 
                   row.names = paste0("Within Corr. (",if(cte) "constant" else "average",")"))
  
  blk <- paste0(paste0(rep(" ",digits-1), collapse=""), "NA", collapse ="")
  
  if(fixed_eff) { 
    
    out <- roundi(dplyr::bind_rows(res, d6), digits = digits)
    out[out== blk] <- blank_sign
    
    if(sig){ 
      
      p.values <- as.numeric(out$"p-value")
      
      Signif <- symnum(p.values, corr = FALSE, 
                       na = FALSE, cutpoints = 
                         c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       symbols = c("***", "**", "*", ".", " ")) 
      
      out <- add_column(out, Sig. = Signif, .after = "p-value")
    }
    
    if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
    if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
    if(!is.null(drop_rows)) out <- out[-drop_rows, ]
    
    out <- dplyr::select(out, -tidyselect::all_of(drop_cols))
    
    if(tidy) out <- cbind(Terms = rownames(out), set_rownames_(out, NULL))
    
    return(out)
  }
  
  res <- rbind(res, "|RANDOM|" = NA)
  
  #on.exit(Sys.setlocale("LC_ALL"))                                    
  #Sys.setlocale(locale = "Greek")
  
  if(fit$withS){
    
    cr <- cr[names(cr) %in% fit$s.names]
    d1 <- data.frame(Sigma = sqrt(fit$sigma2), 
                     row.names = paste0(names(cr), ifelse(cr," [crossed]"," [nested]"))) 
    
    #d1 <- setNames(d1, intToUtf8(963))
  } else { d1 <- NULL}
  
  if(fit$withG){
    
    out_nm <- tail(fit$g.names,1)
    
    outer_nm <- paste(out_nm, if(out_nm %in% names(cr)) "[crossed]" else "[nested]") 
    
    h <- paste(fit$struct[1], "Corr.")
    is_un <- fit$struct[1] == "UN"
    is_gen <- fit$struct[1] == "GEN"
    is_diag <- fit$struct[1] == "DIAG"
    is_simple <- length(fit$tau2) == 1
    
    rnm <- paste("Outer:", outer_nm)
    clnm <- clean_GH_names(fit)
    
    d2 <- data.frame(Tau = sqrt(fit$tau2), 
                     row.names = paste0(if(!is_simple) clnm else fit$g.names[1],
                                        paste0(if(is_diag)" (Uncor. " 
                                               else " (Cor. ",if(!is_simple & !is_gen) paste0(" ", fit$g.names[1]),")")))
    
    #d2 <- setNames(d2, intToUtf8(964))
    
    d2 <- rbind(NA, d2)
    rownames(d2)[1] <- rnm
    
    d3 <- data.frame(Rho = fit$rho, 
                     row.names = if(is_un || is_gen) apply(combn(clnm,2),2,paste0, collapse = "~") 
                     else paste0(h," (",shorten_(clnm, cat_shown),")")) 
    
    #d3 <- setNames(d3, intToUtf8(961))
    
  } else { d2 <- NULL; d3 <- NULL}
  
  if(fit$withH){
    
    out_nm <- tail(fit$h.names,1)
    
    outer_nm <- paste(out_nm, if(out_nm %in% names(cr)) "[crossed]" else "[nested]") 
    
    h <- paste(fit$struct[2], "Corr.")
    is_un <- fit$struct[2] == "UN"
    is_gen <- fit$struct[2] == "GEN"
    is_diag <- fit$struct[2] == "DIAG"
    is_simple <- length(fit$gamma2) == 1
    
    rnm <- paste("Outer:", paste0(outer_nm," "))
    
    clnm <- clean_GH_names(fit, G=FALSE)
    
    d4 <- data.frame(Gamma = sqrt(fit$gamma2), 
                     row.names = paste0(if(!is_simple) clnm else fit$h.names[1],
                                        paste0(if(is_diag)" (Uncor. " 
                                               else " (Cor. ",if(!is_simple) paste0(" ",if(!is_gen)fit$h.names[1]),") "))) 
    
    #d4 <- setNames(d4, intToUtf8(933))
    
    d4 <- rbind(NA, d4)
    rownames(d4)[1] <- rnm
    
    d5 <- data.frame(Phi = fit$phi, 
                     row.names = if(is_un || is_gen) apply(combn(clnm,2),2,paste0, collapse = "~ ")
                     else paste0(h," (",shorten_(clnm, cat_shown),") "))
    
   # d5 <- setNames(d5, intToUtf8(966))
    
  } else { d4 <- NULL; d5 <- NULL}
  
  out <- roundi(dplyr::bind_rows(res, d1, d2, d3, d4, d5, d6), digits = digits)

  out[out== blk] <- blank_sign
  
  if(sig){ 
    
    p.values <- as.numeric(out$"p-value")
    
    Signif <- symnum(p.values, corr = FALSE, 
                     na = FALSE, cutpoints = 
                       c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("***", "**", "*", ".", " ")) 
    
    out <- add_column(out, Sig. = Signif, .after = "p-value")
  }
  
  if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
  if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
  if(!is.null(drop_rows)) out <- out[-drop_rows, ]
  if(!is.null(drop_cols)) out <- dplyr::select(out, -tidyselect::all_of(drop_cols))
  
  if(tidy) out <- cbind(Terms = rownames(out), set_rownames_(out, NULL))
  
  return(out)
}                                          

# H=================================================================================================================================================                     

roundi <- function(x, digits = 7, except = NULL){
  
  if(!inherits(x, c("data.frame","tibble","matrix"))) stop("'x' must be a 'data.frame' or a 'matrix'.", call. = FALSE)
  if(inherits(x,"matrix")) x <- as.data.frame(x)
  
  num <- names(Filter(function(i) is.numeric(i), x[setdiff(names(x), except)]))
  
  x[num] <- lapply(x[num], function(i) formatC(round(i, digits), digits, format = "f"))
  
  return(x)
}                       

# H=================================================================================================================================================

random_GH_form <- function(fit, G = TRUE){
  
  fm <- fit$random
  if(G) fm[[1]] else fm[[2]]
  
}

# H=================================================================================================================================================

clean_GH_names <- function(fit, G = TRUE) {
  
  fmla <- random_left(random_GH_form(fit, G = G))
  vec <- if(G) rownames(fit$G) else rownames(fit$H)
  clean_reg(fmla, vec)
  
}                                     

# H=================================================================================================================================================

clean_reg_names <- function(fit) {
  
  fmla <- fixed_form_rma(fit)
  vec <- rownames(fit$b)
  
  vec <- clean_reg(fmla, vec)
  if(fit$int.only) vec[vec=="(Intercept)"||vec==""] <- "Overall Effect"
  rownames(fit$b) <- vec
  rownames(fit$beta) <- vec
  rownames(fit$vb) <- colnames(fit$vb) <- vec
  return(fit)
}

# H=================================================================================================================================================                       

set_rownames_ <- function (object = nm, nm) 
{
  rownames(object) <- nm
  object
}       

# M=================================================================================================================================================

smooth_vi <- function(data, study, vi, digits = 8, fun = sd, ylab = "Studies", xlab = NULL, max_bar = FALSE, return_list = TRUE,
                      breaks = "Sturges", main = "Closeness of Sampling Variances\n (in each study)"){
  
  dot_cols <- rlang::ensyms(study, vi, fun)
  str_cols <- purrr::map_chr(dot_cols, rlang::as_string)
  
  d_clean <- full_clean(data)
  
  out <- map_dfr(group_split(d_clean, study), 
                 ~data.frame(study = unique(.[[str_cols[[1]]]]), sd = round(fun(.[[str_cols[[2]]]]),digits)))
  
  names(out)[2] <- str_cols[[3]]
  
  res <- hist(out[[2]], plot = FALSE, breaks = breaks)
  graphics.off()
  
  hist(out[[2]], breaks = breaks, xaxt = "n", yaxt = "n", cex.axis = .8, mgp = c(1.5,.4,0), xlab = if(is.null(xlab)) toupper(str_cols[[3]]) else xlab,
       main = main, cex.lab = .9, ylab = ylab, cex.main = .8)
  
  axis(1, at = res$breaks, las=1, cex.axis=.8, mgp=c(1.5,.4,0))
  axis(2, at = c(axTicks(2), max(res$counts)), las=1, cex.axis=.8, mgp=c(1.5,.55,0))
  text(res$mids[res$counts!=0], res$counts[res$counts!=0], res$counts[res$counts!=0], pos = 3, col = 4, font = 2, xpd = NA, cex = .75)
  
  whichbar <- findInterval(out[[2]], res$breaks)
  mm <- sort(unique(whichbar[!is.na(whichbar)]))
  
  out <- if(!max_bar) setNames(lapply(mm, function(i) na.omit(out[whichbar==i,])), 
                               paste0("Bar",seq_along(mm))) else na.omit(out[whichbar==which.max(res$counts),])
  
  return(invisible(if(return_list) out else dplyr::bind_rows(out, .id = "Bar")))
}

# M=================================================================================================================================================

post_rma <- function(fit, specs = NULL, cont_var = NULL, by = NULL,p_value = TRUE, ci = TRUE, mutos_vars_null = NULL, 
                     mutos_vars_contrast = NULL, block = FALSE, adjust = "none", compare = FALSE, plot_pairwise = FALSE, 
                     reverse = FALSE, digits = 3, xlab = "Estimated Effect", shift_up = NULL, shift_down = NULL, 
                     drop_rows = NULL, drop_cols = NULL, contrast_contrasts=FALSE, 
                     na.rm = TRUE, robust = FALSE, cluster, show0df = FALSE, sig = TRUE, contr, horiz = TRUE, at=NULL, 
                     at_vals=NA, get_rows = NULL, get_cols = NULL, df = NULL, tran = NULL, sigma=NULL, data=NULL, 
                     round_except=NULL, var=NULL,
                     ...)
{
  
  if(!inherits(fit, c("rma.uni", "rma.mv"))) stop("Model is not 'rma()/rma.mv()'.", call. = FALSE)
  
  dot_args <- list(...)
  dot_args_nm <- names(dot_args)
  
  cl <- match.call()
  
  if(!is.null(mutos_vars_null) & !is.null(mutos_vars_contrast)) { 
    
    message("Note: Only 'mutos_vars_null' results are shown.")
    mutos_vars_contrast <- NULL
    
  }
  
  if(!is.null(mutos_vars_contrast)) { 
    
    mutos_vars_null <- NULL 
    
    specs <- mutos_vars_contrast <- if(is_bare_formula(mutos_vars_contrast, lhs=FALSE)) 
      .all.vars(mutos_vars_contrast) else 
        if(is.character(mutos_vars_contrast)) mutos_vars_contrast
  }
  
  if(!is.null(mutos_vars_null)) { 
    
    mutos_vars_contrast <- NULL 
    
    specs <- mutos_vars_null <- if(is_bare_formula(mutos_vars_null, lhs=FALSE)) .all.vars(mutos_vars_null) else 
      if(is.character(mutos_vars_null)) mutos_vars_null
    
  }
  
  
  data_. <- if(is.null(data)) get_data_(fit) else data
  
  infer <- c(ci, p_value)
  
  if(robust) { 
    
    if(inherits(fit, "rma.mv")){
      
      fixed_eff <- is.null(fit$random)
      cr <- if(!fixed_eff) is_crossed(fit) else FALSE
      
      if(any(cr) || fixed_eff) { 
        
        robust <- FALSE
        message("Robust estimation not available for models with", if(any(cr))" crossed random-" else " only fixed-", "effects.")
      }
      
      vcov_. <- try(as.matrix(vcovCR(fit, type = "CR2", cluster = cluster)), silent=TRUE)
      
      if(inherits(vcov_., "try-error")) { 
        robust <- FALSE
        message("Robust tests unavailable (likely having <2 clusters for some moderators).\nResults are model-based.")
      }
      
    } else {
      
      if(missing(cluster)) { cluster <- 1:nrow(data_.) ;
      message("The missing 'cluster=' was set to corrospond to each row.")}
      
      vcov_. <- try(as.matrix(vcovCR(fit, type = "CR2", cluster = cluster)), silent=TRUE)
      
      if(inherits(vcov_., "try-error")) { 
        robust <- FALSE
        message("Robust tests unavailable. Results are model-based.")
      }
    }
  }    
  
  rma.mv_fit <- fit
  
  type. <- if('type' %in% dot_args_nm) dot_args$type else FALSE
  
  if(!is.null(at)) {
    
    at <- if(is_bare_formula(at, lhs=FALSE) || is.character(at)) { 
      
      lo_ave_up(at, data_., at_vals)
      
    } else at
    
  }
  
  df. <- if(is.null(df)) {
    
    df_detect(rma.mv_fit)
    
  } else df
  
  tran. <- if(is.null(tran)) {
    
    tran_detect(rma.mv_fit)
    
  } else { tran }
  
  
  sigma. <- if (is.null(sigma)) {
    
    sigma_detect(rma.mv_fit)
    
  } else {
    
    sigma
    
  }
  
  lookup <- c(Contrast="contrast",Estimate="estimate",Mean="emmean",Response="response",t="t.ratio",
              Df="df","p-value"="p.value",Lower="lower.CL",Upper="upper.CL",
              Df1="df1", Df2="df2","F"="F.ratio",Term="model term",
              Lower="asymp.LCL", Upper="asymp.UCL", z="z.ratio", Ratio = "ratio")
  
  if(!is.null(mutos_vars_contrast) & block || !is.null(mutos_vars_null) & block) names(lookup)[13] <- "Block Term"
  if(!is.null(mutos_vars_contrast) & !block || !is.null(mutos_vars_null) & !block) names(lookup)[13] <- "(M)UTOS Term"
  
  fit <- rma2gls(fit)
  
  if(!is.null(specs) & !is.character(specs) & !is_bare_formula(specs)) {stop("The 'specs' must be either a character or a formula containing '~'.", call.=FALSE)}
  
  if(is.null(specs)) specs <- as.formula(bquote(~.(terms(fit)[[3]])))
  
  
  if(robust) { 
    
    rownames(vcov_.)[rownames(vcov_.) %in% "intrcpt"] <- "(Intercept)"
    colnames(vcov_.)[colnames(vcov_.) %in% "intrcpt"] <- "(Intercept)"
    
    fit$varBeta <- vcov_. 
  }
  
  
  is_contr <- !missing(contr)            
  
  ems <- suppressWarnings(suppressMessages(try(if(is.null(cont_var) & is.null(var)){
    
    if(!isFALSE(tran.)) emmeans(object = fit, specs = specs, infer = infer, adjust = adjust, contr = contr, data = data_., tran = tran., sigma = sigma., df = df., at=at, ...)
    else emmeans(object = fit, specs = specs, infer = infer, adjust = adjust, contr = contr, data = data_., sigma = sigma., df = df., at=at, ...)
 
     } else {
       
       cont_var <- if(!is.null(cont_var)) cont_var else var
       
    if(!is_contr){ 
      
      emtrends(object = fit, specs = specs, var = cont_var, infer = infer, adjust = adjust, data = data_., tran = tran., sigma = sigma., df = df., at=at, ...)
      
    } else {
      
      emtrends(object = fit, specs = specs, var = cont_var, infer = infer, adjust = adjust, contr = contr, data = data_., tran = tran., sigma = sigma., df = df., at=at, ...)
      
    }
    
  }, silent = TRUE)))
  
  
  if(inherits(ems,"try-error")) stop("Wrong specification OR no relavant data found for the comparisons.", call.=FALSE)
  
  
  con_methods <- c("pairwise","revpairwise","tukey","consec",
                   "poly","trt.vs.ctrl","trt.vs.ctrlk","trt.vs.ctrl1",
                   "dunnett","mean_chg","eff","del.eff","identity")
  
  is_pair <- any(con_methods %in% as.character(specs))
  
  if(!is.null(cont_var) & is_pair) names(lookup)[2] <- paste0(cont_var,".dif")
  
  out <- if(is_pair){
    
    methd <- as.character(specs[2])
    
    ems <- ems[[if(!contrast_contrasts) 1 else 2]]
    
    if(plot_pairwise) print(plot(ems, by = by, comparisons = compare, horizontal = horiz, adjust = adjust, xlab = xlab)) 
    
    pp <- if(isFALSE(tran.)) contrast(ems, method = methd, each="simple", infer=infer, reverse=reverse, adjust=adjust,...) 
    else contrast(ems, method = methd, each="simple", infer=infer, reverse=reverse, adjust=adjust,tran = tran.,...)
    
    if(!is_contr) pp else if(!isFALSE(tran.)) contrast(pp, contr, tran = tran.,...) else contrast(pp, contr, ...)
    
  }
  
  else {
    
    if(!is.null(mutos_vars_contrast)) {
      
      mutos_vars_contrast <- if(is_bare_formula(mutos_vars_contrast, lhs=FALSE)) .all.vars(mutos_vars_contrast) else 
        if(is.character(mutos_vars_contrast)) mutos_vars_contrast      
      
      if(!block){
        
        message("Testing jointly if EMMs for levels of each categorical moderator are equal to each other.")
        
        joint_tests(ems, by = by, adjust = adjust, show0df = show0df, tran = tran., ...)
        
      } else {
        
        if(length(mutos_vars_contrast) < 2) stop("A block needs at least two categorical moderators.", call. = FALSE)
        
        com <- comb_facs(ems, mutos_vars_contrast)
        
        message("Testing jointly if the EMMs *across* multiple
categorical moderators (a block of them) are equal to each other.")
        
        joint_tests(com)
        
      }
      
      
    } else if (!is.null(mutos_vars_null)){
      
      is_fm <- is_bare_formula(mutos_vars_null, lhs=FALSE)    
      
      mutos_vars_null <- if(is_fm) .all.vars(mutos_vars_null) else 
        if(is.character(mutos_vars_null)) mutos_vars_null      
      
      if(!block){  
        
        zz <- cbind(mod=mutos_vars_null, as.data.frame(map_dfr(mutos_vars_null,~emmeans::test(emmeans(ems,.),joint=TRUE))))
        names(zz)[1] <- "(M)UTOS Term"
        message("Testing jointly if EMMs for levels of each categorical moderator are equal to their null (e.g., 0).")
        
        zz 
        
      } else {
        
        if(length(mutos_vars_null) < 2) stop("A block needs at least two categorical moderators.", call. = FALSE)
        
        zz <- cbind(mod=paste0(mutos_vars_null, collapse="."), as.data.frame(emmeans::test(ems, joint=TRUE)))
        names(zz)[1] <- "Block Term"
        message("Testing jointly if the EMMs *across* multiple
categorical moderators (a block of them) are equal to their null (e.g., 0).")
        
        zz
        
      }     
    }
    
    else {
      
      ems
    }
  }
  
  out <- as.data.frame(out, adjust = adjust, infer = infer, tran = tran., ...) %>%
    dplyr::rename(tidyselect::any_of(lookup)) %>% 
    dplyr::select(-tidyselect::any_of("note"))
  
  
  out <- set_rownames_(out,NULL)
  
  if(p_value){
    
    p.values <- as.numeric(out$"p-value")
    
    if(all(is.na(p.values))) { 
      stop("Comparison(s)/moderator adjustments are non-estimable,\nlikely some combination of moderating or control variables are missing.\nTake those variables out of the model one by one and re-run.",
           call. = FALSE)
    }
    
    if(sig){
      Signif <- symnum(p.values, corr = FALSE, 
                       na = FALSE, cutpoints = 
                         c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       symbols = c("***", "**", "*", ".", " "))
      
      out <- tibble::add_column(out, Sig. = Signif, .after = "p-value")
    }
  }  
  
  out0 <- roundi(out, digits = digits, except=round_except)
  if(na.rm) out <- na.omit(out)
  
  out <- roundi(out, digits = digits, except = round_except)
  
  if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
  if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
  if(!is.null(drop_rows)) out <- out[-drop_rows, ]
  if(!is.null(get_rows)) out <- out[get_rows, ]
  
  if(!is.null(drop_cols)) out <- dplyr::select(out, -tidyselect::all_of(drop_cols))
  if(!is.null(get_cols)) out <- dplyr::select(out, tidyselect::all_of(get_cols))
  
  out <- list(table = out, table0 = out0, specs = specs, call = cl, fit = fit, rma.mv_fit = rma.mv_fit, ems = ems,
              tran. = tran., type. = type., df. = df., sigma. = sigma., digits = digits, is_contr = is_contr)
  
  class(out) <- "post_rma"
  return(out)
}                                 

# M=================================================================================================================================================

R2_rma <- function(..., robust = TRUE, digits = 3, 
                   model_names = NULL, null_model = NULL,
                   level_names = NULL, blank_sign = "", 
                   null_name = "No (M)UTOS", tol_large = 1e4)
{
  
  LL <- list(...)
  if(!all(sapply(LL,inherits,"rma.mv"))) stop("All models must be 'rma.mv()'.", call. = FALSE)
  
  bad <- sapply(LL, function(i) i$withG || i$withH || i$withR || is.null(i$random))
  
  if(any(bad)) stop("Model not supported.", call. = FALSE)
  
  ok <- length(unique(map(LL,~map_chr(strsplit(.$s.names,"/",fixed=TRUE),tail,1))))==1
  
  if(!ok) stop("Models must have the exact same random-effects.", call.=FALSE)
  
  first <- LL[[1]]
  
  Model <- if(is.null(model_names)) as.character(substitute(...())) else model_names
  
  .zolqui_. <- as.formula(paste0(as.character(fixed_form_rma(first))[2],"~1"))
  
  null_fit <- if(is.null(null_model)) update.rma(first, yi = .zolqui_.) else null_model
  
  lvl_names <- if(is.null(level_names)) sapply(strsplit(null_fit$s.names,"/",fixed=TRUE),tail,1) else level_names
  
  sigmasn <- setNames(sqrt(null_fit$sigma2), lvl_names)
  
  sigma_totaln <- sqrt(sum(sigmasn^2))
  
  null_res <- data.frame(Model = null_name,._A_.= sigma_totaln,._D_.=NA,R2=NA)
  
  null_res <- add_column(null_res, as.data.frame(t(sigmasn)), .after = "._A_.")
  
  z <- function(nm) paste0("Sigma(",nm,")") 
  
  #  on.exit(Sys.setlocale("LC_ALL"))               
  #  Sys.setlocale(locale = "Greek")
  
  f <- function(fit){
    
    if(robust){
      
      mc <- try(clubSandwich::Wald_test(fit, constrain_zero(fit$btt), "CR2"), silent=TRUE)   
      
      bad <- inherits(mc,"try-error")
      
      if(bad || !bad && is.na(mc$p_val)) { 
        robust <- FALSE
        message("Note: Robust QM unavailable,likely: \n1- Some moderators in <2 clusters OR/AND \n2- High # of coefficients vs. # of highest clusters.\n3- Some combination of variables in the interactive model are missing.\nQM results are model-based.\n")
      }
      
      if(!bad && mc$Fstat>tol_large) { message("Note: Robust estimation seems unfit for the model (use 'robust=FALSE').") } 
    }
    
    p <- if(robust) mc$p_val else fit$QMp
    
    sigmas <- setNames(sqrt(fit$sigma2), lvl_names) 
    
    sigma_total <- sqrt(sum(sigmas^2)) 
    
    R2 <- (sigma_totaln - sigma_total) / sigma_totaln*1e2
    
    if(R2<0) message("Negative R2 was set to 0.")
    
    R2 <- if(R2>0) R2 else 0
    
    c(._A_.=sigma_total,sigmas,._D_.=p, R2=R2)
  }
  
  out <- map_dfr(LL, f) %>% as.data.frame() %>%
    add_column(Model = Model, .before = "._A_.")
  
  res <- roundi(bind_rows(null_res, out), digits = digits)
  
  res[-1,]$R2 <- paste0(res[-1,]$R2,"%")
  
  names(res)[names(res) %in% c("._A_.",lvl_names,"._D_.")] <- c(z(c("total",lvl_names)),"p-value")
  
  blk <- paste0(paste0(rep(" ",digits-1), collapse=""), "NA", collapse ="")
  
  res[res == blk] <- blank_sign
  
  return(res)
}                                                                         


# H=================================================================================================================================================

sizetree2 <- function (x, left = 0, top, right = 1, lastcenter = NA,
                       showval = TRUE, showcount = FALSE, stacklabels = TRUE, firstcall = TRUE,
                       col = NULL, border = NA, toplab = NULL, base.cex = 1, cex_top = 1, ...) 
{
  dimx <- dim(x)
  colname <- names(x)[1]
  if (firstcall) {
    x <- x[do.call(order, x), ]
    oldmar <- par("mar")
    par(mar = c(1, 2, 2, 1))
    top <- sum(!is.na(x[, 1]))
    if (top < dimx[1]) 
      cat(dimx[1] - top, "NA values dropped from first stack.\n")
    plot(0, xlim = c(0, dimx[2]), ylim = c(0, top), type = "n", 
         axes = FALSE, xlab = "", ylab = "", ...)
  }
  xfreq <- table(x[, 1])
  lenxf <- length(xfreq)
  if (firstcall) {
    if (is.null(col)) {
      col <- list()
      for (index in 1:dimx[2]) col[[index]] <- rainbow(length(table(x[, 
                                                                      index])))
    }
    for (index in 1:dimx[2]) if (is.null(names(col[[index]]))) 
      names(col[[index]]) <- names(table(x[, index]))
  }
  if (lenxf) {
    if (is.list(col)) {
      barcol <- col[[1]]
      barcol <- barcol[names(col[[1]]) %in% names(xfreq)]
    }
    else barcol <- col[names(col) %in% names(xfreq)]
    labels <- names(xfreq)
    squeeze <- (right - left)/10
    for (bar in 1:lenxf) {
      if (length(xfreq[bar])) {
        if (!is.na(xfreq[bar])) {
          if (xfreq[bar] > 0) {
            rect(left + squeeze, top - xfreq[bar], right - 
                   squeeze, top, col = barcol[bar], border = border)
            labelheight <- strheight(labels[bar])
            cex <- ifelse((1.5 * labelheight) > xfreq[bar], 
                          base.cex * 0.75 * xfreq[bar]/labelheight, 
                          base.cex)
            if (showval) {
              textcol <- ifelse(colSums(col2rgb(unlist(barcol[bar])) * 
                                          c(1.4, 1.4, 0.5)) < 350, "white", "black")
              bartext <- ifelse(showcount, paste(labels[bar], 
                                                 " (", xfreq[bar], ")", sep = ""), labels[bar])
              text((left + right)/2, top - xfreq[bar]/2, 
                   bartext, cex = cex, col = textcol)
            }
            if (!is.na(lastcenter)) 
              segments(left + squeeze, top - xfreq[bar]/2, 
                       left - squeeze, lastcenter)
            xvalue <- ifelse(is.numeric(x[, 1]), as.numeric(labels[bar]), 
                             labels[bar])
            if (dimx[2] > 1) {
              newcol <- col
              newcol[[1]] <- NULL
              nextx <- subset(x, x[, 1] == xvalue, 2:dimx[2])
              sizetree2(nextx, right, top, right + 1, 
                        lastcenter = top - xfreq[bar]/2, showval = showval,
                        showcount = showcount, stacklabels = stacklabels, firstcall = FALSE,
                        col = newcol, border = border, base.cex = base.cex)
            }
          }
        }
      }
      top <- top - xfreq[bar]
    }
    if (stacklabels) 
      mtext(colname, side = 1, at = (left + right)/2, line = -0.4)
  }
  if (firstcall) {
    if (!is.null(toplab)) {
      par(xpd = TRUE)
      top <- sum(!is.na(x[, 1]))
      text(seq(0.5,dimx[2] - 0.5,by=1), 1.01 * top, toplab, adj = c(0.5,
                                                                    0), cex = cex_top)
      par(xpd = FALSE)
    }
    par(mar = oldmar)
  }
}

#================================================================================================================================================

contr_rma <- function(post_rma_fit, contr_index){
  
  if(!inherits(post_rma_fit, "post_rma")) stop("post_rma_fit is not 'post_rma()'.", call. = FALSE)
  
  is_contr <- post_rma_fit$is_contr
  
  post_rma_fit <- if(!is_contr) post_rma_fit$table0 else filter(post_rma_fit$table0, Contrast == ".")
  
  ind <- rep(0, nrow(post_rma_fit))
  
  ind[abs(contr_index)] <- sign(contr_index)
  
  return(ind)
}                
#===================================================================================================================================================
## From clubSandwich package kept as the function is get depricated.
                
unblock <- function (A, block = attr(A, "groups")) 
{
  if (is.null(block)) 
    block <- factor(rep(names(A), times = sapply(A, function(x) dim(x)[1])))
  n <- length(block)
  mat <- matrix(0, n, n)
  for (i in levels(block)) {
    index <- i == block
    mat[index, index] <- A[[i]]
  }
  return(mat)
}

isPosDef <- function (x) 
{
  x_na <- is.na(x)
  mis_rows <- apply(x_na, 1, all)
  mis_cols <- apply(x_na, 2, all)
  if (all(mis_rows) | all(mis_cols)) 
    return(TRUE)
  x_nomiss <- x[!mis_rows, !mis_cols]
  x_eig <- eigen(x_nomiss)
  all(x_eig$values > 0)
}


check_PD <- function (vcov_list) 
  {
    PD <- sapply(vcov_list, isPosDef)
    if (!all(PD)) {
      NPD_clusters <- names(vcov_list)[!PD]
      warn_text <- paste(c("The following clusters have non-positive definite covariance matrices:", 
                           NPD_clusters), collapse = "\n")
      warning(warn_text)
    }
    else {
      NULL
    }
  }


impute_covariance_matrix <- function (vi, cluster, r, ti, ar1, smooth_vi = FALSE, subgroup = NULL, 
            return_list = identical(as.factor(cluster), sort(as.factor(cluster))), 
            check_PD = TRUE) 
  {
    cluster <- droplevels(as.factor(cluster))
    vi_list <- split(vi, cluster)
    
    if (smooth_vi) 
      vi_list <- lapply(vi_list, function(x) rep(mean(x, na.rm = TRUE), 
                                                 length(x)))
    if (missing(r) & missing(ar1)) 
      stop("You must specify a value for r or for ar1.")
    if (!missing(r)) {
      r_list <- rep_len(r, length(vi_list))
      if (missing(ar1)) {
        vcov_list <- Map(function(V, rho) (rho + diag(1 - 
                                                        rho, nrow = length(V))) * tcrossprod(sqrt(V)), 
                         V = vi_list, rho = r_list)
      }
    }
    if (!missing(ar1)) {
      if (missing(ti)) 
        stop("If you specify a value for ar1, you must provide a vector for ti.")
      ti_list <- split(ti, cluster)
      ar_list <- rep_len(ar1, length(vi_list))
      if (missing(r)) {
        vcov_list <- Map(function(V, time, phi) (phi^as.matrix(stats::dist(time))) * 
                           tcrossprod(sqrt(V)), V = vi_list, time = ti_list, 
                         phi = ar_list)
      }
      else {
        vcov_list <- Map(function(V, rho, time, phi) (rho + 
                                                        (1 - rho) * phi^as.matrix(stats::dist(time))) * 
                           tcrossprod(sqrt(V)), V = vi_list, rho = r_list, 
                         time = ti_list, phi = ar_list)
      }
      vcov_list <- lapply(vcov_list, function(x) {
        attr(x, "dimnames") <- NULL
        x
      })
    }
    if (!is.null(subgroup)) {
      si_list <- split(subgroup, cluster)
      subgroup_list <- lapply(si_list, function(x) sapply(x, 
                                                          function(y) y == x))
      vcov_list <- Map(function(V, S) V * S, V = vcov_list, 
                       S = subgroup_list)
    }
    if (check_PD) 
      check_PD(vcov_list)
    if (return_list) {
      return(vcov_list)
    }
    else {
      vcov_mat <- unblock(vcov_list)
      cluster_index <- order(order(cluster))
      return(vcov_mat[cluster_index, cluster_index])
    }
  }
                
#===================================================================================================================================================

prob_rma <- function(post_rma_fit, target_effect = 0, condition = c("or larger", "or smaller"), 
                     gain = FALSE, ..., sep = get_emm_option("sep")){
  
  Term <- term_names_(post_rma_fit=post_rma_fit, sep=sep, na.rm = TRUE)
  
  fit <- post_rma_fit$rma.mv_fit
  
  if(fit$withG || fit$withH || fit$withR || is.null(fit$random)) stop("Model not supported.", call. = FALSE)
  
  specs <- post_rma_fit$specs
  
  digits <- post_rma_fit$digits
  
  ems <- if(inherits(post_rma_fit, "post_rma")) post_rma_fit$ems else 
    post_rma_fit$con
  
  post_rma_fit <- type.convert(post_rma_fit$table, as.is=TRUE)
  
  cond <- match.arg(condition)  
  
  lower.tail <- switch(cond, 
                       "or larger" = FALSE, 
                       "or smaller" = TRUE)
  
  lower_ave_eff <- post_rma_fit$Lower
  
  upper_ave_eff <- post_rma_fit$Upper
  
  ave_eff <- na.omit(round(predict(ems), digits = digits))
  
  ci <- as.data.frame(confint.rma.mv(fit, ...))
  
  all_lvls <- ci[odds_.(seq_len(nrow(ci))), , drop=FALSE]
  
  #total_sd: estimate, lower, upper
  total_sd <- if(!gain) sqrt(colSums(all_lvls)) else sqrt(2)*sqrt(ci[nrow(ci),]) #sqrt(2) * sqrt(colSums(all_lvls[-1, ,drop=FALSE]))
  
  # Probability at the estimates
  Probability <- paste0(formatC(round(pnorm(target_effect, ave_eff, total_sd[1], lower.tail=lower.tail), 4)*1e2,digits = 2, format = "f"),"%")
  
  # Probability at the lowest ave_eff and highest variability
  min_Probability <- paste0(formatC(round(pnorm(target_effect, lower_ave_eff, total_sd[3], lower.tail=lower.tail), 4)*1e2,digits = 2, format = "f"),"%")
  
  # Probability at the highest ave_eff and lowest variability
  max_Probability <- paste0(formatC(round(pnorm(target_effect, upper_ave_eff, total_sd[2], lower.tail=lower.tail), 4)*1e2, digits = 2, format = "f"),"%")
  
  data.frame(Term=Term, Target_Effect = paste(target_effect, cond, collapse = " "), Probability = Probability, 
             Min = min_Probability, Max = max_Probability)
}

#M==============================================================================================================================================

sense_rma <- function(post_rma_fit = NULL, var_name, fit = NULL, 
                      r = (3:7)*.1, cluster = NULL, clean_names = NULL,
                      regression = NULL, label_lines = TRUE,
                      cex_labels = .55, plot_coef = TRUE, plot_hetro = TRUE, digits = 3, 
                      ..., sep = get_emm_option("sep")){
  
  
  if(is.null(fit) & is.null(post_rma_fit)) stop("Provide either 'fit=' or 'post_rma_fit='.", call. = TRUE)
  if(!is.null(fit) & !inherits(fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  if(is.null(fit)) fit <- post_rma_fit$rma.mv_fit
  if(!is.null(post_rma_fit) & !inherits(post_rma_fit, "post_rma")) stop("post_rma_fit is not 'post_rma()'.", call. = FALSE) 
  if(fit$withG || fit$withH || fit$withR || is.null(fit$random)) stop("Model not supported.", call. = FALSE)
  
  
  dat <- get_data_(fit)
  
  regression <- if(is.null(regression) & !is.null(post_rma_fit)) {
    
    if("cont_var" %in% as.character(post_rma_fit$call)) TRUE else FALSE
    
  } else if(is.null(regression) & is.null(post_rma_fit)) TRUE else regression 
  
  
  if(regression){
    
    lm_fit <- lm(fixed_form_rma(fit), data = dat, na.action = "na.omit")
    
    cl <- clean_reg(lm_fit, names(coef(lm_fit)))
    
    if(is.null(clean_names)){
      
      if(any_num_vec(cl)) {
        
        clean_names <- FALSE
        
      } else { 
        
        clean_names <- TRUE
        
      }
    } 
    
    if(clean_names) names(lm_fit$coefficients) <- cl
    if(clean_names) fit <- clean_reg_names(fit)
  }   
  
  if(!is.null(post_rma_fit)){
    
    tran. <- post_rma_fit$tran.
    type. <- post_rma_fit$type.
    specs <- post_rma_fit$specs
    ems <- post_rma_fit$ems
    Term <- term_names_(post_rma_fit=post_rma_fit, sep=sep, na.rm=FALSE)      
  }
  
  cluster_name <- if(is.null(cluster)) strsplit(fit$s.names,"/",fixed=TRUE)[[1]] else cluster
  
  V_list <- lapply(r, function(i) impute_covariance_matrix(vi = dat[[var_name]], cluster=dat[[cluster_name]], r=i))
  
  model_list <- lapply(V_list, function(i) suppressWarnings(update.rma(fit, V = i)))
  
  xaxis_lab <- paste0("r=",r)
  
  total_hetros <- sapply(model_list, function(i) sqrt(sum(i$sigma2)))
  

  if(plot_coef & plot_hetro){
    
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  
    par(mfrow = c(2,1), mgp = c(1.5, 0.14, 0), mar = c(1.5, 2.6, 1.5, .5), 
        tck = -0.02)
    
  }  
  
  output <- if(regression){
    
    fixed_eff_list <- lapply(model_list, function(i) setNames(coef(summary(i))$estimate, rownames(fit$b)))
    
    if(plot_coef){  
      
      matplot(t(as.data.frame(fixed_eff_list)), type = "l", xaxt = "n", ylab = "Estimates", ...)
      
      axis(1, at = axTicks(1), labels = xaxis_lab,...)
      
      mn <- mean(seq_len(length(fixed_eff_list)))
      
      if(label_lines) text(mn, as.data.frame(fixed_eff_list)[,mn], rownames(fit$b),
                           cex = cex_labels)
    }
    
    output <- as.data.frame(t(do.call(rbind, fixed_eff_list)))
    
    setNames(output, xaxis_lab)
    
  } else {
    
    if(!is.null(post_rma_fit)){
      
      post_rma_list <- lapply(model_list, function(i) 
        setNames(predict(post_rma(i, specs, tran = tran., type = type.)$ems),Term))
      
    } else {
      
      stop("Please provide a 'post_rma_fit' or use 'regression=TRUE'.", call. = FALSE)
      
    }
    
    if(plot_coef){    
      
      matplot(t(as.data.frame(post_rma_list)), type = 'l', xaxt = "n", ylab = "Mean Effect", xlab = NA,...)
      
      axis(1, at = axTicks(1), labels = xaxis_lab, ...)
      
      mn <- mean(seq_len(length(post_rma_list)))
      
      if(label_lines) text(mn, as.data.frame(post_rma_list)[,mn], Term,
                           cex = cex_labels)
    }    
    output <- as.data.frame(t(do.call(rbind, post_rma_list)))
    
    setNames(output, xaxis_lab)
    
  }
  
  if(plot_hetro){
    
    rng <- range(total_hetros)
    mrng <- mean(rng)
    plot(total_hetros, type = "l", ylim = rng+c(-mrng,mrng), xaxt = "n", xlab = NA, ylab = "Total Variation (SD)",...)
    axis(1, at = axTicks(1), labels = xaxis_lab, ...)
  }  
  
  out <- na.omit(rbind(output, Total_variation_in_SD = total_hetros))
    
  out <- cbind(out, sd =apply(out, 1, sd), change_rate =apply(out, 1, function(i) coef(lm(i~r))[2]))

  roundi(rownames_to_column(out, "Term"), digits = digits)
} 
                                              

#M================================================================================================================================================

con_rma <- function(post_rma_fit, method, type,
                    digits = 3, ci = TRUE, 
                    p_value = TRUE, adjust = "none",
                    na.rm = TRUE, sig = TRUE, round_except=NULL, ...){
  
  if(!inherits(post_rma_fit, "post_rma")) stop("post_rma_fit is not 'post_rma()'.", call. = FALSE)
  
  infer <- c(ci, p_value)
  
  lookup <- c(Contrast="contrast",Estimate="estimate",Mean="emmean",Response="response",t="t.ratio",
              Df="df","p-value"="p.value",Lower="lower.CL",Upper="upper.CL",
              Df1="df1", Df2="df2","F"="F.ratio",Term="model term",
              Lower="asymp.LCL", Upper="asymp.UCL", z="z.ratio", Ratio = "ratio")

  tran. <- post_rma_fit$tran.
  
  con <- if(!isFALSE(tran.)) contrast(post_rma_fit$ems, method = method, type = type, infer = infer, tran = tran., ...)
  else contrast(post_rma_fit$ems, method = method, infer = infer, ...)
  
  out <- as.data.frame(con, adjust = adjust, infer=infer, ...) %>% 
    dplyr::rename(tidyselect::any_of(lookup)) %>% 
    dplyr::select(-tidyselect::any_of("note"))
  
  if(p_value){
    p.values <- as.numeric(out$"p-value")
    
    if(all(is.na(p.values))) { 
      return(message("Error: Comparison(s) are non-estimable,\nlikely some combination of moderating or control variables are missing.\nTake variables out of model one by one and re-run."))
    }
    
    if(sig){
      Signif <- symnum(p.values, corr = FALSE, 
                       na = FALSE, cutpoints = 
                         c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       symbols = c("***", "**", "*", ".", " "))
      
      out <- tibble::add_column(out, Sig. = Signif, .after = "p-value")
    }
  } 

  out0 <- roundi(out, digits = digits, except=round_except)
  if(na.rm) out <- na.omit(out)
  
  out <- roundi(out, digits = digits, except=round_except)
  
  out <- list(table = out, table0 = out0, specs = post_rma_fit$specs, call = post_rma_fit$call, fit = post_rma_fit$fit, rma.mv_fit = post_rma_fit$rma.mv_fit, ems = post_rma_fit$ems,
              tran. = post_rma_fit$tran., type. = post_rma_fit$type., con = con, digits = digits)
  
  class(out) <- "contrast_rma"
  
  return(out)
}      


# H================================================================================================================================================

add_signs <- function(post_rma_fit, con_index, sep = get_emm_option("sep")) 
{
  
  term_names <- term_names_(post_rma_fit=post_rma_fit, sep=sep)
  
  uniq <- unique(con_index)
  ind_abs <- abs(con_index)
  uniq_abs <- unique(abs(uniq))
  same_val <- length(uniq_abs)==1
  dup <- duplicated(ind_abs)
  rps <- any(dup)
  last_sgn <- if(rps) tail(sign(con_index[dup|duplicated(ind_abs, fromLast=TRUE)]), 1) else FALSE
  
  if(same_val & last_sgn > 0) return({ 
    message("Note: No contrast made for ", toString(dQuote(unique(term_names[ind_abs[dup]]))) ,", just returned from the inputted EMM table.\n")
    term_names[uniq_abs] }) 
  
  if(same_val & last_sgn < 0) return({ 
    message("Note: No contrast made for ", toString(dQuote(unique(term_names[ind_abs[dup]]))), ", just returned with opposite sign from the inputted EMM table.\n")
    paste("-", term_names[uniq_abs]) })
  
  if(rps & last_sgn > 0) { 
    message("Note: No contrast made for ", toString(dQuote(unique(term_names[ind_abs[dup]]))) ,", just used from the inputted EMM table.\n")
  } 
  
  if(rps & last_sgn < 0) { 
    message("Note: No contrast made for ", toString(dQuote(unique(term_names[ind_abs[dup]]))), ", just used with opposite sign from the inputted EMM table.\n")
  }
  
  merged <- paste(ifelse(con_index < 0, '-', '+'), 
                  term_names[abs(con_index)], 
                  collapse = sep)
  sub('^\\+ ', '', merged)
}              
                                                              
# H===============================================================================================================================================
                                                              
#not.integer_ <- function(x) { 
  
#  f <- function(y) (abs((y) - floor((y) + .5)) > 1e-7)
  
#  if(is.list(x)) lapply(x, f) else f(x)
#}
                                                              
# M================================================================================================================================================       
                              
contrast_rma <- function(post_rma_fit, con_list, ..., 
                         auto_name=FALSE, sep=get_emm_option("sep"),
                         pretest_name = "baseline", brief=FALSE,
                         posttest_base_name = "posttest", gain_dif = FALSE, 
                         gain_dif_type = c("all","same","different"),
                         emm_style=FALSE)
{
  
  
 # if(!missing(con_list) & !emm_style) {
    
 #   emm_style <- any(sapply(not.integer_(con_list), any))
    
 # }
  
  if(!emm_style){  
    
    if(missing(con_list)) {
      
      con_list <- con_gain_list(post_rma_fit, pretest_name, 
                                brief, posttest_base_name, gain_dif, 
                                gain_dif_type, sep)$con_list
      
      
    }
    
    if(is.numeric(con_list)) con_list <- list(con_list)
    
    con_methods <- c("pairwise","revpairwise","tukey","consec",
                     "poly","trt.vs.ctrl","trt.vs.ctrlk","trt.vs.ctrl1",
                     "dunnett","mean_chg","eff","del.eff","identity")
    
    if(is.character(con_list) && any(!con_list %in% con_methods))
      stop("If not a list, 'con_list' can be one of: ", toString(dQuote(con_methods)),
           ".", call. = FALSE)
    
    con_indx <- if(is.list(con_list)) {
      
      con_list <- if(auto_name) {
        
        sapply(seq_along(con_list), function(i) {
          
          nm <- if(is.null(names(con_list)[i]) || names(con_list)[i] == "") {
            
            add_signs(post_rma_fit, con_list[[i]], sep=sep)
          }
          else { names(con_list)[i] }
          
          return(setNames(con_list[i], nm))
        })
        
      } else
      { con_list }
      
      lapply(con_list, contr_rma, post_rma_fit = post_rma_fit)
      
    } else
    { con_list }
    
    con_rma(post_rma_fit, con_indx, ...)
    
  } else {
    
    con_rma(post_rma_fit, con_list, ...)
    
  }
}

# M======================================================================================================================================================  

plot_rma <- function(fit, formula, ylab, CIs=TRUE, PIs=FALSE, 
                     CIarg = list(lwd = .5, alpha = 1),
                     PIarg = list(lwd = 1.25, alpha = 0.33),
                     linearg = list(lwd = .5, linetype = "solid"),
                     cov.reduce = TRUE, tran = NULL, 
                     sigma = NULL, df = NULL, at = NULL, 
                     at_vals = NA, interpolate_length = 150,
                     dodge=.15, ...){
  
  if(!inherits(fit, c("post_rma", "rma.mv", "rma.uni"))) stop("fit is not 'post_rma()','rma.mv()' or 'rma.uni()'.", call. = FALSE)
  
  dots <- list(...)
  dot_nms <- names(dots)
  cl <- match.call()
  
  no_ciarg <- is.null(cl$CIarg)
  no_piarg <- is.null(cl$PIarg)
  no_linearg <- is.null(cl$linearg)
  
  is_post_rma <- inherits(fit, "post_rma")
  
  df. <- if(is.null(df)) {
    if(!is_post_rma) df_detect(fit) else fit$df.
  } else { df }
  
  
  data_. <- if(is_post_rma) get_data_(fit$rma.mv_fit) else get_data_(fit)
  
  
  if(!is.null(at)) {
    
    at <- if(is_bare_formula(at, lhs=FALSE) || is.character(at)) { 
      
      lo_ave_up(at, data_., at_vals)
      
    } else at
    
  }
  
  
  tran. <- if(is.null(tran)) {
    
    if(!is_post_rma) tran_detect(fit) else fit$tran.
    
  } else { tran }
  
  
  sigma. <-  if (is.null(sigma)) {
    
    if(!is_post_rma) sigma_detect(fit) else fit$sigma.
    
  } else {
    
    sigma
    
  }
  
  is_var1 <- is_post_rma & "cont_var" %in% names(as.list(fit$call)) || is_post_rma & "var" %in% names(as.list(fit$call))
 
   if(is_var1) { 
    
    var <- if("cont_var" %in% names(as.list(fit$call))) fit$call$cont_var else fit$call$var
    
    if(length(var)>1) stop("Only one continous variable can be supplied to 'var='/'cont_var='.", call. = FALSE)
    
    is_post_rma <- FALSE
    fit <- fit$rma.mv_fit
    num_var <- data_.[[var]]
    
    at_add <- setNames(list(seq(min(num_var, na.rm=TRUE),max(num_var, na.rm=TRUE),length.out=interpolate_length)), var)
    
    at <- c(at, at_add)
  }
  
    
  is_var2 <- any(c("cont_var","var") %in% dot_nms)
  
  if(is_var2) {
    
    var <- if(!is.null(cl$cont_var)) cl$cont_var else cl$var
      
    var <- if(is.character(var)) var else 
      stop("'var='/'cont_var=' must be a single-element character vector (ex. 'year').", call. = FALSE)
    
    if(length(var)>1) stop("Only one continous variable can be supplied to 'var='/'cont_var='.", call. = FALSE)
    
    num_var <- data_.[[var]]
    
    is_num_var <- is.numeric(num_var)
    
    at_add <- if(is_num_var)  
      setNames(list(seq(min(num_var, na.rm=TRUE),max(num_var, na.rm=TRUE),length.out=interpolate_length)), var)
    else stop("'var='/'cont_var=' must be a continous variable.")
    
    at <- c(at, at_add)
  } 
  
  
  if(missing(ylab)) ylab <- paste0("Effect Size (",as.character(fixed_form_rma(if(is_post_rma) fit$rma.mv_fit else fit))[2],")")
  
  is_var <- any(is_var1, is_var2)
  
  if(no_ciarg & CIs & is_var) CIarg <- list(lwd = 2, alpha = 0.08)
  if(no_piarg & PIs & is_var) PIarg <- list(lwd = 5, alpha = 0.08)
  if(no_linearg & CIs & is_var || no_linearg & PIs & is_var) linearg <- list(size=1, alpha=1, linetype="solid")
   
  fit <- if(!is_post_rma) { 
    
    ref_grid(rma2gls(fit), cov.reduce = cov.reduce, df = df., sigma = sigma., tran = tran., at = at, ...)
    
  } else fit
  
  emmip(object=if(is_post_rma) fit$ems else fit, formula=formula, 
        ylab=ylab, CIs=CIs, PIs=PIs, CIarg=CIarg, PIarg=PIarg,
        linearg=linearg, cov.reduce=cov.reduce, tran=tran., at=at, 
        dodge=dodge, ...)
  
}
                                
# M===============================================================================================================================================

pct_dif_tran <- list(
  linkfun = function(mu) log(mu/100 + 1),
  linkinv = function(eta) 100 * (exp(eta) - 1),
  mu.eta = function(eta) 100 * exp(eta),
  name = "log(pct_dif)"
)        

# M==============================================================================================================================================

coef.post_rma <- function(post_rma_fit, ..., sep = get_emm_option("sep")){
  
  if(!inherits(post_rma_fit, "post_rma"))
    stop("post_rma_fit can be either from 'post_rma()' or 'contrast_rma()'.", call. = FALSE)
  
  ems <- post_rma_fit$ems
  Term <- term_names_(post_rma_fit=post_rma_fit, sep=sep)
  
  setNames(predict(ems), Term) 
} 

# M================================================================================================================================================

coef.contrast_rma <- function(post_rma_fit, ..., sep = get_emm_option("sep")){
  
  if(!inherits(post_rma_fit, "contrast_rma"))
    stop("post_rma_fit can be either from 'post_rma()' or 'contrast_rma()'.", call. = FALSE)

  if(inherits(post_rma_fit, "post_rma") & post_rma_fit$is_contr)
    stop("Use contrast_rma() instead of 'contr=' in post_rma().", call. = FALSE)
  
  ems <- post_rma_fit$con 
 Term <- term_names_(post_rma_fit=post_rma_fit, sep=sep)
  
  setNames(predict(ems), Term) 
} 

# M================================================================================================================================================

vcov.post_rma <- function(post_rma_fit, ..., sep = get_emm_option("sep")){
  
  if(!inherits(post_rma_fit, "post_rma"))
      stop("post_rma_fit can be either from 'post_rma()' or 'contrast_rma()'.", call. = FALSE)

  if(inherits(post_rma_fit, "post_rma") & post_rma_fit$is_contr)
    stop("Use contrast_rma() instead of 'contr=' in post_rma().", call. = FALSE)
  
  object <- post_rma_fit$ems
  
  vcov.emmGrid = function(object, ..., sep) {
    tol = get_emm_option("estble.tol")
    if (!is.null(hook <- object@misc$vcovHook)) {
      if (is.character(hook)) 
        hook = get(hook)
      hook(object, tol = tol, ...)
    }
    else {
      if(is.null(disp <- object@misc$display))
        disp = seq_len(nrow(object@linfct))
      X = object@linfct[disp, , drop = FALSE]
      estble = estimability::is.estble(X, object@nbasis, tol) 
      X[!estble, ] = NA
      X = X[, !is.na(object@bhat), drop = FALSE]
      rtn = X %*% tcrossprod(object@V, X)
      largs = as.list(object@grid[disp, seq_along(object@levels), drop = FALSE])
      largs$sep = sep
      rownames(rtn) = colnames(rtn) = do.call(paste, largs)
      return(rtn)
    }
  }
  
  vcov.emmGrid(object = object, ..., sep = sep)
}        

# M================================================================================================================================================

vcov.contrast_rma <- function(post_rma_fit, ..., sep = get_emm_option("sep")){
  
  if(!inherits(post_rma_fit, "contrast_rma"))
    stop("post_rma_fit can be either from 'post_rma()' or 'contrast_rma()'.", call. = FALSE)
  
  object <- post_rma_fit$con
  
  vcov.emmGrid = function(object, ..., sep) {
    tol = get_emm_option("estble.tol")
    if (!is.null(hook <- object@misc$vcovHook)) {
      if (is.character(hook)) 
        hook = get(hook)
      hook(object, tol = tol, ...)
    }
    else {
      if(is.null(disp <- object@misc$display))
        disp = seq_len(nrow(object@linfct))
      X = object@linfct[disp, , drop = FALSE]
      estble = estimability::is.estble(X, object@nbasis, tol) 
      X[!estble, ] = NA
      X = X[, !is.na(object@bhat), drop = FALSE]
      rtn = X %*% tcrossprod(object@V, X)
      largs = as.list(object@grid[disp, seq_along(object@levels), drop = FALSE])
      largs$sep = sep
      rownames(rtn) = colnames(rtn) = do.call(paste, largs)
      return(rtn)
    }
  }
  
  vcov.emmGrid(object = object, ..., sep = sep)
}        

# M================================================================================================================================================

AHW_tran <- list(
  linkfun = function(mu) transf.ahw(mu),
  linkinv = function(eta) transf.iahw(eta),
  mu.eta = function(eta) 3*(1-eta)^2,
  valideta = function (eta) 
    all(is.finite(eta)) && all(eta <= 1) && all(eta >= 0),
  name = "AHW"
)

# M================================================================================================================================================

ZR2_tran <- list(
  linkfun = function(mu) transf.r2toz(mu),
  linkinv = function(eta) transf.ztor2(eta),
  mu.eta = function(eta) 2*sinh(eta)/cosh(eta)^3,
  valideta = function (eta) 
    all(is.finite(eta)) && all(eta <= 1) && all(eta >= 0),
  name = "ZR2"
)

# M================================================================================================================================================

ABT_tran <- list(
  linkfun = function(mu) transf.abt(mu),
  linkinv = function(eta) transf.iabt(eta),
  mu.eta = function(eta) 1/(1-eta),
  valideta = function (eta) 
    all(is.finite(eta)) && all(eta <= 1) && all(eta >= 0),
  name = "ABT"
)

# M===============================================================================================================================================
        
r2z_tran <- list(
  linkfun = function(mu) atanh(mu),
  linkinv = function(eta) tanh(eta),
  mu.eta = function(eta) 1/cosh(eta)^2,
  valideta = function (eta) 
    all(is.finite(eta)) && all(abs(eta) <= 1),
  name = "r2z"
)                                        
# M=================================================================================================================================================
# Requires posttest_base_name followed by a number with no space
                                                                               
con_gain_list <- function(post_rma_fit, pretest_name = "baseline", brief=FALSE,
                          posttest_base_name = "posttest", gain_dif = FALSE, 
                          gain_dif_type = c("all","same","different"),
                          sep = get_emm_option("sep")){
  
  tms <- unlist(strsplit(term_names_(post_rma_fit), split = sep)) 
  
  specs <- post_rma_fit$specs 
  
  cond <- match.arg(gain_dif_type) 
  
  posttest_suffix <- "\\d+"
  
  if(!pretest_name %in% tms) stop("Wrong 'pretest_name='.", call. = FALSE)
  if(!posttest_base_name %in% str_remove(tms, posttest_suffix)) stop("Wrong 'posttest_base_name=' or 'posttest_suffix'.", call. = FALSE)
  
  varis <- if(is_bare_formula(specs, lhs=FALSE)) 
    .all.vars(specs) else 
      if(is.character(specs)) specs else stop("Only models with '~moderator * Time' are acceptable.", call. = FALSE)
  
  if(length(varis) > 2) stop("Only models with no more than TWO interacting variables (~moderator * Time) are acceptable.", call. = FALSE)
  
  DAT <- get_data_(post_rma_fit$rma.mv_fit)
  
  sec <- names(which(sapply(DAT, function(i) pretest_name %in% i)))
  
  varis <- varis[order(varis == sec)]
  
  tab <- post_rma_fit$table0[varis]
  
  DATA <- mutate(tab, Variables = apply(tab, 1, paste, collapse="_"),
                 Row = 1:nrow(tab))
  
  gain <- DATA %>%
    filter(!grepl(pretest_name, Variables)) %>%
    mutate(Variables2 = sub(paste0(posttest_base_name,posttest_suffix), pretest_name, Variables),
           Variables = paste0("(", Variables, " - ", Variables2, ")")) %>%
    right_join(filter(DATA, grepl(pretest_name, Variables)), 
               by = c("Variables2" = "Variables"), suffix = c("_post", "_pre")) %>%
    group_by(Variables) %>%
    summarize(Variables, Row = list(c(Row_post, -Row_pre))) %>%
    tibble::deframe()
  
  
  if(brief){  
    nms <- sub(paste0("^.([^_]+)\\D+(", posttest_suffix, ").*"), "Gain\\2(\\1)", names(gain))
    gain <- setNames(gain, nms)
  }
  
  
out <- if(gain_dif){
    
    gd_all <- do.call(c, combn(length(gain), 2, 
                               FUN = function(i)
                                 setNames(list(c(gain[[i[1]]], -gain[[i[2]]])),
                                          paste(names(gain)[i], collapse = " - ")), 
                               simplify = FALSE))
    
    if(brief){
      
      typ <- strcapture("Gain([0-9]+).*Gain([0-9]+)", names(gd_all), list(g1=0L, g2=0L)) 
      
      if(cond=="all") { gd_all
      } else if(cond=="same") { gd_all[with(typ, g1==g2)] 
      } else { gd_all[with(typ, g1!=g2)] }
      
    }
    
    
} else { gain }


res <- list(post_rma_fit = post_rma_fit, con_list = out)
class(res) <- "gain_list"
return(res)
  
}
#================================================================================================================================================

effect_count <- function(data, cluster, ..., arrange_by = NULL, show0 = TRUE,
                         na.rm = FALSE, subset){
  
  data <- full_clean(data)
  
  if(!missing(subset)) {
    
    s <- substitute(subset)
    data <- filter(data, eval(s))
  }
  
  cluster <- rlang::ensym(cluster)
  cat_mod <- rlang::ensyms(...)
  cat_nms <- purrr::map_chr(cat_mod, rlang::as_string)  
  c_nm <- as_string(cluster)
  clus_nm <- paste("n", c_nm)
  nms <- c(c_nm, cat_nms)
  
  idx <- nms %in% names(data)
  if(!all(idx)) stop(toString(dQuote(nms[!idx]))," not found in the 'data'.", call. = FALSE)
  
  dat <-  data %>% 
    count(!!cluster, !!!rlang::syms(cat_nms)) %>% 
    group_by(!!!rlang::syms(cat_nms)) %>% 
    summarize(!!clus_nm := n(),
              `n effect` = sum(n)) %>% 
    ungroup() 
  
  if(na.rm) dat <- drop_na(dat)
  
  if(show0) {
    dat <- dat %>% tidyr::complete(!!!rlang::syms(cat_nms), fill = list2(!!clus_nm := 0, `n effect` = 0))
  }
  
  dat %>% arrange(across(if(is.null(arrange_by)) all_of(cat_nms) else all_of(arrange_by)))    
}                              

#=================================================================================================================================================  
# To get SDs for each group, the best choice is the first formula (https://doi.org/10.1186/1471-2288-5-13):

range_iqr_n2sd <- function(min, max, q1, q3, n, dat) { ( (max-min) / (4*qnorm( (n-.375)/(n+.25)))) + 
    ( (q3-q1) / (4*qnorm( (n*.75-.125)/(n+.25))) ) } 

range_n2sd <- function(min, max, n) ( (max-min) / (2*qnorm( (n-.375)/(n+.25))))

iqr_n2sd <- function(q1, q3, n) ( (q3-q1) / (2*qnorm( (n*.75-.125)/(n+.25))) )

iqr_sd_cochrane <- function(q1, q3) (q3-q1)/1.35

range_f2sd <- function(min, max, f) f * (max-min)

#=================================================================================================================================================

# To get Means for each group, the best choice is the first formula:

med_range_iqr_n2mean <- function(min, max, q1, q3, n, median) {
  
  ((2.2/(2.2+n^.75))*((min+max)/2))+((.7-(.72/n^.55))*((q1+q3)/2))+((.3+(.72/n^.55)-(2.2/(2.2+n^.75)))*median)
  
}  

#=================================================================================================================================================

# Missing standard deviation and no other measure of variability: The Cochrane Handbook (6.5.2.3)
# Note that this SD is the average of the SDs of the two groups and so it this same SD should be 
# inputted into the meta-analysis for both groups.

mdifSE_n2sd <- function(mdifSE, n1, n2)  { mdifSE / ( sqrt ( (1/n1) + (1/n2) ) ) }  

#=================================================================================================================================================

cfactor <- function(df) exp(lgamma(df/2)-log(sqrt(df/2)) - lgamma((df-1)/2))

#=================================================================================================================================================

t2d <- function(t, n1, n2 = NA, g = TRUE){
  
  df <- ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
 d <- t/sqrt(N)
 ifelse(g==TRUE, cfactor(df)*d, d)
}

#=================================================================================================================================================

v_d <- function(d, n1, n2 = NA, g = FALSE, r = .5, cont.grp = FALSE){
  
  df <- ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)
  
  v <- if(is.na(n2) & !cont.grp) (1/n1) + ((d^2)/(2*n1)) 
  else if(is.na(n2) & cont.grp) ((2*(1-r))/n1) + ((d^2)/(2*n1)) 
  else ((n1 + n2)/(n1 * n2)) + ((d^2)/(2 * (n1+n2)))
  
  
  ifelse(g == TRUE, cfactor(df)^2 * v, v)
}

#=================================================================================================================================================

t2smd <- function(t, n1, n2 = NA, g = TRUE, r = .5, cont_sd = FALSE, d_given=NA, g_given=NA){
  
  df <- ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)
  
  d <- ifelse(is.na(d_given) & is.na(g_given), t2d(t, n1, n2, g), 
              ifelse(!is.na(g_given), g_given, 
                     ifelse(!is.na(d_given) & g==TRUE, cfactor(df)*d_given, d_given)))
  
  v <- ifelse(!is.na(g_given), v_d(d, n1, n2, g=TRUE, r, cont_sd), 
              v_d(d, n1, n2, g, r, cont_sd))
  data.frame(yi = d, vi = v)
}
                         
#======================== WCF Meta Dataset ======================================================================================================                

wcf <- read.csv("https://raw.githubusercontent.com/hkil/m/master/wcf.csv")

dat_design <- structure(list(Participant = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 51L, 51L, 51L, 51L, 51L, 51L, 51L, 51L, 51L, 
51L, 51L, 51L, 51L, 51L, 51L, 51L, 51L, 51L, 51L, 51L, 51L, 51L, 
51L, 51L, 51L, 51L, 51L, 51L), Order = c("s2c", "s2c", "s2c", 
"s2c", "s2c", "s2c", "s2c", "s2c", "s2c", "s2c", "s2c", "s2c", 
"s2c", "s2c", "s2c", "s2c", "s2c", "s2c", "s2c", "s2c", "s2c", 
"s2c", "s2c", "s2c", "s2c", "s2c", "s2c", "s2c", "c2s", "c2s", 
"c2s", "c2s", "c2s", "c2s", "c2s", "c2s", "c2s", "c2s", "c2s", 
"c2s", "c2s", "c2s", "c2s", "c2s", "c2s", "c2s", "c2s", "c2s", 
"c2s", "c2s", "c2s", "c2s", "c2s", "c2s", "c2s", "c2s"), Time = c(1L, 
1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 
3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L), Task = c("simple", "simple", "simple", 
"simple", "simple", "simple", "simple", "complex", "complex", 
"complex", "complex", "complex", "complex", "complex", "simple", 
"simple", "simple", "simple", "simple", "simple", "simple", "complex", 
"complex", "complex", "complex", "complex", "complex", "complex", 
"complex", "complex", "complex", "complex", "complex", "complex", 
"complex", "simple", "simple", "simple", "simple", "simple", 
"simple", "simple", "complex", "complex", "complex", "complex", 
"complex", "complex", "complex", "simple", "simple", "simple", 
"simple", "simple", "simple", "simple"), Questionnaire = c("Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", 
"Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng.", "Motiv./Eng."
), Outcomes = c("com_MLTU", "com_DC/T", "com_CN/T", "com_CN/C", 
"ac_EFC/C", "lex_Vocd", "lex_WRDFRQmc", "com_MLTU", "com_DC/T", 
"com_CN/T", "com_CN/C", "ac_EFC/C", "lex_Vocd", "lex_WRDFRQmc", 
"com_MLTU", "com_DC/T", "com_CN/T", "com_CN/C", "ac_EFC/C", "lex_Vocd", 
"lex_WRDFRQmc", "com_MLTU", "com_DC/T", "com_CN/T", "com_CN/C", 
"ac_EFC/C", "lex_Vocd", "lex_WRDFRQmc", "com_MLTU", "com_DC/T", 
"com_CN/T", "com_CN/C", "ac_EFC/C", "lex_Vocd", "lex_WRDFRQmc", 
"com_MLTU", "com_DC/T", "com_CN/T", "com_CN/C", "ac_EFC/C", "lex_Vocd", 
"lex_WRDFRQmc", "com_MLTU", "com_DC/T", "com_CN/T", "com_CN/C", 
"ac_EFC/C", "lex_Vocd", "lex_WRDFRQmc", "com_MLTU", "com_DC/T", 
"com_CN/T", "com_CN/C", "ac_EFC/C", "lex_Vocd", "lex_WRDFRQmc"
)), class = "data.frame", row.names = c(NA, -56L))

#================================================================================================================================================
dat_quasi <- read.table(header=TRUE, text=
"Group time          n  mpre sdpre mpost sdpost  r
 C     pre-post1    25  1.44  1.08  1.08   1.12 .5
 C     pre-post2    25  1.44  1.08  1.48   1.08 .5
 T     pre-post1    25  4.68  1.07  7.4    2.22 .5
 T     pre-post2    25  4.68  1.07  8.08   1.75 .5")        

#================================================================================================================================================

 dat_note <- read.table(header=TRUE, text = 
 "study                 year  g     v_g      assign_type  n_class Nt    Nc
  'Hayati & Jalilifar'  2009  0.213 0.101    student      NA      20    20
  'Hayati & Jalilifar'  2009  0.785 0.108    student      NA      20    20
  'Piri & Shirkhani'    2021  0.319 0.083    class        4       25    24
  'Piri & Shirkhani'    2021  0.986 0.0912   class        4       25    24
  'Piri & Shirkhani'    2021  0.723 0.087    class        4       25    24")

#=================================================================================================================================================

dat_ROM <- structure(list(mT = c(0.896193771626298, 0.496153846153846, 0.855882352941176, 
0.855203619909502, 0.932291666666667, 0.795833333333333, 0.787634408602151, 
0.830357142857143, 0.772972972972973, 0.661111111111111, 0.816993464052288, 
0.894736842105263, 0.824175824175824, 0.765, 0.846153846153846, 
0.818181818181818, 0.690598290598291, 0.688311688311688, 0.787654320987654, 
0.697788697788698, 0.626811594202899, 0.774553571428571, 0.889952153110048, 
0.828282828282828, 0.792016806722689, 0.717948717948718, 0.72280701754386, 
0.682051282051282, 0.632575757575758, 0.838095238095238, 0.884848484848485, 
0.810457516339869, 0.823529411764706, 0.6996336996337, 0.819444444444444, 
0.633333333333333, 0.890909090909091, 0.736111111111111, 0.92962962962963, 
0.660899653979239), mC = c(0.543252595155709, 0.319230769230769, 
0.594117647058824, 0.588235294117647, 0.671875, 0.583333333333333, 
0.502688172043011, 0.5625, 0.627027027027027, 0.466666666666667, 
0.630718954248366, 0.586466165413534, 0.681318681318681, 0.485, 
0.559440559440559, 0.56969696969697, 0.483760683760684, 0.471861471861472, 
0.607407407407407, 0.533169533169533, 0.445652173913043, 0.551339285714286, 
0.669856459330143, 0.646464646464646, 0.489495798319328, 0.602564102564103, 
0.53859649122807, 0.461538461538462, 0.526515151515151, 0.485714285714286, 
0.515151515151515, 0.516339869281046, 0.529411764705882, 0.406593406593407, 
0.578703703703704, 0.45, 0.618181818181818, 0.472222222222222, 
0.518518518518518, 0.480968858131488), sdT = c(0.0870742528177275, 
0.123468506620365, 0.091274943174622, 0.0744661299326225, 0.067726541840677, 
0.166613239299797, 0.117533899123826, 0.0913776106322973, 0.183558420383522, 
0.0962250448649376, 0.0829167278228901, 0.0933394709249929, 0.118864327762549, 
0.186120997313644, 0.113694623759974, 0.153000297088244, 0.11289614280979, 
0.153424121437801, 0.103009479336491, 0.161850601187354, 0.168305513035224, 
0.096190395393192, 0.0834239032256588, 0.0763946278835921, 0.0792152702654301, 
0.161721508012528, 0.120585474440829, 0.156208961925541, 0.157411615126825, 
0.174574312188794, 0.158353268282526, 0.0943194097363517, 0.131866499503455, 
0.113848133830916, 0.142112351884836, 0.121179596002999, 0.147709789175199, 
0.156478742031686, 0.0581736064053768, 0.135628442049199), sdC = c(0.125697727041304, 
0.130339852347371, 0.126933629176639, 0.129322607434048, 0.100159673657978, 
0.206724557648681, 0.146698025321454, 0.111544294473674, 0.189514949056574, 
0.177525072919719, 0.100488437827474, 0.206990841002065, 0.194777442800665, 
0.221359436211787, 0.169451203984782, 0.212845085221694, 0.10453791244664, 
0.169350682548532, 0.147147618575435, 0.154756484207011, 0.149934109375253, 
0.140808320821747, 0.114232884079553, 0.100224215528118, 0.108646527935511, 
0.205895998973177, 0.0984950771744254, 0.151073924009533, 0.17574678536683, 
0.286854866240254, 0.223776132721427, 0.161690416690889, 0.166930200177313, 
0.140039844589661, 0.176380884725896, 0.128803840296624, 0.22175998100985, 
0.149071198499986, 0.0985000064511608, 0.138321772258178), n = c(17L, 
20L, 40L, 13L, 12L, 40L, 31L, 28L, 37L, 12L, 18L, 19L, 13L, 40L, 
13L, 33L, 39L, 33L, 27L, 37L, 23L, 32L, 19L, 11L, 28L, 26L, 38L, 
15L, 24L, 21L, 33L, 17L, 34L, 21L, 24L, 30L, 22L, 16L, 18L, 17L
)), class = "data.frame", row.names = c(NA, -40L))

#=================================================================================================================================================
                         
dat_redundent <- structure(list(study = 1:6, Feedback_Type = c("Explicit", "Explicit", 
"Implicit", "Implicit", "Mixed", "Mixed"), Review_Allowed = c("Yes", 
"Yes", "No", "No", "No", "No"), No_Treats = c(2L, 1L, 4L, 2L, 
5L, 7L), Study_Length = c(3L, 3L, 3L, 2L, 6L, 9L), g = c(0.186773094628834, 
0.591821662111041, 0.0821856937949764, 1.2976404010689, 0.66475388590768, 
0.0897658079409924), v_g = c(0.153053426998667, 0.107615557732061, 
0.165476212999783, 0.124654886312783, 0.157642776239663, 0.198785914224572
)), class = "data.frame", row.names = c(NA, -6L)) 
                         
#=================================================================================================================================================

needzzsf <- c('metafor', 'clubSandwich', 'nlme', 'effects', 'lexicon', 'plotrix', 'rlang', 'emmeans','tidyverse','fastDummies')      

not.have23_. <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23_.)) install.packages(not.have23_.)

suppressWarnings(
  suppressMessages({ 
    
    invisible(lapply(needzzsf, base::require, character.only = TRUE))
    
  }))

options(dplyr.summarise.inform = FALSE)                        
