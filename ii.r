
peta.ci <- function(peta, f = NA, df1, df2, N, conf.level = .9, digits = 9){

ci <- Vectorize(function(peta, f, N, df1, df2, conf.level){
    
  q <- ifelse(is.na(f), peta2F(peta, df1, df2), f) 
  
    alpha <- (1 - conf.level)/2
    
    u <- function (ncp, alpha, q, df1, df2) {
      suppressWarnings(pf(q = q, df1 = df1, df2 = df2, ncp, lower.tail = FALSE)) - alpha
    }
    
    g <- try(uniroot(u, c(0, q+1e7), alpha = alpha, q = q, df1 = df1, df2 = df2)[[1]], silent = TRUE)
    if(inherits(g, "try-error")) g <- 0
    h <- try(uniroot(u, c(0, q+1e7), alpha = 1-alpha, q = q, df1 = df1, df2 = df2)[[1]], silent = TRUE)
    if(inherits(h, "try-error")) h <- 0
    I <- c(g, h)
    
    I <- I / (I + N)
    
    P.eta.sq <- if(is.na(f)) peta else F2peta(f, df1, df2)
    
    return(c(P.eta.sq = P.eta.sq, lower = I[1], upper = I[2], conf.level = conf.level, ncp = peta2ncp(P.eta.sq, N), F.value = q))
})

peta <- if(missing(peta)) NA else peta

round(data.frame(t(ci(peta = peta, f = f, N = N, df1 = df1, df2 = df2, conf.level = conf.level))), digits = digits)
}


R2.ci <- function(R2, n.pred, N, f = NA, df1 = NA, df2 = NA, conf.level = .9, digits = 9){ 
  
  if(is.na(df1) & is.na(df2)){
    df1 <- n.pred
    df2 <- N - n.pred - 1
  }    
  a <- if(!missing(R2)){ peta.ci(peta = R2, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = digits)
    } else { peta.ci(f = f, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = digits) }
  
  names(a)[1] <- "R2"
  a
}
