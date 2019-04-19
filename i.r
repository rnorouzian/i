
Break = "\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "   \"bayesL2\", a suite of R functions for Bayesian estimation in L2 research.
    Copyright (C) 2017-present  Reza Norouzian, rnorouzian@gmail.com

    This set of programs is free software: you can redistribute it under the 
    terms of the GNU General Public License as published by the Free 
    Software Foundation, either version 3 of the License, or any later 
    version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>."

message(Break, notice, Break)

Break = "\n*****************************************************************************\n"

cite = "To cite the package use:\n\nNorouzian, R., de Miranda, M. A., & Plonsky, L. (2018). The Bayesian \nrevolution in second language research: An applied approach. Language Learning, 64, 1032-1075.

\nNorouzian, R., de Miranda, M. A., & Plonsky, L. (2019). A Bayesian approach to measuring evidence \nin L2 research: An empirical investigation. Modern Language Journal, 103, 248-263."

cat(Break, cite, Break)

#==================================================================================================================

HDI <- function(fun, lower = 0, upper = 1, level = .95, eps = 1e-3)
{
  UseMethod("HDI")
}

HDI.default <- function(fun, lower = 0, upper = 1, level = .95, eps = 1e-3){
  
  if(!is.function(fun)) stop("Error: 'fun' must be a function.")
  if(length(formals(fun)) > 1) stop("Error: 'fun' must be a 'single-argument' function.")
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  x <- formals(fun)
  FUN <- function(x) fun(x)
  
  posterior <- function(x) FUN(x)/integrate(FUN, lower, upper)[[1]]
  mode <- optimize(posterior, c(lower, upper), maximum = TRUE)[[1]]
  inverse.posterior <- function(x, side = "left") {
    target <- function(y) posterior(y) - x
    ur <- switch(side,
                 left = try(uniroot(target, interval = c(lower, mode))),
                 right = try(uniroot(target, interval = c(mode, upper))))
    if(inherits(ur, "try-error")) stop("You may change prior parameters or 'lower' & 'upper'.", call. = FALSE)
    return(ur[[1]])
  }
  areafun <- function(h) {
    i1 <- inverse.posterior(h, "left")
    i2 <- inverse.posterior(h, "right")
    return(integrate(posterior, i1, i2)[[1]])
  }
  post.area <- 1
  find.lims <- function(a) {
    ur <- uniroot(function(h) areafun(h) / post.area - a,
                  c(eps, posterior(mode) - eps))
    return(ur[[1]])
  }
  f <- find.lims(level)
  return(c(inverse.posterior(f, "left"),
           inverse.posterior(f, "right")))
}


#==================================================================================================================

hdi <- function(x, y, level = .95)
{
  UseMethod("hdi")
}

hdi.default <- function(x, y, level = .95){
  
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  areas <- diff(x) * .5 * (head(y, -1) + tail(y, -1))
  peak <- which.max(areas)
  range <- c(peak, peak)
  found <- areas[peak]
  while(found < level) {
    if(areas[range[1]-1] > areas[range[2]+1]) {
      range[1] <- range[1]-1
      found <- found + areas[range[1]-1]
    } else {
      range[2] <- range[2]+1
      found <- found + areas[range[2]+1]
    }
  }
  val <- x[range]
  attr(val, "indexes") <- range
  attr(val, "area") <- found
  return(val)
}

#==================================================================================================================

hdir <- function(sample, level = .95)
{
  UseMethod("hdir")
}

hdir.default <- function(sample, level = .95){
  
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  #if(length(sample) < 1e3) message("Warning: \n\tInsufficient sample to produce reliable 'interval' estimates.")  
  sorted <- sort(sample)
  index <- ceiling(level*length(sorted))
  n <- length(sorted)- index
  width <- numeric(n)
  for(i in 1:n){
    width[i] <- sorted[i+ index]- sorted[i]
  }
  lower <- sorted[which.min(width)]
  upper <- sorted[which.min(width)+ index]
  return(c(lower, upper))
}

#==================================================================================================================

cip <- function(fun, lower = 0, upper = 1, level = .95){
    
  if(!is.function(fun)) stop("Error: 'fun' must be a function.")
  if(length(formals(fun)) > 1) stop("Error: 'fun' must be a 'single-argument' function.")
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  
  p <- (1 - level) / 2
  
  inv.cdf(c(p, 1-p), fun, lower, upper)
}

#====================================================================================================================


hdiq <- function(qdist, level = .95, ...)
{
  UseMethod("hdiq")
}


hdiq.default <- function(qdist, level = .95, ...)
  {

  alpha <-  1L - level
  width <- function(lower, qdist, level, ...){
    qdist(level + lower, ...) - qdist(lower, ...)
  }

  low <- optimize(width, c(0, alpha), qdist = qdist, level = level, ...)[[1]]
  
  return(c(qdist(low, ...), qdist(level + low, ...)))
}

#==================================================================================================================


cdf <- Vectorize(function(q, fun, lower, upper){
    
  if(!is.function(fun)) stop("Error: 'fun' must be a function.")
  if(length(formals(fun)) > 1) stop("Error: 'fun' must be a 'single-argument' function.")
  
  x <- formals(fun)
  f <- function(x) fun(x)/integrate(fun, lower, upper)[[1]]
  
  integrate(f, lower, q)[[1]]
  
}, c("q", "lower", "upper"))


#==================================================================================================================


inv.cdf <- Vectorize(function(p, fun, lower, upper){
    
  if(!is.function(fun)) stop("Error: 'fun' must be a function.")
  if(length(formals(fun)) > 1) stop("Error: 'fun' must be a 'single-argument' function.")
  
  uniroot(function(q) cdf(q, fun, lower, upper) - p, c(lower, upper), extendInt = "yes")[[1]]
  
}, c("p", "lower", "upper"))


#===================================================================================================================


eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
                            
                            
cor1 <- function(...)cor(cbind(...))
                            
                            
cor2 <- function(...) {
  vapply(combn(list(...), 2, simplify = FALSE), 
         function(y) cor(y[[1]], y[[2]]),
         numeric(1))
}


#==================================================================================================================                            
                            
                            
set.margin <- function() 
{
  par.mf <- par("mfrow", "mfcol")
  if (all(unlist(par.mf) == 1)) {
    par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1, 
        tck = -0.02)
  }
}
                            
#==================================================================================================================


prop.ci <- function(k, n, conf.level = .95, digits = 1e2)
{
  UseMethod("prop.ci")
}

prop.ci.default <- function(k, n, conf.level = .95, digits = 1e2){
  
  ci <- Vectorize(function(k, n, conf.level){
    
    I = as.numeric(binom.test(k, n, conf.level = conf.level)[[4]])
    return(c(Prop = k/n, lower = I[1], upper = I[2], conf.level = conf.level))
  })
  round(data.frame(t(ci(k = k, n = n, conf.level = conf.level))), digits = digits)
}

                               
#=================================================================================================================================                  
                  
d.cic <- function(d, n1, n2 = NA, conf.level = .95, digits = 1e2){
  
  ci <- Vectorize(function(d, n1, n2, conf.level){
    
    options(warn = -1)  
    alpha = (1 - conf.level)/2

    f <- function(dbase, alpha){
      alpha - suppressWarnings(pcohen(d, dbase, n1, n2, lower.tail = FALSE))
    }
    
    CI <- sapply(c(alpha, 1-alpha),
          function(x) uniroot(f, interval = c(-1e7, 1e7), alpha = x, extendInt = "yes")[[1]])
    
    return(c(Cohen.d = d, lower = CI[1], upper = CI[2], conf.level = conf.level))
  })
  
  round(data.frame(t(ci(d = d, n1 = n1, n2 = n2, conf.level = conf.level))), digits = digits)
} 
                                    

#=================================================================================================================================                
                
cor.ci <- function(r, n, conf.level = .95, digits = 1e2)
{
  UseMethod("cor.ci")
}

cor.ci.default <- function(r, n, conf.level = .95, digits = 1e2){
  
  ci <- Vectorize(function(r, n, conf.level){
    p = (1 - conf.level) / 2 
    g = tanh(atanh(r) + qnorm(c(p, 1-p))*1/sqrt(n - 3))
    g <- if(r == -1 || r == 1) rep(r, 2) else g 
    return(c(r = r, lower = g[1], upper = g[2], conf.level = conf.level))
  }) 
  round(data.frame(t(ci(r = r, n = n, conf.level = conf.level))), digits = digits)
}                            

#==================================================================================================================

beta.id <- function(Low, High, Cover = NA, digits = 1e2)
{
  UseMethod("beta.id")
}

beta.id.default <- function(Low, High, Cover = NA, digits = 1e2){

bet <- Vectorize(function(Low, High, Cover){
  
  options(warn = -1)
  L <- if(is.character(Low)) as.numeric(substr(Low, 1, nchar(Low)-1)) / 1e2 else Low
  U <- if(is.character(High)) as.numeric(substr(High, 1, nchar(High)-1)) / 1e2 else High
  
  if(L <= 0 || U >= 1) stop("NOTE: The smallest LOWER value that you can choose is \".000001\"AND the largest UPPER value is \".999999\".")
  if(L >= U) stop("Put the smaller value for Low, and the larger value for High")
  
  coverage  <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 1e2 else if(is.na(Cover)) .95 else Cover
  
  p1 = (1 - coverage) / 2 
  p2 = 1 - p1
  
  if( p1 <= 0 || p2 >= 1 || Low > High || p1 > p2 || coverage >= 1 ){
    stop("Error: \n\tUnable to find such a prior, make sure you have selected the correct values.") 
  } else {
    
    f.beta <- function(alpha, beta, x, lower = 0, upper = 1){
      p <- pbeta((x-lower)/(upper-lower), alpha, beta)
      log(p/(1-p))
    }
    
    delta <- function(fit, actual) sum((fit-actual)^2)
    
    objective <- function(theta, x, prob, ...){
      ab <- exp(theta)
      fit <- f.beta(ab[1], ab[2], x, ...)
      return (delta(fit, prob))
    }
    
    x.p <- (function(p) log(p/(1-p)))(c(p1, p2))
    
    sol <- nlm(objective, log(c(1e1, 1e1)), x = c(L, U), prob = x.p, lower = 0, upper = 1, typsize = c(1, 1), 
               fscale = 1e-12, gradtol = 1e-12)
    
    parm <- as.numeric(exp(sol$estimate))
    
    q <- qbeta(p = c(p1, p2), parm[1], parm[2])
    
    is.df <- function(a, b, sig = 3) round(a, sig) != round(b, sig)
    
    if(is.df(L, q[1]) || is.df(U, q[2])){
      
      stop("Error: \n\tUnable to find such a prior, make sure you have selected the correct values.")
      
    }else{
      
      return(c(a = parm[1], b = parm[2], Cover = Cover))   
    }
  } 
})
    round(data.frame(t(bet(Low = Low, High = High, Cover = Cover)), row.names = NULL), digits = digits)
}

#===============================================================================================

cauchy.id <- function(Low, High, Cover = NA)
{
  UseMethod("cauchy.id")
}

cauchy.id.default <- Vectorize(function(Low, High, Cover = NA){
  
  options(warn = -1)
  
  coverage  <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 1e2 else if(is.na(Cover)) .95 else Cover
  
  p1 = (1 - coverage) / 2
  p2 = 1 - p1
  
  if(p1 <= 0 || p2 >= 1 || Low > High || p1 > p2 || coverage >= 1) {
    
    stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
    
  } else {
    
    f <- function(x) {   
      y <- c(Low, High) - qcauchy(c(p1, p2), location = x[1],  scale = x[2])
    }
    
    parm <- optim(c(1, 1), function(x) sum(f(x)^2), control = list(reltol = (.Machine$double.eps)))[[1]]
  }
  
  q <- qcauchy(c(p1, p2), parm[[1]], parm[[2]])
  
  is.df = function(a, b, sig = 4) round(a, sig) != round(b, sig)
  
  if(is.df(Low, q[1]) || is.df(High, q[2])) {
    
    stop("\n\tUnable to find such a prior, make sure you have selected the correct values")
    
  } else { 
    
    return(c(mode = round(parm[[1]], 6), scale = round(parm[[2]], 6))) 
  }
})    

#===============================================================================================
      
logis.id <- function(Low, High, Cover = NA)
{
  UseMethod("logis.id")
}

logis.id.default <- Vectorize(function(Low, High, Cover = NA){
  
  options(warn = -1)
  
  coverage  <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 1e2 else if(is.na(Cover)) .95 else Cover
  
  p1 = (1 - coverage) / 2
  p2 = 1 - p1
  
  if(p1 <= 0 || p2 >= 1 || Low > High || p1 > p2 || coverage >= 1) {
    
    stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
    
  } else {
    
    f <- function(x) {   
      y <- c(Low, High) - qlogis(c(p1, p2), location = x[1],  scale = x[2])
    }
    
    parm <- optim(c(1, 1), function(x) sum(f(x)^2), control = list(reltol = (.Machine$double.eps)))[[1]]
  }
  
  q <- qlogis(c(p1, p2), parm[[1]], parm[[2]])
  
  is.df = function(a, b, sig = 4) round(a, sig) != round(b, sig)
  
  if(is.df(Low, q[1]) || is.df(High, q[2])) {
    
    stop("\n\tUnable to find such a prior, make sure you have selected the correct values")
    
  } else { 
    
    return(c(mode = round(parm[[1]], 6), scale = round(parm[[2]], 6))) 
  }
})
      
#============================================================================================================

tdist.id <- function(Low, High, Cover = NA)
{
  UseMethod("tdist.id")
}

tdist.id.default <- Vectorize(function(Low, High, Cover = NA){
  
  options(warn = -1)
  
  coverage  <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 1e2 else if(is.na(Cover)) .95 else Cover
  
  p1 = (1 - coverage) / 2
  p2 = 1 - p1
  
  if(p1 <= 0 || p2 >= 1 || Low > High || p1 > p2 || coverage >= 1) {
    
    stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
    
  } else {
    
    f <- function(x){   
      y <- c(Low, High) - qt(c(p1, p2), df = x[1], ncp = x[2])
    }
    
    parm <- optim(c(1, 0), function(x) sum(f(x)^2), control = list(reltol = (.Machine$double.eps)))[[1]]
  }
  
  q <- qt(c(p1, p2), parm[[1]], parm[[2]])
  
  is.df = function(a, b, sig = 4) round(a, sig) != round(b, sig)
  
  if(is.df(Low, q[1]) || is.df(High, q[2])) {
    
    stop("\n\tUnable to find such a prior, make sure you have selected the correct values")
    
  } else { 
    
    return(c(df = round(parm[[1]], 6), ncp = round(parm[[2]], 6))) 
  }
}) 
      
#============================================================================================================      

norm.id <- function(Low, High, Cover = NA, digits = 6)
{
  UseMethod("norm.id")
}

norm.id.default <- Vectorize(function(Low, High, Cover = NA, digits = 6){
  
  options(warn = -1)
  
  coverage <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 1e2 else if(is.na(Cover)) .95 else Cover
  
  p1 <- (1 - coverage) / 2 
  p2 <- 1 - p1
  
  q <- c(Low, High)  
  alpha <- c(p1, p2)
  
  is.df <- function(a, b, sig = 4) (round(a, sig) != round(b, sig))
  
  if( p1 <= 0 || p2 >= 1 || q[1] >= q[2] || p1 >= p2 ) {
    
    stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
    
  } else {
    
    beta <- qnorm(alpha)
    
    parm <- solve(cbind(1, beta), q)
    
    q <- qnorm(c(p1, p2), parm[[1]], parm[[2]])
  }
  
  if(is.df(Low, q[[1]]) || is.df(High, q[[2]])) {
    
    stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
  } else {
    
    return(c(mean = round(parm[[1]], digits = digits), sd = round(parm[[2]], digits = digits)))
    
  }
})

      
#===============================================================================================

      
prop.bayes <- function(a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", yes = 55, n = 1e2, level = .95, scale = .1, top = 1.1, 
                       show.prior = FALSE, bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .1, p.h0 = .5, digits = 6, 
                       col.depth = .55, labels = NULL, cex.lab = .8, xlab = NULL, ylab = NULL, col.hump = NULL, ...)
{
  UseMethod("prop.bayes")
}

prop.bayes.default <- function(a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", yes = 55, n = 1e2, 
                               level = .95, scale = .1, top = 1.1, show.prior = FALSE, bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .1,
                               p.h0 = .5, digits = 6, col.depth = .55, labels = NULL, cex.lab = .8, xlab = NULL, ylab = NULL, col.hump = NULL, ...){
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
  
  pr = show.prior
  I = eq(a, b, d, lo, hi, yes, n)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]] ; yes = I[[6]] ; n = I[[7]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)     
  
  if(!pr){   
    Bi = round(yes)
    n = round(n)                          
    loop = length(d)
    CI = matrix(NA, loop, 2)
    mode = numeric(loop)
    peak = numeric(loop)
    h = list()
    k = numeric(loop)
    eq.prob = numeric(loop)
    BF10 = numeric(loop)
    estimate = numeric(loop)
    
    if(any(yes > n)) stop("Error: 'yes' cannot be larger than 'n'.")
    for(i in 1:loop){
      p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) dbinom(Bi[i], n[i], x)
      k[i] = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k[i]    
      h[[i]] = list(x = x <- seq(0, 1, length.out = 5e2), y = posterior(x))
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior, level = level)
      peak[i] = posterior(mode[i])
      eq.prob[i] = integrate(posterior, lo[i], eq.hi)[[1]] - integrate(posterior, lo[i], eq.lo)[[1]]
      BF10[i] = k[i]/dbinom(yes[i], n[i], p.h0)
      estimate[i] <- yes[i]/n[i]     
    }
    graphics.off()
    lab <- if(is.null(labels)) substring(d, 2) else labels
    xlab <- if(is.null(xlab)) "Proportion" else xlab
    ylab <- if(is.null(ylab)) NA else ylab
    
    plot(CI, rep(1:loop, 2), type = "n", xlim = 0:1, ylim = c(bottom*1, top*loop), ylab = ylab, yaxt = "n", xaxt = "n", xlab = xlab, font.lab = 2, mgp = c(2, .3, 0), ...)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), labels = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
    axis(2, at = 1:loop, labels = lab, font = 2, las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .3, 0))
    
    for(i in 1:loop){
      col <- if(is.null(col.hump)) i else col.hump[i]    
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(col, col.depth), border = NA, xpd = NA)
    }
    col <- if(is.null(col.hump)) 1:loop else col.hump
    
    legend(x = leg, legend = rev(paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")")), 
           pch = 22, title = "Priors", pt.bg = rev(col), col = rev(col), cex = .7, pt.cex = .6, bg = 0, 
           box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4, adj = c(0, .3))
    box()
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = col, xpd = NA)
    m = scale*peak + 1:loop
    segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0(I[,1], "%", "    ", o, "%", "    ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
    
    rownames <- if(is.null(labels)) {paste0("Prop ", 1:loop, " posterior: ")} else {paste0(1:loop, " ", labels, " posterior:")}
    return(round(data.frame(estimate = estimate, mode = mode, lower = CI[,1], upper = CI[,2], eq.prob = eq.prob, BF10 = BF10, row.names = rownames), digits = digits))    
    
  }else{
      
    xlab <- if(is.null(xlab)) "Proportion" else xlab  
    p = function(x) get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1])
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = xlab, bty = "n", font.lab = 2, lwd = 2, n = 1e3, yaxs = "i", main = bquote(.(xlab)*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}
    
    
#================================================================================================      
      
      
prop.priors <- function(a, ...)
{
  UseMethod("prop.priors")
}  

prop.priors.default <- function(a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", yes = 55, n = 1e2, scale = .1, top = 1.5, show.prior = FALSE, bottom = 1, legend = "topleft"){
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  pr = show.prior
  is.v <- function(...) lengths(list(...)) > 1  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(a, b, d, lo, hi)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)     
  
  if(!pr){   
    Bi = round(yes)
    n = round(n)                          
    loop = length(d)
    CI = matrix(NA, loop, 2)
    mode = numeric(loop)
    peak = numeric(loop)
    h = list()
    
    
    if(any(is.v(yes, n))) stop("Error: 'yes' & 'n' must each have a length of '1'.")  
    if(yes > n) stop("Error: 'yes' cannot be larger than 'n'.")
    for(i in 1:loop){
      p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      h[[i]] = list(x = x <- seq(0, 1, length.out = 5e2), y = posterior(x))
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
      peak[i] = posterior(mode[i])
    }
    plot(CI, rep(1:loop, 2), type = "n", xlim = 0:1, ylim = c(bottom*1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = "Credible Interval (Proportion)", font.lab = 2, mgp = c(2, .3, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
    axis(2, at = 1:loop, lab = substring(d, 2), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
    legend(x = leg, legend = rev(paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")")), pch = 22, title = "Priors", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
    box()
    for(i in 1:loop){
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
    }
    m = scale*peak + 1:loop
    segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0(I[,1], "%", "    ", o, "%", "    ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
  }else{
    p = function(x) get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1])
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, yaxs = "i", main = bquote(Proportion*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}


#==========================================================================================================================


prop.hyper <- function(a, ...)
{
  UseMethod("prop.hyper")
}

prop.hyper.default <- function(a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", yes = 55, n = 1e2, show.prior = FALSE, pos = 3, top = 1.01){
  
  is.v <- function(...) lengths(list(...)) > 1
  
  pr = show.prior
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))   
  eq <- function(...) { lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(a, b, d, lo, hi)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop = length(a)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  
  if(!pr){  
    if(any(is.v(yes, n))) stop("Error: 'yes' & 'n' must each have a length of '1'.")
    if(yes > n) stop("Error: 'yes' cannot be larger than 'n'.")  
    Bi = round(yes)
    n = round(n)   
    for(i in 1:loop){
      p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
    }
    
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    
    par(mgp = c(2.2, .3, 0), mar = c(5.1, 4.1, 4.1, 3))   
    plot(CI, rep(1:loop, 2), type = "n", xlim = c(0, 1), ylim = c(1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = "Credible Interval (Proportion)", font.lab = 2)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4", xpd = NA)
    points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
    axis(2, at = 1:length(a), lab = deci(a), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    axis(4, at = 1:length(b), lab = deci(b), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    text(par('usr')[1:2], par('usr')[4], c("A", "B"), pos = 3, cex = 1.5, xpd = NA, font = 2)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0("[", I[,1], "%", ",  ", o, "%", ",  ", I[,2], "%", "]"), cex = .75, pos = pos, xpd = NA)
  }else{
    p = function(x) get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1])
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, yaxs = "i", main = bquote(Proportion*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}


#===================================================================================================================


prop.hyper.ab <- function(a, ...)
{
  UseMethod("prop.hyper.ab")
}

prop.hyper.ab.default <- function(a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", add = FALSE, 
                                  yes = 55, n = 1e2, col = 1, show.prior = FALSE){
  
  is.v <- function(...) lengths(list(...)) > 1
  pr = show.prior    
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(a, b, d, lo, hi)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  if(is.v(a) & pr || is.v(b) & pr) message("\tNote: You can see only '1 prior' at a time.")
  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop = length(d)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  
  if(!pr){   
    if(any(is.v(yes, n))) stop("Error: 'yes' & 'n' must each have a length of '1'.")
    if(yes > n) stop("Error: 'yes' cannot be larger than 'n'.")  
    Bi = round(yes)
    n = round(n)   
    for(i in 1:loop){
      p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]] 
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
    }
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI, type = "n", ylim = 0:1, xlim = c(1, loop), xlab = "Prior Parameter 'A'", xaxt = "n", yaxt = "n", ylab = "Credible Interval (Proportion)", font.lab = 2, mgp = c(2.3, .3, 0), cex.lab = 1.2)
    abline(v = 1:loop, col = 8, lty = 3)
    axis(2, at = axTicks(2), lab = paste0(axTicks(2)*1e2, "%"), mgp = c(2, .4, 0), las = 1)
    axis(3, at = 1:length(b), lab = deci(b), font = 2, las = 1, cex.axis = .8, mgp = c(2, .2, 0))
    text(mean(par('usr')[1:2]), 1.06*par('usr')[4], "Prior Parameter 'B'", pos = 3, cex = 1.2, xpd = NA, font = 2)
    axis(1, at = 1:length(a), lab = round(a, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .3, 0))
  }
  
  if(!pr){
    segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
    lines(1:loop, mode, col = col, lty = 3)
    points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
  }
  
  if(!add & pr){ 
    p = function(x) get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1])
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, yaxs = "i", main = bquote(Proportion*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}


#====================================================================================================================


prop.diff <- function(yes, n, a = 1.2, b = a, how = c("two.one", "one.two"), level = .95, top = 1.1, bottom = 1, scale = .1, margin = 6, legend = "topleft", eq.level = "2.5%", digits = 6, labels = NA, cex.lab = .8)
{
  UseMethod("prop.diff")
}

prop.diff.default <- function(yes, n, a = 1.2, b = a, how = c("two.one", "one.two"), level = .95, top = 1.1, bottom = 1, scale = .1, margin = 6, legend = "topleft", eq.level = "2.5%", digits = 6, labels = NA, cex.lab = .8){
  
  n <- round(n)
  yes <- round(yes)  
  loop <- length(n)
  
  is.s <- function(...)lengths(list(...)) < 2 
  
  if(any(yes > n)) stop("Error: 'yes' cannot be larger than 'n'.") 
  if(any(is.s(n, yes))) stop("Error: 'yes' & 'n' must each have a length of '2' or larger.")
  
  eq.b <- if(is.character(eq.level)) as.numeric(substr(eq.level, 1, nchar(eq.level)-1)) / 1e2 else eq.level
  legn <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  I = eq(n, yes)   
  n = I[[1]] ; yes = I[[2]] 
  
  comp <- ncol(combn(loop, 2))
  eq <- function(x) c(x, rep(rev(x)[1], ifelse(comp == 1, 1, comp - length(x))))
  
  a = eq(a)
  b = eq(b)[-length(b)]
  
  message(paste0("\n CAUTION: Check to see if you have chosen your desired ", "\"", 2*comp, "\"", " pairs of 'a' and 'b'."))
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  how <- match.arg(how)
  
  delta <- switch(how,
                  one.two = function(x) x[[1]] - x[[2]], 
                  two.one = function(x) x[[2]] - x[[1]])
  
  p <- list()
  for(i in 1:loop){
    p[[i]] <- rbeta(1e6, a[i] + yes[i], b[i] + (n[i] - yes[i]))
  }
  
  pi <- yes/n
  estimate <- combn(pi, 2, FUN = delta) 
  
  ps <- combn(p, 2, FUN = delta)
  
  loop <- ncol(ps)
  
  yeses <- combn(yes, 2, FUN = sum)
  ns <- combn(n, 2, FUN = sum)
  
  H0 <- numeric(loop)
  A0 <- rbeta(1e5, 1, 1)
  B0 <- list()                
  for(i in 1:loop){
    B0[[i]] <- dbinom(yeses[i], ns[i], A0)
    H0[i] <- mean(B0[[i]])
  }
  
  yes1 <- combn(yes, 2, simplify = FALSE)
  n1 <- combn(n, 2, simplify = FALSE)
  
  H1 <- numeric(loop)
  B1 <- list()
  B2 <- list()                
  for(i in 1:loop){
    B1[[i]] <- dbinom(yes1[[i]][1], n1[[i]][1], A0)
    B2[[i]] <- dbinom(yes1[[i]][2], n1[[i]][2], A0)
    H1[i] <- mean(B1[[i]])*mean(B2[[i]])
  }
  
  BF01 <- H0/H1                 
  BF10 <- 1/BF01
  
  CI <- matrix(NA, loop, 2)
  den <- list()
  mode <- numeric(loop)
  peak <- numeric(loop)
  mean <- numeric(loop)
  median <- numeric(loop)                  
  sd <- numeric(loop)
  from <- numeric(loop)                  
  to <- numeric(loop)
  BB <- numeric(loop)                 
  
  for(i in 1:loop){
    CI[i,] <- hdir(ps[, i], level = level)
    den[[i]] <- density(ps[, i], adjust = 2, n = 1e3)
    BB[i] <- mean(abs(ps[, i]) <= eq.b)
    mode[i] <- den[[i]]$x[which.max(den[[i]]$y)]
    peak[i] <- den[[i]]$y[which.max(den[[i]]$y)]
    mean[i] <- mean(ps[, i])
    median[i] <- median(ps[, i])
    sd[i] <- sd(ps[, i])
    from[i] <- mean[i] - margin * sd[i]
    to[i] <- mean[i] + margin * sd[i]
  }
  
  np <- combn(seq_along(p), 2, FUN = function(x){if(how == "one.two") paste0('p', x[1], ' - p', x[2]) else paste0('p', x[2], ' - p', x[1])})
  
  leg <- if(comp == 1) loop else 2
  lab <- if(is.na(labels)) np else labels
  plot(CI, rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*loop), ylab = NA, xaxt = "n", yaxt = "n", xlab = "Credible Interval (Proportion Differences)", font.lab = 2, mgp = c(2, .3, 0))
  axis(1, at = axTicks(1), labels = paste0(round(axTicks(1), 2)*1e2, "%"), mgp = c(2, .3, 0))
  abline(h = 1:loop, col = 8, lty = 3)
  axis(2, at = 1:loop, labels = lab, font = 2, las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .3, 0))
  legend(x = legn, legend = rep(rev(paste0("beta", "(", round(a, 2), ", ", round(b, 2), ")")), leg), pch = 22, title = "Priors", pt.bg = rep(loop:1, each = leg), col = rep(loop:1, each = leg), cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
  segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
  box()
  
  for(i in 1:loop){
    polygon(x = den[[i]]$x, y = scale*den[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
  }
  
  m = scale*peak + 1:loop
  segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
  points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
  I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
  text(mode, 1:loop, paste0(I[,1], "%", "      ", o, "%", "      ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
  
  return(round(data.frame(estimate = estimate, mean = mean, mode = mode, median = median, sd = sd, lower = CI[,1], upper = CI[,2], eq.prob = BB, BF01 = BF01, BF10 = BF10, row.names = paste0(np, ":")), digits = digits))                                                
}     

              
#====================================================================================================================

              
prop.diff.eq <- function(n1, n2, yes1, yes2, a1 = 1.2, b1 = 1.2, a2 = a1, b2 = b1, how = c("two.one", "one.two"), pL = -.025, pU = .025, level = .95, scale = .1)
{
  UseMethod("prop.diff.eq")
}

prop.diff.eq.default <- function(n1, n2, yes1, yes2, a1 = 1.2, b1 = 1.2, a2 = a1, b2 = b1, how = c("two.one", "one.two"), pL = -.025, pU = .025, level = .95, scale = .1){
  
  ro <- function(...){ lapply(list(...), function(x) round(x))}
  I <- ro(n1, n2, yes1, yes2)
  n1 <- I[[1]] ; n2 <- I[[2]] ; yes1 <- I[[3]] ; yes2 <- I[[4]]
  
  if(any(lengths(list(get(formalArgs(prop.diff.eq))))) > 1) stop("Error: Only 'one' comparison is allowed at a time.")
  if(yes1 > n1 || yes2 > n2) stop("Error: 'yes' cannot be larger than 'n'.")
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  p1 <- rbeta(1e6, a1 + yes1, b1 + (n1 - yes1))
  p2 <- rbeta(1e6, a2 + yes2, b2 + (n2 - yes2))
  
  how <- match.arg(how)
  
  delta <- switch(how, 
                  one.two = p1 - p2, 
                  two.one = p2 - p1) 
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(xpd = NA)
  
  d <- density(delta, adjust = 2, n = 1e5)
  
  plot(d, las = 1, type = "n", col = 0, main = NA, bty = "n", zero.line = FALSE,
       xlab = if(how == "one.two") bquote(Delta[~(p[1]-p[2])]) else bquote(Delta[~(p[2]-p[1])]), 
       cex.lab = 2, ylab = NA, axes = FALSE, yaxs = "i")
  
  axis(1, at = seq(min(d$x), max(d$x), length.out = 7), labels = paste0(deci(seq(min(d$x)*1e2, max(d$x)*1e2, length.out = 7), 2), "%"), mgp = c(2, .5, 0))
  
  polygon(x = d$x, y = scale*d$y, border = NA, col = rgb(1, 1, 0, .5)) # adjustcolor(4, .3)
  
  lines(d$x, scale*d$y, lwd = 2)     
  
  legend("topleft", c(paste0("group 1: ", "beta", "(", round(a1, 2), ", ", round(b1, 2), ")"), paste0("group 2: ", "beta", "(", round(a2, 2), ", ", round(b2, 2), ")")), title = "Priors", 
         pch = 22, col = 2, cex = .7, pt.cex = .6, pt.bg = 2, bty = "n", x.intersp = .5, title.adj = .4)
  
  mode <- d$x[which.max(d$y)]
  peak <- d$y[which.max(d$y)]*scale
  
  CI <- hdir(delta, level = level)
  segments(CI[1], 0, CI[2], 0, lend = 1, lwd = 4)
  segments(mode, 0, mode, peak, lend = 1, lty = 3)
  points(mode, 0, pch = 21, cex = 1.5, bg = "cyan")
  
  axis(side = 1, at = 0, mgp = c(3, 1.1, 0), col = 0, col.axis = "magenta", tick = FALSE, line = - 1.4, cex.axis = 1.4, font = 2)
  
  text(c(CI[1], mode, CI[2]), 0, paste0(c(deci(CI[1]*1e2, 2), deci(mode*1e2, 2), deci(CI[2]*1e2, 2)), "%"), pos = 3, 
       font = 2, col = "magenta", cex = .85)
  
  f <- approxfun(d$x, d$y, yleft = 0, yright = 0)
  
  cdf <- Vectorize(function(q){
    integrate(f, -1, q)[[1]]
  })
  
  # invcdf <- function(p){
  #  uniroot(function(q)cdf(q) - p, range(delta))[[1]]  # Not implemented # 
  # }
  
  y1 = y2 = 1.02*peak
  x.text = (pL+pU)/2
  y.text = 1.05*peak
  low.extreme <- par('usr')[3]
  
  segments(c(pL, pU), rep(low.extreme, 2), c(pL, pU), c(y1, y2), col = 'green2', lend = 1, lty = 2)
  
  segments(pL, 0, pU, 0, col = adjustcolor(3, .5), lend = 1, lwd = 80) 
  
  segments(c(pL, pU), c(y1, y2), rep(x.text, 2), rep(y.text*1.015, 2), lwd = 2, col = 'magenta')
  
  text(x.text, y.text, "Practically Equivalent to ZERO", font = 2, pos = 3, col = 'darkgreen', cex = .65, xpd = TRUE)
  
  points(c(pL, pU), c(y1, y2), pch = 21, col = 'green3', bg = 'green3', cex = 1.1)
  
  ## How much is it probable that the equivalence be true in population:
  
  a = cdf(pL)
  b = cdf(pU)
  
  Post.in.ROPE.Y = (b - a)
  Post.in.ROPE.X = (pU - pL) / 2
  
  BB = deci(Post.in.ROPE.Y*1e2, 2)
  
  title(main = paste0("There is ", "''", BB, "%", "''", " probability that TRUE diff. is equivalent to ZERO"), cex.main = .8)
  
  if(CI[1] > pU || CI[2] < pL) {
    
    legend("topright", "NOT Practically equivalent to \"0\" ", bty = 'n', cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
    
  } else
    
    if(CI[1] > pL & CI[2] < pU) {
      
      legend("topright", "Practically equivalent to \"0\" ", bty = 'n', cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
      
    } else {
      
      legend("topright", "No decision can be made ", bty = 'n', cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
    }
}             

       
#====================================================================================================================              

       
d.bayes <- function(t, n1, n2 = NA, m = 0, s = 1, level = .95, lo = -Inf, hi = Inf, dist.name = "dnorm", scale = .1, margin = 7, top = 1.1,
                    show.prior = FALSE, LL = -3, UL = 3, bottom = 1, prior.left = -6, prior.right = 6, legend = "topleft", eq.level = .1, 
                    d.h0 = 0, digits = 6, col.depth = .55, labels = NULL, cex.lab = .8, xlab = NULL, ylab = NULL, col.hump = NULL, ...){
  UseMethod("d.bayes")
}
       

d.bayes.default <- function(t, n1, n2 = NA, m = 0, s = 1, level = .95, lo = -Inf, hi = Inf, dist.name = "dnorm", scale = .1, margin = 7, top = 1.1,
                            show.prior = FALSE, LL = -3, UL = 3, bottom = 1, prior.left = -6, prior.right = 6, legend = "topleft", eq.level = .1,
                            d.h0 = 0, digits = 6, col.depth = .55, labels = NULL, cex.lab = .8, xlab = NULL, ylab = NULL, col.hump = NULL, ...){
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  pr <- show.prior
  
  if(!pr){    
    I = eq(m, s, d, lo, hi, t, n1, n2)
    m = I[[1]] ; s = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]] ; t = I[[6]] ; n1 = I[[7]] ; n2 = I[[8]]
    
    loop = length(d) 
    CI = matrix(NA, loop, 2)
    mode = numeric(loop)
    peak = numeric(loop)
    mean = numeric(loop)
    sd = numeric(loop)
    from = numeric(loop)
    to = numeric(loop) 
    h = list()
    k = numeric(loop)
    eq.prob = numeric(loop)
    BF10 = numeric(loop)
    estimate = numeric(loop)
    
    N = ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)    
    options(warn = -1)
    
    for(i in 1:loop){
      p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]        
      likelihood = function(x) dt(t[i], df[i], x*sqrt(N[i]))
      k[i] = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k[i]
      mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
      sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean[i]^2)
      from[i] = mean[i] - margin * sd[i]
      to[i] = mean[i] + margin * sd[i]
      mode[i] = optimize(posterior, c(from[i], to[i]), maximum = TRUE)[[1]]
      peak[i] = posterior(mode[i])
      CI[i,] = HDI(posterior, LL, UL, level = level)
      h[[i]] = list(x = x <- seq(from[i], to[i], length.out = 5e2), y = posterior(x))
      BF10[i] =  k[i] / dt(t[i], df[i], d.h0*sqrt(N[i]))
      eq.prob[i] = integrate(posterior, lo[i], eq.level)[[1]] - integrate(posterior, lo[i], -eq.level)[[1]]
      estimate[i] <- t[i]/sqrt(N[i])
    }    
    graphics.off()  
    lab <- if(is.null(labels)) substring(d, 2) else labels
    xlab <- if(is.null(xlab)) bquote(bold(delta)) else xlab
    ylab <- if(is.null(ylab)) NA else ylab    
    f = peak + 1:loop
    plot(CI, rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*max(f)), ylab = ylab, yaxt = "n", xlab = xlab, font.lab = 2, mgp = c(2, .5, 0), ...)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(2, at = 1:loop, labels = lab, font = 2, las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .3, 0))
    
    for(i in 1:loop){
      col <- if(is.null(col.hump)) i else col.hump[i]    
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(col, col.depth), border = NA, xpd = NA)
    }
    a = scale*(f-1:loop)+1:loop
    
    col <- if(is.null(col.hump)) 1:loop else col.hump 
    legend(x = leg, legend = rev(paste0(substring(d, 2), "(", round(m, 2), ", ", round(s, 2), ")")), pch = 22, title = "Priors", pt.bg = rev(col), col = rev(col), cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4, adj = c(0, .3))
    box()  
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = col, xpd = NA)                         
    segments(mode, 1:loop, mode, a, lty = 3, xpd = NA, lend = 1)
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = deci(CI) ; o = deci(mode)
    text(c(CI[,1], o, CI[,2]), 1:loop, c(I[,1], o, I[,2]), pos = 3, font = 2, cex = .8, xpd = NA)
    
    eq.decision <- if(CI[1] > -eq.level & CI[2] < eq.level) TRUE else if(CI[1] > eq.level || CI[2] < -eq.level) FALSE else NA
    
    rownames <- if(is.null(labels)) paste0("Cohen's d ", 1:loop, " posterior:") else paste0(1:loop, " ", labels, " posterior:")
    h <- round(data.frame(estimate = estimate, mode = mode, lower = CI[,1], upper = CI[,2], eq.prob = eq.prob, BF10 = BF10, row.names = rownames), digits = digits)
    h$equal <- eq.decision
    return(h)
    
  }else{
    xlab <- if(is.null(xlab)) bquote(bold("Effect Size "(delta))) else xlab  
    p = function(x) { get(d[1])(x, m[1], s[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, prior.left, prior.right, yaxt = "n", ylab = NA, xlab = xlab, bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(.(xlab)*" ~ "*.(if(lo[1] > -Inf || hi[1] < Inf) "truncated-")*.(substring(d[1], 2))(.(round(m[1], 2)), .(round(s[1], 2)))), mgp = c(2, .5, 0), yaxs = "i")
  }
}      
           
   
   
#====================================================================================================================
       
       
d.priors <- function(t, ...)
{
  UseMethod("d.priors")
}

d.priors.default <- function(t, n1, n2 = NA, m = 0, s = 1, lo = -Inf, hi = Inf, dist.name = "dnorm", scale = .1, margin = 7, top = 1.1, show.prior = FALSE, LL = -3, UL = 3, bottom = 1, prior.left = -6, prior.right = 6, legend = "topleft"){
  
  is.v <- function(...) lengths(list(...)) > 1
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))  
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi)
  m = I[[1]] 
  s = I[[2]] 
  d = I[[3]] 
  lo = I[[4]] 
  hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)                           
  loop = length(d) 
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  peak = numeric(loop)
  mean = numeric(loop)
  sd = numeric(loop)
  from = numeric(loop)
  to = numeric(loop) 
  h = list()
  
  if(!pr){    
    if(any(is.v(t, n1, n2))) stop("Error: 't' & 'n1' & 'n2' must each have a length of '1'.") 
    N = ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)   
    
    options(warn = -1)
    
    for(i in 1:loop){
      p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]        
      likelihood = function(x) dt(t, df, x*sqrt(N))
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
      sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean[i]^2)
      from[i] = mean[i] - margin * sd[i]
      to[i] = mean[i] + margin * sd[i]
      mode[i] = optimize(posterior, c(from[i], to[i]), maximum = TRUE)[[1]]
      peak[i] = posterior(mode[i])
      CI[i,] = HDI(posterior, LL, UL)
      h[[i]] = list(x = x <- seq(from[i], to[i], length.out = 5e2), y = posterior(x))
    }
    
    f = peak + 1:loop
    plot(CI, rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*max(f)), ylab = NA, yaxt = "n", xlab = bquote(bold("Credible Interval "(delta))), font.lab = 2, mgp = c(2, .5, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    legend(x = leg, legend = rev(paste0(substring(d, 2), "(", round(m, 2), ", ", round(s, 2), ")")), pch = 22, title = "Priors", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
    box()
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop)
    axis(2, at = 1:loop, lab = substring(d, 2), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
    
    for(i in 1:loop){
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
    }
    a = scale*(f-1:loop)+1:loop
    segments(mode, 1:loop, mode, a, lty = 3, xpd = NA, lend = 1)
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.1, col = 4, xpd = NA)
    I = deci(CI) ; o = deci(mode)
    text(c(CI[,1], o, CI[,2]), 1:loop, c(I[,1], o, I[,2]), pos = 3, font = 2, cex = .8, xpd = NA)
  }else{
    p = function(x) { get(d[1])(x, m[1], s[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, prior.left, prior.right, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo[1] > -Inf || hi[1] < Inf) "truncated-")*.(substring(d[1], 2))(.(round(m[1], 2)), .(round(s[1], 2)))), mgp = c(2, .5, 0))
  }
}


#========================================================================================================================


d.hyper <- function(t, ...)
{
  UseMethod("d.hyper")
}

d.hyper.default <- function(t, n1, n2 = NA, m = 0, s = 1, lo = -Inf, hi = Inf, dist.name = "dnorm", LL = -3, UL = 3, pos = 3, show.prior = FALSE, top = 1.01, margin = 6, prior.left = -6, prior.right = 6){
  
  is.v <- function(...) lengths(list(...)) > 1
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi)
  m = I[[1]] ; s = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]] 
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)  
  
  options(warn = -1)
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  loop = length(m)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  mean = numeric(loop)
  sd = numeric(loop)
  from = numeric(loop)
  to = numeric(loop)
  
  if(!pr){
    if(any(is.v(t, n1, n2))) stop("Error: 't' & 'n1' & 'n2' must each have a length of '1'.")
    N = ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2) 
    
    for(i in 1:loop){
      p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]  
      likelihood = function(x) dt(t, df, x*sqrt(N))
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(LL, UL), maximum = TRUE)[[1]]
      mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
      sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean[i]^2)
      CI[i,] = HDI(posterior, LL, UL)
      from[i] = mean[i] - margin * sd[i]
      to[i] = mean[i] + margin * sd[i] 
    }
    
    par(mgp = c(2, .5, 0), mar = c(5.1, 4.1, 4.1, 3))   
    plot(CI, rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(1, top*loop), ylab = NA, yaxt = "n", xlab = bquote(bold("Credible Interval "(delta))), font.lab = 2)
    abline(h = 1:loop, col = 8, lty = 3)
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4")  
    points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
    axis(2, at = 1:length(m), lab = deci(m), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    axis(4, at = 1:length(s), lab = deci(s), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    text(par('usr')[1:2], par('usr')[4], c("M", "S"), pos = 3, cex = 1.5, xpd = NA, font = 2)
    I = deci(CI) ; o = deci(mode)
    text(mode, 1:loop, paste0("[", I[,1], ",  ", o, ",  ", I[,2], "]"), pos = pos, cex = .8, xpd = NA)
  }else{
    p = function(x) { get(d[1])(x, m[1], s[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, prior.left, prior.right, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo[1] > -Inf || hi[1] < Inf) "truncated-")*.(substring(d[1], 2))(.(round(m[1], 2)), .(round(s[1], 2)))), mgp = c(2, .5, 0))
  }  
}


#===================================================================================================================


d.hyper.ms <- function(t, ...)
{
  UseMethod("d.hyper.ms")
}

d.hyper.ms.default <- function(t, n1, n2 = NA, m = 0, s = 1, lo = -Inf, hi = Inf, dist.name = "dnorm", add = FALSE, 
                               col = 1, top = 6, margin = 1.01, LL = -3, UL = 3, show.prior = FALSE, prior.left = -6, prior.right = 6){
  
  is.v <- function(...) lengths(list(...)) > 1
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))  
  pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi)
  m = I[[1]] ; s = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")
  
  options(warn = -1)
  loop = length(m) 
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  mean = numeric(loop)
  sd = numeric(loop)
  from = numeric(loop)
  to = numeric(loop)
  
  if(!pr){   
    if(any(is.v(t, n1, n2))) stop("Error: 't' & 'n1' & 'n2' must each have a length of '1'.")
    N = ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2) 
    for(i in 1:loop){
      p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) dt(t, df, x*sqrt(N))
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(LL, UL), maximum = TRUE)[[1]]
      mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
      sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean[i]^2)
      CI[i,] = HDI(posterior, LL, UL)
      from[i] = mean[i] - top * sd[i]
      to[i] = mean[i] + top * sd[i]
    }
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI, type = "n", ylim = c(min(from), max(to)), xlim = c(1, margin*loop), xlab = "Prior Parameter 'M'", xaxt = "n", ylab = bquote(bold("Credible Interval "(delta))), font.lab = 2)
    abline(v = 1:loop, col = 8, lty = 3, mgp = c(2, .5, 0))
    axis(3, at = 1:length(s), lab = round(s, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .4, 0))
    text(mean(par('usr')[1:2]), 1.06*par('usr')[4], "Prior Parameter 'S'", pos = 3, cex = 1, xpd = NA, font = 2)
    axis(1, at = 1:length(m), lab = round(m, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .3, 0))
  }
  
  if(!pr){
    segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
    lines(1:loop, mode, col = col, lty = 3)
    points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
  }
  
  if(!add & pr){
    p = function(x){ get(d[1])(x, m[1], s[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, prior.left, prior.right, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo[1] > -Inf || hi[1] < Inf) "truncated-")*.(substring(d[1], 2))(.(round(m[1], 2)), .(round(s[1], 2)))), mgp = c(2, .5, 0))
  }
}


#==================================================================================================================


peta.bayes <- function(f, N, df1, df2, a = 1.2, b = 1.2, level = .95, lo = 0, hi = 1, dist.name = "dbeta", scale = .1, top = 1.1, show.prior = FALSE, 
                       bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .05, peta.h0 = 0, digits = 6, col.depth = .55, labels = NULL, cex.lab = .8, 
                       xlab = NULL, ylab = NULL, col.hump = NULL, ...){ 
  UseMethod("peta.bayes")
}

peta.bayes.default <- function(f, N, df1, df2, a = 1.2, b = 1.2, level = .95, lo = 0, hi = 1, dist.name = "dbeta", scale = .1, top = 1.1, show.prior = FALSE, 
                               bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .05, peta.h0 = 0, digits = 6, col.depth = .55, labels = NULL, cex.lab = .8, 
                               xlab = NULL, ylab = NULL, col.hump = NULL, ...){
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)      
  pr <- show.prior
    
  if(!pr){    
  I <- eq(a, b, d, lo, hi, f, N, df1, df2)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]] ; f = I[[6]] ; N = I[[7]] ; df1 = I[[8]] ; df2 = I[[9]] 
                                                                                                                          
  loop <- length(a)
  CI <- matrix(NA, loop, 2)
  mode <- numeric(loop)
  peak <- numeric(loop)
  h <- list()
  k = numeric(loop)
  eq.prob = numeric(loop)
  BF10 = numeric(loop)
  estimate = numeric(loop)
 
    options(warn = -1)
    
    for(i in 1:loop){
      p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]  
      likelihood = function(x) df(f[i], df1[i], df2[i], (x * N[i]) / (1 - x) )
      k[i] = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k[i]
      h[[i]] = list(x = x <- seq(0, 1, length.out = 5e2), y = posterior(x))
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      peak[i] = posterior(mode[i])
      CI[i,] = HDI(posterior, 0, .9999999, level = level)
      BF10[i] =  k[i] / df(f[i], df1[i], df2[i], (peta.h0 * N[i]) / (1 - peta.h0))
      eq.prob[i] = integrate(posterior, lo[i], eq.hi)[[1]] - integrate(posterior, lo[i], eq.lo)[[1]]
      estimate[i] <- (f[i]*df1[i]) / ((f[i]*df1[i]) + df2[i])
    } 
    graphics.off()
    lab <- if(is.null(labels)) substring(d, 2) else labels
    xlab <- if(is.null(xlab)) bquote(bold((eta[p]^2))) else xlab
    ylab <- if(is.null(ylab)) NA else ylab  
        
    plot(CI, rep(1:loop, 2), type = "n", xlim = 0:1, ylim = c(bottom*1, top*loop), ylab = ylab, yaxt = "n", xaxt = "n", xlab = xlab, font.lab = 2, mgp = c(2, .5, 0), ...)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), labels = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0)) 
    axis(2, at = 1:loop, labels = lab, font = 2, las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .3, 0))
                      
    for(i in 1:loop){
    col <- if(is.null(col.hump)) i else col.hump[i]     
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(col, col.depth), border = NA, xpd = NA)
    }
    m = scale*peak + 1:loop
    col <- if(is.null(col.hump)) 1:loop else col.hump
    legend(x = leg, legend = rev(paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")")), pch = 22, title = "Priors", pt.bg = rev(col), col = rev(col), cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4, adj = c(0, .3))
    box()
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = col, xpd = NA)                   
    segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0(I[,1], "%", "    ", o, "%", "    ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
    
    rownames <- if(is.null(labels)) paste0("P.eta.sq ", 1:loop, " posterior:") else paste0(1:loop, " ", labels, " posterior:")                   
    return(round(data.frame(estimate = estimate, mode = mode, lower = CI[,1], upper = CI[,2], eq.prob = eq.prob, BF10 = BF10, row.names = rownames), digits = digits))  
    
}else{
    xlab <- if(is.null(xlab)) bquote(bold("Partial Eta.Sq"~(eta[p]^2))) else xlab
    p = function(x) { get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = xlab, bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(.(xlab)*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))), yaxs = "i")
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}


#===================================================================================================================


peta.priors <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", scale = .1, top = 1.1, show.prior = FALSE, bottom = 1, legend = "topleft")
{
  UseMethod("peta.priors")
}

peta.priors.default <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", scale = .1, top = 1.1, show.prior = FALSE, bottom = 1, legend = "topleft"){
  
  is.v <- function(...) lengths(list(...)) > 1
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  pr <- show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I <- eq(a, b, d, lo, hi)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)                                                                                                                           
  
  loop <- length(a)
  CI <- matrix(NA, loop, 2)
  mode <- numeric(loop)
  peak <- numeric(loop)
  h <- list()
  
  if(!pr){  
    
    if(any(is.v(f, N, df1, df2))) stop("Error: 'f' & 'N' & 'df1' & 'df2' must each have a length of '1'.")  
    options(warn = -1)
    
    for(i in 1:loop){
      p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]  
      likelihood = function(x) df(f, df1, df2, (x * N) / (1 - x) )
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      h[[i]] = list(x = x <- seq(0, 1, length.out = 5e2), y = posterior(x))
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      peak[i] = posterior(mode[i])
      CI[i,] = HDI(posterior, 0, .9999999)
    } 
    
    plot(CI, rep(1:loop, 2), type = "n", xlim = 0:1, ylim = c(bottom*1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = bquote(bold("Credible Interval"~(eta[p]^2))), font.lab = 2, mgp = c(2, .5, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0)) 
    axis(2, at = 1:loop, lab = substring(d, 2), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
    legend(x = leg, legend = rev(paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")")), pch = 22, title = "Priors", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
    box()
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
    for(i in 1:loop){
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
    }
    m = scale*peak + 1:loop
    segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0(I[,1], "%", "    ", o, "%", "    ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
  }else{
    p = function(x) { get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = bquote(bold("Partial Eta.Sq"~(eta[p]^2))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(eta[p]^2*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}


#===================================================================================================================


peta.hyper <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", show.prior = FALSE, pos = 3, top = 1.01)
{
  UseMethod("peta.hyper")
}

peta.hyper.default <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", show.prior = FALSE, pos = 3, top = 1.01){
  
  is.v <- function(...) lengths(list(...)) > 1
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  pr <- show.prior
  
  eq <- function(...) { lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I <- eq(a, b, d, lo, hi)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop <- length(a)
  CI <- matrix(NA, loop, 2)
  mode <- numeric(loop)
  
  if(!pr){  
    if(any(is.v(f, N, df1, df2))) stop("Error: 'f' & 'N' & 'df1' & 'df2' must each have a length of '1'.")
    options(warn = -1)
    
    for(i in 1:loop){
      p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) df(f, df1, df2, (x * N) / (1 - x) )
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior, 0, .9999999)
    }
    
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    
    par(mgp = c(2.2, .3, 0), mar = c(5.1, 4.1, 4.1, 3))   
    plot(CI, rep(1:loop, 2), type = "n", xlim = 0:1, ylim = c(1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = bquote(bold("Credible Interval"~(eta[p]^2))), font.lab = 2)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4", xpd = NA)
    points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
    axis(2, at = 1:length(a), lab = deci(a), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    axis(4, at = 1:length(b), lab = deci(b), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    text(par('usr')[1:2], par('usr')[4], c("A", "B"), pos = 3, cex = 1.5, xpd = NA, font = 2)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0("[", I[,1], "%", ",  ", o, "%", ",  ", I[,2], "%", "]"), cex = .75, pos = pos, xpd = NA)
  }else{
    p = function(x) { get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = bquote(bold("Partial Eta.Sq"~(eta[p]^2))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(eta[p]^2*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}


#===================================================================================================================


peta.hyper.ab <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", add = FALSE, 
                                  col = 1, show.prior = FALSE)
{
  UseMethod("peta.hyper.ab")
}

peta.hyper.ab.default <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", add = FALSE, 
                                  col = 1, show.prior = FALSE){
  
  is.v <- function(...) lengths(list(...)) > 1
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  pr <- show.prior    
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I <- eq(a, b, d, lo, hi)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  if(is.v(a) & pr || is.v(b) & pr) message("\tNote: You can see only '1 prior' at a time.")
  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop <- length(a)
  CI <- matrix(NA, loop, 2)
  mode <- numeric(loop)
  
  options(warn = -1)
  
  if(!pr){    
    if(any(is.v(f, N, df1, df2))) stop("Error: 'f' & 'N' & 'df1' & 'df2' must each have a length of '1'.")                       
    for(i in 1:loop){
      p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) df(f, df1, df2, (x * N) / (1 - x) )
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior, 0, .9999999)
    }
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI, type = "n", ylim = 0:1, xlim = c(1, loop), xlab = "Prior Parameter 'A'", xaxt = "n", yaxt = "n", ylab = bquote(bold("Credible Interval"~(eta[p]^2))), font.lab = 2, mgp = c(2.3, .3, 0), cex.lab = 1.2)
    abline(v = 1:loop, col = 8, lty = 3)
    axis(2, at = axTicks(2), lab = paste0(axTicks(2)*1e2, "%"), mgp = c(2, .4, 0), las = 1)
    axis(3, at = 1:length(b), lab = round(b, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .2, 0))
    text(mean(par('usr')[1:2]), 1.06*par('usr')[4], "Prior Parameter 'B'", pos = 3, cex = 1.2, xpd = NA, font = 2)
    axis(1, at = 1:length(a), lab = round(a, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .3, 0))
  }
  
  if(!pr){
    segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
    lines(1:loop, mode, col = col, lty = 3)
    points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
  }
  
  if(!add & pr){
    p = function(x) { get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = bquote(bold("Partial Eta.Sq"~(eta[p]^2))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(eta[p]^2*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}


#=================================================================================================================


cor.bayes <- function(r, n, prior.mean = 0, prior.sd = .707, eq.bound = .05, level = .95, top = 1.1, bottom = 1, scale = .1, margin = 7,
                      legend = "topleft", show.prior = FALSE, digits = 6, col.depth = .55, labels = NULL, cex.lab = .8, xlab = NULL, 
                      ylab = NULL, col.hump = NULL, ...){
  UseMethod("cor.bayes")
}

cor.bayes.default <- function(r, n, prior.mean = 0, prior.sd = .707, eq.bound = .05, level = .95, top = 1.1, bottom = 1, scale = .1, 
                              margin = 5, legend = "topleft", show.prior = FALSE, digits = 6, col.depth = .55, labels = NULL, 
                              cex.lab = .8, xlab = NULL, ylab = NULL, col.hump = NULL, ...){ 
  
  pr <- show.prior    
  mu <- prior.mean
  lambda <- 1/(prior.sd^2)
  
  if(!pr){   
    I = eq(n, r, prior.mean, prior.sd)   
    n = I[[1]] ; r = I[[2]] ; prior.mean = I[[3]] ; prior.sd = I[[4]] ;  
    
    deci <- function(x, k = 3) format(round(x, k), nsmall = k)
    leg <- if(is.character(legend)) legend else deparse(substitute(legend))
    
    lambda.post <- (lambda + (n - 3))
    mu.post <- (lambda*mu + (n - 3)*atanh(r))/lambda.post
    
    loop <- length(r)
    p <- list()
    CI <- matrix(NA, loop, 2)
    den <- list()
    mode <- numeric(loop)
    peak <- numeric(loop)
    mean <- numeric(loop)
    median <- numeric(loop)                  
    sd <- numeric(loop)
    from <- numeric(loop)                  
    to <- numeric(loop)
    eq.prob <- numeric(loop)          
    
    for(i in 1:loop){
      p[[i]] <- tanh(rnorm(1e6, mu.post[i], sqrt(1/lambda.post[i])))
      CI[i,] <- hdir(p[[i]], level = level)
      den[[i]] <- density(p[[i]], adjust = 2, n = 1e3)
      eq.prob[i] <- mean(abs(p[[i]]) <= eq.bound)
      mode[i] <- den[[i]]$x[which.max(den[[i]]$y)]
      peak[i] <- den[[i]]$y[which.max(den[[i]]$y)]
      mean[i] <- mean(p[[i]])
      median[i] <- median(p[[i]])
      sd[i] <- sd(p[[i]])
      from[i] <- mean[i] - margin * sd[i]
      to[i] <- mean[i] + margin * sd[i]
    }
    graphics.off()
    lab <- if(is.null(labels)) paste0("r", 1:loop) else labels
    xlab <- if(is.null(xlab)) "Credible Interval (Pearson correlation)" else xlab
    ylab <- if(is.null(ylab)) NA else ylab    
    plot(CI, rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*loop), ylab = ylab, yaxt = "n", xlab = xlab, font.lab = 2, mgp = c(2, .3, 0), ...)
    axis(2, at = 1:loop, labels = lab, font = 2, las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .3, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    
    for(i in 1:loop){
      col <- if(is.null(col.hump)) i else col.hump[i]    
      polygon(x = den[[i]]$x, y = scale*den[[i]]$y +i, col = adjustcolor(col, col.depth), border = NA, xpd = NA)
    }
    
    m = scale*peak + 1:loop
    col <- if(is.null(col.hump)) 1:loop else col.hump
    legend(x = leg, legend = rev(paste0("s.norm", "(", round(prior.mean, 2), ", ", round(prior.sd, 2), ")")), pch = 22, title = "Priors", pt.bg = rev(col), col = rev(col), cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4, adj = c(0, .3))
    box()
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = col, xpd = NA)                            
    segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = deci(CI, 2); o = deci(mode, 2)
    text(mode, 1:loop, paste0(I[,1], "        ", o, "         ", I[,2]), cex = .75, pos = 3, font = 2, xpd = NA)
    
    rownames <- if(is.null(labels)) paste0("r", 1:loop, " posterior:") else paste0(1:loop, " r ", labels, " posterior:")                            
    return(round(data.frame(mean = mean, mode = mode, median = median, sd = sd, lower = CI[,1], upper = CI[,2], eq.prob = eq.prob, row.names = rownames), digits = digits))
    
  }else{
    xlab <- if(is.null(xlab)) bquote(rho[~("Pearson correlation")]) else xlab
    p <- function(x) dnorm(atanh(x), prior.mean[1], prior.sd[1])*1/(1-x^2)
    curve(p, -1, 1, yaxt = "n", ylab = NA, xlab = xlab, bty = "n", font.lab = 2, cex.lab = 1.5, lwd = 2, n = 1e4, yaxs = "i", main = bquote(rho*" ~ "*"scaled.norm"(.(round(prior.mean[1], 3)), .(round(prior.sd[1], 3)))), xpd = NA) 
  }  
}

                              
#==================================================================================================================                              

                              
cor.diff <- function(r, n, prior.mean = 0, prior.sd = .707, how = c("two.one", "one.two"), eq.bound = .05, level = .95, top = 1.1, bottom = 1, scale = .1, margin = 5, legend = "topleft", digits = 6, labels = NA, cex.lab = .8, xlab = NA, ylab = NA, ...)
{
  UseMethod("cor.diff")
}

cor.diff.default <- function(r, n, prior.mean = 0, prior.sd = .707, how = c("two.one", "one.two"), eq.bound = .05, level = .95, top = 1.1, bottom = 1, scale = .1, margin = 5, legend = "topleft", digits = 6, labels = NA, cex.lab = .8, xlab = NA, ylab = NA, ...){ 
  
  is.s <- function(...)lengths(list(...)) < 2
  if(any(is.s(n, r))) stop("Error: 'r' & 'n' must each have a length of '2' or larger.")
  
  I = eq(n, r, prior.mean, prior.sd)   
  n = I[[1]] ; r = I[[2]] ; prior.mean = I[[3]] ; prior.sd = I[[4]]
  
  legn <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  loop <- length(r)
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  mu     <- prior.mean
  lambda <- 1/(prior.sd^2)
  
  how <- match.arg(how)
  
  delta <- switch(how,
                  one.two = function(x) x[[1]] - x[[2]], 
                  two.one = function(x) x[[2]] - x[[1]])
  
  lambda.post <- (lambda + (n - 3))
  mu.post     <- (lambda*mu + (n - 3)*atanh(r))/lambda.post
  
  p <- list()
  for(i in 1:loop){
    p[[i]] <- tanh(rnorm(1e6, mu.post[i], sqrt(1/lambda.post[i])))
  }
  
  ps <- combn(p, 2, FUN = delta)
  
  loop <- ncol(ps)
  
  CI <- matrix(NA, loop, 2)
  den <- list()
  mode <- numeric(loop)
  peak <- numeric(loop)
  mean <- numeric(loop)
  median <- numeric(loop)                  
  sd <- numeric(loop)
  from <- numeric(loop)                  
  to <- numeric(loop)
  BB <- numeric(loop)
  
  for(i in 1:loop){
    CI[i,] <- hdir(ps[, i], level = level)
    den[[i]] <- density(ps[, i], adjust = 2, n = 1e3)
    BB[i] <- mean(abs(ps[, i]) <= eq.bound)
    mode[i] <- den[[i]]$x[which.max(den[[i]]$y)]
    peak[i] <- den[[i]]$y[which.max(den[[i]]$y)]
    mean[i] <- mean(ps[, i])
    median[i] <- median(ps[, i])
    sd[i] <- sd(ps[, i])
    from[i] <- mean[i] - margin * sd[i]
    to[i] <- mean[i] + margin * sd[i]
  }
  
  np <- combn(seq_along(p), 2, FUN = function(x){if(how == "one.two") paste0('r', x[1], ' - r', x[2]) else paste0('r', x[2], ' - r', x[1])})
  
  leg <- if(loop == 1) 1 else 2
  
  lab <- if(is.na(labels)) np else labels
  xlab <- if(is.na(xlab)) "Credible Interval (Correlation Differences)" else xlab
  ylab <- if(is.na(ylab)) NA else ylab
      
  plot(CI, rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*loop), ylab = ylab, yaxt = "n", xlab = xlab, font.lab = 2, mgp = c(2, .3, 0), ...)
  axis(2, at = 1:loop, labels = lab, font = 2, las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .3, 0))
  abline(h = 1:loop, col = 8, lty = 3)
  legend(x = legn, legend = rev(paste0("s.norm", "(", round(rep(prior.mean[1], loop), 2), ", ", round(rep(prior.sd[1], loop), 2), ")")), pch = 22, title = "Priors", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
  segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
  box()
  
  for(i in 1:loop){
    polygon(x = den[[i]]$x, y = scale*den[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
  }
  
  m = scale*peak + 1:loop
  segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
  points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
  I = deci(CI, 2); o = deci(mode, 2)
  text(mode, 1:loop, paste0(I[,1], "        ", o, "         ", I[,2]), cex = .75, pos = 3, font = 2, xpd = NA)
  
  rownames <- if(is.na(labels)) paste0(np, ":") else paste0(labels, ":")                                               
  return(round(data.frame(mean = mean, mode = mode, median = median, sd = sd, lower = CI[,1], upper = CI[,2], eq.prob = BB, row.names = rownames), digits = digits))
}

              
#===================================================================================================================

              
prop.update <- function(n = 100, yes = 55, top = 5, scale = .1, lo = 0, hi = 1, a = 1.2, b = 1.2, dist.name = "dbeta", prior.scale = 1,
                        level = .95, show.prior = FALSE, tol = 1e5, col.depth = .55, labels = NULL, cex.lab = .9, xlab = NULL, ylab = NULL)
{
  UseMethod("prop.update")
}

prop.update.default <- function(n = 100, yes = 55, top = 5, scale = .1, lo = 0, hi = 1, a = 1.2, b = 1.2, dist.name = "dbeta", prior.scale = 1, 
                                level = .95, show.prior = FALSE, tol = 1e5, col.depth = .55, labels = NULL, cex.lab = .9, xlab = NULL, ylab = NULL){
  
  pri <- show.prior
  s <- round(yes)
  n <- round(n)  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))  
  is.v <- function(...) lengths(list(...)) > 1
  if(any(is.v(a, b, d))) stop("Error: Choose only 'one' prior knowledge base at a time.")
  if(any(yes > n)) stop("Error: 'yes' cannot be larger than 'n'.")  
  if(d == "dunif" & a == 0 & b == 1) {d = 'dbeta'; a <- b <- 1.0000001}
  if(d == "dbeta" & a == 1 & b == 1) a <- b <- 1.0000001;    
  if(tol < 1e4) stop("'tol' must be '10,000' or larger.")
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
  I <- eq(n, s) ; n <- I[[1]] ; s <- I[[2]]
  loop <- length(n) 
  
  props <- seq(0, 1, 1/tol)
  prx <- get(d)(props, a, b)*as.integer(props >= lo)*as.integer(props <= hi)
  pr <- tol * prx / sum(prx)
  
  graphics.off()                            
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  xlab <- if(is.null(xlab)) "Proportion" else xlab
  ylab <- if(is.null(ylab)) NA else ylab                            
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~props, ylim = c(0, top*loop), type = "n", yaxs = "i", ylab = ylab, xlab = xlab, font.lab = 2, axes = FALSE, mgp = c(2, .4, 0), main = if(pri) bquote(.(xlab)*" ~ "*.(if(lo > 0 || hi < 1) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))) else NA)
  axis(1, at = axTicks(1), labels = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
  
  if(!pri){  
    labels <- if(is.null(labels)) paste0("Study ", 1:loop) else labels  
    abline(h = 1:loop+1, col = 8, lty = 3)
    axis(2, at = 0:loop+1, labels = c("Base knowledge", labels), las = 1, font = 2, cex.axis = cex.lab, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
  }  
  polygon(x = c(lo, props, hi), y = prior.scale*c(0, pr, 0), col = adjustcolor(8, .8))
  
  I = hdi(x = props, y = pr, level = level)
  
  m = props[which.max(pr)]
  y = prior.scale*(pr[which.max(pr)])
  segments(I[1], 0, I[2], 0, lend = 1, lwd = 4, xpd = NA)
  points(m, 0, pch = 19, xpd = NA, cex = 1.4)
  segments(m, 0, m, y, lty = 3)
  q <- deci(I*1e2, 2) 
  o <- deci(m*1e2, 2)
  text(c(I[1], m, I[2]), 0, c(paste0(q[1], "%"), paste0(o, "%"), paste0(q[2], "%") ), pos = 3, cex = .8, font = 2, xpd = NA)
  
  if(!pri){
    for(i in 1:loop) {
      ps <- dbinom(s[i], n[i], props) * pr
      ps <- tol * ps / sum(ps)
      polygon(y = scale*ps+i+1, x = props, col = adjustcolor(i+1, col.depth), border = NA, xpd = NA)
      I = hdi(x = props, y = ps, level = level)
      m = props[which.max(ps)]
      q = deci(I*1e2 , 2); o = deci(m*1e2, 2)
      y = ps[which.max(ps)]*scale + (i+1)
      segments(I[1], i+1, I[2], i+1, lend = 1, lwd = 3, col = i +1)
      segments(m, i+1, m, y, lty = 3, xpd = NA)
      text(m, i+1, paste0(q[1], "%", "     ", o, "%", "     ", q[2], "%"), pos = 3, cex = .7, font = 2, xpd = NA)
      points(m, i+1, pch = 21, bg = "cyan", col = "magenta")
      
      pr <- ps
    }
  }
}

#=======================================================================================================================

d.update <- function(t = 3.35, n1 = 30, n2 = NA, top = 5, scale = .1, m = 0, s = 1, dist.name = "dnorm", prior.scale = 1, level = .95, show.prior = FALSE, lo = -2, hi = 2, tol = 1e4, margin = hi, col.depth = .55, labels = NULL, cex.lab = .9, xlab = NULL, ylab = NULL)
{
  UseMethod("d.update")
}

d.update.default <- function(t = 3.35, n1 = 30, n2 = NA, top = 5, scale = .1, m = 0, s = 1, dist.name = "dnorm", 
                             prior.scale = 1, level = .95, show.prior = FALSE, lo = -2, hi = 2, tol = 1e4, 
                             margin = hi, col.depth = .55, labels = NULL, cex.lab = .9, xlab = NULL, ylab = NULL){
  
  pri <- show.prior
  d <- dist.name
  if(is.infinite(lo)) lo <- -6
  if(is.infinite(hi)) hi <- 6
  if(tol < 1e4) stop("'tol' must be '10,000' or larger.")
  is.v <- function(...) lengths(list(...)) > 1
  if(any(is.v(m, s, d))) stop("Choose only 'one' prior knowledge base at a time.")
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
  I <- eq(t, n1, n2) 
  t <- I[[1]]  
  n1 <- I[[2]]  
  n2 <- I[[3]] 
  loop <- length(t) 
  
  ds <- seq(-6, 6, 1/tol)
  prx <- get(d)(ds, m, s)*as.integer(ds >= lo)*as.integer(ds <= hi)
  pr <- tol * prx / sum(prx)
  
  graphics.off()                            
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  xlab <- if(is.null(xlab)) bquote(bold("Effect Size"~ (delta))) else xlab
  ylab <- if(is.null(ylab)) NA else ylab                            
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~ds, ylim = c(0, top*loop), xlim = c(-margin, margin), type = "n", xaxs = "i", yaxs = "i", ylab = ylab, xlab = xlab, font.lab = 2, mgp = c(2, .4, 0), main = if(pri) bquote(.(xlab)*" ~ "*.(if(lo > -Inf || hi < Inf) "truncated-")*.(substring(d, 2))(.(round(m, 2)), .(round(s, 2)))) else NA, yaxt = "n", bty = "n")
  
  if(!pri){
    labels <- if(is.null(labels)) paste0("Study ", 1:loop) else labels  
    abline(h = 1:loop+1, col = 8, lty = 3)
    axis(2, at = 0:loop+1, labels = c("Base knowledge", labels), las = 1, font = 2, cex.axis = cex.lab, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
  }  
  
  polygon(x = c(-margin, ds, margin), y = prior.scale*c(0, pr, 0), col = adjustcolor(8, .8))
  
  # I = hdi(x = ds, y = pr, level = level)
  
  if(d != "dunif"){
    mode = ds[which.max(pr)]
    y = prior.scale*(pr[which.max(pr)])
    # segments(I[1], 0, I[2], 0, lend = 1, lwd = 4, xpd = NA)
    points(mode, 0, pch = 19, xpd = NA, cex = 1.4)
    segments(mode, 0, mode, y, lty = 3)
    # text(c(.85*I[1], mode, I[2]), 0, paste0(round(c(I[1], mode, I[2]), 3)), pos = 3, cex = .8, font = 2, xpd = NA)
  }
  
  N = ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2) 
  
  options(warn = -1)
  if(!pri){
    for(i in 1:loop) {
      
      ps <- dt(t[i], df[i], ds*sqrt(N[i])) * pr
      ps <- tol * ps / sum(ps)
      polygon(y = scale*ps+i+1, x = ds, col = adjustcolor(i+1, col.depth), border = NA, xpd = NA)
      I = hdi(x = ds, y = ps, level = level)
      mode = ds[which.max(ps)]
      q = deci(I, 3); o = deci(mode, 3)
      y = ps[which.max(ps)]*scale + (i+1)
      segments(I[1], i+1, I[2], i+1, lend = 1, lwd = 3, col = i +1)
      segments(mode, i+1, mode, y, lty = 3, xpd = NA)
      text(mode, i+1, paste0(q[1], "     ", o, "     ", q[2]), pos = 3, cex = .7, font = 2, xpd = NA)
      points(mode, i+1, pch = 21, bg = "cyan", col = "magenta")
      
      pr <- ps
    }
  }
}


#==================================================================================================================

peta.update <- function(f = 50, N = 120, df1 = 3, df2 = 116, top = 5, scale = .1, a = 2, b = 2, lo = 0, hi = 1, dist.name = "dbeta", prior.scale = 1, level = .95, show.prior = FALSE, tol = 1e5, col.depth = .55, labels = NULL, cex.lab = .9, xlab = NULL, ylab = NULL)
{
  UseMethod("peta.update")
}

peta.update.default <- function(f = 50, N = 120, df1 = 3, df2 = 116, top = 5, scale = .1, a = 2, b = 2, lo = 0, hi = 1, dist.name = "dbeta", prior.scale = 1, level = .95, show.prior = FALSE, tol = 1e5, col.depth = .55, labels = NULL, cex.lab = .9, xlab = NULL, ylab = NULL){
  
  pri <- show.prior
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  if(hi == 1) hi <- .9999999 ;
  if(tol < 1e4) stop("'tol' must be '10,000' or larger.")
  is.v <- function(...) lengths(list(...)) > 1
  if(any(is.v(a, b, d))) stop("Error: Choose only 'one' prior knowledge base at a time.")
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
  I <- eq(f, N, df1, df2) ; f <- I[[1]] ; N <- I[[2]] ; df1 <- I[[3]] ; df2 <- I[[4]] ;
  loop <- length(f) 
  
  peta <- seq(0, .9999999, 1/tol)
  prx <- get(d)(peta, a, b)*as.integer(peta >= lo)*as.integer(peta <= hi)
  pr <- tol * prx / sum(prx)
  
  graphics.off()                            
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  xlab <- if(is.null(xlab)) bquote(bold((eta[p]^2))) else xlab
  ylab <- if(is.null(ylab)) NA else ylab                            
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~peta, ylim = c(0, top*loop), type = "n", yaxs = "i", ylab = ylab, xlab = xlab, font.lab = 2, axes = FALSE, mgp = c(2, .4, 0), main = if(pri) bquote(.(xlab)*" ~ "*.(if(lo > 0 || hi < .9999999) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))) else NA)
  axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
  
  if(!pri){
    labels <- if(is.null(labels)) paste0("Study ", 1:loop) else labels  
    abline(h = 1:loop+1, col = 8, lty = 3)
    axis(2, at = 0:loop+1, labels = c("Base knowledge", labels), las = 1, font = 2, cex.axis = cex.lab, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
  }  
  polygon(x = c(lo, peta, hi), y = prior.scale*c(0, pr, 0), col = adjustcolor(8, .8))
  
  I = hdi(x = peta, y = pr, level = level)
  
  m = peta[which.max(pr)]
  y = prior.scale*(pr[which.max(pr)])
  segments(I[1], 0, I[2], 0, lend = 1, lwd = 4, xpd = NA)
  points(m, 0, pch = 19, xpd = NA, cex = 1.4)
  segments(m, 0, m, y, lty = 3)
  text(c(.85*I[1], m, I[2]), 0, paste0(round(c(I[1], m, I[2])*1e2, 4), "%"), pos = 3, cex = .8, font = 2, xpd = NA)
  
  if(!pri){
    for(i in 1:loop) {
      
      ps <- df(f[i], df1[i], df2[i], (peta * N[i]) / (1 - peta) ) * pr
      ps <- tol * ps / sum(ps)
      m = peta[which.max(ps)]
      polygon(y = c(i+1, scale*ps+i+1, i+1), x = c(lo, peta, hi), col = adjustcolor(i+1, col.depth), border = NA, xpd = NA)
      I = hdi(x = peta, y = ps, level = level)
      
      q = deci(I*1e2 , 2); 
      o = deci(m*1e2, 2)
      y = ps[which.max(ps)]*scale + (i+1)
      segments(I[1], i+1, I[2], i+1, lend = 1, lwd = 3, col = i +1)
      segments(m, i+1, m, y, lty = 3, xpd = NA)
      text(m, i+1, paste0(q[1], "%", "     ", o, "%", "     ", q[2], "%"), pos = 3, cex = .7, font = 2)
      points(m, i+1, pch = 21, bg = "cyan", col = "magenta")
      
      pr <- ps
    }
  }
}

#===================================================================================================================

d.eq.test <- function(t, n1, n2 = NA, m = 0, s = 1, dist.name = "dnorm", dL = -.1, dU = .1, lo = -Inf, hi = Inf)
{
  UseMethod("d.eq.test")
}

d.eq.test.default <- function(t, n1, n2 = NA, m = 0, s = 1, dist.name = "dnorm", dL = -.1, dU = .1, lo = -Inf, hi = Inf){
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  
  if(any(lengths(list(get(formalArgs(d.eq.test))))) > 1) stop("Error: Only 'one' equivalence testing at a time is allowed.")
  if(dL >= dU) stop("Your Upper value must be larger than your Lower value")
  
  if(abs(dL) != abs(dU)) message("\n\tYou have an \"Unequal Equivalence Bound\", thus we can't provide an extra\n\t function showing the effect of choosing various Unequal bounds.")
  
  decimal <- function(x, k){   
    if(is.character(x)){
      return(x)
    }else{
      format(round(x, k), nsmall = k, scientific =
               ifelse(x >= 1e5 || x <= -1e5 || x <= 1e-5 & x >= -1e-5, TRUE, FALSE) )
    }
  }
  
  options(warn = -1)
  
  N <- ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)
  
  p <- function(x) get(d)(x, m, s)*as.integer(x >= lo)*as.integer(x <= hi)
  prior <- function(x) p(x)/integrate(p, lo, hi)[[1]]  
  likelihood <- function(x) dt(t, df, x*sqrt(N))
  k <- integrate(function(x) prior(x)*likelihood(x), lo, hi)[[1]]
  posterior <- function(x) prior(x)*likelihood(x)/k
  
  mean <- integrate(function(x) x*posterior(x), lo, hi)[[1]]
  sd <- sqrt(integrate(function(x) x^2*posterior(x), lo, hi)[[1]] - mean^2)
  
  x.min.1 <- mean - 9 * sd
  x.max.1 <- mean + 9 * sd
  
  ## The dL and dU may be different from x.min.1 and x.max.1 respectively, if so, adjust accordingly.
  x.min <- if(dL < x.min.1) { 1.05*dL } else { x.min.1 }
  x.max <- if(dU > x.max.1) { 1.05*dU } else { x.max.1 }
  
  CI <- HDI(posterior, x.min, x.max)
  
  cdf <- Vectorize(function(q){
    integrate(posterior, lo, q)[[1]]
  })
  
  #  inv.cdf <- Vectorize(function(p){
  #    uniroot(function(q)cdf(q) - p, c(x.min, x.max))[[1]]  # Not implemented #
  #  })
  
  mode <- optimize(posterior, c(x.min, x.max), maximum = TRUE)[[1]]
  peak <- posterior(mode)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mar = c(.1, 4.1, 3.1, 2.1), mfcol = c(2, 1))
  
  h = curve(posterior, from = x.min, to = x.max, las = 1, type = "n",
            xlab = NA, ylab = NA, bty = "n", ylim = c(0, 1.1*peak), 
            xaxt = "n", yaxt = "n", mgp = c(2, .5, 0), n = 1e3)
  
  X <- h$x >= CI[1] &  h$x <= CI[2]
  
  low.extreme <- par('usr')[3]
  
  polygon(c(CI[1], h$x[X], CI[2]), c(low.extreme, h$y[X], low.extreme), col = rgb(1, 1, 0, .5), border = NA)
  
  segments(mode, low.extreme, mode, peak, lty = 3)
  
  text(mode, peak/2, decimal(mode, 2), srt = 90, pos = 3, font = 2)
  
  lines(h, lwd = 2)
  
  segments(CI[1], low.extreme, CI[2], low.extreme, col = 2, lend = 1, lwd = 40)
  
  segments(dL, low.extreme, dU, low.extreme, col = adjustcolor(3, .5), lend = 1, lwd = 40)
  
  points(mode, low.extreme/5, pch = 21, col = 0, bg = 0, cex = 1.5)
  
  axis(side = 1, at = decimal(seq(x.min, x.max, length.out = 7), 2), mgp = c(2, .5, 0))
  axis(side = 1, at = 0, mgp = c(3, 1.1, 0), col = 0, col.axis = "magenta", tick = FALSE, line = - 1.4, cex.axis = 1.4, font = 2)
  
  mtext(side = 1, bquote(bold("Population Effect Size"~(delta))), line = 2, cex = .95)
  
  y1 = y2 = 1.02*peak
  x.text = (dL+dU)/2
  y.text = 1.05*peak
  
  segments(c(dL, dU), rep(low.extreme, 2), c(dL, dU), c(y1, y2), col = 'green2', lend = 1, lty = 2)
  
  segments(c(dL, dU), c(y1, y2), rep(x.text, 2), rep(y.text*1.023, 2), lwd = 2, col = 'magenta')
  
  text(x.text, y.text, "Practically Equivalent to ZERO", font = 2, pos = 3, col = 'darkgreen', cex = .65, xpd = TRUE)
  
  points(c(dL, dU), c(y1, y2), pch = 21, col = 'green3', bg = 'green3', cex = 1.1)
  
  ## How much is it probable that the equivalence be true in population:
  
  a = cdf(dL)
  b = cdf(dU)
  
  Post.in.ROPE.Y = (b - a)
  Post.in.ROPE.X = (dU - dL) / 2
  
  BB = decimal(Post.in.ROPE.Y*1e2, 2)
  
  title(main = paste0("There is ", "''", BB, "%", "''", " probability that TRUE effect size is equivalent to ZERO"), cex.main = .8)
  
  legend("topleft", legend = paste0("95% HDI: [",decimal(CI[1], 2), ", ", decimal(CI[2], 2),"]"),
         bty = "n", inset = c(-.035,.1), text.font = 2, text.col = 'red4', cex = .8)
  
  if(CI[1] > dU || CI[2] < dL) {
    
    legend("topright", "NOT Practically equivalent to \"0\" ", bty = 'n', inset = c(-.01, .1), cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
    
  } else
    
    if(CI[1] > dL & CI[2] < dU) {
      
      legend("topright", "Practically equivalent to \"0\" ", bty = 'n', inset = c(-.01, .1), cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
      
    } else  {
      
      legend("topright", "No decision can be made ", bty = 'n', inset = c(-.01, .1), cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
    }
  
  #########################################################################
  ## How choice of ROPE can affect porortion of posterior that ROPE covers:
  #########################################################################
  
  par(mar = c(3.1, 4.1, 6.1, 2.1), mgp = c(2.5, .5, 0))
  
  eq.low = ifelse(abs(dL) <= .3, 4, 2)*( - ((dU - dL) / 2) )
  eq.upp = ifelse(abs(dL) <= .3, 4, 2)*(   ((dU - dL) / 2) )
  
  L = seq(eq.low, 0, length.out = 1e2)
  U = seq(eq.upp, 0, length.out = 1e2)
  
  aa = cdf(L)
  bb = cdf(U)
  
  Eq = (bb - aa)  # porortion of posterior that ROPE covers
  half = (U - L)/2
  
  plot(half, Eq, type = ifelse(abs(dL) == abs(dU), 'l' ,'n'), lwd = 3, col = 'red4', axes = FALSE,
       xlab = NA, ylab = paste0("%",'Posterior in ROPE'), font.lab = 2, cex.lab = .8)
  
  mtext(side = 1, "Half of ROPE", font = 2, line = 1.5, cex = .9)
  
  axis(1, at = decimal(seq(0, eq.upp[1], length.out = 7), 2), las = 1)
  axis(2, at = seq(0, Eq[1], length.out = 5),
       labels = paste0(1e2*round(seq(0, Eq[1], length.out = 5),
                                 2), "%"), las = 1)
  
  pars = par('usr')
  
  rect(pars[1], pars[3], pars[2], pars[4], col = adjustcolor("grey", .1), border = NA)
  
  rect(pars[1], pars[3], Post.in.ROPE.X, Post.in.ROPE.Y,
       col = adjustcolor("yellow",
                         alpha.f = ifelse(Post.in.ROPE.Y  <= .2, .1,
                                          ifelse(Post.in.ROPE.Y > .2 & Post.in.ROPE.Y <= .3, .15,
                                                 ifelse(Post.in.ROPE.Y > .3 & Post.in.ROPE.Y <= .4, .2, .3)))),
       lty = 2 )
  
  points(Post.in.ROPE.X, Post.in.ROPE.Y,
         pch = 21, cex = 2, bg = 'green')
  
  box()  
}

                       
#======================================================================================================================                       
                      
                       
R <- function(fit)
{
  UseMethod("R")
}                        
                       
R.default <- function(fit){
  
if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")
    
      y <- rstanarm::get_y(fit)
  ypred <- rstanarm::posterior_linpred(fit, transform = TRUE)
 
  if(family(fit)$family == "binomial" && ncol(y) == 2) {
    trials <- rowSums(y)
    y <- y[, 1]
    ypred <- ypred %*% diag(trials)
  }
  e <- -1 * sweep(ypred, 2, y)
  var.ypred <- apply(ypred, 1, var)
  var.e <- apply(e, 1, var)
  var.ypred / (var.ypred + var.e)

}                       
                       
#==============================================================================================================                       


R2.bayes <- function(..., scale = .02, bottom = 1, top = 1.1, margin = 5, legend = "topleft", level = .95, eq.lo = 0, eq.hi = .1, digits = 6, labels = NA, cex.lab = .8, xlab = NA, ylab = NA)
{
  UseMethod("R2.bayes")
}                        
                                           
                       
R2.bayes.default <- function(..., scale = .02, bottom = 1, top = 1.1, margin = 5, legend = "topleft", level = .95, eq.lo = 0, eq.hi = .1, digits = 6, labels = NA, cex.lab = .8, xlab = NA, ylab = NA)
{

if(!(all(sapply(list(...), inherits, "stanreg")))) stop("Error: all '...' must be fitted models from package 'rstanarm's 'stan_glm()'.")  
    
leg <- if(is.character(legend)) legend else deparse(substitute(legend))   
Rs <- lapply(list(...), R)
loop <- length(Rs)

deci <- function(x, k = 3) format(round(x, k), nsmall = k) 

d <- list()
I <- matrix(NA, loop, 2)
mean <- numeric(loop)
median <- numeric(loop)
sd <- numeric(loop)
mad <- numeric(loop)
mode <- numeric(loop)
peak <- numeric(loop)
from <- numeric(loop)                  
to <- numeric(loop)
eq.prob <- numeric(loop)

for(i in 1:loop){
    
d[[i]] <- density(Rs[[i]], adjust = 2, n = 1e3)
I[i,] <- hdir(Rs[[i]], level = level)
mean[i] <- mean(Rs[[i]])
median[i] <- median(Rs[[i]])
sd[i] <-  sd(Rs[[i]])
mad[i] <- mad(Rs[[i]])
mode[i] <- d[[i]]$x[which.max(d[[i]]$y)]
peak[i] <- d[[i]]$y[which.max(d[[i]]$y)]
eq.prob[i] <- mean(eq.lo <= Rs[[i]] & Rs[[i]] <= eq.hi)
from[i] <- mean[i] - margin * sd[i]
to[i] <- mean[i] + margin * sd[i]
}

a = if(min(from) >= 0) min(from) else 0
b = if(max(to) <= 1) max(to) else 1

xlab <- if(is.na(xlab)) bquote(bold("Regression Coefficient " (R^2))) else xlab
ylab <- if(is.na(ylab)) NA else ylab    
graphics.off()    
plot(1, loop, type = "n", xlim = c(a, b), ylim = c(bottom*1, top*loop), ylab = ylab, xaxt = "n", yaxt = "n", xlab = xlab, mgp = c(2, .3, 0))
abline(h = 1:loop, col = 8, lty = 3)

for(i in 1:loop){
polygon(x = d[[i]]$x, y = scale*d[[i]]$y + i, col = adjustcolor(i, .55), border = NA, xpd = NA)
}
lab <- if(any(is.na(labels))) paste0("Model ", 1:loop) else labels
axis(1, at = seq(a, b, length.out = 4), labels = paste0(round(seq(a, b, length.out = 4), 4)*1e2, "%"), mgp = c(2, .5, 0))
axis(2, at = 1:loop, labels = lab, font = 2, las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .3, 0))

legend(x = leg, legend = rev(paste0("Model ", loop:1)), pch = 22, title = "Models ", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, adj = c(0, .3))
segments(I[, 1], 1:loop, I[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
box()

m = scale*peak + 1:loop
segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.5, col = "magenta", xpd = NA)

q = deci(I*1e2 , 2); o = deci(mode*1e2, 2)
text(mode, 1:loop, paste0(q[,1], "%", "    ", o, "%", "    ", q[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)

rownames <- if(is.na(labels)) paste0("Model-", 1:loop, " posterior:") else paste0(1:loop, " ", labels, " posterior:")    
round(data.frame(mode = mode, mean = mean, sd = sd, mad = mad, lower = I[,1], upper = I[,2], coverage = level, eq.prob = eq.prob, row.names = rownames), digits = digits)
}


#=======================================================================
                         
                       
lm.sample2 <- function(fit, n = 1e4, no.names = TRUE)
{
  UseMethod("lm.sample2")
}                       
                
                       
lm.sample2.default <- function(fit, n = 1e4, no.names = TRUE){
 
if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")
    
  output <- as.data.frame(rmvnorm(n = n, mean = c(coef(fit), sigma(fit)), sigma = cov(as.matrix(fit))))
  
  if(no.names){
    for(i in 1:ncol(output)){
     if(colnames(output)[i] == "(Intercept)"){
        colnames(output)[i] <- "Intercept"
      }
      if(colnames(output)[i] == paste0("V", ncol(output))){
        colnames(output)[i] <- "Sigma"
      }
    }
  }
output
}
                     
                       
#======================================================================================
                       
                       
lm.sample <- function(fit, no.names = TRUE)
{
  UseMethod("lm.sample")
}                       
                
                       
lm.sample.default <- function(fit, no.names = TRUE){
 
if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")
    
  output <- as.data.frame(fit)
  
  if(no.names){
    for(i in 1:ncol(output)){
      if(colnames(output)[i] == "(Intercept)"){
        colnames(output)[i] <- "Intercept"
      }
    }
  }
  output
}
                       
                       
#======================================================================================

 
predict.mean <- function(fit, predi, scale = .5, level = .95, col.hump = "cyan", integer = FALSE, xlab = NA, ...)
{
  UseMethod("predict.mean")
} 
       
                       
predict.mean.default <- function(fit, predi, scale = .5, level = .95, col.hump = "cyan", integer = FALSE, xlab = NA, ...){
  
if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")  
if(length(coef(fit)) > 2) stop("Error: 'fit' must contain only 'one' predictor.")  
if(!(predi %in% fit$model[, 2])) message("WARNING: \nYou have no data on ", dQuote(names(fit$model)[2]), " = ", predi, ". The prediction is not highly reliable.")

mus_at_xi <- rstanarm::posterior_linpred(fit, newdata = setNames(data.frame(tmp = predi), names(fit$model)[2]))

d <- density(mus_at_xi, adjust = 2, n = 1e3)

xlab = if(is.na(xlab)) bquote(bold(bolditalic(p)*(Ave.*.(names(fit$model)[1])[i] *" | "* .(names(fit$model)[2])[i] == .(predi)))) else xlab   

graphics.off()    
plot(d, type = "n", ylab = NA, main = NA, yaxt = "n", bty = "n", las = 1, zero.line = FALSE, yaxs = "i",
     xlab = xlab, ...)

  I <- hdir(mus_at_xi, level = level)
med <- mean(mus_at_xi)
peak <- approx(d$x, d$y, xout = med)[[2]]*scale

polygon(d$x, scale*d$y, col = adjustcolor(col.hump, .35), border = NA)
segments(med, 0, med, peak, lty = 3)

segments(I[1], 0, I[2], 0, lend = 1, lwd = 6, col = 'magenta', xpd = NA)
points(med, 0, pch = 21, bg = "cyan", col = 'magenta', cex = 2, xpd = NA)
text(c(I, med), 0, if(integer) round(c(I, med)) else round(c(I, med), 2), pos = 3, font = 2)
}                       
                       
                      
#======================================================================================                       

                       
predict.case <- function(fit, predi, scale = .5, level = .95, col.hump = "gray", integer = FALSE, xlab = NA, ...)
{
  UseMethod("predict.case")
} 
                       
                       
predict.case.default <- function(fit, predi, scale = .5, level = .95, col.hump = "gray", integer = FALSE, xlab = NA, ...){
  
  if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")  
  if(length(coef(fit)) > 2) stop("Error: 'fit' must contain only 'one' predictor.")  
  if(!(predi %in% fit$model[, 2])) message("WARNING: \nYou have no data on ", dQuote(names(fit$model)[2]), " = ", predi, ". The prediction is not highly reliable.")
  
case_at_xi <- rstanarm::posterior_predict(fit, newdata = setNames(data.frame(tmp = predi), names(fit$model)[2]))
  
  d <- density(case_at_xi, adjust = 2, n = 1e3)
    
  xlab = if(is.na(xlab)) bquote(bold(bolditalic(p)*(.(names(fit$model)[1])[i] *" | "* .(names(fit$model)[2])[i] == .(predi)))) else xlab  
  graphics.off()
  plot(d, type = "n", ylab = NA, main = NA, yaxt = "n", bty = "n", las = 1, zero.line = FALSE, yaxs = "i",
       xlab = xlab, ...)
  
  I <- hdir(case_at_xi, level = level)
  med <- mean(case_at_xi)
  peak <- approx(d$x, d$y, xout = med)[[2]]*scale
  
  polygon(d$x, scale*d$y, col = col.hump, border = NA)
  segments(med, 0, med, peak, lty = 3)
  
  segments(I[1], 0, I[2], 0, lend = 1, lwd = 6, xpd = NA)
  points(med, 0, pch = 21, bg = "cyan", col = 'magenta', cex = 2, xpd = NA)
  text(c(I, med), 0, if(integer) round(c(I, med)) else round(c(I, med), 2), pos = 3, font = 2)
}  
                       
                       
#======================================================================================                       

                       
predict.bayes <- function(fit, xlab = NA, ylab = NA, level = .95, line.int = TRUE, pred.int = TRUE, pt.cex = 1, pt.col = 4, col.depth = .3, col.line = "cyan", col.pred = "gray", col.reg = "cyan", ...)
{
  UseMethod("predict.bayes")
} 

                       
predict.bayes.default <- function(fit, xlab = NA, ylab = NA, level = .95, line.int = TRUE, pred.int = TRUE, pt.cex = 1, pt.col = 4, col.depth = .3, col.line = "cyan", col.pred = "gray", col.reg = "cyan", ...){
  
  if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")  
  if(length(coef(fit)) > 2) stop("Error: 'fit' must contain only 'one' predictor. Consider using 'counter.plot()' for multiple predictors.")
  
  pred <- fit$model[, 2]
  dep <- fit$model[, 1]  
  
  pred_lin <- rstanarm::posterior_predict(fit, transform = TRUE)
  pred_lin2 <- rstanarm::posterior_linpred(fit, transform = TRUE)
  graphics.off()
  plot(dep ~ pred, xlab = ifelse(is.na(xlab), names(fit$model)[2], xlab), ylab = ifelse(is.na(ylab), names(fit$model)[1], ylab), type = "n", las = 1, ...)
  
  loop <- length(pred)
  
  I <- matrix(NA, loop, 2)
  for(i in 1:loop){
    I[i,] = hdir(pred_lin[,i], level = level)
  }
  
  OK <- I[,1] < dep & dep < I[,2]
  
  points(dep ~ pred, pch = 19, col = ifelse(OK, adjustcolor(pt.col, .55), 2), cex = pt.cex)
  
  x <- sort(pred)
  o <- order(pred)
  
  if(pred.int){   
    
    y <- I[,1][o]
    z <- I[,2][o]
    
    polygon(c(rev(x), x), c(rev(z), y), col = adjustcolor(col.pred, col.depth), border = NA)
  }
  
  I2 <- matrix(NA, loop, 2)
  for(i in 1:loop){
    I2[i,] = hdir(pred_lin2[,i], level = level)
  }
  
  if(line.int){
    
    y <- I2[,1][o]
    z <- I2[,2][o]
    
    polygon(c(rev(x), x), c(rev(z), y), col = adjustcolor(col.line, col.depth), border = NA)
  }
  
  abline(fit, col = col.reg, lwd = 2, lend = 1)
  box()    
}                                              
                       
#==================================================================================
                  
                       
compare.R2 <- function(..., how = c("two.one", "one.two"), scale = .02, bottom = 1, top = 1.1, margin = 5, legend = "topleft", level = .95, eq.level = "2.5%", digits = 6, labels = NA, cex.lab = .7, xlab = NA, ylab = NA, main = NA)
{
  UseMethod("compare.R2")
} 

                       
compare.R2.default <- function(..., how = c("two.one", "one.two"), scale = .02, bottom = 1, top = 1.1, margin = 5, legend = "topleft", level = .95, eq.level = "2.5%", digits = 6, labels = NA, cex.lab = .7, xlab = NA, ylab = NA, main = NA){
  
if(!(all(sapply(list(...), inherits, "stanreg")))) stop("Error: all '...' must be models from package 'rstanarm's 'stan_glm()'.")  
  
  eq.bound <- if(is.character(eq.level)) as.numeric(substr(eq.level, 1, nchar(eq.level)-1)) / 1e2 else eq.level
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  Rs <- lapply(list(...), R)
  
  if(length(Rs) < 2) stop("Error: You need to have a least '2' fitted models (from 'rstanarm' package) for comparison!")
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  how <- match.arg(how)
  
  delta <- switch(how,
                  one.two = function(x) x[[1]] - x[[2]], 
                  two.one = function(x) x[[2]] - x[[1]])
  
  ps <- combn(Rs, 2, FUN = delta)
  
  loop <- ncol(ps)
   
  CI <- matrix(NA, loop, 2)
  den <- list()
  mode <- numeric(loop)
  peak <- numeric(loop)
  mean <- numeric(loop)
  median <- numeric(loop)                  
  sd <- numeric(loop)
  eq.prob <- numeric(loop)
  from <- numeric(loop)                  
  to <- numeric(loop)
  
                  
  for(i in 1:loop){
      
    CI[i,] <- hdir(ps[, i], level = level)
    den[[i]] <- density(ps[, i], adjust = 2, n = 1e3)
    mode[i] <- den[[i]]$x[which.max(den[[i]]$y)]
    peak[i] <- den[[i]]$y[which.max(den[[i]]$y)]
    mean[i] <- mean(ps[, i])
    median[i] <- median(ps[, i])
    eq.prob[i] <- mean(abs(ps[, i]) <= eq.bound)
    sd[i] <- sd(ps[, i])
    from[i] <- mean[i] - margin * sd[i]
    to[i] <- mean[i] + margin * sd[i]
  }
  
  np <- combn(seq_along(Rs), 2, FUN = function(x){if(how == "one.two") paste0('Model ', x[1], ' - Model ', x[2]) else paste0('Model ', x[2], ' - Model ', x[1])})
  
  
  a <- if(min(from) >= -1) min(from) else -1
  b <- if(max(to) <= 1) max(to) else 1
  
  graphics.off()    
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  xlab <- if(is.na(xlab)) bquote(bold(Delta~R^2~("Model Comparison"))) else xlab
  ylab <- if(is.na(ylab)) NA else ylab
  main <- if(is.na(main)) NA else main  
      
  par(mar = c(5.1, 6.1, 4.1, 2.1))
  labels <- if(is.na(labels)) np else labels
  plot(1, loop, type = "n", xlim = c(a, b), ylim = c(bottom*1, top*loop), ylab = ylab, xaxt = "n", yaxt = "n", xlab = xlab, mgp = c(2, .3, 0), main = main)
  abline(h = 1:loop, col = 8, lty = 3)
  axis(1, at = seq(a, b, length.out = 4), labels = paste0(round(seq(a, b, length.out = 4), 4)*1e2, "%"), mgp = c(2, .5, 0))
  axis(2, at = 1:loop, labels = labels, font = 2, las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .3, 0))
  legend(x = leg, legend = rev(paste0(np)), pch = 22, title = "Comparisons", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5)
  segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
  box()
  
  for(i in 1:loop){
    polygon(x = den[[i]]$x, y = scale*den[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
  }
  
  m = scale*peak + 1:loop
  segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
  points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.4, col = "magenta", xpd = NA)
  q = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
  text(mode, 1:loop, paste0(q[,1], "%", "         ", o, "%", "         ", q[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
  
  return(round(data.frame(mean = mean, mode = mode, median = median, sd = sd, lower = CI[,1], upper = CI[,2], coverage = level, eq.prob = eq.prob, row.names = paste0(np, ":")), digits = digits))
}

                       
#=================================================================================== 

reg <- function(x, y, col) abline(lm(y~x), col = "cyan", lwd = 2) 

              
panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                      cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...){
                      points(x, y, pch = pch, col = col, bg = bg, cex = cex)
                      ok <- is.finite(x) & is.finite(y)
                      if (any(ok)) reg(x[ok], y[ok], col.smooth)
                      }
  
              
get.asp <- function(){
  uy <- diff(grconvertY(1:2, "user", "inches"))
  ux <- diff(grconvertX(1:2, "user", "inches"))
  uy/ux
}

              
panel.cor <- function(x, y, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  asp <- get.asp()
  r <- (cor(x, y))
  txt <- paste0("r = ", round(r, 2))
  text(.5, .5, txt, cex = 1.3, font = 4, srt = 180/pi*atan(r*asp))
}

              
panel.hist <- function(x, col.hist = "khaki", ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = col.hist, ...)
}
              

panel.dens <- function(x, col.dens = adjustcolor(2, .2), ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  hist(x, plot = FALSE)
  d <- density(x); y <- d$y/max(d$y)
  polygon(x = d$x, y = y, col = col.dens, border = NA, ...)
}
              
              
#=====================================================================================
              
              
lm.check.plot <- function(data, pch = 19, cex = 1.6, col = adjustcolor("magenta", .05), 
                          gap = .15, panel = panel.lm, lower.panel = panel.cor, 
                          cex.labels = 1.3, font.labels = 2, font = 2, diag.panel = panel.hist, 
                          mgp = c(2, .6, 0), las = 1, ...)
{
  UseMethod("lm.check.plot")
}


lm.check.plot.default <- function(data, pch = 19, cex = 1.6, col = adjustcolor("magenta", .05), 
                          gap = .15, panel = panel.lm, lower.panel = panel.cor, 
                          cex.labels = 1.3, font.labels = 2, font = 2, diag.panel = panel.hist, 
                          mgp = c(2, .6, 0), las = 1, ...){ 

pairs(data, pch = pch, cex = cex, col = col, gap = gap, panel = panel, 
      cex.labels = cex.labels, font.labels = font.labels, lower.panel = lower.panel, font = font, 
      diag.panel = diag.panel, mgp = mgp, las = las, ...)
 
}
              
             
#===================================================================================
     
              
lm.post.plot <- function(fit, pch = 19, cex = 1.6, col = adjustcolor(4, .3), 
                          gap = .15, panel = panel.lm, lower.panel = panel.cor, diag.panel = panel.dens,
                          cex.labels = 1.3, font.labels = 2, font = 2, mgp = c(2, .6, 0), las = 1, ...)
{
  UseMethod("lm.post.plot")
}


lm.post.plot.default <- function(fit, pch = 19, cex = 1.6, col = adjustcolor(4, .3), 
                         gap = .15, panel = panel.lm, lower.panel = panel.cor, diag.panel = panel.dens,
                         cex.labels = 1.3, font.labels = 2, font = 2, mgp = c(2, .6, 0), las = 1, ...){

if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")

post <- lm.sample(fit)

pairs(post, pch = pch, cex = cex, col = col, gap = gap, panel = panel, 
      cex.labels = cex.labels, font.labels = font.labels, lower.panel = lower.panel, 
      diag.panel = diag.panel, font = font, mgp = mgp, las = las, ...)
}

                           
#===================================================================================
              
              
type.sm <- function(d = .1, obs.d = .6, n1 = 20, n2 = NA, digits = 6)
{
  UseMethod("type.sm")
}

type.sm.default <- function(d = .1, obs.d = .6, n1 = 20, n2 = NA, digits = 6){
  
  graphics.off()  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mfrow = c(2, 1), mgp = c(2, .5, 0), mar = c(4, 4, 3, 2), xpd = TRUE)  
  
     N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  d.SE = 1/sqrt(N) ; ncp = d*sqrt(N)
  
 min.d = qt(1e-4, df)*d.SE  ;  max.d = qt(0.9999, df)*d.SE
  
`d|H0` = curve( dt(x/d.SE, df)/d.SE, min.d, max.d, n = 1e4, xlab = "Effect Size", 
                  ylab = NA, font = 2, font.lab = 2, type = "n", yaxt = "n", bty = "n",
                  cex.axis = 1, cex.lab = 1, yaxs = "i")
  
    CI = qt(c(.025, .975), df)*d.SE
  
     x = seq(min.d, CI[1], l = 1e4) ;  y = dt(x /d.SE, df)/d.SE
    xx = seq(max.d, CI[2], l = 1e4) ; yy = dt(xx/d.SE, df)/d.SE
  
  polygon(c(min.d,  x, CI[1]), c( 0,  y, 0), col = 2, border = NA)  
  polygon(c(max.d, xx, CI[2]), c(0, yy, 0), col = 2, border = NA)   
  
  lines(`d|H0`, lwd = 2)
  
  points(obs.d, 0, pch = 23, bg = 3, cex = 1.4, xpd = TRUE)
  
  legend("topright", "Observed \nS.S. Effect", pch = 23, pt.bg = 3, pt.cex = 1.2, bty = "n", text.font = 2, adj = c(0, .3))   
  abline(v = 0, col = 2, xpd = FALSE) 
  
  par(mar = c(5, 4, 1, 2))
  
`d|H1` = curve( dt(x/d.SE, df, ncp)/d.SE, min.d, max.d, n = 1e4, xlab = "Effect Size", 
                  ylab = NA, font = 2, font.lab = 2, yaxt = "n", bty = "n",
                  cex.axis = 1, cex.lab = 1, yaxs = "i", ty = "n")
  
     x = seq(min.d, CI[1], l = 1e4)   ;  y = dt(x /d.SE, df, ncp)/d.SE
    xx = seq(max.d, CI[2], l = 1e4)   ; yy = dt(xx/d.SE, df, ncp)/d.SE 
  
  polygon(c(min.d,  x, CI[1]), c( 0,  y, 0), col = 2, border = NA)  
  polygon(c(max.d, xx, CI[2]), c(0, yy, 0), col = 2, border = NA)  
  
  lines(`d|H1`, lwd = 2)
  
  axis(1, at = d, col = 4, col.axis = 4, font = 2)
  points(obs.d, 0, pch = 23, bg = 3, cex = 1.4, xpd = TRUE)
  abline(v = d, col = 4, xpd = FALSE)
  
  segments(c(CI[1], CI[2]), 0, c(CI[1], CI[2]), 20, lty = 2, col = 2, xpd = NA)
  
  type.s.area = pt(ifelse(d > 0, CI[1]/d.SE, CI[2]/d.SE), df, ncp, lower.tail = ifelse(d > 0, TRUE, FALSE))
        power = type.s.area + pt(ifelse(d > 0, CI[2]/d.SE, CI[1]/d.SE), df, ncp, lower.tail = ifelse(d > 0, FALSE, TRUE))
       type.s = type.s.area / power
      p.value = 2*pt(abs(obs.d)/d.SE, df, lower.tail = FALSE)
     random.d = rt(n = 1e6, df, ncp)*d.SE
          sig = if(d > 0) abs(random.d) > CI[2] else -abs(random.d) < CI[1]
  exaggration = if(d > 0) mean(abs(random.d)[sig])/ d else mean(-abs(random.d)[sig])/ d
  
  round(data.frame(exaggration = exaggration, type.s = type.s, power = power, Crit.d = CI[2], p.value = p.value, row.names = "Results:"), digits = digits)
}


#=======================================================================


type.sm.fun <- function(n1, n2 = NA, d.min = 0, d.max = 1.4, alpha = .05)
{
  UseMethod("type.sm.fun")
}

type.sm.fun.default <- function(n1, n2 = NA, d.min = 0, d.max = 1.4, alpha = .05){
  
   type.sm <- function(n1, n2, d, alpha){
    
          N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
         df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
       d.SE = 1/sqrt(N) ; ncp = d*sqrt(N)
         CI = qt(c(alpha/2, 1-(alpha/2)), df)*d.SE
    
type.s.area = pt(ifelse(d > 0, CI[1]/d.SE, CI[2]/d.SE), df, ncp, lower.tail = ifelse(d > 0, TRUE, FALSE))
      power = type.s.area + pt(ifelse(d > 0, CI[2]/d.SE, CI[1]/d.SE), df, ncp, lower.tail = ifelse(d > 0, FALSE, TRUE))
     type.s = type.s.area / power
   random.d = rt(1e4, df, ncp)*d.SE
        sig = if(d > 0) abs(random.d) > CI[2] else -abs(random.d) < CI[1]
exaggration = if(d > 0) mean(abs(random.d)[sig])/ d else mean(-abs(random.d)[sig])/ d
    
    list(exaggration = exaggration, type.s = type.s, power = power)
  }
  
  graphics.off()  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mfrow = c(2, 1), mgp = c(2, .5, 0), mar = c(4, 4, 3, 2), las = 1)  
  
       d_range = seq(d.min, d.max, by = 5e-3)
             n = length(d_range)
         power = numeric(n)
        type.s = numeric(n)
   exaggration = numeric(n)
  
  for(i in 1L:n){
             a = type.sm(d = d_range[i], n1 = n1, n2 = n2, alpha = alpha)
      power[i] = a$power
     type.s[i] = a$type.s
exaggration[i] = a$exaggration
  }
  
  plot(power, type.s, type = "l", xaxt = "n", lwd = 2, font.lab = 2, col = 2)
  axis(1, at = c(alpha, seq(.2, 1, by = .2)))
  abline(v = alpha, col = 8)
  plot(power, exaggration, type = "l", ylim = c(1, 10), xaxt = "n", yaxt = "n", lwd = 2, font.lab = 2, col = 4)
  axis(1, at = c(alpha, seq(.2, 1, by = .2)))
  axis(2, at = seq(1, 10, by = 2))
  abline(h = 1, v = alpha, col = 8)
  
}

                      
#=================================================================================================================

              
index <- function(...)
{
  UseMethod("index")
}
              
index.default <- function(...){
  
  L <- list(...)
  if(is.list(L[[1]]) && length(L) == 1L) L <- L[[1L]]
  if(length(L) == 1L){
  return(as.integer(as.factor(as.character(L[[1L]]))))
    
  }else{
    
    var.names <- match.call()
    var.names <- as.character(var.names)[2L:(length(L) + 1L)]
    
    M <- L
    for(i in 1L:length(L)) M[[i]] <- as.character(L[[i]])
    Mall <- M[[1L]]
    for(i in 2L:length(L)) Mall <- c(Mall, M[[i]])
    Mall <- unique(Mall)
    new.levels <- levels(as.factor(Mall))
    for(i in 1L:length(L)){
      M[[i]] <- factor(M[[i]], levels = new.levels)
      M[[i]] <- as.integer(M[[i]])
    }
    names(M) <- paste(var.names, ".idx", sep = "")
    return(M)
  } 
}                       


#=================================================================================================================
              

standard <- function(data = mtcars, center = TRUE, scale = TRUE, na.rm = TRUE)
{
  UseMethod("standard")
}


standard.default <- function(data = mtcars, scale = TRUE, center = TRUE, na.rm = TRUE){
    
 
message("\nNote: You now have new column(s) in your 'data' with suffix '.s' ('.s' for standardized).
      Also, 'NA's are removed by default. Use 'na.rm = FALSE' for otherwise.") 
    
  if(inherits(data, "data.frame") && ncol(data) > 1){ 
    
    data <- if(na.rm) data[complete.cases(data), ]    
    
    ind <- sapply(data, is.numeric)
    
    data[paste0(names(data)[ind], ".s")] <- lapply(data[ind], scale)
    
    return(data)
  }
  
  if(inherits(data, "data.frame") && ncol(data) == 1){
    
    data <- if(na.rm) data[complete.cases(data), , drop = FALSE]    
    
    d <- scale(data, center = center, scale = scale)
    
    data[, paste0(names(data), ".s") ] <- c(d)
      
    return(data)
  }  
  
  if(!inherits(data, "data.frame")){
    
    data <- if(na.rm) data[complete.cases(data), drop = FALSE]
    
    data <- as.data.frame(data)
      
    names(data) <- "V1"
      
    d <- as.data.frame(scale(data, center = center, scale = scale))  
      
    data[, paste0(names(d), ".s") ] <- d
      
    return(data)
  }
}
   

#=================================================================================================================              
              
              
model.standard <- function(..., level = .95, digits = 6)
{
  UseMethod("model.standard")
}


model.standard.default <- function(..., level = .95, digits = 6){
  
if(!(all(sapply(list(...), inherits, "stanreg")))) stop("Error: all '...' must be models from package 'rstanarm's 'stan_glm()'.")  
  
stand.fit <- function(fit, level, digits){
  
  X <- model.matrix(fit)
  sd_X <- apply(X, MARGIN = 2, FUN = sd)[-1]
  sd_Y <- apply(rstanarm::posterior_predict(fit), MARGIN = 1, FUN = sd)
  beta <- as.matrix(fit)[ , 2:ncol(X), drop = FALSE]
  b <- sweep(sweep(beta, MARGIN = 2, STATS = sd_X, FUN = `*`), 
             MARGIN = 1, STATS = sd_Y, FUN = `/`)
  
  loop <- ncol(b)
  mean <- numeric(loop)
  sd <- numeric(loop)
  
  I <- matrix(NA, loop, 2)
  for(i in 1:loop){
    I[i,] <- hdir(b[, i], level = level)
    mean[i] <- mean(b[, i])
    sd[i] <- sd(b[, i])
  }
  
return(round(data.frame(standard.coef = mean, sd = sd, lower = I[,1], upper = I[,2], level = level, row.names = colnames(b)), digits = digits))
 }

b <- lapply(list(...), stand.fit, level = level, digits = digits)  

return(b)
}
              

#==========================================================================================================
              
              
newdata <- function(fit.data, focus.var, n = 1e2, FUN = mean, hold.at = NA)
{
  UseMethod("newdata")
}

newdata.default <- function(fit.data, focus.var, n = 1e2, FUN = mean, hold.at = NA){
  
if(!inherits(fit.data, "data.frame") || ncol(fit.data) < 2) stop("Error: 'fit.data' must be 'data.frame' with '>= 2' columns.")
    
  foc <- if(is.character(focus.var)) focus.var else deparse(substitute(focus.var))
      
  tgt <- fit.data[, foc]
  focus.var.new <- seq(min(tgt), max(tgt), length.out = n)
  fit.data2 <- data.frame(focus.var.new)
  names(fit.data2) <- foc
  
  for(i in names(fit.data)[!names(fit.data) %in% foc]){
    
    if(is.na(hold.at)){
      
    fit.data2[[i]] <- FUN(fit.data[[i]]) 
    
    } else { 
      
    fit.data2[[i]] <- hold.at
      
      }
  }
  
  fit.data2 <- fit.data2[, names(fit.data)]
  
  return(fit.data2[,-1])
}
              
              
#=========================================================================================================
              
 
counter.plot <- function(fit, xlab = NA, ylab = NA, line.int = TRUE, pred.int = TRUE, level = .95,
                       focus.pred, n = 2e2, FUN = mean, hold.at = NA, legend = "topleft", ...)
{
  UseMethod("counter.plot")
}


counter.plot.default <- function(fit, xlab = NA, ylab = NA, line.int = TRUE, pred.int = TRUE, level = .95,
                               focus.pred, n = 2e2, FUN = mean, hold.at = NA, legend = "topleft", ...){
  
if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")  
if(length(coef(fit)) < 3) stop("Error: 'fit' must contain at least 'two' predictors.")
  
  m <- stats::model.frame(fit)
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
  foc <- if(is.character(focus.pred)) focus.pred else deparse(substitute(focus.pred))
    
if(names(m)[1] == foc) message("\nYou're looking at effect of the changing 'dep.var' against itself on the original model's prediction! This leads to a horizontal line showing 'No Change'!")    
  
    pred <- range(m[foc])
   dep <- range(m[names(m)[1]])
  
  pred <- seq(pred[1], pred[2], length.out = n)
   dep <- seq(dep[1], dep[2], length.out = n)
  
nd <- newdata(m, foc, n = n, FUN = FUN, hold.at = hold.at)   
   
  pred_lin <- rstanarm::posterior_predict(fit, newdata = nd)

xlab <- ifelse(is.na(xlab), foc, xlab)
ylab <- ifelse(is.na(ylab), names(m)[1], ylab)

loop <- n

v1 <- deparse(substitute(FUN))
main <- if(is.na(hold.at)) v1 else hold.at
graphics.off()
plot(dep ~ pred, xlab = xlab, ylab = ylab, type = "n", las = 1, main = paste0("Other predictor(s) held at: ", dQuote(main)), ...)


I <- matrix(NA, loop, 2)
for(i in 1:loop){
  I[i,] = hdir(pred_lin[,i], level = level)
}

x <- sort(pred)
o <- order(pred)
    
if(pred.int){   
  
  y <- I[,1][o]
  z <- I[,2][o]
  
  polygon(c(rev(x), x), c(rev(z), y), col = adjustcolor('gray', .5), border = NA)
}

pred_lin2 <- rstanarm::posterior_linpred(fit, newdata = nd)

I2 <- matrix(NA, loop, 2)
for(i in 1:loop){
  I2[i,] = hdir(pred_lin2[,i], level = level)
}

if(line.int){
  
  y <- I2[,1][o]
  z <- I2[,2][o]
  
  polygon(c(rev(x), x), c(rev(z), y), col = adjustcolor('magenta', .35), border = NA)
}

E.mu <- apply(pred_lin2, 2, mean)

lines(pred, E.mu, col = "cyan", lwd = 2)

legend(x = leg, legend = "Counterfactual Plot\n(prediction analysis)", text.font = 4,
      cex = .9, bg = 0, box.col = 0)    
box()
}              
              
              
#=======================================================================================================
              
              
case.fit.plot <- function(fit, level = .95, legend = "topleft", lwd = 2, fit.tol = 1, pt.cex = 1, cex.lab = .7, col.corridor = "yellow", col.depth = .4)
{
  UseMethod("case.fit.plot")
}  


case.fit.plot.default <- function(fit, level = .95, legend = "topleft", lwd = 2, fit.tol = 1, pt.cex = 1, cex.lab = .7, col.corridor = "yellow", col.depth = .4){
  
if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")  

leg <- if(is.character(legend)) legend else deparse(substitute(legend))    
m <- stats::model.frame(fit)
y <- m[, 1]

loop <- nrow(m)

mus <- rstanarm::posterior_linpred(fit, transform = TRUE)
ys <-rstanarm::posterior_predict(fit, transform = TRUE)

E.mus <- apply(mus, 2, mean)  
CI.reg <- apply(mus, 2, hdir, level = level)
CI.y <- apply(ys, 2, hdir, level = level)

e <- y - E.mus
o <- order(e)

CI.e <- matrix(NA, loop, 2)
PI.e <- matrix(NA, loop, 2)

for(i in 1:loop){
  j <- o[i]
  CI.e[i,] <- c(y[j] - c(CI.reg[1,j], CI.reg[2,j]))
  PI.e[i,] <- c(y[j] - c(CI.y[1,j], CI.y[2,j]))
}

ok <- min(e[o]) < e[o] & e[o] < max(e[o])

unit <- fit.tol*sd(e)

graphics.off()    
original.par = par(no.readonly = TRUE)
on.exit(par(original.par))    
par(mar = c(2.6, 4, 2.2, 4))
    
plot(e[o], 1:loop, cex = .6, xlim = range(PI.e), pch = 19, ylab = NA, yaxt = "n", mgp = c(1.3, .4, 0), type = "n", xlab = "Credible Interval (Residuals)", font.lab = 2)
rect(-unit, par('usr')[3], unit, par('usr')[4], border = NA, col = adjustcolor(col.corridor, col.depth), lend = 1)

abline(v = 0, h = 1:loop, lty = 3, col = 8)

good <- -unit < e[o] & e[o] < unit
    
pos <- (1:loop)[o]

axis(2, at = (1:loop)[-range(1:loop)], labels = paste0("subj ", pos[-range(pos)]), las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .2, 0))
axis(2, at = range(1:loop), labels = paste0("subj ", c(pos[1], rev(pos)[1])), las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .2, 0), col.axis = 2)

segments(PI.e[, 1], 1:loop, PI.e[, 2], 1:loop, lend = 1, col = 8, lwd = lwd)

segments(CI.e[, 1], 1:loop, CI.e[, 2], 1:loop, lend = 1, lwd = lwd, col = ifelse(good, "green3", 1))

points(e[o], 1:loop, pch = 21, bg = ifelse(ok, "cyan", 2), col = ifelse(ok, "magenta", 2), cex = pt.cex)

a <- par('usr')[1:2] ; b <- par('usr')[4]
    
text(.1*mean(a), b, "Fit Corridor", pos = 3, cex = 1, xpd = NA, font = 2, col = "green3")
    
text(.8*a, b, rep("Misfit", 2), pos = 3, cex = 1.5, xpd = NA, font = 2, col = 2)

legend(x = leg, legend = c("Worst fit", "Perfect-Good fit", "Fair-Bad fit"), pch = 22, title = "Person Fit", 
       pt.bg = c(2, "green3", 1), col = c(2, "green3", 1), cex = .7, pt.cex = .6, bg = 0, 
       box.col = 0, xpd = NA, x.intersp = .5, title.adj = .2)
box()
}
        
              
#====================================================================================================================
      
              
model.cor <- function(..., cor = TRUE, digits = 6)
{
  UseMethod("model.cor")
}
              
              
model.cor.default <- function(..., cor = TRUE, digits = 6){
  
if(!(all(sapply(list(...), inherits, "stanreg")))) stop("Error: all '...' must be models from package 'rstanarm's 'stan_glm()'.")  
  
cor.fit <- function(fit, cor, digits){
  
m <- if(cor) cov2cor(cov(as.matrix(fit))) else cov(as.matrix(fit))
  
return(round(m, digits = digits))

}

b <- lapply(list(...), cor.fit, cor = cor, digits = digits)  

return(b)

}
              
              
#============================================================================================================================
  
              
cor.norm <- function(x, r, m = 95, s = 4){
  e  = rnorm(length(x), m, s)
  return(y = round(r*x + e))
}

prof <- {
set.seed(1170)  
data.frame(LAA = LAA <- round(rnorm(60, 32, 4)), TOEFL = cor.norm(LAA, .3)  )

}        
  
              
#=================================================================================================================================
     
                  
              
model.hdi <- function(..., level = .95, digits = 6)
{
  UseMethod("model.hdi")
}


model.hdi.default <- function(..., level = .95, digits = 6){
  
if(!(all(sapply(list(...), inherits, "stanreg")))) stop("Error: all '...' must be models from package 'rstanarm's 'stan_glm()'.")  

hdi.fit <- function(fit, level, digits){
  
m <- round(data.frame(t(apply(as.matrix(fit), 2, hdir, level = level)), level), digits = digits)
  
colnames(m) <- c("lower", "upper", "level")
  
return(m)
}  

m <- lapply(list(...), hdi.fit, level = level, digits = digits)  
  
return(m)
}

              
#==================================================================================================================
      
              
model.info <- function(fit, cor = TRUE, level = .95, digits = 6)
{
  UseMethod("model.info")
}  

              
model.info.default <- function(fit, cor = TRUE, level = .95, digits = 6){
  
 if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be a fitted model from package 'rstanarm's 'stan_glm()'.") 
    
 hdi.fit <- function(fit, level, digits){
    
    m <- round(data.frame(t(apply(as.matrix(fit), 2, hdir, level = level)), level), digits = digits)
    
    colnames(m) <- c("lower", "upper", "level")
    
    return(m)
  }  

cor.fit <- function(fit, cor, digits){
  
  m <- if(cor) cov2cor(cov(as.matrix(fit))) else cov(as.matrix(fit))
  
  return(round(m, digits = digits))
}

stand.fit <- function(fit, level, digits){
  
  X <- model.matrix(fit)
  sd_X <- apply(X, MARGIN = 2, FUN = sd)[-1]
  sd_Y <- apply(rstanarm::posterior_predict(fit), MARGIN = 1, FUN = sd)
  beta <- as.matrix(fit)[ , 2:ncol(X), drop = FALSE]
  b <- sweep(sweep(beta, MARGIN = 2, STATS = sd_X, FUN = `*`), 
             MARGIN = 1, STATS = sd_Y, FUN = `/`)
  
  loop <- ncol(b)
  mean <- numeric(loop)
  sd <- numeric(loop)
  I <- matrix(NA, loop, 2)
  
  for(i in 1:loop){
    I[i,] <- hdir(b[, i], level = level)
    mean[i] <- mean(b[, i])
    sd[i] <- sd(b[, i])
  }
  
  return(round(data.frame(standard.coef = mean, sd = sd, lower = I[,1], upper = I[,2], level = level, row.names = colnames(b)), digits = digits))
}

return(list(model.hdi = hdi.fit(fit = fit, level = level, digits = digits), 
           model.cor = cor.fit(fit = fit, cor = cor, digits = digits), 
           model.standard = stand.fit(fit = fit, level = level, digits = digits)))
}
              
#============================================================================================================================
              

models.info <- function(..., cor = TRUE, level = .95, digits = 6)
{
  UseMethod("models.info")
}  

models.info.default <- function(..., cor = TRUE, level = .95, digits = 6){

return(list(model.hdi = model.hdi(... = ..., level = level, digits = digits), 
           model.cor = model.cor(... = ..., cor = cor, digits = digits), 
           model.standard = model.standard(... = ..., level = level, digits = digits)))
}
              
              
#=================================================================================================================================
              
 
icc <- function(a = 1, b = 0, c = 0, from = -4, to = 4, legend = "topleft")
{
  UseMethod("icc")
}

              
icc.default <- function(a = 1, b = 0, c = 0, from = -4, to = 4, legend = "topleft"){
  
  if(any(c > 1) || any(c < 0)) stop("Error: 'Gussing' parameter 'c' must be between '0' and '1'.")
  
  I <- eq(a, b, c)
  
  a = I[[1]] ; b = I[[2]] ; c = I[[3]]
  
  loop <- seq_along(a)
  
  h <- list()
graphics.off()  
for(i in loop){
p <- function(x) c[i] + ((1 - c[i])/(1 + exp(-a[i]*(x - b[i]))))  
h[[i]] <- curve(p, from, to, add = i!= 1, n = 1e3, las = 1, ylim = 0:1,
          font.lab = 2, xlab = bquote(bold("Person Ability"~(theta))),
          type = "n", ylab = bquote(p(theta)), mgp = c(2, .5, 0), tck = -.015)
}  
                              
  u <- par('usr')
  abline(h = c(0, 1, .5), col = 8, lty = c(3, 3, 2))
  axis(2, at = .5, col = 2, col.axis = 2, mgp = c(2, .5, 0), tck = -.015, las = 1)
  
for(i in loop){
  lines(h[[i]], col = i, lwd = 2)
  segments(b[i], u[3], b[i], .5, col = i, lty = 3)  
  points(b[i], .5, pch = 21, col = i, font = 2, cex = 1.5, bg = "cyan", lwd = 2)
}
  
  legend(legend, paste0("a = ", round(a, 2), ", b = ", round(b, 2), ", c = ", round(c, 2)), 
         pch = 22, title = "IRT Parameters", pt.bg = loop, col = loop, cex = .7, pt.cex = .6, 
         bg = 0, box.col = 0, x.intersp = .5, title.adj = .5, title.col = "red4")
  box()
}
                              
                            
#===========================================================================================================================
                              
                              
logit <- function(x){ 
    
 return(stats::qlogis(x)) 
    
}


#===========================================================================================================================
                              
                              
inv.logit <- function(x, percent = FALSE, digits = 4){
  
  p <- stats::plogis(x)
  return(if(percent) 
  noquote(paste0(round(p*1e2, digits = digits), "%")) else p)
  
}            
                              
#==============================================================================================================================
         
log.sum.exp <- function (x) 
{
  xmax <- max(x)
  xsum <- sum(exp(x - xmax))
  xmax + log(xsum)
}         
         
#==============================================================================================================================
     
         
dzpois <- function (x, p, lambda, log = FALSE) 
{
  ll <- rep(0, length(x))
  p_i <- p[1]
  lambda_i <- lambda[1]
  for (i in 1:length(x)) {
    if (length(p) > 1) 
      p_i <- p[i]
    if (length(lambda) > 1) 
      lambda_i <- lambda[i]
    if (x[i] == 0) {
      ll[i] <- log.sum.exp(c(log(p_i), log(1 - p_i) + dpois(x[i], 
                          lambda_i, TRUE)))
    }
    else {
      ll[i] <- log(1 - p_i) + dpois(x[i], lambda_i, log = TRUE)
    }
  }
  if (log == FALSE) 
    ll <- exp(ll)
  return(ll)
}         
         
         
#==============================================================================================================================
         
         
multilogit <- function(...)
{
  UseMethod("multilogit")
}


multilogit.default <- function (...){
  
  X <- list(...)
  K <- length(X)
  X <- as.data.frame(X)
  N <- nrow(X)
  if(N == 1){
    f <- exp(X[1, ])
    below <- sum(f)
    as.numeric(f/below)
  } else {
    f <- lapply(1:N, function(i) exp(X[i, ]))
    below <- sapply(1:N, function(i) sum(f[[i]]))
    p <- sapply(1:N, function(i) unlist(f[[i]])/below[i])
    p <- t(as.matrix(p))
    colnames(p) <- NULL
    p
  }
}
                
                
#====================================================================================================================
             
                
inv.multilogit <- function(x, lambda = 1, diff = TRUE, log = FALSE){
  
  x <- round(x)
  
  if(length(x) == 1){ x <- 0:x  ;
  message("Note: ", length(x), " categories were assumed.")
  }
  
  if(diff){ 
    x <- x - min(x)
    f <- exp(lambda * x)
  }
  if(!log){
    output <- f/sum(f)
  } else {
    output <- log(f) - log(sum(f))
  }
  output
}
         
                
#====================================================================================================================

                
anova.es <- function(fit = NULL, f, df1, df2, N, conf.level = .9, digits = 6)
{
  UseMethod("anova.es")
}
  
                
anova.es.default <- function(fit = NULL, f, df1, df2, N, conf.level = .9, digits = 6){
  
  if(!is.null(fit)){
    
    if(class(fit)[1] != "lm" && class(fit)[1] != "aov") { stop("Error: 'fit' must be a fitted model from base R's 'aov()' or 'lm()' commands.") }        
    N <- nobs(fit)
    sit <- anova(fit)
    #if(!("F value" %in% names(sit))) { stop("Error: Fitted model does not include any 'F value'.") } 
    f <- head(sit[,4], -1)
    df1 <- head(sit$Df, -1)
    df2 <- tail(sit$Df, 1)
  }
  
  if(length(f) != length(df1)){message("Warning: The length of 'f' and 'df1' must be equal. Check your inputted values.\n")}
  I <- eq(f, df1) 
  f = I[[1]] ; df1 = I[[2]]
  
  omega <- (df1 * (f - 1)) / as.numeric(crossprod(df1, f) + df2 + 1)
  eta <- (df1 * f) / as.numeric(crossprod(df1, f) + df2)
  pomega <- (df1 * (f - 1)) / ((df1 * (f - 1)) + N)
  peta <- peta.ci(f = f, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = digits)
  
  result <- round(data.frame(F.value = f, eta.sq = eta, P.eta.sq = peta[,1], lower.P.eta.sq = peta[,2], 
                             upper.P.eta.sq = peta[,3], conf.level = conf.level, omega.sq = omega, 
                             P.omega.sq = pomega, row.names = paste0("effect ", 1:length(f), ":")), digits = digits)
  
  message("Note: If analysis includes random-effects, manually input the right 'df2' to obtain correct 'P.eta- or P.omega-sq.'")
  
  if(is.null(fit)){  
    
    return(result)
    
  }else{
    
    rownames(result) <- head(rownames(sit), -1)
    
    return(result)
  } 
}
                
                
#===========================================================================================================================
                
                
dens.plot <- function(x, adjust = 1, na.rm = TRUE, n = 1e3, from = min(x), to = max(x), add = FALSE, hdi = FALSE, ci = FALSE, level = .95, xlab = deparse(substitute(x)), main = NA, lwd = 2, lty = 1, ...){
  
  UseMethod("dens.plot")
}


dens.plot.default <- function(x, adjust = 1, na.rm = TRUE, n = 1e3, from = min(x), to = max(x), add = FALSE, hdi = FALSE, ci = FALSE, level = .95, xlab = deparse(substitute(x)), main = NA, lwd = 2, lty = 1, ...){
  
  d <- density(x, adjust = adjust, na.rm = na.rm, n = n, from = from, to = to)
  
  if(!add){
    
    graphics.off()                            
    
    plot(d, zero.line = FALSE, xlab = xlab, main = main, lwd = lwd, lty = lty, ...)
    
  } else {
    
    lines(d, lwd = lwd, lty = lty, ...)
    
  }
  
  
  alpha <- (1 - level)/2
  q <- if(ci) quantile(x, probs = c(alpha, 1 - alpha), na.rm = TRUE) else c(NA, NA)
  i <- if(hdi) hdir(x, level = level) else c(NA, NA)
  
  
  mode <- d$x[which.max(d$y)]
  mean <- mean(x)
  median <- median(x)
  sd <- sd(x)
  mad <- mad(x)
  
  if(hdi){
    h <- min(d$y)
    lines(i, rep(h, 2), lend = 1, lwd = 6, lty = 1, xpd = NA, ...)
    text(i, h, round(i, 3), pos = 3, cex = .8, font = 2, xpd = NA)
    points(mode, h, pch = 21, bg = "cyan", col = "magenta", cex = 1.7, xpd = NA)
  }
  
  invisible(list(lower = i[1], upper = i[2], level = level, mean = mean, mode = mode, median = median, 
                 mad = mad, sd = sd, q1 = q[[1]], q2 = q[[2]], x = d$x, y = d$y))
}

                
#===================================================================================================================
                
                
count.plot <- function(x, freq = FALSE, type = "h", lwd = 2, lend = 2, xlab = "Trials", ylab = NA, xaxt = "s", add = FALSE, ...)
{
  UseMethod("count.plot")
}


count.plot.default <- function(x, freq = FALSE, type = "h", lwd = 2, lend = 2, xlab = "Outcomes", ylab = NA, xaxt = "s", add = FALSE, ...)
{  
  x <- sapply(x, round)
  ylab <- if(is.na(ylab) & freq) "Frequency" else if(is.na(ylab) & !freq) "Probability" else ylab
  z <- if(freq) table(x) else table(x)/length(x)
  x <- as.numeric(names(z))
  y <- as.numeric(z)
  graph(x, y, type = type, lwd = lwd, lend = lend, xlab = xlab, ylab = ylab, xaxt = "n", add = add, ...)
  if(xaxt != "n") axis(1, at = min(x):max(x), labels = if(add) FALSE else TRUE, tick = if(add) FALSE else TRUE)
  invisible(list(x = x, y = y))
}
                
#=========================================================================================================================

likert <- function(x){
  
  x <- sapply(x, round)
  fq <- table(x)
  prop <- fq/length(x)
  cumprop <- cumsum(prop)
  logcumodd <- log(cumprop/(1-cumprop))
  x <- as.numeric(names(fq))
  dif <- cumprop[x[2]:tail(x, 1)] - cumprop[x[1]:tail(x, 2)[-2]]
  
  list(x = x, prop = prop, cumprop = cumprop, logcumodd = logcumodd, ordlike = c(cumprop[1], dif))
}
    
    
#=========================================================================================================================

not.integer <- function(x) (abs((x) - floor((x) + .5)) > 1e-7)  

#=========================================================================================================================
                
dbetabinom <- function (x, size, mu.p, disp, shape1 = NULL, shape2 = NULL, log = FALSE) 
{
  if(missing(mu.p) && !is.null(shape1) && !is.null(shape2)){
    mu.p <- shape1/(shape1 + shape2)
    disp <- shape1 + shape2
  }
    
if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
mu.p[mu.p < 0] <- 0 ;
mu.p[mu.p > 1] <- 1 }   
  
t <- disp * mu.p
u <- disp * (1 - mu.p)
h <- lbeta(x + t, size - x + u) - lbeta(t, u) + lchoose(size, x)
  
  if(any(g <- not.integer(x))){
    message("Warning: For non-integer 'x' (successes), probability of \"ZERO\" is returned.")
    h[g] <- -Inf
  }
  if(log) h else exp(h)
}


#==============================================================================================================================


dbetab <- function (x, mu.p, disp, log = FALSE){

if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
mu.p[mu.p < 0] <- 0 ;
mu.p[mu.p > 1] <- 1 }
    
  shape1 <- mu.p * disp
  shape2 <- (1 - mu.p) * disp
  dbeta(x, shape1 = shape1, shape2 = shape2, log = log)
}

    
#=================================================================================================================================

    
qbetab <- function(p, mu.p, disp, lower.tail = TRUE, log.p = FALSE){

if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
mu.p[mu.p < 0] <- 0 ;
mu.p[mu.p > 1] <- 1 }
    
  shape1 <- mu.p * disp
  shape2 <- (1 - mu.p) * disp
  qbeta(p, shape1 = shape1, shape2 = shape2, lower.tail = lower.tail, log.p = log.p)
}
    

#===========================================================================================================


pbetab <- function(q, mu.p, disp, lower.tail = TRUE, log.p = FALSE){
    
if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
mu.p[mu.p < 0] <- 0 ;
mu.p[mu.p > 1] <- 1 }
    
  shape1 <- mu.p * disp
  shape2 <- (1 - mu.p) * disp
  pbeta(q, shape1 = shape1, shape2 = shape2, lower.tail = lower.tail, log.p = log.p)
}


#==================================================================================================================    

pbetabinom <- function(q, size, mu.p, disp){
  
  k <- eq(q, size, mu.p, disp)
  q <- k[[1]] ; size <- k[[2]] ; mu.p <- k[[3]] ; disp <- k[[4]]
 
if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
mu.p[mu.p < 0] <- 0 ;
mu.p[mu.p > 1] <- 1 }
    
  if(any(g <- not.integer(q))){
    message("Warning: For non-integer 'q' (successes),  'q' is rounded.")
    q[g] <- round(q)
  }
  
  t <- disp * mu.p
  u <- disp * (1 - mu.p)
  
  prob <- numeric(length(q))
  for (i in 1:length(q)){
    h <- 0:q[i]
    prob[i] <- sum(exp(lbeta(h + t[i], size[i] - h + u[i]) - lbeta(t[i], u[i]) + lchoose(size[i], h)))
  }
  prob
}


#==========================================================================================================


qbetabinom <- function(p, size, mu.p, disp){

k <- eq(p, size, mu.p, disp)
p <- k[[1]] ; size <- k[[2]] ; mu.p <- k[[3]] ; disp <- k[[4]]

if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
mu.p[mu.p < 0] <- 0 ;
mu.p[mu.p > 1] <- 1 }
    
p[p < 0] <- 0
p[p > 1] <- 1
    
h <- function(g) {
  t <- disp[i] * mu.p[i]
  u <- disp[i] * (1 - mu.p[i])
  d <- 0:g
  sum(exp(lbeta(d + t, size[i] - d + u) - lbeta(t, u) + 
            lchoose(size[i], d))) - p[i]
}

qs <- numeric(length(p))

for(i in 1:length(p)){
  interval <- c(0, size[i])
  qs[i] <- if(h(interval[1]) * h(interval[2]) > 0) 
    0
  else uniroot(h, interval)[[1]]
}
round(qs)
}
    
    
#====================================================================================================================================

    
rbetab <- function(n, mu.p, disp){
  
if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
mu.p[mu.p < 0] <- 0 ;
mu.p[mu.p > 1] <- 1 }
    
  shape1 <- mu.p * disp
  shape2 <- (1 - mu.p) * disp
  rbeta(n, shape1 = shape1, shape2 = shape2)
}


#====================================================================================================================================


rbetabinom <- function(n, size, mu.p, disp, shape1 = NULL, shape2 = NULL){
    
  if(missing(mu.p) && !is.null(shape1) && !is.null(shape2)){
    mu.p <- shape1/(shape1 + shape2)
    disp <- shape1 + shape2
  } 
    
if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
mu.p[mu.p < 0] <- 0 ;
mu.p[mu.p > 1] <- 1 }

  rbinom(n, size, rbetab(n, mu.p, disp))
}
    
    
#=======================================================================================================================================
    
  
 dcohen <- function(x, dbase = 0, n1, n2 = NA, log = FALSE){
  
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  ncp <- dbase*sqrt(N)
  
  dt(x*sqrt(N), df, ncp, log = log)*sqrt(N)
}

#=======================================================================================================================================


qcohen <- function(p, dbase = 0, n1, n2 = NA, lower.tail = TRUE, log.p = FALSE){
  
  q <- Vectorize(function(p, dbase, n1, n2, lower.tail, log.p){
    
    N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    ncp <- dbase*sqrt(N)
    
    qt(p, df, ncp, lower.tail = lower.tail, log.p = log.p)/sqrt(N)
  })
  q(p = p, dbase = dbase, n1 = n1, n2 = n2, lower.tail = lower.tail, log.p = log.p)
}

    
#=======================================================================================================================================


pcohen <- function(q, dbase = 0, n1, n2 = NA, lower.tail = TRUE, log.p = FALSE){
  
  p <- Vectorize(function(q, dbase, n1, n2, lower.tail, log.p){
  
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  ncp <- dbase*sqrt(N)
  
  pt(q*sqrt(N), df, ncp, lower.tail = lower.tail, log.p = log.p)
})
p(q = q, dbase = dbase, n1 = n1, n2 = n2, lower.tail = lower.tail, log.p = log.p)
}

    
#=======================================================================================================================================


rcohen <- function(n, dbase = 0, n1, n2 = NA){
  
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  ncp <- dbase*sqrt(N)
  
  rt(n, df, ncp)/sqrt(N)
}
    
    
#=========================================================================================================================================
    
    
dpeta <- function(x, df1, df2, pbase = 0, N, log = FALSE){
  x[x > .9999999] <- .9999999
  x[x < 0] <- 0
  pbase[pbase > .9999999] <- .9999999
  pbase[pbase < 0] <- 0
  ncp <- (pbase * N) / (1 - pbase)
  d <- df2 / df1
  f <- x / (1 - x) * d
  df(f, df1, df2, ncp, log = log) * d * ( 1 / (1 - x) + x / (1 - x)^2 )
}


#=========================================================================================================================================


ppeta <- function(q, df1, df2, pbase = 0, N, lower.tail = TRUE, log.p = FALSE){
  
  p <- Vectorize(function(q, df1, df2, pbase, N, lower.tail, log.p){
    
  q[q > .9999999] <- .9999999
  q[q < 0] <- 0
  pbase[pbase > .9999999] <- .9999999
  pbase[pbase < 0] <- 0
  ncp <- (pbase * N) / (1 - pbase)
  d <- df2 / df1
  f <- q / (1 - q) * d
  pf(f, df1, df2, ncp, lower.tail = lower.tail, log.p = log.p)
})
p(q = q, df1 = df1, df2 = df2, pbase = pbase, N = N, lower.tail = lower.tail, log.p = log.p)
}


#=========================================================================================================================================


qpeta <- function(p, df1, df2, pbase = 0, N, lower.tail = TRUE, log.p = FALSE){
  
  q <- Vectorize(function(p, df1, df2, pbase, N, lower.tail, log.p){
    
  p[p > 1] <- 1
  p[p < 0] <- 0
  pbase[pbase > .9999999] <- .9999999
  pbase[pbase < 0] <- 0
  ncp <- (pbase * N) / (1 - pbase)
  d <- df2 / df1
  f <- qf(p, df1, df2, ncp, lower.tail = lower.tail, log.p = log.p)
  f / (f + d)
})
q(p = p, df1 = df1, df2 = df2, pbase = pbase, N = N, lower.tail = lower.tail, log.p = log.p)
}


#=========================================================================================================================================


rpeta <- function(n, df1, df2, pbase = 0, N){
  pbase[pbase > .9999999] <- .9999999
  pbase[pbase < 0] <- 0
  ncp <- (pbase * N) / (1 - pbase)
  d <- df2 / df1
  f <- rf(n, df1, df2, ncp)
  f / (f + d)
}
    
    
#==================================================================================================================

dpetab <- function(x, df1, df2, ncp = 0, log = FALSE){
  x[x > .9999999] <- .9999999
  x[x < 0] <- 0
  d <- df2 / df1
  f <- x / (1 - x) * d
  df(f, df1, df2, ncp, log = log) * d * ( 1 / (1 - x) + x / (1 - x)^2 )
}


#=========================================================================================================================================


ppetab <- function(q, df1, df2, ncp = 0, lower.tail = TRUE, log.p = FALSE){
  q[q > .9999999] <- .9999999
  q[q < 0] <- 0
  d <- df2 / df1
  f <- q / (1 - q) * d
  pf(f, df1, df2, ncp, lower.tail = lower.tail, log.p = log.p)
}


#=========================================================================================================================================


qpetab <- function(p, df1, df2, ncp = 0, lower.tail = TRUE, log.p = FALSE){
  p[p > 1] <- 1
  p[p < 0] <- 0
  d <- df2 / df1
  f <- qf(p, df1, df2, ncp, lower.tail = lower.tail, log.p = log.p)
  f / (f + d)
}


#=========================================================================================================================================


rpetab <- function(n, df1, df2, ncp = 0){
  d <- df2 / df1
  f <- rf(n, df1, df2, ncp)
  f / (f + d)
}
    
#=========================================================================================================================================
   
dbern <- function(x, prob, log = FALSE){
  dbinom(x, 1, prob, log = log)
}

    
#==================================================================================================================


pbern <- function(q, prob, lower.tail = TRUE, log.p = FALSE){
 pbinom(q, 1, prob, lower.tail = lower.tail, log.p = log.p)
}

#==================================================================================================================


qbern <- function(p, prob, lower.tail = TRUE, log.p = FALSE){
qbinom(p, 1, prob, lower.tail = lower.tail, log.p = log.p)
}

#==================================================================================================================


rbern <- function(n, prob){
rbinom(n, 1, prob)
}
    
#====================================================================================================================
    
    
likert <- function(x){
  
  x <- sapply(x, round)
  fq <- table(x)
  len <- length(x)
  prop <- fq/len
  cumprop <- cumsum(prop)
  logcumodd <- log(cumprop/(1-cumprop))
  x <- as.numeric(names(fq))
  dif <- cumprop[x[2]:tail(x, 1)] - cumprop[x[1]:tail(x, 2)[-2]]
  
  list(x = x, prop = prop, cumprop = cumprop, logcumodd = logcumodd, ordlike = c(cumprop[1], dif))
}
    

#==========================================================================================================================
    
    
typem.anova <- function(peta.h1, df1, df2, N, alpha = .05, peta.h0 = 0, peta.obs = .1, xlab = bquote(eta[p]^2), from = 0, to = .2){

graphics.off()  
original.par = par(no.readonly = TRUE)
on.exit(par(original.par))
  
par(mfrow = c(2, 1), mgp = c(2, .5, 0), mar = c(4, 4, 3, 2), xpd = TRUE)
  
h0 = curve(dpeta(x, df1, df2, peta.h0, N), from = from, to = to, lwd = 2, n = 1e4, xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i")
a <- qpeta(alpha, df1, df2, peta.h0, N, lower.tail = FALSE)
x = seq(a, 1, l = 1e3) ; y = dpeta(x, df1, df2, peta.h0, N)
polygon(c(a, x, 1), c(0, y, 0), col = 2, border = NA)
lines(h0, lwd = 2)
abline(v = peta.h0, col = 2, xpd = FALSE) 

h1 = curve(dpeta(x, df1, df2, peta.h1, N), from = from, to = to, lwd = 2, n = 1e4, xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i")
x = seq(a, 1, l = 1e3) ; y = dpeta(x, df1, df2, peta.h1, N)
polygon(c(a, x, 1), c(0, y, 0), col = 2, border = NA)
lines(h1, lwd = 2)
abline(v = peta.h1, col = 4, xpd = FALSE)

points(peta.obs, 0, pch = 23, bg = "cyan", col = "magenta", cex = 1.5)

abline(v = a, col = 2, lty = 2, xpd = NA)

power <- ppeta(a, df1, df2, peta.h1, N, lower.tail = FALSE)
p.value <- ppeta(peta.obs, df1, df2, peta.h0, N, lower.tail = FALSE)

random.p <- rpeta(1e6, df1, df2, peta.h1, N)
sig <- random.p > a
exaggeration <- mean(random.p[sig]) / peta.h1

data.frame(power = power, p.value = p.value, exaggeration = exaggeration, row.names = "Result:")
}
    
#=====================================================================================================================================
    
    
typem.anova.fun <- function(df1, df2, N, peta.h0 = 0, peta.min = 0, peta.max = .5, alpha = .05){

peta <- function(df1, df2, peta.h1, peta.h0, N, alpha){
  
a <- qpeta(alpha, df1, df2, peta.h0, N, lower.tail = FALSE)
power <- ppeta(a, df1, df2, peta.h1, N, lower.tail = FALSE)
random.p <- rpeta(1e4, df1, df2, peta.h1, N)
sig <- random.p > a
exaggration <- mean(random.p[sig]) / peta.h1

list(power = power, exaggration = exaggration)
}

peta_range = seq(peta.min, peta.max, by = .001)
n = length(peta_range)
power = numeric(n)
exaggration = numeric(n)

for(i in 1L:n){
  g = peta(peta.h1 = peta_range[i], df1 = df1, df2 = df2, N = N, alpha = alpha, peta.h0 = peta.h0)
  power[i] = g$power
  exaggration[i] = g$exaggration
}

plot(power, exaggration, type = "l", ylim = c(1, 10), xaxt = "n", yaxt = "n", lwd = 2, font.lab = 2, col = 4)
axis(1, at = c(alpha, seq(.2, 1, by = .2)))
axis(2, at = seq(1, 10, by = 2), las = 1)
abline(h = 1, v = alpha, col = 8)
}
    
#========================================================================================================================
    

plan.t.tests <- function(d = .1, sig.level = .05, power = .8, base.rate = 1, paired = FALSE, d.range = seq(.1, .5, .05),
                          two.tailed = TRUE, xlab = "Cohen's d", xlim = c(NULL, NULL), ylim = NULL)
{
  
  UseMethod("plan.t.tests")
}
    
    
plan.t.tests.default <- function(d = .1, sig.level = .05, power = .8, base.rate = 1, paired = FALSE, d.range = seq(.1, .5, .05),
                                  two.tailed = TRUE, xlab = "Cohen's d", xlim = c(NULL, NULL), ylim = NULL){
  
  graphics.off()  
  original.par <- par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  d[d == 0] <- 1e-4
  if(d == -Inf) d <- -6  
  if(d == Inf) d <- 6
  if(power == 0) power <- sig.level
  
  from <- xlim[1]
  to <- xlim[2]
  
  d2 <- d
  d3 <- d.range
  d <- abs(d)
  sig.level <- if(two.tailed) sig.level/2 else sig.level
  k <- base.rate / (1 + base.rate)
  options(warn = -1)
  
  f <- if(two.tailed){ function(x){
    
    power - (pt(qt(sig.level, df = x), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2))) + pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE))
  }
    
  } else {
    
    function(x){
      
      power - pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE)
      
    }
  }
  
  df <- ceiling(uniroot(f, c(1e-8, 1e6), extendInt = "downX")[[1]])
  
  n1 <- df + 1
  n2 <- if(paired) NA else round(base.rate*n1)
  
  loop <- length(d3)
  
  dfb <- numeric(loop)
  n1b <- numeric(loop)
  n2b <- numeric(loop)
  
  for(i in 1:loop){
    
    f <- if(two.tailed){ function(x){
      
      power - (pt(qt(sig.level, df = x), df = x, ncp = d3[i]*sqrt(if(paired) x + 1 else k*(x + 2))) + pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d3[i]*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE))
    }
      
    } else {
      
      function(x){
        
        power - pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d3[i]*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE)
        
      }
    }
    
    dfb[i] <- ceiling(uniroot(f, c(1e-8, 1e6), extendInt = "downX")[[1]])
    
    n1b[i] <- dfb[i] + 1
    n2b[i] <- if(paired) NA else round(base.rate*n1b[i])
  }
  
  base.rate <- if(paired) NA else base.rate
  method <- paste(if(paired) "One- or Paired sample" else "Two-sample", "t test power analysis")
  note <- paste("Use 'base.rate' to specify how many times one group might be larger than the other e.g.,'base.rate = 1.1' ")
  
  a <- qcohen(sig.level, 0, n1, n2)
  b <- -a
  
  est.power <- if(two.tailed) pcohen(a, d, n1, n2) + pcohen(b, d, n1, n2, lower.tail = FALSE) else pcohen(b, d, n1, n2, lower.tail = FALSE)
  
  from <- if(is.null(from)) min(qcohen(1e-5, 0, n1, n2), qcohen(1e-5, d2, n1, n2), na.rm = TRUE) else from
  to <- if(is.null(to)) max(qcohen(.99999, 0, n1, n2), qcohen(.99999, d2, n1, n2), na.rm = TRUE) else to
  
  x <- seq(from, to, 1e-4)
  ylimb <- c(0, max(c(dcohen(x, 0, n1, n2), dcohen(x, d2, n1, n2)), na.rm = TRUE) )
  
  ylim <- if(is.infinite(ylimb[2]) & is.null(ylim)) NULL else if(is.null(ylim)) ylimb else ylim
  
  par(mfrow = c(2, 1), mgp = c(2.5, .5, 0), mar = c(4, 4, 2, 2), tck = -.02)
      
  h0 <- curve(dcohen(x, 0, n1, n2), from, to, n = 1e4, xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i", ylim = ylim, font.lab = 2)
  
  x1 <- seq(from, a, length.out = 1e3) ; y1 <- dcohen(x1, 0, n1, n2) 
  x2 <- seq(b, to, length.out = 1e3) ; y2 <- dcohen(x2, 0, n1, n2)
  
  if(d2 < 0 & !two.tailed || two.tailed) polygon(c(from, x1, a), c(0, y1, 0), col = adjustcolor(2, .25), border = NA)
  if(d2 > 0 & !two.tailed || two.tailed) polygon(c(b, x2, to), c(0, y2, 0), col = adjustcolor(2, .25), border = NA) 
  lines(h0, lwd = 2, col = 2, xpd = TRUE)
  
  g <- if(d2 < 0 & !two.tailed) a else if(d2 > 0 & !two.tailed) b else c(a, b)
  p <- if(two.tailed) rep(par('usr')[4], 2) else par('usr')[4]
  abline(v = g, col = 2, lty = 2)
  
  crit <- round(g, 4) 
  
  points(g, p, pch = 19, col = 2, xpd = NA)  
  
  #text(g, par('usr')[4], bquote(bold(critical~ bolditalic(d) == .(crit))), pos = 3, cex = .7, font = 2, xpd = TRUE) 
  text(g, par('usr')[4], paste("critical d =", crit), pos = 3, cex = .7, font = 2, xpd = TRUE)
  
  h1 <- curve(dcohen(x, d2, n1, n2), from, to, n = 1e4, add = TRUE) 
  x1 <- seq(from, a, length.out = 1e3) ; y1 <- dcohen(x1, d2, n1, n2) 
  x2 <- seq(b, to, length.out = 1e3) ; y2 <- dcohen(x2, d2, n1, n2)
  if(d2 < 0 & !two.tailed || two.tailed) polygon(c(from, x1, a), c(0, y1, 0), border = NA, density = 15, col = 4, xpd = TRUE)
  if(d2 > 0 & !two.tailed || two.tailed) polygon(c(b, x2, to), c(0, y2, 0), border = NA, density = 15, col = 4, xpd = TRUE) 
  lines(h1, lwd = 2, col = 4, xpd = TRUE)
  
  legend("topleft", legend = c("Sig. Area(s)", "Power"), inset = c(-.15, 0), density = c(NA, 35), x.intersp = c(.3, .3),
         bty = "n", xpd = NA, cex = .7, text.font = 2, angle = c(NA, 45), fill = c(adjustcolor(2, .4), 4), border = c(2, 4), adj = c(0, .4))
  
  plot(d3, n1b, type = "b", pch = 19, lwd = 2, xlab = xlab, las = 1, col = 4, font.lab = 2, ylab = "Group Sample Size", ylim = c(0, max(n1b, n2b, na.rm = TRUE)), xaxt = "n")
  axis(1, at = d3)
  if(!paired)lines(d3, n2b, col = 2, lty = 3, lwd = 2, type = "b")
  
  points(if(!paired)rep(d, 2) else d, if(!paired) c(n1, n2) else n1, col = "magenta", bg = "cyan", pch = 21, cex = 1.5)
  
  text(d3, n1b, n1b, pos = 1, font = 2, col = 4, cex = .7, xpd = TRUE)
  if(!paired) text(d3, n2b, n2b, pos = 3, font = 2, col = 2, cex = .7, xpd = TRUE)
  
  if(paired){
    
    legend("topright", legend = "Group 1", col = 4, pch = 19, cex = .7, text.font = 2, lwd = 1,
           pt.cex = 1, bty = "n")
  } else {
    
    legend("topright", paste("Group", 1:2), col = c(4, 2), pch = c(19, 1), cex = .7, text.font = 2, x.intersp = c(.6, .6),
           lwd = 1, adj = c(0, .4), pt.cex = 1, pt.lwd = c(1, 2), lty = c(1, 3), bty = "n")
  }
  
  box()
  sig.level <- if(two.tailed) sig.level*2 else sig.level
  two.tailed <- if(two.tailed) "Yes" else "No"
  
  structure(list(n1 = n1, n2 = n2, base.rate = base.rate, d = d, est.power = est.power, sig.level = sig.level, 
                 two.tailed = two.tailed, method = method, note = note), class = "power.htest")
}
       
                  
#=============================================================================================================================
                  
                  
gpower.peta <- function(spss, df2, N, design){
    
(spss * df2) / (N - (spss * design))
    
}
                  
gpower.peta.bw <- function(peta, rho = .5, N, m, n.group){
    
((1 - rho)*peta*(N - n.group)*(m -1)) / ((1 - peta)*(m*N) + (1-rho)*peta*(N - n.group)*(m -1))
    
}



gpower.peta.b <- function(peta, rho = .5, N, m, n.group){
    
(peta*(N - n.group)*(1 + (m-1)*rho)) / ((1 - peta)*(m*N) + peta*(N - n.group)*(1 + (m-1)*rho))
    
}
   
#===============================================================================================================================

                  
plan.f.tests <- function(pov, n.level, design, sig.level = .05, n.covar = 0, n.pred = NA, power = .8, peta.range = seq(1e-1, .9, 1e-1),
                         xlab = NULL, ylim = NULL, to = NULL, regress = FALSE, d = NA)
{
  
  UseMethod("plan.f.tests")
}


plan.f.tests.default <- function(pov, n.level, design, sig.level = .05, n.pred = NA, n.covar = 0, power = .8, peta.range = seq(1e-1, .9, 1e-1),
                                 xlab = NULL, ylim = NULL, to = NULL, regress = FALSE, d = NA){
  
  graphics.off()  
  original.par <- par(no.readonly = TRUE)
  on.exit(par(original.par))
  options(warn = -1)
  if(regress & !is.na(n.pred)) n.level <- n.pred
  
  peta2 <- peta.range
  
  peta2[peta2 == 0] <- 1e-2
  peta2[peta2 == 1] <- .99
  
  if(!is.na(d) & !regress) { message("\nNote: You are doing reseach planning for 'pairwise' comparisons.") ;  n.level <- design <- 2 }
  if(!is.na(d)) pov <- d2peta(d = d, n1 = 300, n2 = 300) 
  peta <- pov
  
  if(n.level <= 1 & !regress) stop("Error: You must have at least '2 levels'.")
  if(n.level < 1) stop("Error: You must have at least '2 levels' or '1 predictor' for regression.")
  xlab <- if(is.null(xlab) && !regress) bquote(eta[p]^2) else if(is.null(xlab) && regress) bquote(bold(R^2)) else xlab
  if(!regress && missing(design)) stop("Error: 'design' must be numerically specified e.g., 'design = 2 * 4'.", call. = FALSE)
  if(regress){ n.level <- n.level + 1 ; design <- n.level }
  df1 <- n.level - 1
  if(n.covar < 0) n.covar <- 0
  x <- sapply(list(n.level, design, n.covar), round)
  n.level <- x[1] ; design <- x[2] ; n.covar <- x[3]
  
  f <- function(x){
    
    power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = (peta * (x + design + n.covar) ) /(1 - peta), lower.tail = FALSE))
  }
  
  df2 <- ceiling(uniroot(f, c(1e-8, 1e6), extendInt = "downX")[[1]])
  
  N <- df2 + design + n.covar
  
  loop <- length(peta2)
  
  Nb <- numeric(loop)
  df2b <- numeric(loop)
  
  for(i in 1:loop){
    
    f <- function(x){
      
      power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = (peta2[i] * (x + design + n.covar) ) /(1 - peta2[i]), lower.tail = FALSE))
    }
    
    df2b[i] <- ceiling(uniroot(f, c(1e-8, 1e6), extendInt = "downX")[[1]])
    
    
    Nb[i] <- df2b[i] + design + n.covar
    
  }
  
  a <- qpeta(sig.level, df1, df2, 0, N, lower.tail = FALSE)
  
  to <- if(is.null(to)) max(qpeta(.999999, df1, df2, 0, N), qpeta(.999999, df1, df2, peta, N), na.rm = TRUE) else to
  x <- seq(0, to, 1e-4)
  ylimb <- c(0, max(dpeta(x, df1, df2, 0, N), dpeta(x, df1, df2, peta, N), na.rm = TRUE))
  
  ylim <- if(is.infinite(ylimb[2]) & is.null(ylim)) NULL else if(is.null(ylim)) ylimb else ylim
  
  est.power <- ppeta(a, df1, df2, peta, N, lower.tail = FALSE)
  
par(mfrow = c(2, 1), mgp = c(2.5, .5, 0), mar = c(4, 4, 2, 2), tck = -.02)
    
  h0 <- curve(dpeta(x, df1, df2, 0, N), from = 0, to = to, n = 1e4, xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i", ylim = ylim, font.lab = 2)
  
  x = seq(a, to, length.out = 1e3) ; y = dpeta(x, df1, df2, 0, N)
  polygon(c(a, x, to), c(0, y, 0), col = adjustcolor(2, .25), border = NA)
  lines(h0, lwd = 2, col = 2, xpd = TRUE)
  abline(v = a, col = 2, lty = 2) ; crit <- round(a, 4) ; points(a, par('usr')[4], pch = 19, col = 2, xpd = NA)
  
  es <- if(regress) bquote(R^2) else bquote(eta[p]^2)
  
  text(a, par('usr')[4], bquote(bold("critical"~ .(es) == .(crit)~"or"~.(crit*1e2)*"%")), pos = 3, cex = .7, font = 4, xpd = TRUE)
  
  h1 <- curve(dpeta(x, df1, df2, peta, N), from = 0, to = to, n = 1e4, add = TRUE) # xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i", main = bquote(bolditalic(H[1]))
  x <- seq(a, to, length.out = 1e3) ; y <- dpeta(x, df1, df2, peta, N)
  polygon(c(a, x, to), c(0, y, 0), border = NA, density = 15, col = 4, xpd = TRUE)
  lines(h1, lwd = 2, col = 4, xpd = TRUE)
  
  legend("topleft", legend = c("Sig. Area", "Power"), inset = c(-.15, 0), density = c(NA, 35), x.intersp = c(.3, .3),
         bty = "n", xpd = NA, cex = .7, text.font = 2, angle = c(NA, 45), fill = c(adjustcolor(2, .4), 4), border = c(2, 4), adj = c(0, .4))
  
  plot(peta2, Nb, type = "b", lwd = 2, xlab = xlab, las = 1, col = "green4", font.lab = 2, xaxt = "n", ylab = "Total Sample Size")
  axis(1, at = peta2)
  
  legend("topright", legend = bquote(bold("Current required"~ bolditalic("\"N\""))), col = "magenta", pt.bg = "cyan", pch = 21, cex = .7,
         pt.cex = 1.2, bty = "n")
  box()
  points(peta, N, col = "magenta", bg = "cyan", pch = 21, cex = 1.5)
  
  text(peta2, Nb, Nb, pos = 3, font = 2, col = "gray40", cex = .8, xpd = TRUE)
  
  method <- paste("fixed-effects", if(regress) "Regression" else if(n.covar == 0) "ANOVA" else "ANCOVA", "power analysis") 
  
  balannced.N <- if(!regress) ceiling(N/design) * design else NA
  
  n.level <- if(regress) n.level-1 else n.level
  design <- if(regress) n.level else design
  
  r  <- structure(list(method, pov, est.power, a, sig.level, n.covar, design, n.level, df1, df2, N, balannced.N), class = "power.htest")
  
  setNames(r, c("method", ifelse(regress, "R-squared", "peta squared"), "est.power", ifelse(regress, "crit.Rsq", "crit.peta"), 
                "sig.level", "n.covar", "design", ifelse(regress, "n.pred", "n.level"), "df1", "df2", "total.N", "balanced.N"))
}

         
                                 
#=====================================================================================================================================
                  
                  
harmonic <- function(x, na.rm = TRUE, zero.rm = TRUE){

  if(zero.rm) {
    x[x == 0] <- NA
  } 
  
  if(is.null(nrow(x))){
  1 / mean(1/x, na.rm = na.rm)
} 
  else {
  1/(apply(1/x, 2, mean, na.rm = na.rm))
  }
}

#=====================================================================================================================================

                  
geometric <- function (x, na.rm = TRUE){
  
    if (is.null(nrow(x))) {
      exp(mean(log(x), na.rm = TRUE))
    }
    else {
      exp(apply(log(x), 2, mean, na.rm = na.rm))
    }
  }
                  
 
#======================================================================================================================================
                  
                  
#cell.makeup <- function(N, design)
#  {
#  y <- arrangements::partitions(N, design)
#  y <- y[nrow(y):1, ncol(y):1]
#  rownames(y) <- paste("form", 1:nrow(y))
#  colnames(y) <- paste0("group.", 1:ncol(y))
#  y
#}                  
         
                  
#================================================================================================================================
                  
d2f <- function(d, n1, n2) sqrt(d2peta(d, n1, n2) / (1 - d2peta(d, n1, n2) ))
f2d <- function(f, n1, n2) peta2d(f2peta(f), n1, n2)                  
peta2f <- function(peta) sqrt(peta / (1 - peta))
f2peta <- function(f) (f^2) / (1 + f^2)
peta2F <- function(peta, df1, df2) (peta / df1) / ((1 - peta)/df2)
F2peta <- function(F.value, df1, df2) (F.value*df1) / ((F.value*df1) + df2)
d2r <- function(d, n1 = 300, n2 = 300) sqrt((d^2) / ((d^2) + (((n1 + n2)^2) - (2*(n1 + n2))) / (n1 * n2)))
r2d <- function(r, n1 = 300, n2 = 300) sqrt((r^2)*(((n1 + n2)^2)-(2*(n1 + n2)))/(n1 * n2)/(1-(r^2)))
d2t <- function(d, n1, n2 = NA){
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  d*sqrt(N)
}

t2d <- function(t, n1, n2 = NA){
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  t/sqrt(N)
}
                  
ncp2peta <- function(ncp, N) { ncp / (ncp + N) }

   
peta2ncp <- function(peta, N) { (peta*N) / (1 - peta) }
                  
peta2N <- function(peta, ncp) { (ncp - (peta * ncp)) / peta }
                  
d2peta <- function(d, n1 = 300, n2 = 300) (d^2) / ((d^2) + (((n1 + n2)^2) - (2*(n1 + n2))) / (n1 * n2))
                  
peta2d <- function(peta, n1 = 300, n2 = 300) sqrt((peta)*(((n1 + n2)^2)-(2*(n1 + n2)))/(n1 * n2)/(1-(peta)))
                  
F2pomega <- function(F.value, df1, N){
  (df1 * (F.value - 1)) / ((df1 * (F.value - 1)) + N)
}

pomega2F <- function(pomega, df1, N) {
  1 - ( (N * pomega )/(df1 * (pomega - 1)) )
}


peta2pomega <- function(peta, df1, df2, N){
  f <- peta2F(peta, df1, df2)  
  F2pomega(f, df1, N)
}


pomega2peta <- function(pomega, df1, df2, N){
  f <- pomega2F(pomega, df1, N)  
  F2peta(f, df1, df2)
}
   

exp.pov <- exp.peta <- Vectorize(function(pbase = 0, df1, df2, N){
  
  integrate(function(x, df1, df2, pbase, N){
    
    x * dpeta(x = x, df1 = df1, df2 = df2, pbase = pbase, N = N)
    
  }, 0, 1, df1 = df1, df2 = df2, pbase = pbase, N = N)[[1]]
  
})
                  

exp.d <- Vectorize(function(dbase = 0, n1, n2 = NA){
  
  integrate(function(x, n1, n2, dbase){
    
    x * dcohen(x = x, dbase = dbase, n1 = n1, n2 = n2)
    
  }, -Inf, Inf, n1 = n1, n2 = n2, dbase = dbase)[[1]]
  
})
                  
                  
exp2peta <- Vectorize(function(exp.val, df1, df2, N){
  
  optimize(function(x){
    
    abs(exp.val - exp.peta(pbase = x, df1 = df1, df2 = df2, N = N))
    
  }, 0:1, tol = 1e-9)[[1]]
  
})


exp2d <- Vectorize(function(exp.val, n1, n2 = NA){
  
  uniroot(function(x){
    
    exp.val - exp.d(dbase = x, n1 = n1, n2 = n2)
    
  }, c(-4, 4), extendInt = "yes")[[1]]
})
                  
#==================================================================================================================================
                  
                  
f.balance <- function(F.unbalance, cell.makeup, df1, df2, N, conf.level = .9)
{
  
  fbalance <- F.unbalance * (mean(cell.makeup) / harmonic(cell.makeup))
  
  ci <- peta.ci(f = c(fbalance, F.unbalance), df1 = df1, df2 = df2, N = N, conf.level = conf.level)
  l <- length(F.unbalance)
  rownames(ci) <- paste((2*l):1, c(rep("balanced", l), rep("Unbalace", l)))
  ci[nrow(ci):1,]
}
                  
                  
#====================================================================================================================================
                  
                  
is.whole <- function(x)  abs(x - round(x)) < .Machine$double.eps^.5

                  
#====================================================================================================================================
                  

plan.rep.measure <- function(peta, n.rep, n.group, factor.type = c("between", "within", "bw"), sig.level = .05, n.covar = 0, power = .8, eps = .9,
                             peta.range = seq(1e-1, .9, 1e-1), rho = .5, xlab = NULL, ylim = NULL, to = NULL, d = NA)
{
  
  UseMethod("plan.rep.measure")
}


plan.rep.measure.default <- function(peta, n.rep, n.group, factor.type = c("between", "within", "bw"), sig.level = .05, n.covar = 0, power = .8, eps = .9,
                                     peta.range = seq(1e-1, .9, 1e-1), rho = .5, xlab = NULL, ylim = NULL, to = NULL, d = NA){
  
  graphics.off()  
  original.par <- par(no.readonly = TRUE)
  on.exit(par(original.par))
  options(warn = -1)
  if(!is.na(d)) { peta <- d2peta(d = d, n1 = 300, n2 = 300) ;
  message("\nNote: For 'pairwise' comparisons, 'total.N' is for '2' groups.") }
  
  m <- n.rep
  
  peta2 <- peta.range
  
  peta2[peta2 == 0] <- 1e-2
  peta2[peta2 == 1] <- .99
  
  factor.type <- match.arg(factor.type)
  
  if(rho <= 0) rho <- 1e-7 else if(rho >= 1) rho <-.9999999
  if(eps < .5) eps <- .5 else if(eps > 1) eps <- 1
  if(n.group < 1) stop("Error: You must have at least '1 group' in your design.")
  if(m < 1) stop("Incorrect # of measurements.", call. = FALSE)
  if(factor.type != "between" & m < 2) stop("Error: You must have at least '2 repeated measurements' in your design.", call. = FALSE)
  xlab <- if(is.null(xlab)) bquote(eta[p]^2) else xlab
  if(missing(n.group)) stop("Error: 'n.group' must be numerically specified.", call. = FALSE)
  
  df1 <- switch(factor.type, between = n.group - 1, within = (m - 1)*eps, bw = (n.group - 1)*(m - 1)*eps)
  
  if(n.covar < 0) n.covar <- 0
  g <- sapply(list(n.group, n.covar, m), round)
  n.group <- g[1] ; n.covar <- g[2] ; m <- g[3]
  
  u <- if(factor.type == "between") m / (1 + (m - 1)*rho) else m / (1 - rho)
  
  f <- if(factor.type == "between"){ function(x){
    
    power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = ((peta * ( x + n.group + n.covar) ) /(1 - peta))*u, lower.tail = FALSE))
  } 
    
  } else {
    
    function(x){ 
      power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = ((peta * ( ((x)/(m-1)) + n.group + n.covar) ) /(1 - peta))*eps*u, lower.tail = FALSE))
    }
  }
  
  df2 <- uniroot(f, c(1e-8, 1e6), extendInt = "downX")[[1]]
    
  df2 <- if(factor.type == "between") ceiling(df2) else df2
      
  N <- if(factor.type == "between") ceiling(df2 + n.group) + n.covar else ceiling((df2 / ((m - 1)*eps)) + n.group) + n.covar
  
  
  loop <- length(peta2)
  
  Nb <- numeric(loop)
  df2b <- numeric(loop)
  
  
  for(i in 1:loop){
    
    f <- if(factor.type == "between"){ function(x){
      
      power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = ((peta2[i] * ( x + n.group + n.covar) ) /(1 - peta2[i]))*u, lower.tail = FALSE))
    } 
      
    } else {
      
      function(x){ 
        power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = ((peta2[i] * ( ((x)/(m-1)) + n.group + n.covar) ) /(1 - peta2[i]))*eps*u, lower.tail = FALSE))
      }
    }
    
    df2b[i] <- uniroot(f, c(1e-8, 1e6), extendInt = "downX")[[1]]
        
    Nb[i] <- if(factor.type == "between") ceiling(df2b[i] + n.group) + n.covar else ceiling((df2b[i] / ((m - 1)*eps)) + n.group) + n.covar
    
  }
  
  a <- qpetab(sig.level, df1, df2, 0, lower.tail = FALSE)
  
  ncp <- if(factor.type == "between") (peta2f(peta)^2)*N*u else (peta2f(peta)^2)*N*u*eps
  
  to <- if(is.null(to)) max(qpetab(.999999, df1, df2, 0), qpetab(.999999, df1, df2, ncp), na.rm = TRUE) else to
  x <- seq(0, to, 1e-4)
  ylimb <- c(0, max(dpetab(x, df1, df2, 0), dpetab(x, df1, df2, ncp), na.rm = TRUE))
  
  ylim <- if(is.infinite(ylimb[2]) & is.null(ylim)) NULL else if(is.null(ylim)) ylimb else ylim
  
  est.power <- ppetab(a, df1, df2, ncp, lower.tail = FALSE)
  
  par(mfrow = c(2, 1), mgp = c(2.5, .5, 0), mar = c(4, 4, 2, 2))
  
  h0 <- curve(dpetab(x, df1, df2, 0), from = 0, to = to, n = 1e4, xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i", ylim = ylim, font.lab = 2)
  
  x = seq(a, to, length.out = 1e3) ; y = dpetab(x, df1, df2, 0)
  polygon(c(a, x, to), c(0, y, 0), col = adjustcolor(2, .25), border = NA)
  lines(h0, lwd = 2, col = 2, xpd = TRUE)
  abline(v = a, col = 2, lty = 2) ; crit <- round(a, 4) ; points(a, par('usr')[4], pch = 19, col = 2, xpd = NA)
  
  text(a, par('usr')[4], bquote(bold("critical"~ eta[p]^2 == .(crit)~"or"~.(crit*1e2)*"%")), pos = 3, cex = .7, font = 4, xpd = TRUE)
  
  h1 <- curve(dpetab(x, df1, df2, ncp), from = 0, to = to, n = 1e4, add = TRUE)
  x <- seq(a, to, length.out = 1e3) ; y <- dpetab(x, df1, df2, ncp)
  polygon(c(a, x, to), c(0, y, 0), border = NA, density = 15, col = 4, xpd = TRUE)
  lines(h1, lwd = 2, col = 4, xpd = TRUE)
  
  legend("topleft", legend = c("Sig. Area", "Power"), inset = c(-.15, 0), density = c(NA, 35), x.intersp = c(.3, .3),
         bty = "n", xpd = NA, cex = .7, text.font = 2, angle = c(NA, 45), fill = c(adjustcolor(2, .4), 4), border = c(2, 4), adj = c(0, .4))
  
  plot(peta2, Nb, type = "b", lwd = 2, xlab = xlab, las = 1, col = "green4", font.lab = 2, xaxt = "n", ylab = "Total Sample Size")
  axis(1, at = peta2)
  
  legend("topright", legend = bquote(bold("Current required"~ bolditalic("\"N\""))), col = "magenta", pt.bg = "cyan", pch = 21, cex = .7,
         pt.cex = 1.2, bty = "n")
  box()
  points(peta, N, col = "magenta", bg = "cyan", pch = 21, cex = 1.5)
  
  text(peta2, Nb, Nb, pos = 3, font = 2, col = "gray40", cex = .8, xpd = TRUE)
  
  method <- paste("fixed-effects repeated-measures", if(n.covar == 0) "ANOVA" else "ANCOVA", "power analysis") 
  
  bal <- ceiling(N/n.group) * n.group
  
  note <- if(n.group != 0 & N %% n.group != 0) paste("We suggest recruiting", bal, "subjects (instead of", N, "subjects) to achieve",  bal/n.group, "subjects per group.")
  
  message("\nIMPORTANT: Always pick the factor with largest # of levels to obtain required 'total.N'.")
  
  r  <- structure(list(factor.type, est.power, a, sig.level, n.covar, n.group, m, df1, df2, N, method, note), class = "power.htest")
  
  setNames(r, c("factor.type", "est.power", "crit.peta", 
                "sig.level", "n.covar", "n.group", "n.rep", "df1", "df2", "total.N", "method", "note"))
}
                  

#================================================================================================================================
                  
                  
 peta.rep.bayes <- function(f = NULL, peta, N, df1, df2, n.rep, factor.type = "between", a = 1.2, b = 1.2, level = .95, lo = 0, hi = 1, dist.name = "dbeta", rho = .5, eps = .9, scale = .1, top = 1.1, show.prior = FALSE, 
                                   bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .05, peta.h0 = 0, digits = 6, col.depth = .55, labels = NULL, cex.lab = .8, 
                                   xlab = NULL, ylab = NULL, col.hump = NULL, ...)
  {

  UseMethod("peta.rep.bayes")
}


peta.rep.bayes.default <- function(f = NULL, peta, N, df1, df2, n.rep, factor.type = "between", a = 1.2, b = 1.2, level = .95, lo = 0, hi = 1, dist.name = "dbeta", rho = .5, eps = .9, scale = .1, top = 1.1, show.prior = FALSE, 
                               bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .05, peta.h0 = 0, digits = 6, col.depth = .55, labels = NULL, cex.lab = .8, 
                               xlab = NULL, ylab = NULL, col.hump = NULL, ...){
  
  if(is.null(f) & missing(peta)) stop("Error: Either 'f' or 'peta' must be provided.")
  if(!is.null(f) & !missing(peta)) stop("Error: Only one of 'f' or 'peta' must be provided.")
  f <- if(is.null(f)) peta2F(peta, df1, df2) else f
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)      
  pr <- show.prior
  
  
  if(rho <= 0) rho <- 1e-7 else if(rho >= 1) rho <-.9999999
  if(eps < .5) eps <- .5 else if(eps > 1) eps <- 1 else eps
  
  
  if(!pr){
    m <- n.rep
    I <- eq(a, b, d, lo, hi, f, N, df1, df2, m, factor.type)
    a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]] ; f = I[[6]] ; N = I[[7]] ; df1 = I[[8]] ; df2 = I[[9]] ; m = I[[10]] ; factor.type = I[[11]]
    
    u <- if(factor.type == "between") m / (1 + (m - 1)*rho) else m / (1 - rho)
    
    loop <- length(a)
    CI <- matrix(NA, loop, 2)
    mode <- numeric(loop)
    peak <- numeric(loop)
    h <- list()
    k = numeric(loop)
    eq.prob = numeric(loop)
    BF10 = numeric(loop)
    estimate = numeric(loop)
    
    options(warn = -1)
    
    for(i in 1:loop){
      
      p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
      prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]  
      likelihood = function(x) df(f[i], df1[i], df2[i], ((x * N[i]) / (1 - x))*if(factor.type == "between") u[i] else u[i]*eps[i] )
      k[i] = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k[i]
      h[[i]] = list(x = x <- seq(0, 1, length.out = 5e2), y = posterior(x))
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      peak[i] = posterior(mode[i])
      CI[i,] = HDI(posterior, 0, .9999999, level = level)
      BF10[i] =  k[i] / df(f[i], df1[i], df2[i], (peta.h0 * N[i]) / (1 - peta.h0))
      eq.prob[i] = integrate(posterior, lo[i], eq.hi)[[1]] - integrate(posterior, lo[i], eq.lo)[[1]]
      estimate[i] <- (f[i]*df1[i]) / ((f[i]*df1[i]) + df2[i])
    } 
    graphics.off()
    lab <- if(is.null(labels)) substring(d, 2) else labels
    xlab <- if(is.null(xlab)) bquote(bold("Credible Interval"~(eta[p]^2))) else xlab
    ylab <- if(is.null(ylab)) NA else ylab  
    
    plot(CI, rep(1:loop, 2), type = "n", xlim = 0:1, ylim = c(bottom*1, top*loop), ylab = ylab, yaxt = "n", xaxt = "n", xlab = xlab, font.lab = 2, mgp = c(2, .5, 0), ...)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), labels = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0)) 
    axis(2, at = 1:loop, labels = lab, font = 2, las = 1, cex.axis = cex.lab, tck = -.006, mgp = c(2, .3, 0), padj = rep(.3, loop))
    
    for(i in 1:loop){
      col <- if(is.null(col.hump)) i else col.hump[i]     
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(col, col.depth), border = NA, xpd = NA)
    }
    m <- scale*peak + 1:loop
    col <- if(is.null(col.hump)) 1:loop else col.hump
    legend(x = leg, legend = rev(paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")")), pch = 22, title = "Priors", pt.bg = rev(col), col = rev(col), cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4, adj = c(0, .3))
    box()
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = col, xpd = NA)                   
    segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I <- deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0(I[,1], "%", "    ", o, "%", "    ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
    
    rownames <- if(is.null(labels)) paste0("P.eta.sq ", 1:loop, " posterior:") else paste0(1:loop, " ", labels, " posterior:")                   
    return(round(data.frame(estimate = estimate, mode = mode, lower = CI[,1], upper = CI[,2], eq.prob = eq.prob, BF10 = BF10, row.names = rownames), digits = digits))  
    
  }else{
    xlab <- if(is.null(xlab)) bquote(bold("Partial Eta.Sq"~(eta[p]^2))) else xlab
    p <- function(x) { get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = xlab, bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(eta[p]^2*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))), yaxs = "i", xpd = TRUE)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}
                                                                                                                                             
                                                                                                                                             
#=======================================================================================================================================================
                                                                                                                                             
                                                                                                                                             
plan.d.cib <- function(d, conf.level = .95, width, assure = NULL, paired = FALSE)
{
  UseMethod("plan.d.cib")
}


plan.d.cib.default <- function(d, conf.level = .95, width, assure = NULL, paired = FALSE)
{ 
  
n.d.d <- function(d, conf.level, width, paired){ 
  
  n.d <- Vectorize(function(d, conf.level, width, paired) 
  {
    alpha <- (1 - conf.level)/2
    
    n0 <- if(paired) 1 else 2 * (qnorm(1 - alpha)/(width/2))^2
    n <- if(paired) ceiling(n0) else 2 * ((qt(1 - alpha, 2 * n0 - 2))/(width/2))^2
    if(!paired) dif <- abs(n - n0)
    
    if(!paired){
      
      while(dif > 1e-6){
        np <- n
        n <- 2 * ((qt(1 - alpha, 2 * n - 2))/(width/2))^2
        dif <- abs(n - np)
      }
      n <- ceiling(n)
      n <- max(4, n - 5)
      
    } else {
      
      n <- max(4, n - 5)
    }
    
    ci <- d.ci(d = d, n1 = n, n2 = if(paired) NA else n, conf.level = conf.level)
    dif.full <- abs(ci$upper - ci$lower) - width
    
    while(dif.full > 0){
      n <- n + 1
      cis <- d.ci(d = d, n1 = n, n2 = if(paired) NA else n, conf.level = conf.level)
      width.now <- abs(cis$upper - cis$lower)
      dif.full <- width.now - width
    }
    
    return(c(d = d, n = n, width = width, conf.level = conf.level))
  })
  
  data.frame(t(n.d(d = d, width = width, conf.level = conf.level, paired = paired)), paired = paired, row.names = NULL)
}    
    
if(is.null(assure)){
  
  n.d.d(d = d, conf.level = conf.level, width = width, paired = paired)
  
} else {
  
n.as <- Vectorize(function(d, conf.level, width, assure, paired){
  
    n0 <- n.d.d(d = d, conf.level = conf.level, width = width, paired = paired)$n
    
    a <- d.ci(d = d, n1 = n0, n2 = if(paired) NA else n0, conf.level = assure)$upper
    b <- d.ci(d = d, n1 = n0, n2 = if(paired) NA else n0, conf.level = assure - (1 - assure))$upper
    
    limits <- function(limit.now = limit.now, nn = n0, d = d, assure = assure){
      lower <- pcohen(-limit.now, d, n1 = nn, n2 = if(paired) NA else nn)
      upper <- pcohen( limit.now, d, n1 = nn, n2 = if(paired) NA else nn, lower.tail = FALSE)
      total <- lower + upper
      return((total - (1 - assure))^2)
    }
    
    d.opt <- optimize(limits, c(a, b), d = d, assure = assure)[[1]]
    n <- n.d.d(d = d.opt, conf.level = conf.level, width = width, paired = paired)
    return(n)
    })

  data.frame(t(n.as(d = d, conf.level = conf.level, width = width, assure = assure, paired = paired)), assure = assure, row.names = NULL)
  }
}
                      
                      
#================================================================================================================================
                      
                      
power.t <- function(d = .1, sig.level = .05, power = .8, base.rate = 1, paired = FALSE, two.tailed = TRUE)
{
  
 pwr <- Vectorize(function(d, sig.level, power, base.rate, paired, two.tailed)
   {
   
d <- abs(d)
sig.level <- if(two.tailed) sig.level/2 else sig.level
k <- base.rate / (1 + base.rate)
options(warn = -1)

f <- if(two.tailed){ function(x){
  
  power - (pt(qt(sig.level, df = x), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2))) + pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE))
}
  
} else {
  
  function(x){
    
  power - pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE)
    
  }
}

df <- ceiling(uniroot(f, c(1e-8, 1e7), extendInt = "downX")[[1]])

n1 <- df + 1
n2 <- if(paired) NA else round(base.rate*n1)

return(c(n1 = n1, n2 = n2))
})

 data.frame(t(pwr(d = d, sig.level = sig.level, power = power, base.rate = base.rate, paired = paired, two.tailed = two.tailed)))
}
             

               
#=================================================================================================================
             
plan.t.ci <- function(d, t = NA, n1, n2 = NA, conf.level = .95, width = NA, base.rate = 1, paired = FALSE, assure = .99, expect = FALSE, reduce.by = "0%", increase.by = "0%")
{
  UseMethod("plan.t.ci")
}


plan.t.ci.default <- function(d, t = NA, n1, n2 = NA, conf.level = .95, width = NA, base.rate = 1, paired = FALSE, assure = .99, expect = FALSE, reduce.by = "0%", increase.by = "0%"){
  
  if(any(conf.level >= 1) || any(conf.level <= 0) || any(assure >= 1) || any(assure <= 0)) stop("'conf.level' and 'assure' must be between '0' and '1'.", call. = FALSE)
  if(is.na(width) & missing(n1) || is.na(width) & is.na(t) & missing(d) || is.na(width) & !paired & missing(n2)) stop("Either provide 'width' or provide 't or d', 'n1' and/or 'n2' from prior study.", call. = FALSE)  
  if(!is.na(t)) d <- t2d(t = t, n1 = n1, n2 = n2)
  if(is.na(width)) width <- d.width(d = d, t = t, n1 = n1, n2 = n2, conf.level = conf.level)
  if(expect) assure <- .5
  
  inc <- if(is.character(increase.by)) as.numeric(substr(increase.by, 1, nchar(increase.by)-1))/ 1e2 else increase.by
  red <- if(is.character(reduce.by)) as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2 else reduce.by
  
  fac <- if(inc != 0 & red == 0) { 1 + inc
  } else if(red != 0 & inc == 0) { 1 - red 
  } else { 1 }
  
  
  if(fac <= 0 || inc == 0 & fac > 1) fac <- 1
  
  width <- width * fac
  
  G <- Vectorize(function(d, conf.level, width, base.rate, paired, assure, expect){
    
    n.d <- function(d, conf.level, width, base.rate, paired, assure){
      
      alpha <- (1 - conf.level)/2
      k <- base.rate 
      
      f <- function(ncp, alpha, d, df){
        alpha - suppressWarnings(pt(d*sqrt(if(paired) df + 1 else ((k/(1 + k))^2)*(df + 2)), df, ncp, lower.tail = FALSE))
      }
      
      dbase <- function(df){
        sapply(c(alpha, 1 - alpha),
               function(x) uniroot(f, c(-5e1, 5e1), alpha = x, d = d, df = df, extendInt = "yes")[[1]]/sqrt(if(paired) df + 1 else ((k/(1 + k))^2)*(df + 2)))
      }
      
      m <- function(df, width){
        abs(abs(diff(dbase(df))) - width)
      }
      
      df <- optimize(m, c(1, 1e7), width = width)
      
      if(round(df$objective, 4) != 0) plan.t.ci <- function(d, t = NA, n1, n2 = NA, conf.level = .95, width = NA, base.rate = 1, paired = FALSE, assure = .99, expect = FALSE, reduce.by = "0%", increase.by = "0%")
{
  UseMethod("plan.t.ci")
}


plan.t.ci.default <- function(d, t = NA, n1, n2 = NA, conf.level = .95, width = NA, base.rate = 1, paired = FALSE, assure = .99, expect = FALSE, reduce.by = "0%", increase.by = "0%"){
  
  if(any(conf.level >= 1) || any(conf.level <= 0) || any(assure >= 1) || any(assure <= 0)) stop("'conf.level' and 'assure' must be between '0' and '1'.", call. = FALSE)
  if(is.na(width) & missing(n1) || is.na(width) & is.na(t) & missing(d) || is.na(width) & !paired & missing(n2)) stop("Either provide 'width' or provide 't or d', 'n1' and/or 'n2' from prior study.", call. = FALSE)  
  if(!is.na(t)) d <- t2d(t = t, n1 = n1, n2 = n2)
  if(is.na(width)) width <- d.width(d = d, t = t, n1 = n1, n2 = n2, conf.level = conf.level)
  if(expect) assure <- .5
  
  inc <- if(is.character(increase.by)) as.numeric(substr(increase.by, 1, nchar(increase.by)-1))/ 1e2 else increase.by
  red <- if(is.character(reduce.by)) as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2 else reduce.by
  
  fac <- if(inc != 0 & red == 0) { 1 + inc
  } else if(red != 0 & inc == 0) { 1 - red 
  } else { 1 }
  
  
  if(fac <= 0 || inc == 0 & fac > 1) fac <- 1
  
  width <- width * fac
  
  G <- Vectorize(function(d, conf.level, width, base.rate, paired, assure, expect){
    
    n.d <- function(d, conf.level, width, base.rate, paired, assure){
      
      alpha <- (1 - conf.level)/2
      k <- base.rate 
      
      f <- function(ncp, alpha, d, df){
        alpha - suppressWarnings(pt(d*sqrt(if(paired) df + 1 else ((k/(1 + k))^2)*(df + 2)), df, ncp, lower.tail = FALSE))
      }
      
      dbase <- function(df){
        sapply(c(alpha, 1 - alpha),
               function(x) uniroot(f, c(-5e1, 5e1), alpha = x, d = d, df = df, extendInt = "yes")[[1]]/sqrt(if(paired) df + 1 else ((k/(1 + k))^2)*(df + 2)))
      }
      
      m <- function(df, width){
        abs(abs(diff(dbase(df))) - width)
      }
      
      df <- optimize(m, c(1, 1e7), width = width)
      
      if(round(df$objective, 4) != 0) stop("Impossible planning: change input values.", call. = FALSE)
      
      n1 <- ceiling(if(paired) df[[1]] + 1 else (df[[1]] + 2)/(1 + k))
      n2 <- if(paired) NA else round(k * n1)
      
      list(d = d, n1 = n1, n2 = n2, base.rate = base.rate, width = width, conf.level = conf.level, assure = assure, paired = paired)
    }
    
    n <- n.d(d = d, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure)
    
    a <- d.ci(d = d, n1 = n$n1, n2 = n$n2, conf.level = c(assure, 2*assure - 1))$upper
    
    dnew <- function(dnew = dnew, n1 = n$n1, n2 = n$n2, d = d, assure = assure){
      total <- sum(pcohen(c(-dnew, dnew), d, n1 = n1, n2 = n2, lower.tail = c(TRUE, FALSE)))
      return(abs(total - (1 - assure)))
    }
    
    dnew <- optimize(dnew, a, d = d, assure = assure)[[1]]
    n.d(d = dnew, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure)
  })
  
  if(paired) base.rate <- NA
  a <- data.frame(t(G(d = d, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure, expect = expect)), row.names = NULL)
  a[,1] <- d
  a                                                                                                          
}
      
      n1 <- ceiling(if(paired) df[[1]] + 1 else (df[[1]] + 2)/(1 + k))
      n2 <- if(paired) NA else round(k * n1)
      
      list(d = d, n1 = n1, n2 = n2, base.rate = base.rate, width = width, conf.level = conf.level, assure = assure, paired = paired)
    }
    
    n <- n.d(d = d, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure)
    
    a <- d.ci(d = d, n1 = n$n1, n2 = n$n2, conf.level = c(assure, 2*assure - 1))$upper
    
    dnew <- function(dnew = dnew, n1 = n$n1, n2 = n$n2, d = d, assure = assure){
      total <- sum(pcohen(c(-dnew, dnew), d, n1 = n1, n2 = n2, lower.tail = c(TRUE, FALSE)))
      return(abs(total - (1 - assure)))
    }
    
    dnew <- optimize(dnew, a, d = d, assure = assure)[[1]]
    n.d(d = dnew, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure)
  })
  
  if(paired) base.rate <- NA
  a <- data.frame(t(G(d = d, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure, expect = expect)), row.names = NULL)
  a[,1] <- d
  a                                                                                                          
}
                                                                                                                
                                                                                                                                                                                                          
#===================================================================================================================================================
                                                                                                     
 
plan.f.ci <- function(pov, design = 2 * 2, n.level = 2, n.pred = NULL, n.covar = 0, conf.level = .95, width = NA, assure = .99, expect = FALSE, reduce.by = "0%", d = NA, lower, upper, increase.by = "0%", tol = 1e3)
{
  
  UseMethod("plan.f.ci")
  
}

plan.f.ci.default <- function(pov, design = 2 * 2, f = NA, n.level = 2, n.pred = NULL, n.covar = 0, conf.level = .95, width = NA, assure = .99, expect = FALSE, reduce.by = "0%", d = NA, lower, upper, increase.by = "0%", tol = 1e3){
  
  
  if(any(conf.level >= 1) || any(conf.level <= 0) || any(assure >= 1) || any(assure <= 0)) stop("'conf.level' and 'assure' must be between '0' and '1'.", call. = FALSE)
  if(expect) assure <- .5
  regress <- if(!is.null(n.pred)) TRUE else FALSE
  if(regress) n.level <- n.pred
  if(!is.na(d)) { pov <- d2peta(d = d, n1 = 500, n2 = 500) ; n.level <- 2 ;
  message("\nNote: For 'pairwise' comparisons, 'total.N' is for '2' groups.") }
  if(!is.na(d)) { pov <- d2peta(d = d, n1 = 500, n2 = 500) ; n.level <- 2 }
  if(!is.na(d) & is.na(width)) width <- d.width.meta(lower = lower, upper = upper)
  if(!is.na(d) & width >= .3) width <- .3
  if(!is.na(d) & pov <= .15) pov <- .15
  
  peta <- pov
  
  inc <- if(is.character(increase.by)) as.numeric(substr(increase.by, 1, nchar(increase.by)-1))/ 1e2 else increase.by
  red <- if(is.character(reduce.by)) as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2 else reduce.by
  
  fac <- if(inc != 0 & red == 0) { 1 + inc
  } else if(red != 0 & inc == 0) { 1 - red 
  } else { 1 }
  
  
  if(fac <= 0 || inc == 0 & fac > 1) fac <- 1
  
  width <- width * fac
  
  
  G <- Vectorize(function(peta, conf.level, width, assure, design, n.level, n.covar, regress, expect){
    
    
    n.f <- function(peta, conf.level, width, assure, design, n.level, n.covar, regress){
      
      alpha <- (1 - conf.level)/2
      if(regress){ n.level <- n.level + 1 ; design <- n.level }
      df1 <- n.level - 1
      if(n.covar < 0) n.covar <- 0
      options(warn = -1)
      
      f <- function(alpha, q, df1, df2, ncp){
        alpha - suppressWarnings(pf(peta2F(peta, df1, df2), df1, df2, ncp, lower.tail = FALSE))
      }
      
      pbase <- function(df2){      
        
        b <- sapply(c(alpha, 1 - alpha), function(x) 
          tryCatch(uniroot(f, c(0, 1e7), alpha = x, q = q, df1 = df1, df2 = df2, extendInt = "yes")[[1]], error = function(e) NA))
        if(any(is.na(b))) b <- c(1, tol)     
        ncp2peta(b, df2 + design + n.covar)
      }
      
      m <- function(df2, width){
        abs(diff(pbase(df2))) - width
      }
      
      df2 <- uniroot(m, c(0, 1e3), width = width, extendInt = "yes")
      
      
      if(round(df2$f.root, 3) != 0) stop("\nImpossible planning: You may change your 'width'.", call. = FALSE)
      
      
      df2 <- ceiling(df2[[1]])
      
      N <- ceiling(df2 + design + n.covar)
      n.covar <- if(n.covar == 0) NA else n.covar
      n.level <- if(regress) n.level-1 else n.level
      design <- if(regress) n.level else design
      df1 <- if(regress) n.level else df1
      
      list(peta = peta, total.N = N, width = width, n.level = n.level, conf.level = conf.level, assure = assure, df1 = df1, df2 = df2, design = design)
    }
    
    n <- n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)  
    
    peta <- exp.peta(pbase = n$peta, df1 = n$df1, df2 = n$df2, N = n$total.N)
    
    n <- n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)
    
    peta.max <- root(pov = peta, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = conf.level)$m
    
    a <- peta.ci(peta = peta, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = 2*assure - 1)
    
    nLU <- sapply(c(a$lower, a$upper), function(x) n.f(peta = x, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)$total.N)
    
    NN1 <- max(nLU, na.rm = TRUE)
    
    b <- peta.ci(peta = peta.max, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = 1 - assure)
    
    nLU <- sapply(c(b$lower, b$upper), function(x) n.f(peta = x, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)$total.N)
    
    NN2 <- max(nLU, na.rm = TRUE)
    
    NN3 <- if(!(peta.max %inn% c(a$lower, a$upper))) NN1 else max(NN1, NN2)
    
    max.w <- round(root(pov = peta, df1 = n$df1, N = NN3, df2 = NN3 - design - n.covar, conf.level = conf.level)$w, 3)
    
    return(c(peta = peta, total.N = NN3, width = width, n.level = n.level, design = design, conf.level = conf.level, max.width = max.w))
    
  })
  
  a <- data.frame(t(G(peta = peta, conf.level = conf.level, width = width, design = design, n.level = n.level, n.covar = n.covar, assure = assure, regress = regress, expect = expect)), regress = regress, assure = assure, row.names = NULL)
  names(a)[4] <- if(regress) "n.pred" else "n.level"
  names(a)[1] <- if(regress) "R2" else if(!is.na(d)) "d" else "peta"
  a[, 1] <- if(is.na(d)) pov else d
  a
}                                                                                  
                                                                                                                
                                                                                                     
#=====================================================================================================================================================
                                                                                                     
plan.f.cib <- function(peta = .2, design = 2 * 2, n.level = 2, n.covar = 0, conf.level = .9, width = .2, regress = FALSE, n.groups = 0, assure = .99){

  if(any(conf.level >= 1) || any(conf.level <= 0) || any(assure >= 1) || any(assure <= 0)) stop("'conf.level' and 'assure' must be between '0' and '1'.", call. = FALSE)
  
  G <- Vectorize(function(peta, conf.level, width, assure, design, n.level, n.covar, regress, n.groups){
  
n.f <- function(peta, conf.level, width, assure, design, n.level, n.covar, regress, n.groups){
  
  alpha <- (1 - conf.level)/2
  if(regress){ n.level <- n.level + 1 ; design <- n.level }
  df1 <- n.level - 1
  if(n.covar < 0) n.covar <- 0
  options(warn = -1)
  
  f <- function(alpha, q, df1, df2, ncp){
   abs(alpha - suppressWarnings(pf((peta / df1) / ((1 - peta)/df2), df1, df2, ncp, lower.tail = FALSE))) 
  }
  
  pbase <- function(df2){
  
    g <- try(uniroot(f, c(0, 1e7), alpha = alpha, q = q, df1 = df1, df2 = df2)[[1]], silent = TRUE)
    if(inherits(g, "try-error")) g <- 0
    h <- try(uniroot(f, c(0, 1e7), alpha = 1-alpha, q = q, df1 = df1, df2 = df2)[[1]], silent = TRUE)
    if(inherits(h, "try-error")) h <- 0
    b <- c(g, h)
    
    b / (b + (df2 + design))
  }
  
  m <- function(df2, width){
    abs(abs(diff(pbase(df2))) - width)
  }
  
  df2 <- optimize(m, c(1, 1e4), width = width)
  if(round(df2$objective, 4) != 0) return(c(NaN, message("Error: NaN produced. Are input values correct?")))
  
  df2 <- ceiling(df2[[1]] - n.covar)
  
  N <- ceiling(df2 + design)
  bal <- ceiling(N/design) * design
  if(n.groups != 0){ N <- n.groups * (bal/2) ; message("\nNote: You are doing reseach planning for 'pairwise' comparisons.") }
  N <- if(design != 0 & N %% design != 0) bal else N
  n.covar <- if(n.covar == 0) NA else n.covar
  n.level <- if(regress) n.level-1 else n.level

list(peta = peta, total.N = N, width = width, n.level = n.level, conf.level = conf.level, assure = assure, df1 = df1, df2 = df2)
}

n <- n.f(peta = peta, conf.level = conf.level, width = width, design = design, n.level = n.level, n.covar = n.covar, regress = regress, n.groups = n.groups, assure = assure)

a <- peta.ci(peta = peta, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = c(assure, assure - (1 - assure)))$upper

petanew <- function(petanew = petanew, df1 = n$df1, df2 = n$df2, peta = n$peta, N = n$total.N, assure = assure){
  total <- sum(ppeta(c(1-petanew, petanew), peta, df1 = df1, df2 = df2, N = N, lower.tail = c(TRUE, FALSE)))
  return(abs(total - (1 - assure)))
}

petanew <- optimize(petanew, a, peta = peta, assure = assure)[[1]]
n.f(peta = petanew, conf.level = conf.level, width = width, design = design, n.level = n.level, n.covar = n.covar, regress = regress, n.groups = n.groups, assure = assure)
})
  
data.frame(t(G(peta = peta, conf.level = conf.level, width = width, design = design, n.level = n.level, n.covar = n.covar, regress = regress, n.groups = n.groups, assure = assure)), row.names = NULL)[, 1:6, drop = FALSE]
}
                
#=====================================================================================================================================================================================================

plan.t.cic <- function(d = .4, conf.level = .95, width = .2, base.rate = 1, paired = FALSE, assure = .99){
  
  k <- base.rate
  
G <- Vectorize(function(d, conf.level, width, base.rate, paired, assure){
    
  nd <- function(d, width, conf.level, base.rate){
    
    f <- function(n1, n2){
      as.numeric(d.ci(d = d, n1 = n1, n2 = if(paired) NA else k * n1, conf.level = conf.level)[, 2:3])
    }
    
    m <- function(n1, n2, width){
      abs(abs(diff(f(n1 = n1, n2 = n2))) - width)
    } 
    
    n1 <- optimize(m, c(1, 1e7), width = width)
    
    if(round(n1$objective, 4) != 0) return(c(NaN, message("Warning: NaN produced. Are input values correct?")))
    
    n1 <- if(paired) ceiling(n1[[1]])
    n2 <- if(paired) NA else round(k * n1)
    
    list(d = d, n1 = n1, n2 = n2, base.rate = base.rate, width = width, conf.level = conf.level, assure = assure, paired = paired)
  }
  
 n <- nd(d = d, width = width, base.rate = base.rate, conf.level = conf.level)
 a <- d.ci(d = d, n1 = n$n1, n2 = n$n2, conf.level = c(assure, 2*assure - 1))$upper
 
 dnew <- function(dnew = dnew, n1 = n$n1, n2 = n$n2, d = d, assure = assure){
   total <- sum(pcohen(c(-dnew, dnew), d, n1 = n1, n2 = n2, lower.tail = c(TRUE, FALSE)))
   return(abs(total - (1 - assure)))
 }
 
 dnew <- optimize(dnew, a, d = d, assure = assure)[[1]]
 nd(d = dnew, conf.level = conf.level, width = width, base.rate = base.rate)
})

data.frame(t(G(d = d, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure)), row.names = NULL)
}
                 
#=====================================================================================================================================================================================================                

d.unbias <- function(d, n1, n2 = NA, t = NA){

  df <- ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)
   N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
   d <- if(!is.na(t)) t/sqrt(N) else d
   d * exp(lgamma(df/2)-log(sqrt(df/2)) - lgamma((df-1)/2))
}


#=======================================================================================================================================
                
                
graph <- function(x, y = NULL, type = "p", xlim = NULL, ylim = NULL, 
                   log = "", main = NULL, sub = NULL, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)), 
                   ann = par("ann"), axes = TRUE, frame.plot = axes, panel.first = NULL, 
                   panel.last = NULL, asp = NA, add = FALSE, show = TRUE, ...)
{
  
  if(!add){ 
    
  if(!show){type <- "n"; axes <- FALSE; ann <- FALSE} 
    
    graphics.off()  
    
    plot(x = x, y = y, type = type, xlim = xlim, ylim = ylim, 
         log = log, main = main, sub = sub, xlab = xlab, ylab = ylab, 
         ann = ann, axes = axes, frame.plot = frame.plot, panel.first = panel.first, 
         panel.last = panel.last, asp = asp, ...)
    
  }else{
    
    lines(x = x, y = y, type = type, ...)
  }
}
                
                
#========================================================================================================================================================
                
              
peta2d.fun <- function(peta = seq(.1, .5, .1), n = seq(30, 300, 10), base.rate = 1, xlab = "Group Sample Size", ylab = "Cohen's d", ylim = NA, ...)
{

n <- sort(n)
peta <- sort(peta)  

d <- lapply(peta, function(x) peta2d(x, n1 = n, n2 = base.rate*n))
ylim <- if(is.na(ylim)) range(d) else ylim
for(i in 1:length(peta)){
  graph(n, d[[i]], type = "l", add = i!= 1, ylim = ylim, xlab = xlab, ylab = ylab, ...)
  text(mean(n), mean(d[[i]]), bquote(eta[p]^2 == .(round(peta[i], 3))), col = 2, pos = 3, xpd = NA, cex = .8)
  }
}
            

#================================================================================================================================================================
            
            
#exp.pov <- function(P2, K, N)
#{ 
#  expect <- 1 - ((N - K - 1)/(N - 1)) * (1 - P2) * gsl::hyperg_2F1(1, 1, (N + 1)/2, P2)
#  max(0, expect)
#}
       
#====================================================================================================================================================================
            

plan.f.cic <- function(peta = .2, design = 2 * 2, n.level = 2, n.covar = 0, conf.level = .9, width = .2, regress = FALSE,  pair.design = 0, assure = .99){
  
  if(any(conf.level >= 1) || any(conf.level <= 0) || any(assure >= 1) || any(assure <= 0)) stop("'conf.level' and 'assure' must be between '0' and '1'.", call. = FALSE)
  
  G <- Vectorize(function(peta, conf.level, width, assure, design, n.level, n.covar, regress, pair.design){
    
    n.f <- function(peta, conf.level, width, assure, design, n.level, n.covar, regress, pair.design){
      
      alpha <- (1 - conf.level)/2
      if(regress){ n.level <- n.level + 1 ; design <- n.level }
      if(pair.design != 0) design <- 2
      df1 <- n.level - 1
      if(n.covar < 0) n.covar <- 0
      options(warn = -1)
      
      f <- function(alpha, q, df1, df2, ncp){
        alpha - suppressWarnings(pf(peta2F(peta, df1, df2), df1, df2, ncp, lower.tail = FALSE))
      }
      
      pbase <- function(df2){      
        
        b <- sapply(c(alpha, 1 - alpha), function(x) 
          tryCatch(uniroot(f, c(0, 1e7), alpha = x, q = q, df1 = df1, df2 = df2)[[1]], error = function(e) NA))
        if(any(is.na(b))) b <- c(1, 1e4)     
        ncp2peta(b, df2 + design)
      }
      
      m <- function(df2, width){
        abs(diff(pbase(df2))) - width
      }
      
      df2 <- uniroot(m, c(0, 1e3), width = width, extendInt = "yes")[[1]]
      
      df2 <- if(regress) df2 else df2 - n.covar
      
      N <- ceiling(df2 + design)
      bal <- ceiling(N/design) * design
      if(pair.design != 0){ N <- pair.design * (bal/2) ; message("\nNote: You are doing reseach planning for accurate 'pairwise' comparisons.") }
      N <- if(!regress & design != 0 & N %% design != 0) bal else N
      n.covar <- if(n.covar == 0) NA else n.covar
      n.level <- if(regress) n.level-1 else n.level
      design <- if(regress) n.level else design
      
      list(peta = peta, total.N = N, width = width, n.level = n.level, conf.level = conf.level, assure = assure, df1 = df1, df2 = df2, design = design)
    }
    
    n <- n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar, pair.design = pair.design)  
    
    peta <- exp.peta(pbase = n$peta, df1 = n$df1, df2 = n$df2, N = n$total.N)
         
    n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar, pair.design = pair.design)
    
  })
  
  data.frame(t(G(peta = peta, conf.level = conf.level, width = width, design = design, n.level = n.level, n.covar = n.covar, regress = regress, pair.design = pair.design, assure = assure)), row.names = NULL)
} 
                                       
                    
#==========================================================================================================================================================================================================================                    

                    
d.width.meta <- Vectorize(function(lower, upper, n1 = 50, n2 = 50){
  
  abs(diff(d2peta(c(lower, upper), n1 = n1, n2 = n2)))
})
  
#==========================================================================================================================================================================================================================                    
                    

plan.f.cig <- function(pov, design = 2 * 2, f = NA, n.level = 2, n.pred = NULL, n.covar = 0, conf.level = .95, width = NA, assure = .99, expect = FALSE, reduce.by = "0%", d = NA, lower, upper, increase.by = "0%", tol = 1e3){
  
  
  if(any(conf.level >= 1) || any(conf.level <= 0) || any(assure >= 1) || any(assure <= 0)) stop("'conf.level' and 'assure' must be between '0' and '1'.", call. = FALSE)
  if(expect) assure <- .5
  regress <- if(!is.null(n.pred)) TRUE else FALSE
  if(regress) n.level <- n.pred
  if(!is.na(d)) { pov <- d2peta(d = d, n1 = 300, n2 = 300) ; n.level <- 2 ;
  message("\nNote: For 'pairwise' comparisons, 'total.N' is for '2' groups.") }
  if(!is.na(d) & is.na(width)) width <- d.width.meta(lower = lower, upper = upper)
  
  peta <- pov
  
  fac <- if(increase.by != "0%" & reduce.by == "0%") { 1 + as.numeric(substr(increase.by, 1, nchar(increase.by)-1))/ 1e2 
  } else if(reduce.by != "0%" & increase.by == "0%") { 1 - (as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2) 
  } else { 1 }
  
  
  if(fac <= 0 || increase.by == "0%" & fac > 1) fac <- 1
  
  width <- width * fac
  
  
  G <- Vectorize(function(peta, conf.level, width, assure, design, n.level, n.covar, regress, expect){
    
    
    n.f <- function(peta, conf.level, width, assure, design, n.level, n.covar, regress){
      
      alpha <- (1 - conf.level)/2
      if(regress){ n.level <- n.level + 1 ; design <- n.level }
      df1 <- n.level - 1
      if(n.covar < 0) n.covar <- 0
      options(warn = -1)
      
      f <- function(alpha, q, df1, df2, ncp){
        alpha - suppressWarnings(pf(peta2F(peta, df1, df2), df1, df2, ncp, lower.tail = FALSE))
      }
      
      pbase <- function(df2){      
        
        b <- sapply(c(alpha, 1 - alpha), function(x) 
          tryCatch(uniroot(f, c(0, 1e7), alpha = x, q = q, df1 = df1, df2 = df2, extendInt = "yes")[[1]], error = function(e) NA))
        if(any(is.na(b))) b <- c(1, tol)     
        ncp2peta(b, df2 + design + n.covar)
      }
      
      m <- function(df2, width){
        abs(diff(pbase(df2))) - width
      }
      
      df2 <- uniroot(m, c(0, 1e3), width = width, extendInt = "yes")
      
      
      if(round(df2$f.root, 3) != 0) stop("\n****************************\nImpossible planning: You may change your 'width' or 'lower' & 'upper' or 'tol'.\n****************************\n", call. = FALSE)
      
      df2 <- ceiling(df2[[1]])
      
      N <- ceiling(df2 + design) + n.covar
      # bal <- ceiling(N/design) * design
      # N <- if(!regress & design != 0 & N %% design != 0) bal else N
      n.covar <- if(n.covar == 0) NA else n.covar
      n.level <- if(regress) n.level-1 else n.level
      design <- if(regress) n.level else design
      df1 <- if(regress) n.level else df1
      
      list(peta = peta, total.N = N, width = width, n.level = n.level, conf.level = conf.level, assure = assure, df1 = df1, df2 = df2, design = design)
    }
    
    n <- n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)  
    
    peta <- try(exp.peta(pbase = n$peta, df1 = n$df1, df2 = n$df2, N = n$total.N), silent = TRUE)
    
    if(inherits(peta, "try-error")) stop("\n****************************\nImpossible planning: You may change your 'width' or 'lower' & 'upper' or 'tol'.\n****************************\n", call. = FALSE)
    
    n <- n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)
    
    peta.max <- root(pov = peta, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = conf.level)$m
    
    a <- peta.ci(peta = peta, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = 2*assure - 1)
    
    nLU <- sapply(c(a$lower, a$upper), function(x) n.f(peta = x, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)$total.N)
    
    NN1 <- max(nLU, na.rm = TRUE)
    
    b <- peta.ci(peta = peta.max, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = 1 - assure)
    
    nLU <- sapply(c(b$lower, b$upper), function(x) n.f(peta = x, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)$total.N)
    
    NN2 <- max(nLU, na.rm = TRUE)
    
    NN3 <- if(!(peta.max %inn% c(a$lower, a$upper))) NN1 else max(NN1, NN2)
    
    return(c(peta = peta, total.N = NN3, width = width, n.level = n.level, design = design, conf.level = conf.level))
    
  })
  
  a <- data.frame(t(G(peta = peta, conf.level = conf.level, width = width, design = design, n.level = n.level, n.covar = n.covar, assure = assure, regress = regress, expect = expect)), regress = regress, assure = assure, row.names = NULL)
  names(a)[4] <- if(regress) "n.pred" else "n.level"
  names(a)[1] <- if(regress) "R2" else if(!is.na(d)) "d" else "peta"
  a[, 1] <- if(is.na(d)) pov else d
  a
}                                                                                                                       
                    
#==========================================================================================================================================================================================================================
                    
                    
"%inn%" <- function(x = 3.5, interval = c(3, 5)){
  
  r <- range(interval, na.rm = TRUE)
  
    x >= r[1] & x <= r[2] 
}
                                       
#=============================================================================================================================================================================================================================
                    
                    
root <- function(pov = .6, df1 = 3, df2 = 108, N = 100, conf.level = .95, show = FALSE, ...){
  
  f <- function(x){ 
    
    ci <- peta.ci(peta = x, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = 1e2)
    
    abs(ci$upper - ci$lower)
  }
  
  m <- optimize(f, 0:1, maximum = TRUE)[[1]]
  
  est <- uniroot(function(x) f(pov) - f(x), if(pov >= m) c(0, m) else c(m, 1))[[1]]
  
  if(show) curve(f, panel.f = abline(v = c(pov, est), h = f(pov), col = 2, lty = c(2, 1, 1)), ...) 
  
  list(m = m, est = est, w = f(m))
}

#=======================================================================================================================================================================================================================================
                 
                 
plan.r.ci <- function(rho = .4, width = .4, conf.level = .95, assure = .99, expect = FALSE){

  rho[rho <= -1] <- -.99999  
  rho[rho >= 1] <- .99999
  
  if(expect) assure <- .5
  
  G <- Vectorize(function(rho, width, conf.level, assure){
  
n.r <- function(rho, width, conf.level){  
  
  f <- function(n){
    as.numeric(cor.ci(r = rho, n, conf.level = conf.level)[, 2:3])
  }
  
  m <- function(n, width){
    abs(abs(diff(f(n = n))) - width)
  }
  
  n <- optimize(m, c(2, 1e7), width = width)
  
  n <- if(round(n$objective, 4) != 0) { c(NaN, message("Warning: NaN produced. Are input values correct?"))
  } else { ceiling(n[[1]]) }
  
  return(n)
}
 
n <- n.r(rho = rho, width = width, conf.level = conf.level)

a <- cor.ci(r = rho, n = n, conf.level = 2*assure - 1)
 
nLU <- sapply(c(a$lower, a$upper), function(x) n.r(rho = x, width = width, conf.level = conf.level))

NN1 <- max(nLU, na.rm = TRUE) 
  
b <- cor.ci(r = 0, n = n, conf.level = 1 - assure)

nLU <- sapply(c(b$lower, b$upper), function(x) n.r(rho = x, width = width, conf.level = conf.level))

NN2 <- max(nLU, na.rm = TRUE)

NN3 <- if(!(0 %inn% c(a$lower, a$upper))) NN1 else max(NN1, NN2)

return(c(rho = rho, n = NN3, width = width, conf.level = conf.level, assure = assure))
})

 data.frame(t(G(rho = rho, width = width, conf.level = conf.level, assure = assure)), row.names = NULL)
}
              
          
#=====================================================================================================================
              
              
plan.r.test <- function (rho = .3, sig.level = .05, power = .8, two.tailed = TRUE) {
  
  rho[rho <= -1] <- -.99999  
  rho[rho >= 1] <- .99999
  
  G <- Vectorize(function(rho, sig.level, power, two.tailed)
  {
    r <- abs(rho)
    
    f <- function(x) { 
      
      tcrit <- qt(if(two.tailed) sig.level/2 else sig.level, df = x - 2, lower.tail = FALSE)
      rc <- sqrt(tcrit^2/(tcrit^2 + x - 2))
      zr <- atanh(r) + r/(2 * (x - 1))
      zrc <- atanh(rc)
      power - if(two.tailed) sum(pnorm((c(zr, -zr) - zrc) * sqrt(x - 3))) else pnorm((zr - zrc) * sqrt(x - 3))
}      
    
return(c(rho = rho, n = ceiling(uniroot(f, c(4, 1e7), extendInt = "yes")[[1]]), sig.level = sig.level, power = power))
})

data.frame(t(G(rho = rho, sig.level = sig.level, power = power, two.tailed = two.tailed)), two.tailed = two.tailed)

}

#=================================================================================================================================
      
      
power.f <- Vectorize(function(peta = .2, n.level = 2, design = 2 * 3, sig.level = .05, n.covar = 0, power = .8, regress = FALSE, pair.design = NULL){

if(n.level <= 1) stop("Error: You must have at least '2 levels' or '2 predictors'.")
#if(!regress & missing(design)) stop("Error: 'design' must be numerically specified e.g., 'design = 2 * 4'.")
if(regress){ n.level <- n.level + 1 ; design <- n.level }
df1 <- n.level - 1
if(n.covar < 0) n.covar <- 0
x <- sapply(list(n.level, design, n.covar), round)
n.level <- x[1] ; design <- x[2] ; n.covar <- x[3]

f <- function(x){
  
  power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = (peta * (x + design) ) /(1 - peta), lower.tail = FALSE))
}

df2 <- ceiling(uniroot(f, c(1, 1e6), extendInt = "yes")[[1]]) - n.covar

N <- df2 + design + n.covar

bal <- ceiling(N/design) * design

N <- if(!is.null(pair.design)) pair.design * (bal/2) else N

return(N)

})
      
      
#====================================================================================================================================
      
      
 rcor <- function(n, mean = 0, N){ 
  tanh(rnorm(n = n, mean = atanh(mean), sd = 1/sqrt(N - 3)))
}

#=======================================================================================================================================
      
      
mrnorm <- function(n = 1, mu, sd, tol = 1e-6, random = TRUE)
{
 
  UseMethod("mrnorm") 
  
}


mrnorm.default <- function(n = 1, mu, sd, tol = 1e-6, random = TRUE) 
{
  p <- length(mu)
  if (!all(dim(sd) == c(p, p))) 
    stop("'mu' and 'sd' don't match.", call. = FALSE)
  eS <- eigen(sd, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'sd' is not positive definite.", call. = FALSE)
  X <- matrix(rnorm(p * n), n)
  if (!random) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(sd))) 
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) 
    drop(X)
  else t(X)
}
      

#=============================================================
      
      
late.penalty <- function(due, submit)
{
  UseMethod("late.penalty")
}

late.penalty.default <- function(due, submit)
  {
  due = strptime(due,  format = "%a %b %d %H:%M:%S %Y")
  sub = strptime(submit, format = "%a %b %d %H:%M:%S %Y")
  dayslate = as.numeric( difftime(sub, due, units = "days"))
  halflife = 7 # days until half credit
  expshape = 1 # shape of decay function
  round(exp( log(.5)/halflife^expshape*(dayslate)^expshape ), 2)
}
      
#=========================================================================
      
Anova <- function(eta.sq = .25, n = 5, min.score = 0, max.score = 25, coef = 1.2, sig.level = .05, ...){
  
  beta = qnorm(c(1e-16, .9999999999999999))
  q = c(min.score, max.score)
  musd = solve(cbind(1L, beta), q)  
  m1 = musd[[1]]  
  sdw = musd[[2]] 
  
  x = ((sdw^2) * eta.sq) / (1 - eta.sq)
  m2 = coef * m1
  A = m1 + m2
  
  a = function(m3) abs((m1 - (A + m3)/3)^2 + (m2 - (A + m3)/3)^2 + (m3 - (A + m3)/3)^2 - 3*x)
  
  m3 = optimize(a, c(-3*max.score, 3*max.score))[[1]]
  
  mus = c(m1, m2, m3)
  k = length(mus)
  sigb = var(mus)*(k-1)/k
  eta.p = sigb / (sigb + sdw^2)
  group = gl(k, n, labels = paste0("GROUP ", 1:k))
  y = as.vector(mapply(rnorm, n = rep(n, k), mean = mus, sd = rep(sdw, k)))
  sim = anova(aov(y ~ group))
  eta = sim[1, 2] / sum(sim[, 2])
  msb = sim[1, 3]
  mse = sim[2, 3]
  omega = (sim[1, 2] - sim[1, 1]*mse) / ((sim[1, 2] + sim[2, 2]) + mse)
  eps = (sim[1, 2] - sim[1, 1]*mse) / (sim[1, 2] + sim[2, 2])
  lab = c(paste0("subj #", rev(n)[1]), paste0(rep(".", n - 2)), paste0("subj #", 1L))
  
  par(font.lab = 2, font = 2, mgp = c(2, .2, 0), ...)
  dotchart(y, groups = group, gcol = 2:(k+1), pch = 19, color = (2:(k+1))[group],
           xlab = "Participants' Scores", labels = rep(lab, k)) 
  
  ms = unique(ave(y, group))
  sd = round(tapply(y, group, sd), 3)
  g = rev(cumsum(rev(tapply(group, group, length)) + 2) - 1)
  u = par("usr")
  
  segments(ms, c(u[4], g[2:3]), ms, c(g[2:3], u[3]), col = 2:(k+1), lty = 2)
  arrows(ms[1:2], g[2:3], ms[2:3], g[2:3], code = 3, length = .12, angle = 40, lwd = 2)
  
  ms = round(ms, 3)
  legend("topright", paste0("Mean = ", ms[1], "\n", "sd = ", sd[1]), bty = "n", text.col = "red4")
  legend("left", paste0("Mean = ", ms[2], "\n", "sd = ", sd[2]), bty = "n", text.col = "darkgreen")
  legend("bottomleft", paste0("Mean = ", ms[3], "\n", "sd = ", sd[3]), bty = "n", text.col = 4)
  legend("right",  paste0(" Eta Sq. = ", signif(eta, 3)), bty = "n")
  
  uneq = function(a, b, sig = 4) round(a, sig) != round(b, sig)
  if(uneq(eta.sq, eta.p)) message("\nNote: You may need to change 'coef' to get requested 'eta.sq'.")
  p = power.anova.test(groups = k, n = n, between.var = msb, within.var = mse, sig.level = sig.level)[[6]]
  
  cbind(sim, Pop.Eta.Sq = c(eta.p, NA), Eta.Sq = c(eta, NA), Omega.Sq = c(omega, NA), Epsilon.Sq = c(eps, NA), Power = c(p, NA))
}
      
      
#======================================================================================================================================
      
      
players <- function(Step = 16, n.Players = 16, n.Steps = 16, step.size = .5, adj = 3){
  
  graphics.off()
  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))
  
  par(mar = c(2.2, 1.8, 1.5, 1.8) );
  m = matrix( c(1, 1, 1, 1, 1, 1,   2, 2), nrow = 8, ncol = 1 );
  layout(m)
  
  if(n.Players < 2) { n.Players = 2;
  message("\n\tYou can't have less then \"2 players\" on the field.")}
  
  if(Step > n.Steps) { Step = n.Steps;
  message("\n\tYou can't take more steps than what you planned to.\n\tWe picked the maximum steps that you planned to.")}
  
  if(Step < 0 || n.Steps < 1) { Step = 0; n.Steps = 1;
  message("\n\tYou can't have less than \"0\" steps.\n\tAlso, You can't have less than \"1\" as your total steps.")}
  
  plot(-6:6, -6:6, type = "n", axes = F, ann = F)
  
  axis(1, at = c(-6:-1, 1:6), font.axis = 2, cex.axis = 1.5 )
  axis(1, at = 0, font.axis = 2, cex.axis = 1.9, col.axis = 'red' )
  
  par = par('usr')
  rect(par[1], par[3], par[2], par[4], col = 'darkseagreen1' )
  
  points( 0, 0, cex = 7, pch = 20, col = 0)
  points( 0, 0, cex = 40, lwd = 5, col = 0)
  abline(v = 0, lwd = 10, col = 0)
  
  rect(-6, -6, 6, 6, lwd = 5, border = 0)
  rect(-6.5, -2, -5.5, 2, col = 'darkseagreen1', border = 0, lwd = 5)
  rect(rep(-6.5, 2), rep(-2, 2), rep(-5.5, 2), rep(2, 2), border = 0, density = 10, angle = c(180, 90), col = 0)
  
  rect(6.5, -2, 5.5, 2, col = 'darkseagreen1', border = 0, lwd = 5)
  rect(rep(6.5, 2), rep(-2, 2), rep(5.5, 2), rep(2, 2), border = 0, density = 10, angle = c(180, 90), col = 0)
  
  box(col = 'darkseagreen1')
  
  points( c( rep(par[1]+.01, 2), rep(par[2]-.01, 2) ), rep( c(-2, 2), 2 ), cex = 3, pch = 20, col = 0 )
  
  x <- rep(0, n.Players)                       ## Initial position of players
  y <- seq(from = -6, to = 6, len = n.Players) ## y-position for players
  
  ## Sample movement of players:
  xStepsMx <- matrix(sample(c(-1, 1)*step.size, n.Players*n.Steps, replace = TRUE),
                     nrow = n.Players, ncol = n.Steps)
  
  ## Position of players:
  xPosMx <- t(sapply(1:nrow(xStepsMx), function(ii) cumsum(xStepsMx[ii,]))) + x
  
  positions = if (Step > 0){ xPosMx[,Step] } else {  x  }
  
  segments(positions, y, 0, y, lty = 2, col = 'red')
  points(positions, y, cex = 7, lwd = 3, pch = 21, bg = "white")
  text(positions, y, 1:n.Players, font = 2, cex = 1.5)
  
  if (Step == 0) {
    
    plot(1, 1, ty = 'n', axes = F, ann = F)
    
    text(1, 1, "Let players take at least \"1\" step on the field.", cex = 2.5, font = 2, col = 'red4')
    
  } else {
    par( mar = c(3.5, 4, 2, 2.1), mgp = c(2, .5, 0) )
    
    DENS = density(positions, adjust = adj, n = 1e4)

    plot( table(positions)/length(positions), col = 2, lwd = 3, xaxt = "n", xlab = "Positions", lend = 1,
          ylab = "Probability", font.axis = 2, font.lab = 2, xlim = c(-6, 6), main = NA, bty = 'n', yaxs = "i")
    
    lines(DENS, col = 4, lwd = 2)
    
    polygon(DENS, col = rgb(0, 1, 1, .1), border = NA )
    
    axis(1, at = seq(-6, 6, len = 13), xlab = 'positions', font = 2, cex.axis = 1.2 )
    
    legend("topleft", legend = c(paste("Steps Taken = ", Step), paste("Players =", n.Players) ),
           bty="n", inset = c(-.02, .01), text.font=2, text.col = 'red4', cex = 1.5)
  }
}
                     
#======================================================================================================================
                     
                     
postmode <- function(x, ...) 
{
  d <- density(x, ...)
  d$x[which.max(d$y)]
}
                     
               
#========================================================================================================================
             
                     
d.width.plot <- function(d.range = seq(.15, .92, l = 5), n1, n2 = NA, reduce.by = "20%", conf.level = .95, expect = FALSE,
                         base.rate = 1, assure = .99, xlab = NA, ylim = NULL, xlim = NULL){
  
  fac <-  if(is.character(reduce.by)) (1 - (as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2))  else 1 - reduce.by
  
  paired <- if(is.na(n2)) TRUE else FALSE
      
  ci <- d.ci(d = d.range, n1 = n1, n2 = n2, conf.level = conf.level)
  
  current <- abs(ci$upper - ci$lower)
  
  desired <- current * fac
  
  x <- 1:length(current)
  y <- current
  z <- desired
  xlim <- if(is.null(xlim)) NULL else xlim
  ylim <- if(is.null(ylim)) range(1.1*y, .9*z) else ylim
  par(xpd = NA)
  
  xlab <- if(is.na(xlab)) bquote(bold("Cohen's d common in L2")) else xlab
  
  plot(x, y, ylim = ylim, xlim = xlim, cex = 1.5, xaxt = "n", panel.f = points(x, z, col = 2, cex = 1.5), 
       panel.l = arrows(x, .98*y, x, 1.02*z, len = .1), las = 1, ylab = "CI width", font.lab = 2, xlab = xlab, 
       main = paste0((1- fac)*1e2, "% ", "reduction in CI width"))
  
  axis(1, at = x, labels = round(d.range, 2))
  legend("topleft", c("Current", "Desired"), pch = 1, col = c(1, 2), cex = .7, text.font = 2, pt.cex = 1.1, adj = c(0, .35), x.intersp = c(.8, .8), bty = "n")
  box()
  
  text(x, y, round(y, 3), pos = 3, cex = .6, font = 2)
  text(x, z, round(z, 3), pos = 1, cex = .6, font = 2, col = 2)
  par(xpd = FALSE)
  
  plan.t.ci(d = d.range, conf.level = conf.level, width = desired, base.rate = base.rate, paired = paired, assure = assure, expect = expect)
}

                     
#========================================================================================================================
                     
                     
R2.width.plot <- function(R2.range = seq(.18, .51, l = 5), n.pred, N, reduce.by = "30%", conf.level = .95, xlab = NA, ylim = NULL, xlim = NULL, assure = .99, expect = FALSE){
  
  fac <- if(is.character(reduce.by)) (1 - (as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2))  else 1 - reduce.by
  
  ci <- R2.ci(R2 = R2.range, n.pred = n.pred, N = N, conf.level = conf.level)
  
  current <-  abs(ci$upper - ci$lower)
  
  desired <- current * fac
  
  x <- 1:length(current)
  y <- current
  z <- desired
  ylim <- if(is.null(ylim)) range(1.1*y, .9*z) else ylim
  xlim <- if(is.null(xlim)) NULL else xlim
  xlab <- if(is.na(xlab)) bquote(bold("Common"~R^2~" in L2")) else xlab
  
  par(xpd = NA)
  
  plot(x, y, ylim = ylim, xlim = xlim, cex = 1.5, xaxt = "n", panel.f = points(x, z, col = 2, cex = 1.5), xlab = xlab,
       panel.l = arrows(x, .98*y, x, 1.02*z, len = .1), las = 1, ylab = "CI width", font.lab = 2,  
       main = paste0((1- fac)*1e2, "% ", "reduction in CI width"))
  
  axis(1, at = x, labels = round(R2.range, 2))
  legend("topleft", c("Current", "Desired"), pch = 1, col = c(1, 2), cex = .7, text.font = 2, pt.cex = 1.1, adj = c(0, .4), x.intersp = c(.8, .8), bty = "n")
  box()
  text(x, y, round(y, 3), pos = 3, cex = .6, font = 2)
  text(x, z, round(z, 3), pos = 1, cex = .6, font = 2, col = 2)

  par(xpd = FALSE)
  
  plan.f.ci(pov = R2.range, n.pred = n.pred, conf.level = conf.level, expect = expect, assure = assure, width = desired)
  
}                     
           
                     
#========================================================================================================================    

                     
peta.width.plot <- function(peta.range = seq(.26, .5, l = 5), n.level = 2, design = 2*2, N = 80, assure = .99, n.covar = 0,
                            reduce.by = "30%", conf.level = .95, ylim = NULL, xlim = NULL, xlab = NA, expect = FALSE){
  
  df1 <- n.level - 1
  df2 <- N - design - n.covar
  
  fac <-  if(is.character(reduce.by)) (1 - (as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2))  else 1 - reduce.by
  
  ci <- peta.ci(peta = peta.range, df1 = df1, df2 = df2, N = N, conf.level = conf.level)
  
  current <-  abs(ci$upper - ci$lower)
  
  desired <- current * fac
  
  x <- 1:length(current)
  y <- current
  z <- desired
  ylim <- if(is.null(ylim)) range(1.1*y, .9*z) else ylim
  xlim <- if(is.null(xlim)) NULL else xlim
  xlab <- if(is.na(xlab)) bquote(bold("Common"~eta[p]^2~" in L2")) else xlab
  
  par(xpd = NA)
  
  plot(x, y, ylim = ylim, xlim = xlim, cex = 1.5, xaxt = "n", panel.f = points(x, z, col = 2, cex = 1.5), xlab = xlab,
       panel.l = arrows(x, .98*y, x, 1.02*z, len = .1), las = 1, ylab = "CI width", font.lab = 2,  
       main = paste0((1- fac)*1e2, "% ", "reduction in CI width"))
  
  axis(1, at = x, labels = round(peta.range, 2))
  legend("topleft", c("Current", "Desired"), pch = 1, col = c(1, 2), cex = .7, text.font = 2, pt.cex = 1.1, adj = c(0, .4), x.intersp = c(.8, .8), bty = "n")
  box()
  text(x, y, round(y, 3), pos = 3, cex = .6, font = 2)
  text(x, z, round(z, 3), pos = 1, cex = .6, font = 2, col = 2)
  
  plan.f.ci(pov = peta.range, n.level = n.level, conf.level = conf.level, expect = expect, assure = assure, width = desired, n.covar = n.covar)
  
}                
                     
#=========================================================================================================================

cor.width.plot <- function(r.range = seq(.25, .68, l = 5), n = 20, reduce.by = "30%", conf.level = .95, assure = .99, ylim = NULL, xlim = NULL, xlab = NA, expect = FALSE){
  
  fac <- if(is.character(reduce.by)) (1 - (as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2))  else 1 - reduce.by
  
  ci <- cor.ci(r = r.range, n = n, conf.level = conf.level)
  
  current <- abs(ci$upper - ci$lower)
  
  desired <- current * fac
  
  x <- 1:length(current)
  y <- current
  z <- desired
  xlim <- if(is.null(xlim)) NULL else xlim
  ylim <- if(is.null(ylim)) range(1.1*y, .9*z) else ylim
  xlab <- if(is.na(xlab)) "common r in L2" else xlab

  par(xpd = NA)
  
  plot(x, y, ylim = ylim, xlim = xlim, cex = 1.5, xaxt = "n", panel.f = points(x, z, col = 2, cex = 1.5), 
       panel.l = arrows(x, .98*y, x, 1.02*z, len = .1), las = 1, ylab = "CI width", font.lab = 2, xlab = xlab, 
       main = paste0((1- fac)*1e2, "% ", "reduction in CI width"))
  
  axis(1, at = x, labels = round(r.range, 2))
  legend("top", c("Current", "Desired"), pch = 1, col = c(1, 2), cex = .7, text.font = 2, pt.cex = 1.1, adj = c(0, .35), x.intersp = c(.8, .8), bty = "n")
  box()
  
  text(x, y, round(y, 3), pos = 3, cex = .6, font = 2)
  text(x, z, round(z, 3), pos = 1, cex = .6, font = 2, col = 2)

  par(xpd = FALSE)
  
  plan.r.ci(rho = r.range, conf.level = conf.level, width = desired, assure = assure, expect = expect)
}                     
                     
#=========================================================================================================================

d.width <- Vectorize(function(d, t = NA, n1, n2 = NA, conf.level = .95){
  
  if(!missing(d)) { diff(as.numeric(d.ci(d = d, n1 = n1, n2 = n2, conf.level = conf.level)[2:3]))
  } else { diff(as.numeric(d.ci(t = t, n1 = n1, n2 = n2, conf.level = conf.level)[2:3])) }
})


R2.width <- Vectorize(function(R2, n.pred, N, f = NA, df1 = NA, df2 = NA, conf.level = .95){
  
  if(!missing(R2)) { diff(as.numeric(R2.ci(R2 = R2, n.pred = n.pred, df1 = df1, df2 = df2, N = N, conf.level = conf.level)[2:3]))
  } else { diff(as.numeric(R2.ci(f = f, n.pred = n.pred, df1 = df1, df2 = df2, N = N, conf.level = conf.level)[2:3])) }
})


peta.width <- Vectorize(function(peta, N, f = NA, df1, df2, conf.level = .95){
  
  if(!missing(peta)) { diff(as.numeric(peta.ci(peta = peta, df1 = df1, df2 = df2, N = N, conf.level = conf.level)[2:3]))
  } else { diff(as.numeric(peta.ci(f = f, df1 = df1, df2 = df2, N = N, conf.level = conf.level)[2:3])) }
  
})


cor.width <- Vectorize(function(r, n, conf.level = .95){
  
  diff(as.numeric(cor.ci(r = r, n = n, conf.level = conf.level)[2:3]))
})
                     
#==========================================================================================================================

                     
d2peta.fun <- function(d = seq(.1, 2, l = 5), n = seq(30, 300, 10), base.rate = 1, xlab = "Group Sample Size", ylab = bquote(eta[p]^2), ylim = NA, ...)
{
  
  n <- sort(n)
  d <- sort(d)  
  
  peta <- lapply(d, function(x) d2peta(x, n1 = n, n2 = base.rate*n))
  ylim <- if(is.na(ylim)) range(peta) else ylim
  for(i in 1:length(d)){
    graph(n, peta[[i]], type = "l", add = i!= 1, ylim = ylim, xlab = xlab, ylab = ylab, ...)
    text(mean(n), mean(peta[[i]]), bquote(d == .(round(d[i], 3))), col = 2, pos = 3, xpd = NA, cex = .8)
  }
}                     
 
#===========================================================================================================================
                 
                 
int.plot <- function(y, factor1, factor2, fun = mean, pch = 15:25, col = 1:20, 
                    type = "o", leg.bty = "o", leg.bg = NA, lty = 1:6, 
                    leg.horiz = FALSE, leg.pos = "top", pt.cex = 1.3,
                    ylab = paste("Interaction", "of ", deparse(substitute(y))), ...){

  options(warn = -1)
    a = tapply(y, list(factor1, factor2), fun)   
   mc = colMeans(a)
   mr = rowMeans(a)
   mg = mean(a)
cells = a - outer(mr, mc, "+") + mg

matplot(t(cells), type = type, xaxt = "n", col = col, lty = lty, pch = pch, las = 1, ylab = ylab, ...)

axis(1, at = 1:ncol(cells), labels = colnames(a), ...)

legend(leg.pos, legend = rownames(a), col = col, pch = pch, lty = lty, bty = leg.bty, bg = leg.bg, horiz = leg.horiz, text.font = 2, pt.cex = pt.cex)
}                 
                 
#===========================================================================================================================
                 
                 
boxdens.plot <- function(data = rnorm(1e4), 
                        adjust = 1, 
                        range = 1.5, 
                        hist.bars = "Sturges",
                        random = FALSE,
                        descriptives = FALSE,
                        outliers = FALSE,
                        box.col = 1,
                        hist.col = "yellow",
                        dens.col = "magenta",
                        histogram = TRUE,
                        dens.curve = FALSE,
                        box.plot = FALSE,
                        sampling = FALSE,
                        violin = FALSE,
                        only.curve = FALSE,
                        reflect.col = "lightblue",
                        reflect.fade = .2,
                        y.axis = TRUE,
                        x.axis = TRUE,
                        conf.interval = FALSE,
                        xlab = NA,
                        conf.level
                        
                         ){ 
  
  original_par = par(no.readonly = TRUE)
  on.exit(par(original_par))
  
  options(warn = -1)
  
  par(mar = c(5.1, 5.1, 4.1, 3.1), mgp = c(3.5, 1, 0))
  
  if( hist.col %in% colors() & hist.col != 0 & hist.col > 1 & !is.null(hist.col) || hist.col  <= length( colors() ) & hist.col != 0 & hist.col > 1 & !is.null(hist.col)) {
    
    hist.col } else if (hist.col == 0){ hist.col = "white"} else { hist.col = "yellow" }
  
  
  
  if(reflect.fade >= 1 || reflect.fade <= 0) { reflect.fade = .2 ; message("\n\tPick a number for \"reflect.fade =\" between 0 and 1.")  }
    
  
  
  if(!random) { set.seed(0) } else { set.seed(NULL) }
  
  
  
  data = na.omit(data)

  h1 = hist(data, plot = FALSE, breaks = hist.bars)
  
  
  q = range(h1$breaks)  
  den.x = quantile(data, c(0, 1) ) # range(den$x)
  
  
  xlimit = if ( diff(range(q)) > diff(range(den.x)) ) q else den.x
  
  
  den = density(data, n = 5e4, adjust = adjust, na.rm = TRUE, from = xlimit[1], to = xlimit[2])
  
  
  a = c(-max(h1$density), max(h1$density) )
  b = c(-max(den$y)     , max(den$y)      )
  

  ylimit = if ( diff(range(a)) > diff(range(b)) ) a else b
  
  ylim = c(ylimit[1], ylimit[2])
  xlim = c(xlimit[1], xlimit[2])
  
  
  at.1 = seq(xlim[1], xlim[2], length.out = 7)
  
  
  at.grey = seq(0, -ylim[2], length.out = 5 )
  at.grey.labels = seq(0, ylim[2], length.out = 5  )
  at.2 = seq(0, ylim[2], length.out = 5 )
  
  
  
  xcall <- match.call()$data
  if(class(xcall) == "call" && xcall[[1]] == "rcauchy") {
    
    message("\n\tCauchy distribution will have very long tails and does not have \"sd\" and \"mean\".")
    
    
    quant = quantile(x = data, probs = c(.25, .75) )
    
    
    f <- function(x) {
      
      y = c(quant[1], quant[2]) - qcauchy(c(.25, .75), location = x[1],  scale = x[2]) 
      
    }
    
    ## SOLVE:  
    
    AA = optim(c(1, 1), function(x) sum(f(x)^2), control = list(reltol = (.Machine$double.eps)) ) 
    
    
     Mode = unname(AA$par)[1]
    Scale = unname(AA$par)[2]
  
    
  } else if(dens.curve & b[2] < (.88*a[2]) || dens.curve & (.88*b[2]) > a[2]) {
    
  message("\n\tYou may need to use \"adjust =\" to adjust the density curve (often between .1 to 3) or change \"hist.bars =\".")
    
    }
  
  
  
  if(length(data) <= 1e2) { 
    
  message("\n\tSmall datasets (i.e., <= 100) may not provide a known shape (e.g., normal).
        Thus \"adjust =\" may not help. You may try \"hist.bars =\" for better visualization.")  }
  
  
  
  if(conf.interval){
    
  if(missing(conf.level) & class(xcall) == "call" && xcall[[1]] == "rcauchy"){ conf.level = 90 
  
  }else if(missing(conf.level)) { conf.level = 95 } 
    
    
  coverage  <- if (is.character(conf.level)) { as.numeric(substr(conf.level, 1, nchar(conf.level)-1)) / 100 
    
  } else if (is.numeric(conf.level)) { conf.level / 100 
    } else if(class(xcall) == "call" && xcall[[1]] == "rcauchy"){ .9 } else { .95 }
  
  
  Low.percentile = (1 - coverage) / 2 
  
  
  p1 = Low.percentile
  p2 = Low.percentile + coverage
  
  
  CI = quantile(x = data, probs = c(p1, p2) )
  
  }
  
 
  if(class(xcall) == "call" && xcall[[1]] == "rcauchy") {

    uncut.cauchy.q = quantile(x = data, probs = c(.25, .75))

    cuts <- quantile(x = data, probs = c(.025,.975) )
    data = data[data >= cuts[1] & data <= cuts[2] ]
    
 
    if(missing(hist.bars)) { hist.bars = 40 } 
    
    
    h1 = hist( data , plot = FALSE, breaks = hist.bars)
    
    
    q = cuts
    den.x = range(h1$breaks)

    
    xlimit = if ( diff(range(q)) > diff(range(den.x)) ) q else den.x
    
    
    den = density(data, n = 5e4, adjust = adjust, na.rm = TRUE, from = xlimit[1], to = xlimit[2])
    
    
    a = c(-max(h1$density), max(h1$density) )
    b = c(-max(den$y)     , max(den$y)      )
    
    ylimit = if ( diff(range(a)) > diff(range(b)) ) a else b


    xlim = c(xlimit[1], xlimit[2])
    ylim = c(ylimit[1], ylimit[2])
    
    at.1 = seq(xlim[1], xlim[2], length.out = 7)
    
    at.grey = seq(0, -ylim[2], length.out = 5 )
    at.grey.labels = seq(0, ylim[2], length.out = 5  )
    at.2 = seq(0, ylim[2], length.out = 5 )
    
    
    if( dens.curve & b[2] < (.88*a[2]) || dens.curve & (.88*b[2]) > a[2] ) {
      
      message("\n\tYou may need to use \"adjust =\" to adjust the density curve (often between .1 to 3) or change \"hist.bars =\".")
      
    }
    
}
  
  
  if(!violin){
    ylim = c( at.2[1], rev(at.2)[1] )
  }

  
  if(only.curve){ hist.col = "white"  ;  dens.col = "white" ; reflect.col = "white" }
  
  
  if(histogram) {
  
  par(xpd = TRUE) 
    
  h2 = h1
  h2$counts = -h1$density
  xlab <- if(is.na(xlab)) "Data" else xlab
  
  plot(h1, axes = FALSE, ylim = ylim, freq = FALSE, las = 1, font.axis = 2, 
       font.lab = 2, xlim = xlim, xaxs = "r", main = NA, cex.lab = 1.4, 
       col = adjustcolor(hist.col, .4), xlab = xlab, border = ifelse(only.curve, NA, 1))
  
  if(violin){
  lines(h2, col = adjustcolor(reflect.col, reflect.fade), border = adjustcolor(reflect.col, reflect.fade) ) }
  
  } else {

  plot(1, type = "n", axes = FALSE, ylim = ylim , las = 1, font.axis = 2, 
       font.lab = 2, xlim = xlim, xaxs = "r", main = NA, cex.lab = 1.4, 
       col = adjustcolor(hist.col, .4), ylab = "Density", xlab = xlab)
    
  }
  
 
  if(!histogram){
    
    if(violin) par(xpd = TRUE)
    
    q = quantile(x = data, probs = c(.25, .5, .75, 1), na.rm = TRUE)
    

    x25 = den$x[ den$x >= min(den$x) &  den$x <= q[1] ]
    y25 = den$y[ den$x >= min(den$x) &  den$x <= q[1] ]
    x50 = den$x[ den$x >= q[1]       &  den$x <= q[2] ]
    y50 = den$y[ den$x >= q[1]       &  den$x <= q[2] ]
    x75 = den$x[ den$x >= q[2]       &  den$x <= q[3] ]
    y75 = den$y[ den$x >= q[2]       &  den$x <= q[3] ]
   x100 = den$x[ den$x >= q[3]       &  den$x <= max(den$x) ]
   y100 = den$y[ den$x >= q[3]       &  den$x <= max(den$x) ]
    
    
    polygon( c(x25, rev(x25)), c(y25, rev(-y25)), col = adjustcolor(dens.col, .5), border = NA)
    
    
    polygon( c(x50, rev(x50)), c(y50, rev(-y50)), col = adjustcolor(dens.col, .2), border = NA)
    
    
    polygon( c(x75, rev(x75)), c(y75, rev(-y75)), col = adjustcolor(dens.col, .2), border = NA)
    
    
    polygon( c(x100, rev(x100)), c(y100, rev(-y100)), col = adjustcolor(dens.col, .1), border = NA)
    
    
    corner = par("usr")
    
    if(!violin) rect(corner[1], corner[3], max(den$x), 0, col = 0, border = NA)
    
  }
  

  if(x.axis){axis(1, at = at.1, labels = round( at.1, 2), font = 2 )}
  
  
  if(violin) {
    axis(2, at = at.grey, labels = round( at.grey.labels , 2), 
         las = 1, font.axis = 2, font.lab = 2, col = "grey", col.axis = "grey")
  }
  
 if(y.axis) {axis(2, at = at.2, labels = round( at.2 , 2), 
        las = 1, font.axis = 2, font.lab = 2)}
  
  
  if(sampling) rug(data, col = "red4")
  
  if(dens.curve || only.curve){
    
  lines(den, col = "red", lwd = 2, xpd = TRUE)
    
  if(violin) lines(den$x, -den$y, col = "grey", lwd = 2)
  
  }
  
  bxp.outlie.xs = boxplot.stats(data)$out
  
  bxp = boxplot.stats(data, coef = range)$stats
  
  low.whisk = bxp[1] ; low.hing = bxp[2] ; med = bxp[3] ; upp.hing = bxp[4] ; upp.whisk = bxp[5] ; Mean = mean(data) ; sd = sd(data)
  
  
  if(violin)left.box = 1/8*-max(den$y) else left.box = 1/14*-max(den$y)
  if(violin)right.box = 1/8*max(den$y) else right.box = 1/14*max(den$y)
  
  if(box.plot){
  
  if(!violin) par(xpd = TRUE) 
   
    
  if(!conf.interval){   
  arrows(low.whisk, 0, upp.whisk, 0, lwd = 2, angle = 90, code = 3, lend = 1, col = box.col) }
    

  rect(rep(low.hing, 2), rep(left.box, 2), rep(upp.hing, 2), rep(right.box, 2), col = c(NA, "white" ), lwd = 2, border = box.col)
  
  segments(med, left.box, med, right.box,  col = 'green3', lend = 1, lwd = 4)
  
  segments(Mean, left.box, Mean, right.box, col = "magenta", lend = 1, lty = 2)
  
  rect(low.hing, left.box, upp.hing, right.box, col = NA, lwd = 2, border = box.col)
  
  points(med, 0, pch = 19, cex = 1.5, col = 'red')
  
}
  
  if(conf.interval)  {
    
    arrows(CI[1], 0, CI[2], 0, lwd = 2, angle = 90, code = 3, lend = 1, col = "blue", length = .2)
    
    text(c(CI[1] , CI[2]), rep(0, 2), round(c(CI[1] , CI[2]), 2), cex = 1.3, font = 2, col = "green4", pos = 3)
    
  }
  
  
  if(outliers) {   
                  
    points(bxp.outlie.xs, rep(0, length(bxp.outlie.xs)), col = 'blue' )
    
    if(sampling) points(bxp.outlie.xs, rep(-ylim[2], length(bxp.outlie.xs)), col = 'blue' )
    
    }
  
  
if(descriptives) {
  
  par(xpd = TRUE)

  
  if(class(xcall) == "call" && xcall[[1]] == "rcauchy") {
    
  legend("topright", legend = bquote(bold(sd == "Indet.")), inset = c(-.01, -.08), text.font=2, cex = 1.5, bty = "n")
  
  legend("topleft", legend = bquote(bold(Scale == .(round(Scale, 2) )  )), inset = c(-.08, -.08), text.font=2, cex = 1.5, bty = "n",
           x.intersp = 2)
    
  legend("top", legend = bquote(bold(Mode == .(round(Mode), 2))), inset = c(-.05, -.08), text.font=2, cex = 1.5, bty = "n",
           x.intersp = .3)
    
    } else {
    
      
    legend("topleft", legend = bquote(bold(Mean == .(round(Mean, 2)))), inset = c(-.08, -.08), lwd = 2, lty = 2, col = "magenta", text.font=2, bty = "n",
             x.intersp = .1)
      
    legend("top", legend = bquote(bold(Median == .(round(med, 2)))), inset = c(-.05, -.08), lwd = 4, col = "green3", text.font=2, bty = "n",
             x.intersp = .3)  
      
    legend("topright", legend = bquote(bold(sd == .(round(sd, 2)))), inset = c(-.01, -.08), text.font = 2, bty = "n")
    
  }
  
  if(box.plot){
    
  arrows(c(low.whisk, upp.whisk), rep(max(at.2)/2, 2) , c(low.hing, upp.hing), rep(right.box/2, 2), code = 2, length = .15, angle = 20, lwd = 2, col = ifelse(histogram, 1, "blue2") )

    
    if(class(xcall) == "call" && xcall[[1]] == "rcauchy") {
    
  text(c(low.whisk, upp.whisk), rep(max(at.2)/2, 2), c(round(uncut.cauchy.q[1], 2), round(uncut.cauchy.q[2], 2)), pos = 3, font = 2, cex = 1.2, col = "magenta" )
    
      } else {
      
      text(c(low.whisk, upp.whisk), rep(max(at.2)/2, 2), c(round(low.hing, 2), round(upp.hing, 2)), pos = 3, font = 2, cex = 1.2, col = "magenta" )
      
           }
      
       }
  
    }

}                 
                 
#===========================================================================================================================
               
               
plan.mrm <- function(peta, n.rep, n.group, factor.type = c("between", "within", "bw"), sig.level = .05, n.covar = 0, power = .8, eps = .9,
                     rho = .5, d = NA)
{
  
  UseMethod("plan.mrm")
}


plan.mrm.default <- function(peta, n.rep, n.group, factor.type = c("between", "within", "bw"), sig.level = .05, n.covar = 0, power = .8, eps = .9,
                             rho = .5, d = NA){
  
  if(!is.na(d)) peta <- d2peta(d = d, n1 = 300, n2 = 300) 
  if(!is.na(d) & n.group == 2) message("\nNote: For 'pairwise' comparisons, 'total.N' is for '2' groups.")    
  if(!is.na(d) & n.group == 1) message("\nNote: For 'pairwise' comparisons, 'total.N' is for '1' group.") 
  
  
  options(warn = -1)
  
  G <- Vectorize(function(peta, n.rep, n.group, factor.type = c("between", "within", "bw"), sig.level, n.covar, power, eps, rho, d){
    
    m <- n.rep
    if(rho <= 0) rho <- 1e-7 else if(rho >= 1) rho <-.9999999
    if(eps < .5) eps <- .5 else if(eps > 1) eps <- 1
    if(n.group < 1) stop("You must have at least '1 group' in your design.", call. = FALSE)
    if(m < 1) stop("Incorrect # of measurements, change 'n.rep'.", call. = FALSE)
    if(factor.type != "between" & m < 2) stop("You must have at least '2 repeated measurements' in your design.", call. = FALSE)
    if(missing(n.group)) stop("'n.group' must be numerically specified.", call. = FALSE)
    peta <- if(missing(peta)) NA else peta
    if(n.covar < 0) n.covar <- 0
    g <- sapply(list(n.group, n.covar, m), round)
    n.group <- g[1] ; n.covar <- g[2] ; m <- g[3]
    
    factor.type <- match.arg(factor.type)
    
    df1 <- switch(factor.type, between = n.group - 1, within = (m - 1)*eps, bw = (n.group - 1)*(m - 1)*eps)
    
    u <- if(factor.type == "between") m / (1 + (m - 1)*rho) else m / (1 - rho)
    
    f <- if(factor.type == "between"){ function(x){
      
      power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = ((peta * ( x + n.group + n.covar) ) /(1 - peta))*u, lower.tail = FALSE))
    } 
      
    } else {
      
      function(x){ 
        power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = ((peta * ( ((x)/(m-1)) + n.group + n.covar) ) /(1 - peta))*eps*u, lower.tail = FALSE))
      }
    }
    
    df2 <- uniroot(f, c(1e-8, 1e3), extendInt = "yes")[[1]]
    
    df2 <- if(factor.type == "between") ceiling(df2) else df2
    
    N <- if(factor.type == "between") ceiling(df2 + n.group + n.covar)  else ceiling((df2 / ((m - 1)*eps)) + n.group + n.covar) 
    
    balanced.N <- if(factor.type == "between") ceiling(N/n.group) * n.group else NA
    
    a <- qpetab(sig.level, df1, df2, 0, lower.tail = FALSE)
    
    ncp <- if(factor.type == "between") (peta2f(peta)^2)*N*u else (peta2f(peta)^2)*N*u*eps
    
    est.power <- ppetab(a, df1, df2, ncp, lower.tail = FALSE)
    
    list(peta = peta, total.N = N, balanced.N = balanced.N, factor.type = factor.type, n.group = n.group, n.rep = n.rep, n.covar = n.covar, sig.level = sig.level, crit.peta = a, est.power = est.power)
  })
  
  a <- data.frame(t(G(peta = peta, n.rep = n.rep, n.group = n.group, factor.type = factor.type, sig.level = sig.level, 
                 n.covar = n.covar, power = power, eps = eps, rho = rho, d = d)))
  
  names(a)[1] <- if(!is.na(d)) "d" else "peta"
  a[, 1] <- if(is.na(d)) peta else d
  a
}
                          
#===========================================================================================================================
          
               
rrbinom <- function(n, size, p1, p2, rho = 0){
  
  UseMethod("rrbinom")
  
}
               
rrbinom.default <- function(n, size, p1, p2, rho = 0){

p <- p1
q <- p2

a <- function(rho, p, q) {
  rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
}

a.0 <- a(rho, p, q)

prob <- c(`(0,0)`= a.0, `(1,0)`= 1-q-a.0, `(0,1)`= 1-p-a.0, `(1,1)`= a.0+p+q-1)
if(min(prob) < 0) {
  print(prob)
  stop("A probability is negative.", call. = FALSE)
}

u <- sample.int(4, n * size, replace = TRUE, prob = prob)
y <- floor((u-1)/2)
x <- 1 - u %% 2
x <- colSums(matrix(x, nrow = size)) 
y <- colSums(matrix(y, nrow = size)) 

list(x = x, y = y)

}              
               
#===========================================================================================================================
               
rzpois <- function (n, p, lambda) 
{
  z <- rbinom(n, 1, p)
  (1 - z) * rpois(n, lambda)
}               
               
#===========================================================================================================================
               
            
pre.post.cont <- function(n1 = 12, n2 = 12, min.score = 0, max.score = 25, subjects = TRUE, conf.level = .95,
                          descriptives = TRUE, correlation = .7, effect.size = 1, digits = 6, ...){
  
  decimal <- function(x, k){
    if(is.character(x)){ x 
    }else{
      format(round(x, k), nsmall = k, scientific =
               ifelse(x >= 1e5 || x <= -1e5 || x <= 1e-5 & x >= -1e-5, TRUE, FALSE) )
    }
  }  
  
  if(min.score >= max.score){
    stop("\n\tYour \"min.score\" must be smaller than your \"max.score\".")  }
  
  beta = qnorm(c(1e-10, .9999999999))
  q = c(min.score, max.score)
  
  mu.sigma = solve(cbind(1L, beta), q)
  
  mean = mu.sigma[[1]]
  sd = mu.sigma[[2]]
  
  coeff = effect.size*sd
  
  aa = mean + .5*mean
  bb = mean + .3*mean
  cc = aa - coeff
  
  mean.g2 = min(bb, cc)
  mean.g1 = mean.g2 + coeff
  
  TRUE.d = (mean.g1 - mean.g2) / sd      
  
  cor.pop = correlation
  
  mu <- c(0, 0)
  cov.pop <- matrix(c(1, cor.pop, cor.pop, 1), nrow = 2)
  
  mvnorm.mat <- mrnorm(n1, sd = cov.pop, mu = mu)
  
  a <- mvnorm.mat[ , 1] * sd + mean.g1
  b <- mvnorm.mat[ , 2] * sd + mean.g2
  
  y1 = c(a - b)
  
  mvnorm.mat <- mrnorm(n2, sd = cov.pop, mu = mu)
  
  a <- mvnorm.mat[ , 1] * sd + mean.g2
  b <- mvnorm.mat[ , 2] * sd + mean.g2
  
  y2 = c(a - b)
  
  y = c(y1, y2)
  
  groups = factor(rep(1:2, c(n1, n2)), labels = c("Treatment", "Control"))
  
  mean.g1 = mean(y[groups == "Treatment"])  
  mean.g2 = mean(y[groups == "Control"])    
  
  sd.g1 = sd(y[groups == "Treatment"])
  sd.g2 = sd(y[groups == "Control"])
  
  groups.for.t = factor(rep(1:2, c(n1, n2)))
  
  test = t.test(y ~ groups.for.t, var.equal = TRUE)
  
  t.value = unname(test$statistic) ; p.value = test$p.value
  
  Cohend = t.value / sqrt((n1*n2)/(n1+n2)) 
  
  lab1 = if(n1 < 10 || n2 < 10) paste0("subj #", rev(1L:n1)) else c(paste0("subj #", rev(n1)[1]), paste0(rep(".", n1 - 2)), paste0("subj #", 1L))
  lab2 = if(n1 < 10 || n2 < 10) paste0("subj #", rev(1L:n2)) else c(paste0("subj #", rev(n2)[1]), paste0(rep(".", n2 - 2)), paste0("subj #", 1L))
  
  if(subjects) {
    graphics.off()
    par(font.lab = 2, mar = c(4.2, 1, 2, 1), mpg = c(1, .2, 0), xaxt = "n", ...)
    dotchart(y, groups = groups, color = c(4, 2)[groups], 
             font = 2, pch = 19, gcolor = c(4, 2), xlab = "Participants' Gain Scores",
             pt.cex = ifelse(n1 <= 20 || n2 <= 20, 1.5, .8), labels = c(lab1, lab2), main = NA,
             cex.main = 2) 
  } else {
    graphics.off()
    par(font.lab = 2, mar = c(4.2, 1, 2, 1), mpg = c(1, .2, 0), xaxt = "n", ...)
    dotchart(y, groups = groups, color = c(4, 2)[groups], 
             font = 2, pch = 19, gcolor = c(4, 2), xlab = "Participants' Gain Scores",
             pt.cex = ifelse(n1 <= 20 || n2 <= 20, 1.5, .8), labels = NA, main = NA)
  }
  par(xaxt = "s") ; axis(1, font = 2)
  
  gpos = rev(cumsum(rev(tapply(groups, groups, length)) + 2) - 1)
  
  u = par("usr")  
  
  segments(c(mean.g2, mean.g1), c(u[3], u[4]), c(mean.g2, mean.g1), rep(gpos[[2]], 2), lty = 2,
           col = c(2, 4))
  
  arrows(mean.g2, gpos[[2]], mean.g1, gpos[[2]], code = 3, length = .08, col = "darkgreen")
  
  mean.diff = mean.g1 - mean.g2
  
  text((mean.g1+mean.g2)/2, gpos[[2]], bquote(bold("Mean diff." == .(decimal((mean.diff), 2)))), font = 2, pos = 3, col = "green4", cex = 1.15 )
  
  legend("topright", legend = bquote(bold("Cohen's"~ bolditalic(d) == .(decimal(Cohend, 2)) )), bty = "n", text.col = "red4", cex = 1.15, bg = NA)
  
  if(descriptives) {
    
    legend("topleft", legend = bquote(bold(Mean == .(decimal(mean.g1, 2)))), text.col = 4, bty = "n", bg = NA)
    
    legend("topleft", legend = bquote(bold(sd == .(decimal(sd.g1, 2)))), text.col = 4, bty = "n", bg = NA,
           inset = .03, adj =  c(.2, 0.5) )
    
    legend("bottomleft", legend = bquote(bold(Mean == .(decimal(mean.g2, 2)))), text.col = 2, bty = "n", bg = NA, 
           inset = .03, adj = .1)
    legend("bottomleft", legend = bquote(bold(sd == .(decimal(sd.g2, 2)))), text.col = 2, bty = "n", bg = NA,
           adj =  c(-.1, 0.5))
  }
  m = matrix(c("R", "R", "O1", "O3", "T", "", "O2", "O4", "->", "->", "O2 - O1", "O4 - O3", "->", "->", "GainT", "GainC"), nrow = 2)
  dimnames(m) = list("PRE-POST-CONTROL DESIGN:" = c("", ""), c(rep("", 8)))
  m <- noquote(m)
  
  u = d.ci(t = test[[1]], n1 = n1, n2 = n2, conf.level = conf.level, digits = digits)
  test <- data.frame(test[1:3], Cohend, u[2], u[3], u[4], row.names = "result:")
  colnames(test) <- c("t.value", "df", "p.value", "cohen.d", "d.lower", "d.upper", "conf.level")
  print(list(m, test), digits = digits)
} 

#===========================================================================================================================
               

model.ci <- function(fit, level = .95){
  
  est <- coef(fit)
  se <- summary(fit)$coef[, 2]
  
  n <- length(est)
  mat <- matrix(c(rep(-1, n), rep(1, n)), nrow = n)
  
  p <- (1 - level)/2
  z <- -qnorm(p)
  ci <- est + mat * (se * z)
  rownames(ci) <- names(est)
  col1 <- paste0(format(p * 1e2, nsmall = 1), "%")
  col2 <- paste0(format((1 - p) * 1e2, nsmall = 1), "%")
  colnames(ci) <- c(col1, col2)
  ci
}              


#===========================================================================================================================
               
               
d.curve <- dcurve <- function(d = seq(0,2,.5), n1 = 30, n2 = NA, biased = TRUE, labels = TRUE){
  
  
  if(!biased) d <- exp2d(d, n1, n2)
  options(warn = -1) ; d = sort(d)
  min.d = qcohen(1e-5, min(d), n1, n2)  ;  max.d = qcohen(.99999, max(d), n1, n2)  
  
  for(i in 1:length(d)){      
    H = curve(dcohen(x, d[i], n1, n2), min.d, max.d, n = 1e3, xlab = "Effect Size (d)", 
              ylab = NA, type = "n", add = i!= 1, bty = "n", axes = FALSE, font.lab = 2, yaxs = "i")
    
    polygon(H, col = adjustcolor(i, .7), border = NA, xpd = NA)
    if(labels) text(d[i], max(H$y), bquote(bolditalic(H[.(i-1)])), pos = 3, xpd = NA)
    axis(1, at = round(d[i], 3), col = i, col.axis = i, font = 2)
    segments(d[i], 0, d[i], dcohen(d[i], d[i], n1, n2), lty = 3)
  }
}
               
#===========================================================================================================================
               
               
pov.curve <- povcurve <- function(pov = seq(0, .5, .1), df1 = 3, df2 = 73, N = 78, biased = TRUE, labels = TRUE){
  
  if(!biased) pov <- exp2peta(pov, df1, df2, N)
  
  options(warn = -1) ; p = sort(pov)
  max.p = qpeta(.999999, df1, df2, max(p), N)  
  
  for(i in 1:length(p)){      
    H = curve(dpeta(x, df1, df2, p[i], N), 0, max.p, n = 1e3, xlab = "Effect Size (POV)",
              ylab = NA, type = "n", add = i!= 1, bty = "n", axes = FALSE, font.lab = 2, yaxs = "i")
    
    polygon(H, col = adjustcolor(i, .7), border = NA, xpd = NA)
    if(labels) text(p[i], max(H$y), bquote(bolditalic(H[.(i-1)])), pos = 3, xpd = NA)
    axis(1, at = round(p[i], 3), col = i, col.axis = i, font = 2)
    segments(p[i], 0, p[i], dpeta(p[i], df1, df2, p[i], N), lty = 3)
  }
}               
               
               
#===========================================================================================================================
               
               
random <- function(N, n.groups, keep = FALSE)
{
  UseMethod("random")
}
               
               
random.default <- function(N, n.groups, keep = FALSE){
  
  r <- lapply(list(N, n.groups), round)
  N <- r[[1]]
  n.groups <- r[[2]]
  
  n <- N
  N <- ceiling(N/n.groups) * n.groups
  n.set <- N/n.groups
  
  if(keep) set.seed(9036)
  a <- sample(rep(seq_len(n.groups), n.set))
  b <- table(if(n != N) head(a, -(N - n)) else a, dnn = NULL)
  if(n != N) a[(n+1):N] <- NA
  a <- matrix(a, n.groups, n.set)
  
  a[] <- sprintf("%s%d: %d", "p", seq_len(N), a)
  rownames(a) <- rep(" ", n.groups)
  colnames(a) <- rep(" ", n.set)
  
  list(Assignments = noquote(a), Groups = b)
}

               
#===========================================================================================================================


random.block <- Vectorize(function(N, n.groups, keep = FALSE){
  
  N <- ceiling(N/n.groups) * n.groups
  n.set <- N/n.groups  
  
  if(keep) set.seed(9036)
  
  a <- replicate(n.set, sample(seq_len(n.groups)))
  colnames(a) <- paste0("block", seq_len(n.set))  
  a[] <- sprintf("%s%d: %d", "p", seq_len(N), a)
  rownames(a) <- rep(" ", n.groups)
  noquote(a)
  
}, SIMPLIFY = FALSE)              
               
               
#===========================================================================================================================
             

plot.pr <- function(fun = dbinom(0:5, 5, .1), type = "h", lwd = 2, lend = 2, xlab = "Outcomes", ylab = "Probability", xaxt = "s", add = FALSE, ...)
{
  x <- match.call()$fun
  if(class(x) == "call") x <- as.numeric(as.character(x[[2]])) else stop("Use an appropriate R function for count probability distributions.", call. = FALSE)
  graph(x[2]:x[3], fun, type = type, lwd = lwd, lend = lend, xlab = xlab, ylab = ylab, xaxt = "n", add = add, ...)
  if(xaxt != "n") axis(1, at = x[2]:x[3], labels = if(add) FALSE else TRUE, tick = if(add) FALSE else TRUE)
}              
             

#===========================================================================================================================    
    

dens.curve <- function(..., adjust = 1, na.rm = TRUE, n = 1e3, hdi = FALSE, ci = FALSE, level = .95, xlab = "x", main = NA, lwd = 2, lty = 1, col = FALSE, ylim = NA, labels = TRUE){
  
  L <- if(all(sapply(list(...), inherits, "data.frame"))) as.list(...) else list(...)
  a <- list() 
  m <- if(all(sapply(list(...), inherits, "data.frame"))) names(L) else substitute(...())
  
  y <- max(sapply(L, function(x) max(density(x, adjust = adjust, na.rm = na.rm, n = n)$y)))
    
  for(i in 1:length(L)){
    
    a[[i]] <- dens.plot(L[[i]], add = i!= 1, adjust = adjust, na.rm = na.rm, n = n, from = min(L[[i]]), to = max(L[[i]]), hdi = hdi, ci = ci, level = level, xlab = xlab, main = main, lwd = lwd, lty = lty, col = if(col) col[i] else i, ylim = if(is.na(ylim)) c(0, y) else ylim)
    
    if(labels) text(a[[i]]$mode, max(a[[i]]$y), m[[i]], pos = 3, cex = .8, font = 2, col = if(col) col[i] else i, xpd = NA)
  }
  return(invisible(a))
}                      
                       
                       
#===========================================================================================================================
                       
                       
c2fac <- function(x, breaks = NULL) {
   if(is.null(breaks)) breaks <- unique(quantile(x, 0:10/10))
   x <- cut(x, breaks, include.lowest = TRUE, right = FALSE)
   levels(x) <- paste0(breaks[-length(breaks)], ifelse(diff(breaks) > 1,
   c(paste("-", breaks[-c(1, length(breaks))] - 1, sep = ""), "+"), ""))
   return(x)
   }    

                       
#=====================================================================

dgammab <- function (x, mu, scale, log = FALSE) 
{
  dgamma(x, shape = mu/scale, scale = scale, log = log)
}

#=====================================================================

rgammab <- function (n, mu, scale) 
{
  rgamma(n, shape = mu/scale, scale = scale)
}

#=====================================================================

dgampois <- function (x, mu, scale, log = FALSE) 
{
  shape <- mu/scale
  prob <- 1/(1 + scale)
  dnbinom(x, size = shape, prob = prob, log = log)
}

#=====================================================================

rgampois <- function (n, mu, scale) 
{
  shape <- mu/scale
  prob <- 1/(1 + scale)
  rnbinom(n, size = shape, prob = prob)
}

                    
#=========================================================================
                    
decimal <- function(x, k = 3) format(round(x, k), nsmall = k) 

#==========================================================================
                    
denscurve <- function(..., adjust = 1, na.rm = TRUE, n = 1e3, hdi = FALSE, level = .95, xlab = "x", ylim = NA, xlim = NA, labels = NA, bottom = 1, top = 1, scale = 1){
  
  L <- if(all(sapply(list(...), inherits, "data.frame"))) as.list(...) else list(...)
  lab <- if(all(sapply(list(...), inherits, "data.frame"))) names(L) else substitute(...())
  
  loop <- length(L)
  soop <- seq_len(loop)
  
  a <- lapply(L, function(x) density(x, adjust = adjust, na.rm = na.rm, n = n))

  from <- numeric(loop)
  to <- numeric(loop)
  hi <- numeric(loop)
  if(hdi) CI <- matrix(NA, loop, 2)
  mode <- numeric(loop)
  
  for(i in soop){
    from[i] <- min(a[[i]]$x)
    to[i] <- max(a[[i]]$x)
    hi[i] <- max(a[[i]]$y)
    if(hdi) CI[i,] <- hdir(L[[i]], level = level)
    mode[i] <- a[[i]]$x[which.max(a[[i]]$y)]
  }
  
  f = hi + soop
  m = scale*hi + soop
  
  plot(rep(soop, 2), rep(soop, 2), type = "n", xlim = if(is.na(xlim)) c(min(from), max(to)) else xlim, ylim = if(is.na(ylim)) c(bottom*1, top*max(f)) else ylim, ylab = NA, yaxt = "n", xlab = xlab, mgp = c(2, .3, 0))
  axis(2, at = soop, labels = if(is.na(labels)) lab else labels, font = 2, las = 1, cex.axis = .8, tck = -.012, mgp = c(2, .3, 0), padj = rep(.35, loop))
  abline(h = soop, col = 8, lty = 3)
  
  for(i in soop){
  polygon(x = a[[i]]$x, y = scale*a[[i]]$y +i, col = adjustcolor(i, .4), border = NA, xpd = NA)
  }
 
if(hdi){   
  segments(CI[, 1], soop, CI[, 2], soop, lend = 1, lwd = 4, col = soop, xpd = NA)                            
  segments(mode, soop, mode, m, lty = 3, xpd = NA, lend = 1)  
  points(mode, soop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
  I = decimal(CI, 2); o = decimal(mode, 2)
  text(c(CI[,1], o, CI[,2]), soop, c(I[,1], o, I[,2]), pos = 3, font = 2, cex = .8, xpd = NA)
}  
  return(invisible(a))
}                                 

#===========================================================================================================================
       
plot.count <- function(..., freq = FALSE, type = "h", lwd = 2, lend = 2, col = NA, col.adj = 1, xlab = "Outcomes", ylab = NA, xaxt = "s", labels = NA, cex.lab = .8, yaxt = "s", xaxs = "r", yaxs = "r", x.same = FALSE, y.same = FALSE, digits = 1e2){
  
  L <- if(all(sapply(list(...), inherits, "data.frame"))) as.list(...) else list(...)
  m <- if(all(sapply(list(...), inherits, "data.frame"))) names(L) else substitute(...())
  
  y <- lapply(L, function(x) if(freq)table(x) else table(x)/length(x))
  x <- lapply(y, function(x) as.numeric(names(x)))
  
  xlim <- if(x.same) range(x, finite = TRUE) else NULL
  ylim <- if(y.same) range(y, finite = TRUE) else NULL
  
  ylab <- if(is.na(ylab) & freq) "Frequency" else if(is.na(ylab) & !freq) "Probability" else ylab
  h <- length(L)
  
  graphics.off()             
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  
  if(h > 1L) { par(mfrow = n2mfrow(h)) ; set.margin2()}
  
  for(i in 1:h){
    
    plot(x[[i]], y[[i]], type = type, lend = lend, xlab = xlab, lwd = lwd, ylab = ylab, col = if(is.na(col)) adjustcolor(i, col.adj) else adjustcolor(col[i], col.adj), ylim = ylim, xlim = xlim, xaxt = "n", yaxt = "n", xaxs = xaxs, yaxs = yaxs)
    
    mtext(if(is.na(labels)) m[[i]] else labels[i], cex = cex.lab, font = 2, col = if(is.na(col)) i else col[i], xpd = NA)
    
    if(xaxt != "n") axis(1, at = x[[i]][round(y[[i]], digits) != 0])
    if(yaxt != "n") axis(2)
  }
}

#===========================================================================================================================
                                 
mode.count <- function(x){
  
  z <- table(x)
  x <- as.numeric(names(z))
  y <- as.numeric(z)
  x[which.max(y)]
}       

#===========================================================================================================================                                 
                                 
dzbinom <- function (x, p.zero, size, prob, log = FALSE) 
{
  ll <- numeric(length(x))
  pz_i <- p.zero[1]
  size_i <- size[1]
  prob_i <- prob[1]
  for (i in 1:length(x)) {
    if (length(p.zero) > 1) 
      pz_i <- p.zero[i]
    if (length(size) > 1) 
      size_i <- size[i]
    if (length(prob) > 1) 
      prob_i <- prob[i]
    if (x[i] == 0) {
      ll[i] <- log.sum.exp(c(log(pz_i), log(1 - pz_i) + 
                               dbinom(x[i], size_i, prob_i, TRUE)))
    }
    else {
      ll[i] <- log(1 - pz_i) + dbinom(x[i], size_i, prob_i, 
                                      TRUE)
    }
  }
  if (log == FALSE) 
    ll <- exp(ll)
  return(ll)
}

#===========================================================================================================================

rzbinom <- function (n, p.zero, size, prob) 
{
  z <- rbinom(n, size = 1, prob = p.zero)
  (1 - z) * rbinom(n, size, prob)
  
}                                 
 
#===========================================================================================================================
         
         
mode.find <- function(x, y, na.rm = TRUE, finite = TRUE){
  
x <- x[if(na.rm) !is.na(x) & if(finite) !is.infinite(x)]
y <- y[if(na.rm) !is.na(y) & if(finite) !is.infinite(y)]  

  x[which.max(y)]
}

       
#===========================================================================================================================
       
       
at.set <- function(x, y, axis.tol, simplify = TRUE){
  
  y <- lapply(y, round, digits = axis.tol)
  
  x <- mapply(function(u, v) u[v != 0], x, y, SIMPLIFY = FALSE)
  
  smallest.max <- which.min(sapply(x, max))
  
  x <- c(x[smallest.max], x[-smallest.max])
  
  x[-1] <- mapply(setdiff, x[-1], x[-length(x)], SIMPLIFY = FALSE)
  
  if(simplify) as.numeric(unlist(x)) else x
}

#===========================================================================================================================

set.margin2 <- function() 
{
    par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + .1, 
        tck = -0.02)
  }

#===========================================================================================================================
              
plot.prob <- function(..., type = "h", lwd = 2, lend = 2, xlab = "Outcomes", ylab = "Probability", xaxt = "s", col = NA, col.adj = 1, labels = NA, cex.lab = .8, yaxt = "s", xaxs = "r", yaxs = "r", x.same = FALSE, y.same = FALSE, digits = 1e2){
  
  x <- match.call()[-1]
  y <- lapply(x, eval)
  if(!is.null(names(y))) y <- y[names(y) == ""]
  m <- substitute(x)
  L <- length(y)
  x <- lapply(1:L, function(i) eval(parse(text = as.character(x[[i]])[2])))
  
  xlim <- if(x.same) range(x, finite = TRUE) else NULL
  ylim <- if(y.same) range(y, finite = TRUE) else NULL
  
  graphics.off()             
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  
  if(L > 1L) { par(mfrow = n2mfrow(L)) ; set.margin2()}
  
  for(i in 1:L){
    
    plot(x[[i]], y[[i]], type = type, lwd = lwd, lend = lend, xlab = xlab, ylab = ylab, xaxt = "n", ylim = ylim, xlim = xlim, col = if(is.na(col)) adjustcolor(i, col.adj) else adjustcolor(col[i], col.adj), yaxt = "n", xaxs = xaxs, yaxs = yaxs)
    
    mtext(if(is.na(labels)) m[[i]] else labels[i], cex = cex.lab, font = 2, col = if(is.na(col)) i else col[i], xpd = NA)
    
    if(xaxt != "n") axis(1, at = x[[i]][round(y[[i]], digits) != 0])
    if(yaxt != "n") axis(2)
  }
}     
              


#===========================================================================================================================
                                                                            
                                                                            
order.list <- function(x, decreasing = FALSE, na.rm = TRUE, finite = FALSE){
  
  if(na.rm) x <- Map(Filter, list(Negate(is.na)), x)
  if(finite) x <- Map(Filter, list(Negate(is.infinite)), x)
  
  maxs <- sapply(x, max) 
  result <- list()
  i <- 1
  while(length(maxs) > 0){
    result[[i]] <- x[[which.max(maxs)]]
    x <- x[-which.max(maxs)]
    maxs <- maxs[-which.max(maxs)]
    i <- i+1
  }

if(decreasing) result else rev(result)
}
            
 
#===========================================================================================================================         
 
compare.model <- function(..., digits = 1e2){
  
  L <- list(...)
  if(length(L) < 2) stop("You need to have a least '2' fitted models for a comparison.", call. = FALSE)
  names(L) <- substitute(...())
  combs <- t(combn(x = names(L), m = 2))
  
  p.value <- round(apply(combs, 1, function(i) pchisq(2 * abs(logLik(L[[i[2]]]) - logLik(L[[i[1]]])), df = abs(L[[i[1]]]$df.residual - L[[i[2]]]$df.residual), lower.tail = FALSE)), digits)
  result <- data.frame(combs, p.value)
  names(result) <- c("model.1", "model.2", "p.value")
  
  Sig. <- symnum(result$p.value, cut = c(0, .001, .01, .05, .1, 1), na = FALSE, symbols = c("***", "**", "*", ":-(", ":-(("), corr = FALSE)
  output <- cbind(result, Sig.)
  names(output)[4] <- " "
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 ':-(' 0.1 ':-((' 1\n-------------\n") 
  return(output)
}
    
#===========================================================================================================================
                         
BICz <- function(...) {
  
  a <- AIC(..., k = sapply(list(...), function(x) log(nrow(model.frame(x)))))
  
  if(length(a) >= 2){
  names(a)[2] <- "BIC"
  a[,1] <- NULL
  rownames(a) <- as.character(substitute(...()))
  }
  a
}                         
                                                  
#===========================================================================================================================                   
                   
is.false <- function(x) identical(FALSE, x)

#===========================================================================================================================
                           
test.dist <- fitdistr <- function (x, densfun, start, ...) 
{
  myfn <- function(parm, ...) -sum(log(dens(parm, ...)))
  mylogfn <- function(parm, ...) -sum(dens(parm, ..., log = TRUE))
  mydt <- function(x, m, s, df, log) dt((x - m)/s, df, log = TRUE) - 
    log(s)
  Call <- match.call(expand.dots = TRUE)
  if (missing(start)) 
    start <- NULL
  dots <- names(list(...))
  dots <- dots[!is.element(dots, c("upper", "lower"))]
  if (missing(x) || length(x) == 0L || mode(x) != "numeric") 
    stop("'x' must be a non-empty numeric vector")
  if (any(!is.finite(x))) 
    stop("'x' contains missing or infinite values")
  if (missing(densfun) || !(is.function(densfun) || is.character(densfun))) 
    stop("'densfun' must be supplied as a function or name")
  control <- list()
  n <- length(x)
  if (is.character(densfun)) {
    distname <- tolower(densfun)
    densfun <- switch(distname, beta = dbeta, cauchy = dcauchy, 
      `chi-squared` = dchisq, exponential = dexp, f = df, 
      gamma = dgamma, geometric = dgeom, `log-normal` = dlnorm, 
      lognormal = dlnorm, logistic = dlogis, `negative binomial` = dnbinom, 
      normal = dnorm, poisson = dpois, t = mydt, weibull = dweibull, 
      NULL)
    if (is.null(densfun)) 
      stop("unsupported distribution")
    if (distname %in% c("lognormal", "log-normal")) {
      if (!is.null(start)) 
        stop(gettextf("supplying pars for the %s distribution is not supported", 
          "log-Normal"), domain = NA)
      if (any(x <= 0)) 
        stop("need positive values to fit a log-Normal")
      lx <- log(x)
      sd0 <- sqrt((n - 1)/n) * sd(lx)
      mx <- mean(lx)
      estimate <- c(mx, sd0)
      sds <- c(sd0/sqrt(n), sd0/sqrt(2 * n))
      names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      vc <- matrix(c(sds[1]^2, 0, 0, sds[2]^2), ncol = 2, 
        dimnames = list(names(sds), names(sds)))
      names(estimate) <- names(sds) <- c("meanlog", "sdlog")
      return(structure(list(estimate = estimate, sd = sds, 
        vcov = vc, n = n, loglik = sum(dlnorm(x, mx, 
          sd0, log = TRUE))), class = "fitdistr"))
    }
    if (distname == "normal") {
      if (!is.null(start)) 
        stop(gettextf("supplying pars for the %s distribution is not supported", 
          "Normal"), domain = NA)
      sd0 <- sqrt((n - 1)/n) * sd(x)
      mx <- mean(x)
      estimate <- c(mx, sd0)
      sds <- c(sd0/sqrt(n), sd0/sqrt(2 * n))
      names(estimate) <- names(sds) <- c("mean", "sd")
      vc <- matrix(c(sds[1]^2, 0, 0, sds[2]^2), ncol = 2, 
        dimnames = list(names(sds), names(sds)))
      return(structure(list(estimate = estimate, sd = sds, 
        vcov = vc, n = n, loglik = sum(dnorm(x, mx, 
          sd0, log = TRUE))), class = "fitdistr"))
    }
    if (distname == "poisson") {
      if (!is.null(start)) 
        stop(gettextf("supplying pars for the %s distribution is not supported", 
          "Poisson"), domain = NA)
      estimate <- mean(x)
      sds <- sqrt(estimate/n)
      names(estimate) <- names(sds) <- "lambda"
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("lambda", 
        "lambda"))
      return(structure(list(estimate = estimate, sd = sds, 
        vcov = vc, n = n, loglik = sum(dpois(x, estimate, 
          log = TRUE))), class = "fitdistr"))
    }
    if (distname == "exponential") {
      if (any(x < 0)) 
        stop("Exponential values must be >= 0")
      if (!is.null(start)) 
        stop(gettextf("supplying pars for the %s distribution is not supported", 
          "exponential"), domain = NA)
      estimate <- 1/mean(x)
      sds <- estimate/sqrt(n)
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("rate", 
        "rate"))
      names(estimate) <- names(sds) <- "rate"
      return(structure(list(estimate = estimate, sd = sds, 
        vcov = vc, n = n, loglik = sum(dexp(x, estimate, 
          log = TRUE))), class = "fitdistr"))
    }
    if (distname == "geometric") {
      if (!is.null(start)) 
        stop(gettextf("supplying pars for the %s distribution is not supported", 
          "geometric"), domain = NA)
      estimate <- 1/(1 + mean(x))
      sds <- estimate * sqrt((1 - estimate)/n)
      vc <- matrix(sds^2, ncol = 1, nrow = 1, dimnames = list("prob", 
        "prob"))
      names(estimate) <- names(sds) <- "prob"
      return(structure(list(estimate = estimate, sd = sds, 
        vcov = vc, n = n, loglik = sum(dgeom(x, estimate, 
          log = TRUE))), class = "fitdistr"))
    }
    if (distname == "weibull" && is.null(start)) {
      if (any(x <= 0)) 
        stop("Weibull values must be > 0")
      lx <- log(x)
      m <- mean(lx)
      v <- var(lx)
      shape <- 1.2/sqrt(v)
      scale <- exp(m + 0.572/shape)
      start <- list(shape = shape, scale = scale)
      start <- start[!is.element(names(start), dots)]
    }
    if (distname == "gamma" && is.null(start)) {
      if (any(x < 0)) 
        stop("gamma values must be >= 0")
      m <- mean(x)
      v <- var(x)
      start <- list(shape = m^2/v, rate = m/v)
      start <- start[!is.element(names(start), dots)]
      control <- list(parscale = c(1, start$rate))
    }
    if (distname == "negative binomial" && is.null(start)) {
      m <- mean(x)
      v <- var(x)
      size <- if (v > m) 
        m^2/(v - m)
      else 100
      start <- list(size = size, mu = m)
      start <- start[!is.element(names(start), dots)]
    }
    if (is.element(distname, c("cauchy", "logistic")) && 
      is.null(start)) {
      start <- list(location = median(x), scale = IQR(x)/2)
      start <- start[!is.element(names(start), dots)]
    }
    if (distname == "t" && is.null(start)) {
      start <- list(m = median(x), s = IQR(x)/2, df = 10)
      start <- start[!is.element(names(start), dots)]
    }
  }
  if (is.null(start) || !is.list(start)) 
    stop("'start' must be a named list")
  nm <- names(start)
  f <- formals(densfun)
  args <- names(f)
  m <- match(nm, args)
  if (any(is.na(m))) 
    stop("'start' specifies names which are not arguments to 'densfun'")
  formals(densfun) <- c(f[c(1, m)], f[-c(1, m)])
  dens <- function(parm, x, ...) densfun(x, parm, ...)
  if ((l <- length(nm)) > 1L) 
    body(dens) <- parse(text = paste("densfun(x,", paste("parm[", 
      1L:l, "]", collapse = ", "), ", ...)"))
  Call[[1L]] <- quote(stats::optim)
  Call$densfun <- Call$start <- NULL
  Call$x <- x
  Call$par <- start
  Call$fn <- if ("log" %in% args) 
    mylogfn
  else myfn
  Call$hessian <- TRUE
  if (length(control)) 
    Call$control <- control
  if (is.null(Call$method)) {
    if (any(c("lower", "upper") %in% names(Call))) 
      Call$method <- "L-BFGS-B"
    else if (length(start) > 1L) 
      Call$method <- "BFGS"
    else Call$method <- "Nelder-Mead"
  }
  res <- eval.parent(Call)
  if (res$convergence > 0L) 
    stop("optimization failed")
  vc <- solve(res$hessian)
  sds <- sqrt(diag(vc))
  structure(list(estimate = res$par, sd = sds, vcov = vc, 
    loglik = -res$value, n = n), class = "fitdistr")
}

#===========================================================================================================================
                           
plot.c.model <- function(object, max = max(object$model[,1]), scale = "raw", ...) {
  UseMethod("rootogram")
}

plot.c.model.default <- function(object, fitted, breaks = NULL,
                              style = c("hanging", "standing", "suspended"),
                              scale = c("sqrt", "raw"), plot = TRUE,
                              width = 0, xlab = NULL, ylab = NULL, main = NULL, ...)
{
  ## rectangle style
  scale <- match.arg(scale)
  style <- match.arg(style)
  
  ## default annotation
  if(is.null(xlab)) {
    xlab <- if(is.null(names(dimnames(object)))) {
      deparse(substitute(object))
    } else {
      names(dimnames(object))[1L]
    }
  }
  if(is.null(ylab)) {
    ylab <- if(scale == "raw") "Frequency" else "sqrt(Frequency)" 
  }
  if(is.null(main)) main <- deparse(substitute(fitted))
  
  ## breaks, midpoints, widths
  if(is.null(breaks)) {
    x <- as.numeric(names(object))
    if(length(x) < 1L) x <- 0L:(length(object) - 1L)
    breaks <- (head(x, -1L) + tail(x, -1L))/2
    breaks <- c(2 * head(x, 1L) - head(breaks, 1L), breaks,
                2 * tail(x, 1L) - tail(breaks, 1L))
    if(is.null(width)) width <- 0.9
  } else {
    x <- (head(breaks, -1L) + tail(breaks, -1L))/2
    if(is.null(width)) width <- 1
  }
  
  ## raw vs. sqrt scale
  if(scale == "sqrt") {
    obsrvd <- sqrt(as.vector(object))
    expctd <- sqrt(as.vector(fitted))
  } else {
    obsrvd <- as.vector(object)
    expctd <- as.vector(fitted)
  }
  
  ## height/position of rectangles
  y <- if(style == "hanging") expctd - obsrvd else 0
  height <- if(style == "suspended") expctd - obsrvd else obsrvd
  
  ## collect everything as data.frame
  rval <- data.frame(observed = as.vector(object), expected = as.vector(fitted),
                     x = x, y = y, width = diff(breaks) * width, height = height,
                     line = expctd)
  attr(rval, "style") <- style
  attr(rval, "scale") <- scale
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  class(rval) <- c("rootogram", "data.frame")
  
  ## also plot by default
  if(plot) plot(rval, ...)
  
  ## return invisibly
  invisible(rval)
}


plot.rootogram <- function(x,
                           xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, main = NULL,
                           border = "black", fill = "lightgray", col = "#B61A51",
                           lwd = 2, pch = 19, lty = 1, max = NULL, type = NULL, axes = TRUE, ...)
{
  ## handling of groups
  if(is.null(x$group)) x$group <- 1L
  n <- max(x$group)
  if(is.null(type)) type <- ifelse(any(table(x$group) > 20L), "l", "b")
  
  ## annotation
  if(is.null(xlab)) xlab <- TRUE
  if(is.null(ylab)) ylab <- TRUE
  if(is.null(main)) main <- TRUE
  xlab <- rep(xlab, length.out = n)
  ylab <- rep(ylab, length.out = n)
  main <- rep(main, length.out = n)
  if(is.logical(xlab)) xlab <- ifelse(xlab, attr(x, "xlab"), "")
  if(is.logical(ylab)) ylab <- ifelse(ylab, attr(x, "ylab"), "")
  if(is.logical(main)) main <- ifelse(main, attr(x, "main"), "")
  
  ## plotting function
  rootogram1 <- function(d, ...) {
    ## rect elements
    xleft <- d$x - d$width/2
    xright <- d$x + d$width/2
    ybottom <- d$y
    ytop <- d$y + d$height
    j <- unique(d$group)
    
    ## defaults
    if(is.null(xlim)) xlim <- range(c(xleft, xright))
    if(is.null(ylim)) ylim <- range(c(ybottom, ytop, d$line))
    
    ## draw rootogram
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
         xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...)
    if(axes) {
      axis(1, at = d$x)
      axis(2)
    }
    rect(xleft, ybottom, xright, ytop, border = border, col = fill)
    abline(h = 0, col = border, lty = 3)
    lines(d$x, d$line,
          col = col, pch = pch, type = type, lty = lty, lwd = lwd)
    legend("topright", c("Model", "Data"), pch = c(19, 124), lwd = c(2, NA), col = c(col, 1), bty = "n")
  }
  
  ## draw plots
  if(n > 1L) par(mfrow = n2mfrow(n))
  for(i in 1L:n) rootogram1(x[x$group == i, ], ...)
}

rootogram.numeric <- function(object, fitted, breaks = NULL,
                              start = NULL, width = NULL, xlab = NULL, ylab = NULL, main = NULL, ...)
{
  ## distribution to be fitted
  dist <- fitted
  if(!is.character(fitted)) fitted <- deparse(substitute(fitted))
  if(substr(fitted, 1L, 1L) != "d") {
    fitted <- match.arg(tolower(fitted),
                        c("beta", "cauchy", "chi-squared", "chisquared", "exponential", "f",
                          "gamma", "geometric", "log-normal", "lognormal", "logistic", "logseries", 
                          "negative binomial", "negbin", "normal", "gaussian", "poisson", "t", "weibull"))
    fitted <- switch(fitted,
                     "chisquared" = "chi-squared",
                     "lognormal" = "log-normal",
                     "negbin" = "negative binomial",
                     "gaussian" = "normal",
                     fitted)
    if(fitted == "logseries") {
      dist <- function(x, logit, ...) dlogseries(x, plogis(logit), ...)
      if(is.null(start)) start <- list(logit = 0)
    } else {
      if(is.character(dist)) dist <- fitted
    }
  }
  
  ## labels
  if(is.null(xlab)) xlab <- deparse(substitute(object))
  if(is.null(main)) main <- sprintf('fitdistr(%s, "%s")', deparse(substitute(object)), fitted)
  

  xfit <- suppressWarnings(try(fitdistr(object, dist, start = start), silent = TRUE))
  if(!inherits(xfit, "fitdistr")) stop("could not obtain fitted distribution")
  
  ## fitted probability distribution function
  pdist <- switch(fitted,
                  "beta" = pbeta,
                  "cauchy" = pcauchy, 
                  "chi-squared" = pchisq,
                  "exponential" = pexp,
                  "f" = pf, 
                  "gamma" = pgamma,
                  "geometric" = pgeom,
                  "log-normal" = plnorm, 
                  "logistic" = plogis,
                  "negative binomial" = pnbinom, 
                  "normal" = pnorm,
                  "poisson" = ppois,
                  "t" = function(x, m, s, df) pt((x - m)/s, df),
                  "weibull" = pweibull, 
                  paste("p", substr(fitted, 2L, nchar(fitted)), sep = ""))
  if(fitted == "logseries") pdist <- function(x, logit, ...) plogseries(x, plogis(logit), ...)
  if(is.character(pdist)) pdist <- try(get(pdist), silent = TRUE)
  if(!is.function(pdist)) stop("invalid specification of fitted distribution")
  
  ## second argument should be the full parameter vector
  f <- formals(pdist)
  args <- names(f)
  m <- match(names(xfit$estimate), args)
  formals(pdist) <- c(f[c(1, m)], f[-c(1, m)])
  pfun <- function(x, parm, ...) pdist(x, parm, ...)
  l <- length(xfit$estimate)
  if(l > 1L) {
    body(pfun) <- parse(text = paste("pdist(x, ", paste("parm[", 1L:l, "]", collapse = ", "), ", ...)"))
  }
  
  ## different default breaks for discrete distributions
  if(is.null(breaks)) {
    if(tolower(fitted) %in% c("geometric", "negative binomial", "poisson", "binomial", "logseries")) {
      breaks <- (if(fitted == "logseries") 0L else -1L):max(object) + 0.5
      if(is.null(width)) width <- 0.9
    } else {
      breaks <- "Sturges"
    }
  }
  
  ## observed and expected frequencies
  xhist <- hist(object, plot = FALSE, breaks = breaks)
  expctd <- xfit$n * (pfun(tail(xhist$breaks, -1L), xfit$estimate) -
                        pfun(head(xhist$breaks, -1L), xfit$estimate))
  
  ## call base rootogram function
  plot.c.model.default(xhist$counts, expctd, breaks = xhist$breaks,
                    xlab = xlab, main = main, width = width, ...)
}

rootogram.zeroinfl <- rootogram.hurdle <- function(object, newdata = NULL,
                                                   max = NULL, xlab = NULL, main = NULL, width = 0.9, ...)
{
  ## observed response
  mt <- terms(object)
  mf <- if(is.null(newdata)) {
    model.frame(object)
  } else {
    model.frame(mt, newdata, na.action = na.omit)
  }
  y <- model.response(mf)
  w <- model.weights(mf)
  if(is.null(w)) w <- rep(1, NROW(y))
  
  ## observed and expected frequencies
  max0 <- if(is.null(max)) max(1.5 * max(y[w > 0]), 20L) else max  
  obsrvd <- as.vector(xtabs(w ~ factor(y, levels = 0L:max0)))
  expctd <- if(is.null(newdata)) {
    colSums(predict(object, type = "prob", at = 0L:max0) * w)
  } else {
    colSums(predict(object, newdata = newdata, type = "prob", at = 0L:max0, na.action = na.omit) * w)
  }
  
  ## try to guess a good maximum
  if(is.null(max)) {
    max <- if(all(expctd >= 1L)) max0 else max(ceiling(mean(y)), min(which(expctd < 1L)) - 1L)
    max <- min(max, length(expctd) - 1L)
  }
  
  ## observed and expected frequencies
  obsrvd <- obsrvd[1L:(max + 1L)]
  expctd <- expctd[1L:(max + 1L)]
  
  if(is.null(xlab)) xlab <- as.character(attr(mt, "variables"))[2L]
  if(is.null(main)) main <- deparse(substitute(object))
  plot.c.model.default(obsrvd, expctd, breaks = -1L:max + 0.5,
                    xlab = xlab, main = main, width = width, ...)  
}

rootogram.zerotrunc <- function(object, newdata = NULL,
                                max = NULL, xlab = NULL, main = NULL, width = 0.9, ...)
{
  ## observed response
  mt <- terms(object)
  mf <- if(is.null(newdata)) {
    model.frame(object)
  } else {
    model.frame(mt, newdata, na.action = na.omit)
  }
  y <- model.response(mf)
  w <- model.weights(mf)
  if(is.null(w)) w <- rep(1, NROW(y))
  
  ## observed and expected frequencies
  max0 <- if(is.null(max)) max(1.5 * max(y[w > 0]), 20L) else max  
  obsrvd <- as.vector(xtabs(w ~ factor(y, levels = 1L:max0)))
  expctd <- if(is.null(newdata)) {
    colSums(predict(object, type = "prob", at = 1L:max0) * w)
  } else {
    colSums(predict(object, newdata = newdata, type = "prob", at = 1L:max0, na.action = na.omit) * w)
  }
  
  ## try to guess a good maximum
  if(is.null(max)) {
    max <- if(all(expctd >= 1L)) max0 else max(ceiling(mean(y)), min(which(expctd < 1L)))
    max <- min(max, length(expctd))
  }
  
  ## observed and expected frequencies
  obsrvd <- obsrvd[1L:max]
  expctd <- expctd[1L:max]
  
  if(is.null(xlab)) xlab <- as.character(attr(mt, "variables"))[2L]
  if(is.null(main)) main <- deparse(substitute(object))
  plot.c.model.default(obsrvd, expctd, breaks = 0L:max + 0.5,
                    xlab = xlab, main = main, width = width, ...)  
}

rootogram.glm <- function(object, newdata = NULL, breaks = NULL,
                          max = NULL, xlab = NULL, main = NULL, width = NULL, ...) 
{
  family <- substr(family(object)$family, 1L, 17L)
  if(!(family %in% c("negbin", "Negative Binomial", "poisson", "binomial", "gaussian"))) {
    stop("family currently not supported")
  }
  
  ## observed response
  mt <- terms(object)
  mf <- if(is.null(newdata)) {
    model.frame(object)
  } else {
    model.frame(mt, newdata, na.action = na.omit)
  }
  y <- model.response(mf)
  w <- model.weights(mf)
  if(is.null(w)) w <- rep(1, NROW(y))
  mu <- predict(object, newdata = newdata, type = "response", na.action = na.omit)
  
  if(family == "gaussian") {
    ## estimated standard deviation (ML)
    s <- sqrt(weighted.mean(residuals(object)^2, w))
    
    ## breaks
    if(is.null(breaks)) breaks <- "Sturges"
    breaks <- hist(y[w > 0], plot = FALSE, breaks = breaks)$breaks
    obsrvd <- as.vector(xtabs(w ~ cut(y, breaks, include.lowest = TRUE)))
    
    ## expected frequencies
    p <- matrix(NA, nrow = length(y), ncol = length(breaks) - 1L)
    for(i in 1L:ncol(p)) p[, i] <- pnorm(breaks[i + 1L], mean = mu, sd = s) -
      pnorm(breaks[i], mean = mu, sd = s)
    expctd <- colSums(p * w)
  } else if(family == "binomial") {
    ## successes and failures
    if(NCOL(y) < 2L) y <- cbind(y, 1L - y)
    
    ## number of attempts
    size <- unique(rowSums(y))
    if(length(size) > 1L) stop("rootogram only applicable to binomial distributions with same size")
    at <- 0L:size
    breaks <- -1L:size + 0.5
    
    ## observed and expected
    obsrvd <- as.vector(xtabs(w ~ factor(y[, 1L], levels = at)))
    p <- matrix(NA, length(mu), length(at))
    for(i in at) p[, i + 1L] <- dbinom(i, prob = mu, size = size)
    expctd <- colSums(p * w)
  } else {
    ## observed frequencies
    max0 <- if(is.null(max)) max(1.5 * max(y[w > 0]), 20L) else max  
    obsrvd <- as.vector(xtabs(w ~ factor(y, levels = 0L:max0)))
    
    ## expected frequencies
    at <- 0L:max0
    p <- matrix(NA, length(mu), length(at))
    if(family == "poisson") {
      for(i in at) p[, i + 1L] <- dpois(i, lambda = mu)
    } else {
      theta <- object$theta
      if(is.null(theta)) theta <- get(".Theta", environment(family(object)$variance))
      for(i in at) p[, i + 1L] <- dnbinom(i, mu = mu, size = theta)
    }
    expctd <- colSums(p * w)
    
    ## try to guess a good maximum
    if(is.null(max)) {
      max <- if(all(expctd >= 1L)) max0 else max(ceiling(mean(y)), min(which(expctd < 1L)) - 1L)
      max <- min(max, length(expctd) - 1L)
    }
    breaks <- -1L:max + 0.5
    
    ## observed and expected frequencies
    obsrvd <- obsrvd[1L:(max + 1L)]
    expctd <- expctd[1L:(max + 1L)]
  }
  
  if(is.null(xlab)) xlab <- as.character(attr(mt, "variables"))[2L]
  if(is.null(main)) main <- deparse(substitute(object))
  plot.c.model.default(obsrvd, expctd, breaks = breaks,
                    xlab = xlab, main = main,
                    width = if(family == "gaussian") 1 else width, ...)  
}



rootogram.gam <- function(object, newdata = NULL, breaks = NULL,
                          max = NULL, xlab = NULL, main = NULL, width = NULL, ...) 
{
  family <- substr(family(object)$family, 1L, 17L)
  if(!(family %in% c("Negative Binomial", "poisson", "binomial", "gaussian"))) {
    stop("family currently not supported")
  }
  
  ## observed response
  mt <- terms(object)
  mf <- if(is.null(newdata)) {
    model.frame(object)
  } else {
    model.frame(mt, newdata, na.action = na.omit)
  }
  y <- model.response(mf)
  mu <- predict(object, newdata = object$model, type = "response", na.action = na.omit)
  
  if(family == "gaussian") {
    ## estimated standard deviation (ML)
    s <- sqrt(mean(residuals(object)^2))
    
    ## breaks
    if(is.null(breaks)) breaks <- "Sturges"
    yhist <- hist(y, plot = FALSE, breaks = breaks)
    breaks <- yhist$breaks
    obsrvd <- yhist$count
    
    ## expected frequencies
    p <- matrix(NA, nrow = length(y), ncol = length(breaks) - 1L)
    for(i in 1L:ncol(p)) p[, i] <- pnorm(yhist$breaks[i + 1L], mean = mu, sd = s) -
      pnorm(yhist$breaks[i], mean = mu, sd = s)
    expctd <- colSums(p)
  } else if(family == "binomial") {
    ## successes and failures
    if(NCOL(y) < 2L) y <- cbind(y, 1L - y)
    
    ## number of attempts
    size <- unique(rowSums(y))
    if(length(size) > 1L) stop("rootogram only applicable to binomial distributions with same size")
    at <- 0L:size
    breaks <- -1L:size + 0.5
    
    ## observed and expected
    obsrvd <- table(factor(y[, 1L], levels = at))
    p <- matrix(NA, length(mu), length(at))
    for(i in at) p[, i + 1L] <- dbinom(i, prob = mu, size = size)
    expctd <- colSums(p)
  } else {
    ## observed frequencies
    max0 <- if(is.null(max)) max(1.5 * max(y), 20L) else max  
    obsrvd <- table(factor(y, levels = 0L:max0))
    
    ## expected frequencies
    at <- 0L:max0
    p <- matrix(NA, length(mu), length(at))
    if(family == "poisson") {
      for(i in at) p[, i + 1L] <- dpois(i, lambda = mu)
    } else {
      theta <- object$theta
      if(is.null(theta)) {
        theta <- if (inherits(family(object), "extended.family")) { # family = nb
          ## for nb, theta is on log scale; transform
          family(object)$getTheta(trans = TRUE)
        } else {                                                    # family = negbin
          family(object)$getTheta()
        }
      }
      for(i in at) p[, i + 1L] <- dnbinom(i, mu = mu, size = theta)
    }
    expctd <- colSums(p)
    
    ## try to guess a good maximum
    if(is.null(max)) {
      max <- if(all(expctd >= 1L)) max0 else max(ceiling(mean(y)), min(which(expctd < 1L)) - 1L)
      max <- min(max, length(expctd) - 1L)
    }
    breaks <- -1L:max + 0.5
    
    ## observed and expected frequencies
    obsrvd <- obsrvd[1L:(max + 1L)]
    expctd <- expctd[1L:(max + 1L)]
  }
  
  if(is.null(xlab)) xlab <- as.character(attr(mt, "variables"))[2L]
  if(is.null(main)) main <- deparse(substitute(object))
  plot.c.model.default(obsrvd, expctd, breaks = breaks,
                    xlab = xlab, main = main,
                    width = if(family == "gaussian") 1 else width, ...)  
}                           
 
#===========================================================================================================================
                       
zzz <- function(...){

m <- as.character(substitute(...()))   
  
z <- function(fit){
  
mu <- predict(fit, type = "response")

if(inherits(fit, "glm") && fit$family[[1]] == "poisson") exptd <- round(sum(dpois(x = 0, lambda = mu)))

if(inherits(fit, "glm") && fit$family[[1]] == "binomial") exptd <- round(sum(dbinom(x = 0, size = sum(fit$model[1][1,]), prob = mu)))

if(class(fit) %in% c("zeroinfl", "hurdle")) exptd <- round(sum(predict(fit, type = "prob")[,1]))

obs <- sum(fit$y < 1)

data.frame(obs.zeros = obs, pred.zeros = exptd)
  
 }

output <- lapply(list(...), z)

for(i in 1:length(m)) rownames(output[[i]]) <- m[i]

return(output)
}

#===========================================================================================================================
                       
zero.count <- function(...){
  
  m <- as.character(substitute(...()))   
  
  z <- function(fit){
    
    z <-  plot.c.model(fit, plot = FALSE)
    
    obs <- z$observed[z$x == 0]
    exptd <- round(z$expected[z$x == 0])
    
    data.frame(obs.zeros = obs, pred.zeros = exptd)
    
  }
  
  output <- lapply(list(...), z)
  
  for(i in 1:length(m)) rownames(output[[i]]) <- m[i]
  
  return(output)
}                         
                       
                       
#===========================================================================================================================
                       
                       
plot.c.fit <- function(..., main = TRUE, max = NULL, zero = FALSE, scale = "raw"){
  
  m <-list(...)
  L <- length(m)
  
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  
  if(L > 1L) { par(mfrow = n2mfrow(L)) ; set.margin2() }
  lab <- as.character(substitute(...()))
  
  if(!zero) {
    
    invisible(lapply(1:L, function(i) 
      plot.c.model(m[[i]], width = 0, main = if(main) lab[i] else NA, max = max, scale = scale, type = "b")) )
  }
  else {
    
    zs <- zero.count(...)
    
    for(i in 1:L){
      z <- as.numeric(zs[[i]])
      r <- range(z)
      plot(0:1, z, type = "h", lend = 1, lwd = 8, col = c(2, 4), xaxt = "n", xlab = "Zero Outcome", ylab = "Frequency", xlim = c(-1.5, 2.5), ylim = c(r[1], r[2]*1.1) )
      axis(1, at = 0:1, labels = c("Obs.0s", "Pred.0s"), font = 2, mgp = c(2, .3, 0))
      text(0:1, z, z, col = c(2, 4), pos = 3, xpd = NA)
    }
  }
}
                      

#===========================================================================================================================
                 
normalize <- function(x) 
{
  x <- x - min(x)
  x/max(x)
}                 

                   
#===========================================================================================================================
                   
                   
dens2freq <- function(dens.obj){
  
  if(class(dens.obj)[1] != "density") stop("'dens.obj' must be a 'density' object.", call. = FALSE)
  
  count <- eval(parse(text = as.character(dens.obj$call[2])))
  
  if(!(all(is.whole(count)))) { message("warning: Non-integer (i.e., count) vector detected."); count <- round(count)}
  
  xout <- as.numeric(names(table(count)))
  
  approx(dens.obj$x, dens.obj$y * length(count), xout = xout)

}                  
                   
#===========================================================================================================================
                   
                   
 get.f.ci <- Vectorize(function(pov, N, design = 2 * 2, n.level = 2, n.pred = NULL, n.covar = 0, conf.level = .95){
  
  round(uniroot(function(x){
    N - plan.f.ci(pov = pov, width = x, design = design, n.level = n.level, n.pred = n.pred, n.covar = n.covar, conf.level = conf.level)$to
  }, c(.1, .35), extendInt = "yes")[[1]], 3)
  
})  
                
#===========================================================================================================================
                   
get.t.ci <- Vectorize(function(d, n1, n2 = NA, conf.level = .95){
   
n <- min(n1, n2)  
  
k <- if(!is.na(n2)) max(n1, n2) / n else 1
  
round(uniroot(function(x){
  n - as.numeric(plan.t.ci(d = d, conf.level = conf.level, width = x, paired = if(is.na(n2)) T else F, base.rate = k)$n1)
  }, c(1e-2, 1), extendInt = "yes")[[1]], 3)
  
})      
      
#===========================================================================================================================      
      
CI.d <- function(d, n1, n2 = NA, conf.level = .95, CI){
  
  min.d <- min(qcohen(1e-4, CI[1], n1, n2), qcohen(1e-4, CI[2], n1, n2 ))
  max.d <- max(qcohen(.9999, CI[1], n1, n2), qcohen(.9999, CI[2], n1, n2))
  
  ylim <- c(0, max(dcohen(seq(min.d, max.d, l = 5e2), CI[1], n1, n2), dcohen(seq(min.d, max.d, l = 5e2), CI[2], n1, n2), na.rm = TRUE))
  
  L <- curve( dcohen(x, CI[1], n1, n2), min.d, max.d, n = 5e2, col = 4, lwd = 2, xpd = TRUE, ylab = "Density", xlab = "Cohen's d", font.lab = 2, mgp = c(1.5, .5, 0), ylim = ylim)
  U <- curve( dcohen(x, CI[2], n1, n2), n = 5e2, col = 2, add = TRUE, lwd = 2, xpd = TRUE)
  lines(CI, c(0, 0), lend = 1, lwd = 4) 
  abline(v = c(CI[1], CI[2], d), col = c(4, 2, 1), lty = 2 ) ; points(d, 0, pch = 21, bg = "cyan", col = "magenta", cex = 2)
  text(CI, c(max(L$y)/2, max(U$y)/2), round(CI, 3) , srt = 90, pos = 3, col = c(4, 2), font = 2)
  
}


#===========================================================================================================================
                   
                   
CI.peta <- function(peta, df1, df2, N, conf.level = .95, CI){
  
  min.p <- min(qpeta(1e-5, df1, df2, CI[1], N), qpeta(1e-5, df1, df2, CI[2], N))
  max.p <- max(qpeta(.99999, df1, df2, CI[1], N), qpeta(.99999, df1, df2, CI[2], N))
  
  ylim <- c(0, max(dpeta(seq(0, 1, l = 5e2), df1, df2, CI[1], N), dpeta(seq(0, 1, l = 5e2), df1, df2, CI[2], N), na.rm = TRUE))
  
  L <- curve( dpeta(x, df1, df2, CI[1], N), min.p, max.p, n = 5e2, col = 4, lwd = 2, xpd = TRUE, ylab = "Density", xlab = bquote(eta[p]^2), font.lab = 2, mgp = c(1.8, .5, 0), ylim = ylim)
  U <- curve( dpeta(x, df1, df2, CI[2], N), n = 5e2, col = 2, add = TRUE, lwd = 2, xpd = TRUE)
  lines(CI, c(0, 0), lend = 1, lwd = 4) 
  abline(v = c(CI[1], CI[2], peta), col = c(4, 2, 1), lty = 2 ); points(peta, 0, pch = 21, bg = "cyan", col = "magenta", cex = 2)
  text(CI, c(max(L$y)/2, max(U$y)/2), round(CI, 3) , srt = 90, pos = 3, col = c(4, 2), font = 2)
  
}


#===========================================================================================================================
                   
                   
CI.R2 <- function(R2, df1, df2, N, conf.level = .95, CI){
  
  min.r <- min(qpeta(1e-5, df1, df2, CI[1], N), qpeta(1e-5, df1, df2, CI[2], N))
  max.r <- max(qpeta(.99999, df1, df2, CI[1], N), qpeta(.99999, df1, df2, CI[2], N))
  
  ylim <- c(0, max(dpeta(seq(0, 1, l = 5e2), df1, df2, CI[1], N), dpeta(seq(0, 1, l = 5e2), df1, df2, CI[2], N), na.rm = TRUE))
  
  L <- curve( dpeta(x, df1, df2, CI[1], N), min.r, max.r, n = 5e2, col = 4, lwd = 2, xpd = TRUE, ylab = "Density", xlab = bquote(R^2), font.lab = 2, mgp = c(1.75, .5, 0), ylim = ylim)
  U <- curve( dpeta(x, df1, df2, CI[2], N), n = 5e2, col = 2, add = TRUE, lwd = 2, xpd = TRUE)
  lines(CI, c(0, 0), lend = 1, lwd = 4)
  abline(v = c(CI[1], CI[2], R2), col = c(4, 2, 1), lty = 2 ) ; points(R2, 0, pch = 21, bg = "cyan", col = "magenta", cex = 2)
  text(CI, c(max(L$y)/2, max(U$y)/2), round(CI, 3) , srt = 90, pos = 3, col = c(4, 2), font = 2)
}                       
                   

#===========================================================================================================================
                                 
                   
R2.ci <- function(R2, n.pred, N, f = NA, df1 = NA, df2 = NA, conf.level = .95, digits = 20, show = FALSE){ 
  
  if(is.na(df1)) df1 <- n.pred 
  if(missing(n.pred) & df1) n.pred <- df1
  if(is.na(df2)) df2 <- N - n.pred - 1
  if(missing(N)) N <- df1 + df2 + 1  
  if(missing(df2) & N) df2 <- N - df1 - 1
  
  a <- if(!missing(R2)){ peta.ci(peta = R2, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = digits)
  } else { peta.ci(f = f, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = digits) }
  
  names(a)[1] <- "R2"
  
  if(show){
    
    r <- nrow(a)
    graphics.off()
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    if(r > 1) { par(mfrow = n2mfrow(r)) ; set.margin2() }
    
    I <- eq(a$R2, df1, df2, N, conf.level)
    
    R2 <- I[[1]] ; df1 <- I[[2]] ; df2 <- I[[3]]; N <- I[[4]] ; conf.level <- I[[5]]
    
    for(i in 1:r) CI.R2(R2 = R2[i], df1 = df1[i], df2 = df2[i], N = N[i], conf.level = conf.level[i], CI = c(a$lower[i], a$upper[i]))
    
  }
  
  return(a)
}

#===========================================================================================================================
                   

peta.ci <- function(peta, f = NA, df1, df2, N, conf.level = .95, digits = 1e2, show = FALSE)
{
  UseMethod("peta.ci")
} 

peta.ci.default <- function(peta, f = NA, df1, df2, N, conf.level = .95, digits = 1e2, show = FALSE){
  
  ci <- Vectorize(function(peta, f, N, df1, df2, conf.level){
    
    q <- ifelse(is.na(f), peta2F(peta, df1, df2), f) 
    
    alpha <- (1 - conf.level)/2
    
    u <- function (ncp, alpha, q, df1, df2) {
      suppressWarnings(pf(q = q, df1 = df1, df2 = df2, ncp, lower.tail = FALSE)) - alpha
    }
    
    g <- try(uniroot(u, c(0, 1e7), alpha = alpha, q = q, df1 = df1, df2 = df2, extendInt = "yes")[[1]], silent = TRUE)
    if(inherits(g, "try-error")) g <- 0
    h <- try(uniroot(u, c(0, 1e7), alpha = 1-alpha, q = q, df1 = df1, df2 = df2, extendInt = "yes")[[1]], silent = TRUE)
    if(inherits(h, "try-error")) h <- 0
    I <- c(g, h)
    
    I <- I / (I + N)
    
    P.eta.sq <- if(is.na(f)) peta else F2peta(f, df1, df2)
    
    return(c(P.eta.sq = P.eta.sq, lower = I[1], upper = I[2], conf.level = conf.level, ncp = peta2ncp(P.eta.sq, N), F.value = q))
  })
  
  peta <- if(missing(peta)) NA else peta
  
 a <- round(data.frame(t(ci(peta = peta, f = f, N = N, df1 = df1, df2 = df2, conf.level = conf.level))), digits = digits)
 
 if(show){
   
   r <- nrow(a)
   graphics.off()
   original.par = par(no.readonly = TRUE)
   on.exit(par(original.par))
   if(r > 1) { par(mfrow = n2mfrow(r)) ; set.margin2() }
   
   I <- eq(a$P.eta.sq, df1, df2, N, conf.level)
   
   peta <- I[[1]] ; df1 <- I[[2]] ; df2 <- I[[3]]; N <- I[[4]] ; conf.level <- I[[5]]
   
   for(i in 1:r) CI.peta(peta = peta[i], df1 = df1[i], df2 = df2[i], N = N[i], conf.level = conf.level[i], CI = c(a$lower[i], a$upper[i]))
   
 }
 
 return(a)
 
}

#===========================================================================================================================
                   
                   
d.ci <- function(d, t = NA, n1, n2 = NA, conf.level = .95, digits = 1e2, show = FALSE)
{
  UseMethod("d.ci")
}

d.ci.default <- function(d, t = NA, n1, n2 = NA, conf.level = .95, digits = 1e2, show = FALSE){
  
  ci <- Vectorize(function(d, t, n1, n2, conf.level){
    
    options(warn = -1)  
    alpha = (1 - conf.level)/2
    N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    d.SE = 1/sqrt(N)
    q = ifelse(is.na(t), d/d.SE, t)
    
    f <- function(ncp, alpha, q, df){
      alpha - suppressWarnings(pt(q, df, ncp, lower.tail = FALSE))
    }
    
    CI <- sapply(c(alpha, 1-alpha),
                 function(x) uniroot(f, interval = c(-1e7, 1e7), alpha = x, q = q, df = df, extendInt = "yes")[[1]]*d.SE)
    
    Cohen.d = ifelse(is.na(t), d, t*d.SE)
    
    return(c(Cohen.d = Cohen.d, lower = CI[1], upper = CI[2], conf.level = conf.level, ncp = q))
  })
  
  d <- if(missing(d)) NA else d
  
  a <- round(data.frame(t(ci(d = d, t = t, n1 = n1, n2 = n2, conf.level = conf.level))), digits = digits)
  
  if(show){
    
    r <- nrow(a)
    graphics.off()
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    if(r > 1) { par(mfrow = n2mfrow(r)) ; set.margin2() }
    
    I <- eq(a$Cohen.d, n1, n2, conf.level)
    
    d <- I[[1]] ; n1 <- I[[2]] ; n2 <- I[[3]]; conf.level <- I[[4]]
    
    for(i in 1:r) CI.d(d = d[i], n1 = n1[i], n2 = n2[i], conf.level = conf.level[i], CI = c(a$lower[i], a$upper[i]))
    
  }
  
  return(a)
  
}                               
                   

#===========================================================================================================================
                  
                  
dhnorm <- function(x, scale = 1, log = FALSE) 
{
  if(scale < 0) stop("'scale' must be larger than '0'.", call. = FALSE)
  result <- rep(-Inf, length(x))
  result[x >= 0] <- log(2) + dnorm(x[x >= 0], mean = 0, sd = scale, 
                                   log = TRUE)
  if(!log) 
    result <- exp(result)
  return(result)
}                  
                  
#===========================================================================================================================
                  
inv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X)) 
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
       t(Xsvd$u[, Positive, drop = FALSE]))
}

#===========================================================================================================================
                  
ave.drr <- function(d, r, weight){
  
p <- eq(d, weight)
d <- p[[1]] ; weight <- p[[2]]
as.numeric(crossprod(d, weight)) / sqrt(r*sum(weight)^2 + (1-r)*sum(weight^2))

}                  
                  
#===========================================================================================================================

                  
cov.d <- function(d, n1, n2, r, no.names = FALSE){
  
D <- diag(meta.se(d, n1, n2))

m <- D%*%r%*%D
if(!no.names) rownames(m) <- colnames(m) <- paste0("d", 1:length(d))

return(m)
}


#===========================================================================================================================


autoreg.cov.d <- function(d, n1, n2, r = .5, no.names = FALSE){
 
x <- diag(length(d)) 
R <- r^abs(row(x)-col(x))

cov.d(d, n1, n2, R, no.names)
}

                  
#===========================================================================================================================
                  
cor.mat <- function(r, dim) { 
  
m <- diag(dim) 
m[lower.tri(m)] <- r
m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
m
}
                  
#===========================================================================================================================
                                        
need <- c("rstanarm", "pscl", "glmmTMB")  #, "arrangements")
have <- need %in% rownames(installed.packages())
if(any(!have)){ install.packages( need[!have] ) }
 
options(warn = -1)
suppressMessages({ 
    library("rstanarm")
    library("pscl")
    library("glmmTMB")
  # library("arrangements")
})
                     
