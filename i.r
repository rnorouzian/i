
Break = "\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "   \"bayesL2\", a suite of R functions for Bayesian estimation.
    Copyright (C) 2018  Reza Norouzian, rnorouzian@gmail.com

    This program is free software: you can redistribute it under the 
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

cite = "To cite the package use:\n\nNorouzian, R., de Miranda, M. A., & Plonsky, L. (under review). The Bayesian \nrevolution in L2 research: An applied approach.

\nNorouzian, R., de Miranda, M. A., & Plonsky, L. (under review). A Bayesian \napproach to measuring evidence in L2 research: An empirical investigation."

cat(Break, cite, Break)

#==================================================================================================================

HDI <- function(FUN, lower = 0, upper = 1, level = .95, eps = 1e-3)
{
  UseMethod("HDI")
}

HDI.default <- function(FUN, lower = 0, upper = 1, level = .95, eps = 1e-3){
  
  if(!is.function(FUN)) stop("Error: 'FUN' must be a function.")
  if(length(formals(FUN)) > 1) stop("Error: 'FUN' must be a 'single-argument' function.")
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  x <- formals(FUN)
  fun <- function(x) FUN(x)
  
  posterior <- function(x) fun(x)/integrate(fun, lower, upper)[[1]]
  mode <- optimize(posterior, c(lower, upper), maximum = TRUE)[[1]]
  inverse.posterior <- function(x, side = "left") {
    target <- function(y) posterior(y) - x
    ur <- switch(side,
                 left = try(uniroot(target, interval = c(lower, mode))),
                 right = try(uniroot(target, interval = c(mode, upper))))
    if(inherits(ur, "try-error")) stop("Error: You may change prior hyperparameters or limits.")
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
  if(length(sample) < 1e3) message("Warning: \n\tInsufficient sample to produce reliable 'interval' estimates.")  
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

prop.ci <- function(k, n, conf.level = .95)
{
  UseMethod("prop.ci")
}

prop.ci.default <- Vectorize(function(k, n, conf.level = .95){
  
I = round(as.numeric(binom.test(k, n, conf.level = conf.level)[[4]]), 6)
data.frame(Prop = k/n, lower = I[1], upper = I[2], conf.level = conf.level, row.names = "Prop CI:")
})

#==================================================================================================

d.ci <- function(d, t = NA, n1, n2 = NA, conf.level = .95)
{
  UseMethod("d.ci")
}

d.ci.default <- Vectorize(function(d, t = NA, n1, n2 = NA, conf.level = .95){
  
  options(warn = -1)  
  alpha = (1 - conf.level)/2
  N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  d.SE = 1/sqrt(N)
  q = ifelse(is.na(t), d/d.SE, t)
  
  f <- function(ncp, alpha, q, df){
    abs(suppressWarnings(pt(q = q, df = df, ncp, lower.tail = FALSE)) - alpha)
  }
  
  a = if(is.na(t)){ lapply(14:ifelse(d!= 0, q+2e2, 30), function(x) c(-x, x))
    }else{ lapply(14:ifelse(t!= 0, q+2e2, 30), function(x) c(-x, x)) }
  
  CI = matrix(NA, length(a), 2)
  
  for(i in 1:length(a)){
    CI[i,] = sapply(c(alpha, 1-alpha),
    function(x) optimize(f, interval = a[[i]], alpha = x, q = q, df = df)[[1]]*d.SE)
  }  
  
  I = CI[which.max(ave(1:nrow(CI), do.call(paste, round(data.frame(CI), 3)), FUN = seq_along)), ]  
  
  Cohen.d = ifelse(is.na(t), d, t*d.SE)
  
 round(data.frame(Cohen.d = Cohen.d, lower = I[1], upper = I[2], conf.level = conf.level, ncp = q, row.names = "Cohen's d CI:"), 6)
})

#=================================================================================================================================

peta.ci <- function(peta, F.value = NA, N, df1, df2, conf.level = .9)
{
  UseMethod("peta.ci")
}
                
peta.ci.default <- Vectorize(function(peta, F.value = NA, N, df1, df2, conf.level = .9){

options(warn = -1) 
  
    q = ifelse(is.na(F.value), (-peta * df2) / ((peta * df1) - df1), F.value) 
alpha = (1 - conf.level)/2

f <- function (ncp, alpha, q, df1, df2) {
abs(suppressWarnings(pf(q = q, df1 = df1, df2 = df2, ncp, lower.tail = FALSE)) - alpha)
}

a = if(is.na(F.value)){ lapply(20:ifelse(peta!= 0, q+3e2, 30), function(x) c(-x, x))
  }else{ lapply(20:ifelse(F.value!= 0, q+3e2, 30), function(x) c(-x, x)) }

CI = matrix(NA, length(a), 2)

for(i in 1:length(a)){
CI[i,] <- sapply(c(alpha, 1-alpha), 
function(x) optimize(f, interval = a[[i]], alpha = x, q = q, df1 = df1, df2 = df2)[[1]])
}

I = CI[which.max(ave(1:nrow(CI), do.call(paste, round(data.frame(CI), 3)), FUN = seq_along)), ] 

I <- I[1:2] / (I[1:2] + N)

P.eta.sq <- if(is.na(F.value)) peta else (F.value * df1)/ ((F.value * df1) + df2)

round(data.frame(P.eta.sq = P.eta.sq, lower = I[1], upper = I[2], conf.level = conf.level, ncp = (P.eta.sq * N) / (1 - P.eta.sq), F.value = q, row.names = "P.eta.sq CI:"), 6)
})               

#=================================================================================================================================                
                
cor.ci <- function(r, n, conf.level = .95)
{
  UseMethod("cor.ci")
}
                              
cor.ci.default <- Vectorize(function(r, n, conf.level = .95){
  p = (1 - conf.level) / 2 
  I = round(tanh(atanh(r) + qnorm(c(p, 1-p))*1/sqrt(n - 3)), 6)
  data.frame(lower = I[1], upper = I[2], conf.level = conf.level, row.names = "Correlation CI:")
})                              

#==================================================================================================================

beta.id <- function(Low, High, Cover = NA)
{
  UseMethod("beta.id")
}

beta.id.default <- Vectorize(function(Low, High, Cover = NA){
  
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
    
    q <- qbeta(p = c(p1, p2), parm[[1]], parm[[2]])
    
    is.df <- function(a, b, sig = 3) round(a, sig) != round(b, sig)
    
    if(is.df(L, q[1]) || is.df(U, q[2])){
      
      stop("Error: \n\tUnable to find such a prior, make sure you have selected the correct values.")
      
    }else{
      
      return(c(alpha = round(parm[[1]], 6), beta = round(parm[[2]], 6)))    
    }
  } 
})

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

norm.id <- function(Low, High, Cover = NA)
{
  UseMethod("norm.id")
}

norm.id.default <- Vectorize(function(Low, High, Cover = NA){
  
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
    
    return(c(mean = round(parm[[1]], 6), sd = round(parm[[2]], 6)))
    
  }
})

      
#===============================================================================================

      
prop.bayes <- function(a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", yes = 55, n = 1e2, scale = .1, top = 1.5, show.prior = FALSE, bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .1, p.h0 = .5)
{
  UseMethod("prop.bayes")
}

prop.bayes.default <- function(a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", yes = 55, n = 1e2, scale = .1, top = 1.5, show.prior = FALSE, bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .1, p.h0 = .5){
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
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
      CI[i,] = HDI(posterior)
      peak[i] = posterior(mode[i])
      eq.prob[i] = integrate(posterior, lo[i], eq.hi)[[1]] - integrate(posterior, lo[i], eq.lo)[[1]]
      BF10[i] = k[i]/dbinom(yes[i], n[i], p.h0)
      estimate[i] <- yes[i]/n[i]     
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
  
    return(round(data.frame(estimate = estimate, mode = mode, lower = CI[,1], upper = CI[,2], eq.prob = eq.prob, BF10 = BF10, row.names = paste0("Prop ", 1:loop, " posterior: ")), 6))    

}else{
    p = function(x) get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1])
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, yaxs = "i", main = bquote(Proportion*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
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


prop.diff <- function(yes, n, a = 1.2, b = a, how = c("two.one", "one.two"), level = .95, top = 1, bottom = 1, scale = .1, margin = 6, legend = "topleft", eq.level = "2.5%")
{
  UseMethod("prop.diff")
}

prop.diff.default <- function(yes, n, a = 1.2, b = a, how = c("two.one", "one.two"), level = .95, top = 1, bottom = 1, scale = .1, margin = 6, legend = "topleft", eq.level = "2.5%"){
  
  n <- round(n)
  yes <- round(yes)  
  loop <- length(n)
  
  is.s <- function(...)lengths(list(...)) < 2 
  
  if(any(yes > n)) stop("Error: 'yes' cannot be larger than 'n'.") 
  if(any(is.s(n, yes))) stop("Error: 'yes' & 'n' must each have a length of '2' or larger.")
  
  eq.b <- if(is.character(eq.level)) as.numeric(substr(eq.level, 1, nchar(eq.level)-1)) / 1e2 else eq.level
  legn <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
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
  
  plot(CI, rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*loop), ylab = NA, xaxt = "n", yaxt = "n", xlab = "Credible Interval (Proportion Differences)", font.lab = 2, mgp = c(2, .3, 0))
  axis(1, at = axTicks(1), labels = paste0(round(axTicks(1), 2)*1e2, "%"), mgp = c(2, .3, 0))
  abline(h = 1:loop, col = 8, lty = 3)
  axis(2, at = 1:loop, labels = np, font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
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
  
  return(round(data.frame(estimate = estimate, mean = mean, mode = mode, median = median, sd = sd, lower = CI[,1], upper = CI[,2], eq.prob = BB, BF01 = BF01, BF10 = BF10, row.names = paste0(np, ":")), 6))                                                
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

       
d.bayes <- function(t, n1, n2 = NA, m = 0, s = 1, lo = -Inf, hi = Inf, dist.name = "dnorm", scale = .1, margin = 7, top = .8, show.prior = FALSE, LL = -3, UL = 3, bottom = 1, prior.left = -6, prior.right = 6, legend = "topleft", eq.level = .1, d.h0 = 0)
{
  UseMethod("d.bayes")
}
       
 d.bayes.default <- function(t, n1, n2 = NA, m = 0, s = 1, lo = -Inf, hi = Inf, dist.name = "dnorm", scale = .1, margin = 7, top = .8, show.prior = FALSE, LL = -3, UL = 3, bottom = 1, prior.left = -6, prior.right = 6, legend = "topleft", eq.level = .1, d.h0 = 0){
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi, t, n1, n2)
  m = I[[1]] ; s = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]] ; t = I[[6]] ; n1 = I[[7]] ; n2 = I[[8]]
  
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
  k = numeric(loop)
  eq.prob = numeric(loop)
  BF10 = numeric(loop)
  estimate = numeric(loop)
  
 if(!pr){    
    
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
      CI[i,] = HDI(posterior, LL, UL)
      h[[i]] = list(x = x <- seq(from[i], to[i], length.out = 5e2), y = posterior(x))
   BF10[i] =  k[i] / dt(t[i], df[i], d.h0*sqrt(N[i]))
   eq.prob[i] = integrate(posterior, lo[i], eq.level)[[1]] - integrate(posterior, lo[i], -eq.level)[[1]]
   estimate[i] <- t[i]/sqrt(N[i])
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
  
    return(round(data.frame(estimate = estimate, mode = mode, lower = CI[,1], upper = CI[,2], eq.prob = eq.prob, BF10 = BF10, row.names = paste0("Cohen's d ", 1:loop, " posterior: ")), 6))
        
}else{
    p = function(x) { get(d[1])(x, m[1], s[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, prior.left, prior.right, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo[1] > -Inf || hi[1] < Inf) "truncated-")*.(substring(d[1], 2))(.(round(m[1], 2)), .(round(s[1], 2)))), mgp = c(2, .5, 0))
  }
}       
   
   
#====================================================================================================================
       
       
d.priors <- function(t, ...)
{
  UseMethod("d.priors")
}

d.priors.default <- function(t, n1, n2 = NA, m = 0, s = 1, lo = -Inf, hi = Inf, dist.name = "dnorm", scale = .1, margin = 7, top = .8, show.prior = FALSE, LL = -3, UL = 3, bottom = 1, prior.left = -6, prior.right = 6, legend = "topleft"){
  
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


peta.bayes <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", scale = .1, top = 1.5, show.prior = FALSE, bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .05, peta.h0 = 0)
{
  UseMethod("peta.bayes")
}

peta.bayes.default <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", scale = .1, top = 1.5, show.prior = FALSE, bottom = 1, legend = "topleft", eq.lo = 0, eq.hi = .05, peta.h0 = 0){
  
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name))
  leg <- if(is.character(legend)) legend else deparse(substitute(legend))
      
  pr <- show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I <- eq(a, b, d, lo, hi, f, N, df1, df2)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]] ; f = I[[6]] ; N = I[[7]] ; df1 = I[[8]] ; df2 = I[[9]] 
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)                                                                                                                           
  
  loop <- length(a)
  CI <- matrix(NA, loop, 2)
  mode <- numeric(loop)
  peak <- numeric(loop)
  h <- list()
  k = numeric(loop)
  eq.prob = numeric(loop)
  BF10 = numeric(loop)
  estimate = numeric(loop)
  
  if(!pr){  
    
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
      CI[i,] = HDI(posterior, 0, .9999999)
      BF10[i] =  k[i] / df(f[i], df1[i], df2[i], (peta.h0 * N[i]) / (1 - peta.h0))
      eq.prob[i] = integrate(posterior, lo[i], eq.hi)[[1]] - integrate(posterior, lo[i], eq.lo)[[1]]
      estimate[i] <- (f[i]*df1[i]) / ((f[i]*df1[i]) + df2[i])
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
    return(round(data.frame(estimate = estimate, mode = mode, lower = CI[,1], upper = CI[,2], eq.prob = eq.prob, BF10 = BF10, row.names = paste0("P.eta.sq ", 1:loop, " posterior: ")), 6))  
    
}else{
    p = function(x) { get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = bquote(bold("Partial Eta.Sq"~(eta[p]^2))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(eta[p]^2*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}


#===================================================================================================================


peta.priors <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", scale = .1, top = 1.5, show.prior = FALSE, bottom = 1, legend = "topleft")
{
  UseMethod("peta.priors")
}

peta.priors.default <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", scale = .1, top = 1.5, show.prior = FALSE, bottom = 1, legend = "topleft"){
  
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


cor.bayes <- function(r, n, prior.mean = 0, prior.sd = .707, eq.bound = .05, level = .95, top = 1, bottom = 1, scale = .1, margin = 5, legend = "topleft", show.prior = FALSE)
{
  UseMethod("cor.bayes")
}

cor.bayes.default <- function(r, n, prior.mean = 0, prior.sd = .707, eq.bound = .05, level = .95, top = 1, bottom = 1, scale = .1, margin = 5, legend = "topleft", show.prior = FALSE){ 
  
  pr <- show.prior    
  mu <- prior.mean
  lambda <- 1/(prior.sd^2)
    
  if(!pr){   
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
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
  
  plot(CI, rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*loop), ylab = NA, yaxt = "n", xlab = "Credible Interval (Pearson correlation)", font.lab = 2, mgp = c(2, .3, 0))
  axis(2, at = 1:loop, labels = paste0("r", 1:loop), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
  abline(h = 1:loop, col = 8, lty = 3)
  legend(x = leg, legend = rev(paste0("norm", "(", round(prior.mean, 2), ", ", round(prior.sd, 2), ")")), pch = 22, title = "Priors", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
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
  
  return(round(data.frame(mean = mean, mode = mode, median = median, sd = sd, lower = CI[,1], upper = CI[,2], eq.prob = eq.prob, row.names = paste0("r", 1:loop, " posterior: ")), 6))

  }else{
  
  p <- function(x) dnorm(atanh(x), prior.mean[1], prior.sd[1])*1/(1-x^2)
  curve(p, -1, 1, yaxt = "n", ylab = NA, xlab = bquote(rho[~("Pearson correlation")]), bty = "n", font.lab = 2, cex.lab = 1.5, lwd = 2, n = 1e4, yaxs = "i", main = bquote(rho*" ~ "*"scaled.norm"(.(round(prior.mean[1], 3)), .(round(prior.sd[1], 3)))), xpd = NA) 
  }  
}

                              
#==================================================================================================================                              

                              
cor.diff <- function(r, n, prior.mean = 0, prior.sd = .707, how = c("two.one", "one.two"), eq.bound = .05, level = .95, top = 1, bottom = 1, scale = .1, margin = 5, legend = "topleft")
{
  UseMethod("cor.diff")
}

cor.diff.default <- function(r, n, prior.mean = 0, prior.sd = .707, how = c("two.one", "one.two"), eq.bound = .05, level = .95, top = 1, bottom = 1, scale = .1, margin = 5, legend = "topleft"){ 
  
  is.s <- function(...)lengths(list(...)) < 2
  if(any(is.s(n, r))) stop("Error: 'r' & 'n' must each have a length of '2' or larger.")
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
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
  
  plot(CI, rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*loop), ylab = NA, yaxt = "n", xlab = "Credible Interval (Correlation Differences)", font.lab = 2, mgp = c(2, .3, 0))
  axis(2, at = 1:loop, labels = np, font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
  abline(h = 1:loop, col = 8, lty = 3)
  legend(x = legn, legend = rev(paste0("norm", "(", round(rep(prior.mean[1], loop), 2), ", ", round(rep(prior.sd[1], loop), 2), ")")), pch = 22, title = "Priors", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
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
  
  return(round(data.frame(mean = mean, mode = mode, median = median, sd = sd, lower = CI[,1], upper = CI[,2], eq.prob = BB, row.names = paste0(np, ":")), 6))
}

              
#===================================================================================================================

              
prop.update <- function(n = 100, yes = 55, top = 5, scale = .1, lo = 0, hi = 1, a = 1.2, b = 1.2, dist.name = "dbeta", prior.scale = 1, level = .95, show.prior = FALSE, tol = 1e5)
{
  UseMethod("prop.update")
}

prop.update.default <- function(n = 100, yes = 55, top = 5, scale = .1, lo = 0, hi = 1, a = 1.2, b = 1.2, dist.name = "dbeta", prior.scale = 1, level = .95, show.prior = FALSE, tol = 1e5){
  
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
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
  I <- eq(n, s) ; n <- I[[1]] ; s <- I[[2]]
  loop <- length(n) 
  
  props <- seq(0, 1, 1/tol)
  prx <- get(d)(props, a, b)*as.integer(props >= lo)*as.integer(props <= hi)
  pr <- tol * prx / sum(prx)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~props, ylim = c(0, top*loop), type = "n", yaxs = "i", ylab = NA, xlab = "Proportion", font.lab = 2, axes = FALSE, mgp = c(2, .4, 0),main = if(pri) bquote(Proportion*" ~ "*.(if(lo > 0 || hi < 1) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))) else NA)
  axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
  
  if(!pri){  
    abline(h = 1:loop+1, col = 8, lty = 3)
    axis(2, at = 0:loop+1, lab = c("Base knowledge", paste0("Study ", 1:loop)), las = 1, font = 2, cex.axis = .9, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
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
      polygon(y = scale*ps+i+1, x = props, col = adjustcolor(i+1, .5), border = NA, xpd = NA)
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

d.update <- function(t, n1, n2 = NA, top = 5, scale = .1, m = 0, s = 1, dist.name = "dnorm", prior.scale = 1, level = .95, show.prior = FALSE, lo = -2, hi = 2, tol = 1e4, margin = hi)
{
  UseMethod("d.update")
}

d.update.default <- function(t, n1, n2 = NA, top = 5, scale = .1, m = 0, s = 1, dist.name = "dnorm", prior.scale = 1, level = .95, show.prior = FALSE, lo = -2, hi = 2, tol = 1e4, margin = hi){
  
  pri <- show.prior
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  if(is.infinite(lo)) lo <- -6
  if(is.infinite(hi)) hi <- 6
  if(tol < 1e4) stop("'tol' must be '10,000' or larger.")
  is.v <- function(...) lengths(list(...)) > 1
  if(any(is.v(m, s, d))) stop("Choose only 'one' prior knowledge base at a time.")
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
  I <- eq(t, n1, n2) 
  t <- I[[1]]  
  n1 <- I[[2]]  
  n2 <- I[[3]] 
  loop <- length(t) 
  
  ds <- seq(-6, 6, 1/tol)
  prx <- get(d)(ds, m, s)*as.integer(ds >= lo)*as.integer(ds <= hi)
  pr <- tol * prx / sum(prx)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~ds, ylim = c(0, top*loop), xlim = c(-margin, margin), type = "n", xaxs = "i", yaxs = "i", ylab = NA, xlab = bquote(bold("Effect Size"~ (delta))), font.lab = 2, mgp = c(2, .4, 0),main = if(pri) bquote("Effect Size"*" ~ "*.(if(lo > -Inf || hi < Inf) "truncated-")*.(substring(d, 2))(.(round(m, 2)), .(round(s, 2)))) else NA, yaxt = "n", bty = "n")
  
  if(!pri){  
    abline(h = 1:loop+1, col = 8, lty = 3)
    axis(2, at = 0:loop+1, lab = c("Base knowledge", paste0("Study ", 1:loop)), las = 1, font = 2, cex.axis = .9, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
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
      polygon(y = scale*ps+i+1, x = ds, col = adjustcolor(i+1, .5), border = NA, xpd = NA)
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

peta.update <- function(f, N, df1, df2, top = 5, scale = .1, a = 2, b = 2, lo = 0, hi = 1, dist.name = "dbeta", prior.scale = 1, level = .95, show.prior = FALSE, tol = 1e5)
{
  UseMethod("peta.update")
}

peta.update.default <- function(f, N, df1, df2, top = 5, scale = .1, a = 2, b = 2, lo = 0, hi = 1, dist.name = "dbeta", prior.scale = 1, level = .95, show.prior = FALSE, tol = 1e5){
  
  pri <- show.prior
  d <- if(is.character(dist.name)) dist.name else deparse(substitute(dist.name)) 
  if(hi == 1) hi <- .9999999 ;
  if(tol < 1e4) stop("'tol' must be '10,000' or larger.")
  is.v <- function(...) lengths(list(...)) > 1
  if(any(is.v(a, b, d))) stop("Error: Choose only 'one' prior knowledge base at a time.")
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
  I <- eq(f, N, df1, df2) ; f <- I[[1]] ; N <- I[[2]] ; df1 <- I[[3]] ; df2 <- I[[4]] ;
  loop <- length(f) 
  
  peta <- seq(0, .9999999, 1/tol)
  prx <- get(d)(peta, a, b)*as.integer(peta >= lo)*as.integer(peta <= hi)
  pr <- tol * prx / sum(prx)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~peta, ylim = c(0, top*loop), type = "n", yaxs = "i", ylab = NA, xlab = bquote(bold("Partial Eta.Sq"~(eta[p]^2))), font.lab = 2, axes = FALSE, mgp = c(2, .4, 0),main = if(pri) bquote(eta[p]^2*" ~ "*.(if(lo > 0 || hi < .9999999) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))) else NA)
  axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
  
  if(!pri){  
    abline(h = 1:loop+1, col = 8, lty = 3)
    axis(2, at = 0:loop+1, lab = c("Base knowledge", paste0("Study ", 1:loop)), las = 1, font = 2, cex.axis = .9, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
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
      polygon(y = c(i+1, scale*ps+i+1, i+1), x = c(lo, peta, hi), col = adjustcolor(i+1, .5), border = NA, xpd = NA)
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
   
                       
need <- c("rstanarm") #, "MASS")
have <- need %in% rownames(installed.packages())
if(any(!have)){ install.packages( need[!have] ) }
                       
suppressMessages({ 
    library("rstanarm")
  #  library("MASS")
})
                      
                       
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


R2.bayes <- function(..., scale = .02, bottom = 1, top = 1, margin = 5, legend = "topleft", level = .95, eq.lo = 0, eq.hi = .1)
{
  UseMethod("R2.bayes")
}                        
                                           
                       
R2.bayes.default <- function(..., scale = .02, bottom = 1, top = 1, margin = 5, legend = "topleft", level = .95, eq.lo = 0, eq.hi = .1)
{

if(!(all(sapply(list(...), inherits, "stanreg")))) stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")    
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

plot(1, loop, type = "n", xlim = c(a, b), ylim = c(bottom*1, top*loop), ylab = NA, xaxt = "n", yaxt = "n", xlab = bquote(bold("Regression Coefficient " (R^2))), mgp = c(2, .3, 0))
abline(h = 1:loop, col = 8, lty = 3)

for(i in 1:loop){
polygon(x = d[[i]]$x, y = scale*d[[i]]$y + i, col = adjustcolor(i, .55), border = NA, xpd = NA)
}

axis(1, at = seq(a, b, length.out = 4), labels = paste0(round(seq(a, b, length.out = 4), 4)*1e2, "%"), mgp = c(2, .5, 0))
axis(2, at = 1:loop, labels = paste0("Model ", 1:loop), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))

legend(x = leg, legend = rev(paste0("Model ", loop:1)), pch = 22, title = "Models ", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5)
segments(I[, 1], 1:loop, I[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
box()

m = scale*peak + 1:loop
segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.5, col = "magenta", xpd = NA)

q = deci(I*1e2 , 2); o = deci(mode*1e2, 2)
text(mode, 1:loop, paste0(q[,1], "%", "    ", o, "%", "    ", q[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)

round(data.frame(mode = mode, mean = mean, sd = sd, MAD = mad, lower = I[,1], upper = I[,2], coverage = level, eq.prob = eq.prob, row.names = paste0("Model-", 1:loop, " posterior:")), 6)
}


#=======================================================================
                         
                       
#lm.sample2 <- function(fit, n = 1e4, no.names = TRUE)
#{
#  UseMethod("lm.sample2")
#}                       
                
                       
#lm.sample2.default <- function(fit, n = 1e4, no.names = TRUE){
 
#if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")
    
#  output <- as.data.frame(MASS::mvrnorm(n = n, mu = c(coef(fit), sigma(fit)), Sigma = cov(as.matrix(fit))))
  
#  if(no.names == TRUE){
#    for(i in 1:ncol(output)){
#     if(colnames(output)[i] == "(Intercept)"){
#        colnames(output)[i] <- "Intercept"
#      }
#      if(colnames(output)[i] == paste0("V", ncol(output))){
#        colnames(output)[i] <- "Sigma"
#      }
#    }
#  }
#  output
# }
                     
                       
#======================================================================================
                       
                       
lm.sample <- function(fit, no.names = TRUE)
{
  UseMethod("lm.sample")
}                       
                
                       
lm.sample.default <- function(fit, no.names = TRUE){
 
if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")
    
  output <- as.data.frame(fit)
  
  if(no.names == TRUE){
    for(i in 1:ncol(output)){
      if(colnames(output)[i] == "(Intercept)"){
        colnames(output)[i] <- "Intercept"
      }
    }
  }
  output
}
                       
                       
#======================================================================================

 
predict.mean <- function(fit, predi, scale = .5, level = .95, col.hump = "cyan", ...)
{
  UseMethod("predict.mean")
} 
       
                       
predict.mean.default <- function(fit, predi, scale = .5, level = .95, col.hump = "cyan", ...){
  
if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")  
if(length(coef(fit)) > 2) stop("Error: 'fit' must contain only 'one' predictor.")  

post <- lm.sample(fit)

mus_at_xi = post[,1] + post[,2] * predi

d <- density(mus_at_xi, adjust = 2, n = 1e3)
plot(d, type = "n", ylab = NA, main = NA, yaxt = "n", bty = "n", las = 1, zero.line = FALSE, yaxs = "i",
     xlab = bquote(bold((mu[i] *" | "* .(names(fit$model)[2])[i] == .(predi)))), ...)

  I <- hdir(mus_at_xi, level = level)
med <- mean(mus_at_xi)
peak <- approx(d$x, d$y, xout = med)[[2]]*scale

polygon(d$x, scale*d$y, col = adjustcolor(col.hump, .35), border = NA)
segments(med, 0, med, peak, lty = 3)

segments(I[1], 0, I[2], 0, lend = 1, lwd = 6, col = 'magenta', xpd = NA)
points(med, 0, pch = 21, bg = "cyan", col = 'magenta', cex = 2, xpd = NA)
text(c(I, med), 0, round(c(I, med), 2), pos = 3, font = 2)
}                       
                       
 
                       
#======================================================================================                       

                       
                       
predict.case <- function(fit, predi, scale = .5, level = .95, col.hump = "gray", ...)
{
  UseMethod("predict.case")
} 
                       
                       
predict.case.default <- function(fit, predi, scale = .5, level = .95, col.hump = "gray", ...){
  
  if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")  
  if(length(coef(fit)) > 2) stop("Error: 'fit' must contain only 'one' predictor.")  
  
  post <- lm.sample(fit)
  
  case_at_xi = post[,1] + post[,2] * predi + post[,3]
  
  d <- density(case_at_xi, adjust = 2, n = 1e3)
  plot(d, type = "n", ylab = NA, main = NA, yaxt = "n", bty = "n", las = 1, zero.line = FALSE, yaxs = "i",
       xlab = bquote(bold((subj[i] *" | "* .(names(fit$model)[2])[i] == .(predi)))), ...)
  
  I <- hdir(case_at_xi, level = level)
  med <- mean(case_at_xi)
  peak <- approx(d$x, d$y, xout = med)[[2]]*scale
  
  polygon(d$x, scale*d$y, col = col.hump, border = NA)
  segments(med, 0, med, peak, lty = 3)
  
  segments(I[1], 0, I[2], 0, lend = 1, lwd = 6, xpd = NA)
  points(med, 0, pch = 21, bg = "cyan", col = 'magenta', cex = 2, xpd = NA)
  text(c(I, med), 0, round(c(I, med), 2), pos = 3, font = 2)
} 
                       
                       
#======================================================================================                       

                       
predict.bayes <- function(fit, xlab = NA, ylab = NA, level = .95, line.int = TRUE, pred.int = TRUE, pt.cex = 1, pt.col = 4, col.depth = .3, col.line = "cyan", col.pred = "gray", col.reg = "cyan", ...)
{
  UseMethod("predict.bayes")
} 

                       
predict.bayes.default <- function(fit, xlab = NA, ylab = NA, level = .95, line.int = TRUE, pred.int = TRUE, pt.cex = 1, pt.col = 4, col.depth = .3, col.line = "cyan", col.pred = "gray", col.reg = "cyan", ...){

if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")  
if(length(coef(fit)) > 2) stop("Error: 'fit' must contain only 'one' predictor.")
    
pred <- fit$model[, 2]
 dep <- fit$model[, 1]  

plot(dep ~ pred, xlab = ifelse(is.na(xlab), names(fit$model)[2], xlab), ylab = ifelse(is.na(ylab), names(fit$model)[1], ylab), type = "n", las = 1, ...)

pred_lin <- rstanarm::posterior_predict(fit, transform = TRUE)

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

polygon(c(rev(x), x), c(rev(z), y), col = adjustcolor(col.pred, .5), border = NA)
}
    
pred_lin2 <- rstanarm::posterior_linpred(fit, transform = TRUE)

I2 <- matrix(NA, loop, 2)
for(i in 1:loop){
I2[i,] = hdir(pred_lin2[,i], level = level)
}

if(line.int){
    
y <- I2[,1][o]
z <- I2[,2][o]

polygon(c(rev(x), x), c(rev(z), y), col = adjustcolor(col.line, col.depth), border = NA)
}
    
abline(fit, col = col.reg, lwd = 2)
box()    
}                       
                       
#==================================================================================
                  
                       
compare.R2 <- function(..., how = c("two.one", "one.two"), scale = .02, bottom = 1, top = 1, margin = 5, legend = "topleft", level = .95, eq.level = "2.5%")
{
  UseMethod("compare.R2")
} 

                       
compare.R2.default <- function(..., how = c("two.one", "one.two"), scale = .02, bottom = 1, top = 1, margin = 5, legend = "topleft", level = .95, eq.level = "2.5%"){
  
  if(!(all(sapply(list(...), inherits, "stanreg")))) stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")
  
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
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mar = c(5.1, 6.1, 4.1, 2.1))
  
  plot(1, loop, type = "n", xlim = c(a, b), ylim = c(bottom*1, top*loop), ylab = NA, xaxt = "n", yaxt = "n", xlab = bquote(bold(Delta~R^2~("Model Comparison"))), mgp = c(2, .3, 0))
  abline(h = 1:loop, col = 8, lty = 3)
  axis(1, at = seq(a, b, length.out = 4), labels = paste0(round(seq(a, b, length.out = 4), 4)*1e2, "%"), mgp = c(2, .5, 0))
  axis(2, at = 1:loop, labels = np, font = 2, las = 1, cex.axis = .7, tck = -.006, mgp = c(2, .3, 0))
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
  
  return(round(data.frame(mean = mean, mode = mode, median = median, sd = sd, lower = CI[,1], upper = CI[,2], coverage = level, eq.prob = eq.prob, row.names = paste0(np, ":")), 6))
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
              
              
type.sm <- function(d = .1, obs.d = .6, n1 = 20, n2 = NA)
{
  UseMethod("type.sm")
}

type.sm.default <- function(d = .1, obs.d = .6, n1 = 20, n2 = NA){
  
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
  
  legend("topright", "Observed \nS.S. Effect", pch = 23, pt.bg = 3, pt.cex = 1.2, bty = "n", text.font = 2)   
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
  
  round(data.frame(exaggration = exaggration, type.s = type.s, power = power, Crit.d = CI[2], p.value = p.value, row.names = "Results:"), 6)
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
              
              
standard.fit <- function(fit, level = .95, digit = 6)
{
  UseMethod("standard.fit")
}


standard.fit.default <- function(fit, level = .95, digit = 6){

if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")
  
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

round(data.frame(standard.coef = mean, sd = sd, lower = I[,1], upper = I[,2], coverage = level, row.names = colnames(b)), digit = digit)
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
              
 
count.plot <- function(fit, xlab = NA, ylab = NA, line.int = TRUE, pred.int = TRUE, level = .95,
                       focus.pred, n = 2e2, FUN = mean, hold.at = NA, legend = "topleft", ...)
{
  UseMethod("count.plot")
}


count.plot.default <- function(fit, xlab = NA, ylab = NA, line.int = TRUE, pred.int = TRUE, level = .95,
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
              
              
case.fit.plot <- function(fit, level = .95, legend = "topleft", lwd = 2, fit.tol = 1, pt.cex = 1, cex.axis = .8)
{
  UseMethod("case.fit.plot")
}  


case.fit.plot.default <- function(fit, level = .95, legend = "topleft", lwd = 2, fit.tol = 1, pt.cex = 1, cex.axis = .8){
  
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

original.par = par(no.readonly = TRUE)
on.exit(par(original.par))    
par(mar = c(2.6, 4, 2.2, 4))
    
plot(e[o], 1:loop, cex = .6, xlim = range(PI.e), pch = 19, ylab = NA, yaxt = "n", mgp = c(1.3, .4, 0), type = "n", xlab = "Credible Interval (Residuals)", font.lab = 2)
rect(-unit, par('usr')[3], unit, par('usr')[4], border = NA, col = adjustcolor("yellow", .4), lend = 1)

abline(v = 0, h = 1:loop, lty = 3, col = 8)

good <- -unit < e[o] & e[o] < unit
    
pos <- (1:loop)[o]

axis(2, at = (1:loop)[-range(1:loop)], labels = paste0("subj ", pos[-range(pos)]), las = 1, cex.axis = cex.axis, tck = -.006, mgp = c(2, .2, 0))
axis(2, at = range(1:loop), labels = paste0("subj ", c(pos[1], rev(pos)[1])), las = 1, cex.axis = cex.axis, tck = -.006, mgp = c(2, .2, 0), col.axis = 2)

segments(PI.e[, 1], 1:loop, PI.e[, 2], 1:loop, lend = 1, col = 8, lwd = lwd)

segments(CI.e[, 1], 1:loop, CI.e[, 2], 1:loop, lend = 1, lwd = lwd, col = ifelse(good, "green3", 1))

points(e[o], 1:loop, pch = 21, bg = ifelse(ok, "cyan", 2), col = ifelse(ok, "magenta", 2), cex = pt.cex)

a <- par('usr')[1:2] ; b <- par('usr')[4]
    
text(.1*mean(a), b, "Fit Corridor", pos = 3, cex = 1, xpd = NA, font = 2, col = 3)
    
text(.8*a, b, rep("Misfit", 2), pos = 3, cex = 1.5, xpd = NA, font = 2, col = 2)

legend(x = leg, legend = c("Worst fit", "Perfect-Good fit", "Fair-Bad fit"), pch = 22, title = "Person Fit", 
       pt.bg = c(2, "green3", 1), col = c(2, "green3", 1), cex = .7, pt.cex = .6, bg = 0, 
       box.col = 0, xpd = NA, x.intersp = .5, title.adj = .2)
box()
}
        
              
#====================================================================================================================
      
              
cor.fit <- function(fit, cor = TRUE){
  
  UseMethod("cor.fit")
}

cor.fit.default <- function(fit, cor = TRUE){

if(class(fit)[1] != "stanreg") stop("Error: 'fit' must be from package 'rstanarm's 'stan_glm()'.")    

m <- if(cor) cov2cor(cov(as.matrix(fit))) else cov(as.matrix(fit))

colnames(m)[1] <- "intercept"
rownames(m)[1] <- "intercept"

return(round(data.frame(m), 6))
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
              
