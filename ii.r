


install.packages("distr") 
library(distr)   # load the package


d.interact <- function(dppc, dppt, nc, nt, digits = 6, d.per.study = NA, long = NA, extract = NA){

ll <- d.per.study

if(!is.na(long)) nm <- if(long) "long" else "short"

G <- Vectorize(function(dppc, dppt, nc, nt, digits){

like1 <- function(x) dt(dppc*sqrt(nc), df = nc - 1, ncp = x*sqrt(nc))
like2 <- function(x) dt(dppt*sqrt(nt), df = nt - 1, ncp = x*sqrt(nt))

d1 <- AbscontDistribution(d = like1)
d2 <- AbscontDistribution(d = like2)

like.dif <- function(x) d(d2 - d1)(x)

Mean <- integrate(function(x) x*like.dif(x), -Inf, Inf)[[1]]
  SE <- sqrt(integrate(function(x) x^2*like.dif(x), -Inf, Inf)[[1]] - Mean^2)

  return(round(c(d.interact = dppt-dppc, SE = SE), digits))
})

out <- data.frame(t(G(dppc = dppc, dppt = dppt, nc = nc, nt = nt, digits = digits)))
if(is.na(ll)) out else {

if(sum(ll) != nrow(out)) stop("Incorrect 'd.per.study' detected.", call. = FALSE)
if(!is.na(long))out[nm] <- long
  
h <- split(out, rep(seq_along(ll), ll))
names(h) <- paste0("Study", seq_along(h))
h <- lapply(h, `row.names<-`, NULL)

set <- if(extract == "long") long else if (extract == "short") !long
if(!is.na(long) & !is.na(extract)) lapply(h, function(x) subset(x, subset = set))
 }
}

# Example of Use:
d.interact(1:4, 1:4 , 30, 30, d.per.study = c(1, 2, 1), long = c(T, F, T, T), extract = "long")

