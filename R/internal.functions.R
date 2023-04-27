

.First <- function() {
  options(
    repos = c(CRAN = "http://cran.rstudio.com/"),
    browserNLdisabled = TRUE,
    deparse.max.lines = 2)
}

if (interactive()) {
  suppressMessages(require(devtools))
}


# function to stop the function when exprs is not satisfied
assert_stopifnot <-
  function(exprs) eval.parent(substitute(stopifnot(exprs = exprs)))

# ci.auc


# a function to cap the 95%CI of measures of accuracy (0-1)
capci <- function(x){
  # lower bound
  x[1] = max(x[1],0)
  # upper bound
  x[2] = min(x[2],1)
  return(x)
}

# function to bootstrap for kappa confidence interval
# independent values
kappa2.boot <-
  function(data, x, weight2 = NULL){
    irr::kappa2(data[x,],
           weight = weight2)$value
  }

# function to bootstrap for kappa confidence interval
# clustered data
kappa2.clu.boot <-
  function(data, x, datafull, weight2 = NULL){
    # id indicator is listed in the 3rd column of the data
    # When a cluster id is sampled, include all observations of the cluster id.
    # data1 <- data[data[,3] %in% data[x,3], ]
    data1 <- datafull[datafull[,3] %in% x, ]
    irr::kappa2(data1[,1:2],
           weight = weight2)$value
  }
