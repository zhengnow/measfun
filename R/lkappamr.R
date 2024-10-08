
#'
#' Calculate Lights's kappa with the confidence interval
#'
#' The function estimates Light's kappa for multiple ratings with the confidence interval.
#' Weights can be specified.
#'
#' @param ratings a matrix, categorical ratings
#' @param weight one of c("equal","squared", "unweighted"), default is "equal"
#' #' @param rating.level default is NULL. Recommend specifying when ratings are factors.
#'
#' @return a list including n, the number of ratings, Light's kappa, and a p value
#'
#' @details
#' lkappamr calculates Lights's kappa with 95% bootstrapped confidence intervals.
#' The confidence interval is the bias corrected ("bca") confidence interval. For clustered data, use kappa2r.clu
#' Check the categories in every rating are in the same order.
#'
#' @examples
#' set.seed(234); x1 <- rbinom(n=50, size=2, prob=0.3)
#' set.seed(393); x2 <- rbinom(n=50, size=2, prob=0.4)
#' set.seed(395); x3 <- rbinom(n=50, size=2, prob=0.2)
#' lkappamr(cbind(x1,x2,x3), weight="equal")
#'
#' @export
#'


# --------------------------
# Light's kappa, add weights
# --------------------------
# under construction
lkappamr <-
  function(ratings,
           weight=c("equal", "unweighted", "squared"),
           rating.level = NULL,
           ci.level = 0.95, ...) {

    ratings <- na.omit(ratings)

    if(!is.null(rating.level)) lev <- rating.level
    else {
      ratings <- as.matrix(ratings)
      lev <- levels(as.factor(ratings))
    }
    levlen <- length(levels(as.factor(ratings)))

    if (is.character(weight))
      weight = match.arg(weight)

    ns <- nrow(ratings)
    nr <- ncol(ratings)
    for (i in 1:(nr - 1)) for (j in (i + 1):nr) {
      if ((i == 1) & (j == (i + 1)))
        kappas <- irr::kappa2(ratings[, c(i, j)], weight = weight)$value
      else kappas <- c(kappas, irr::kappa2(ratings[, c(i, j)], weight = "u")$value)
    }
    value <- mean(kappas)

    for (nri in 1:(nr - 1)) for (nrj in (nri + 1):nr) {
      for (i in 1:levlen) for (j in 1:levlen) {
        if (i != j) {
          r1i <- sum(ratings[, nri] == lev[i])
          r2j <- sum(ratings[, nrj] == lev[j])
          if (!exists("dis"))
            dis <- r1i * r2j
          else dis <- c(dis, r1i * r2j)
        }
      }
      if (!exists("disrater"))
        disrater <- sum(dis)
      else disrater <- c(disrater, sum(dis))
      rm(dis)
    }
    B <- length(disrater) * prod(disrater)
    chanceP <- 1 - B/ns^(choose(nr, 2) * 2)
    varkappa <- chanceP/(ns * (1 - chanceP))
    SEkappa <- sqrt(varkappa)
    clvl <- qnorm((1+ci.level)/2)
    upper_ci <- value + clvl*SEkappa
    lower_ci <- value - clvl*SEkappa

    rval <-
      structure(
        list(method = paste0("Light's Kappa for m Ratings (Weights: ",
                             paste(weight, collapse = ","), ")"),
             subjects = ns, raters = nr,

             meas.name = "kappa",
             value = value,
             lower_ci = lower_ci, upper_ci = upper_ci),
        class = "measfunlist")
    return(rval)
  }
