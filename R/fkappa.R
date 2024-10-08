#'
#' Calculate Fleiss kappa with the confidence interval
#'
#' The function estimates Fleiss kappa for multiple ratings with the confidence interval.
#'
#' @param ratings a matrix, categorical ratings
#' @param rater_per_subject integer
#' @param rating.level default is NULL. Recommend specifying when ratings are factors.
#' @param detail T/F, show rating details
#'
#' @return a list including n, the number of ratings, Fleiss kappa with confidence interval
#'
#' @details
#' kappa2r calculates Cohen's kappa with 95% bootstrapped confidence intervals.
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
# Fleiss kappa
# m readers per subject
# --------------------------
fkappa <- function(ratings,
                   rater_per_subject = NULL,
                   rating.level = NULL,
                   ci.level = 0.95,
                   detail = FALSE, ...)
{
  # if number of readers per subject is not specified, quit and print error message
  assert_stopifnot(is.na(rater_per_subject))

  ratings <- na.omit(ratings)
  if(!is.null(rating.level)) lev <- rating.level
    else {
      ratings <- as.matrix(ratings)
      lev <- levels(as.factor(ratings))
    }

  ns <- nrow(ratings)
  nr <- reader_per_subject
  ntolr <- ncol(ratings)

  for (i in 1:ns) {
    frow <- factor(ratings[i, ], levels = lev)
    if (i == 1)
      ttab <- as.numeric(table(frow))
    else ttab <- rbind(ttab, as.numeric(table(frow)))
  }
  ttab <- matrix(ttab, nrow = ns)
  agreeP <- sum((apply(ttab^2, 1, sum) - nr)/(nr * (nr - 1))/ns)

  method <- paste0("Fleiss' Kappa for ", ntolr," Raters, ",
                   rater_per_subject , " Raters per Subject")

  chanceP <- sum(apply(ttab, 2, sum)^2)/(ns * nr)^2

  value <- (agreeP - chanceP)/(1 - chanceP)

    pj <- apply(ttab, 2, sum)/(ns * nr)
    qj <- 1 - pj
    varkappa <- (2/(sum(pj * qj)^2 * (ns * nr * (nr - 1)))) *
      (sum(pj * qj)^2 - sum(pj * qj * (qj - pj)))
    SEkappa <- sqrt(varkappa)
    u <- value/SEkappa
    p.value <- 2 * (1 - pnorm(abs(u)))

    if (detail) {
      pj <- apply(ttab, 2, sum)/(ns * nr)
      pjk <- (apply(ttab^2, 2, sum) - ns * nr * pj)/(ns *
                                                       nr * (nr - 1) * pj)
      kappaK <- (pjk - pj)/(1 - pj)
      varkappaK <- 2/(ns * nr * (nr - 1))
      SEkappaK <- sqrt(varkappaK)
      clvl <- qnorm((1+ci.level)/2)
      upper_ci <- kappaK + clvl*SEkappaK
      lower_ci <- kappaK - clvl*SEkappaK

      #uK <- kappaK/SEkappaK
      #p.valueK <- 2 * (1 - pnorm(abs(uK)))
      tableK <- ttab
    }

    if (!detail) {
      rval <-
        structure(
          list(method = method, subjects = ns, raters = ntolr,
                   meas.name = "kappa",
                   value = value,
                   lower_ci = lower_ci, upper_ci = upper_ci),
          class = "measfunlist")
    }
    else {
      rval <-
        structure(
          list(method = method, subjects = ns, raters = ntolr,
                    value = value,
                   lower_ci = lower_ci, upper_ci = upper_ci,
                   detail = tableK),
          class = "measfunlist")
    }

  return(rval)
}
