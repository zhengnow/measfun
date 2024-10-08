#'
#' Calculate Cohen's kappa with a confidence interval
#'
#' The function estimates Cohen's kappa for 2 ratings per subject with a confidence interval, for independent subjects.
#' For clustered subjects, the confidence interval is calculated using bootstrapping approach.
#'
#' @param ratings dataframe/matrix, 2 columns
#' @param weight one of c("equal","squared", "unweighted"), default is "equal"
#' @param rating.level default is NULL. Recommend specifying when ratings are factors.
#' @param cluster TRUE/FALSE, whether the data is clustered
#' @param clu.id subject id, when data is clustered
#' @param n.boot number of repetitions in bootstrapping, default is 500
#' @param table If True, a detailed cross table will be created, along with the kappa estimates. Default is False.
#'
#' @return a vector for the estimated kappa with the confidence interval
#'
#' @details
#' kappa2r calculates Cohen's kappa with bootstrapped bias corrected confidence intervals for clustered subjects.
#'
#'
#' @examples
#' set.seed(234); x1 <- rbinom(n=50, size=2, prob=0.3)
#' set.seed(393); x2 <- rbinom(n=50, size=2, prob=0.4)
#' id <- c(1:25, 1:25)
#' kappa2r(cbind(x1,x2), weight="squared", n.boot=500, cluster=TRUE, clu.id = id)
#'
#' @export
#'


# -----------------------------
# kappa and bootstrapped 95%CI
# -----------------------------
# agreement on categorical variables

kappa2r <-
  function (ratings,
            weight = c("unweighted", "equal", "squared"),
            rating.level = NULL,
            ci.level = 0.95,
            cluster = FALSE,
            clu.id = NULL,
            n.boot = 500,
            sort.levels = FALSE,
            ci.level = 0.95, ...) {

  ratings <- na.omit(ratings)
  if (is.character(weight)) weight = match.arg(weight)
  ns <- nrow(ratings)
  nr <- ncol(ratings)
  if (nr > 2) {
    stop("Number of raters exeeds 2. Try fleiss.kappa or lights.kappam")
  }

  r1 <- ratings[, 1]
  r2 <- ratings[, 2]

  if ((is.numeric(r1)) | (is.numeric(r2)))
    sort.levels <- TRUE
  if (!is.factor(r1))
    r1 <- factor(r1)
  if (!is.factor(r2))
    r2 <- factor(r2)
  if (length(levels(r1)) >= length(levels(r2))) {
    lev <- c(levels(r1), levels(r2))
  }
  else {
    lev <- c(levels(r2), levels(r1))
  }

  if (sort.levels) lev <- sort(lev)

  lev <- lev[!duplicated(lev)]

  r1 <- factor(ratings[, 1], levels = lev)
  r2 <- factor(ratings[, 2], levels = lev)
  ttab <- table(r1, r2)
  nc <- ncol(ttab)
  if (is.numeric(weight))
    w <- 1 - (weight - min(weight))/(max(weight) - min(weight))
  else if (weight == "equal")
    w <- (nc - 1):0/(nc - 1)
  else if (weight == "squared")
    w <- 1 - (0:(nc - 1))^2/(nc - 1)^2
  else w <- c(1, rep(0, nc - 1))
  wvec <- c(sort(w, decreasing = FALSE), w[2:length(w)])
  nw <- length(w)
  weighttab <- matrix(0, nrow = nw, ncol = nw)
  for (i in 1:nw) {
    weighttab[i, ] <- wvec[(nw - (i - 1)):(2 * nw - i)]
  }
  agreeP <- sum(ttab * weighttab)/ns
  tm1 <- apply(ttab, 1, sum)
  tm2 <- apply(ttab, 2, sum)
  eij <- outer(tm1, tm2)/ns
  chanceP <- sum(eij * weighttab)/ns
  value <- (agreeP - chanceP)/(1 - chanceP)

  N. <- 1 - (1 - ci.level)/2
  z <- qnorm(N., mean = 0, sd = 1)

  if (!cluster){
  w.i <- apply(rep(tm2/ns, nc) * weighttab, 2, sum)
  w.j <- apply(rep(tm1/ns, each = nc) * weighttab, 1, sum)
  var.matrix <- (eij/ns) * (weighttab - outer(w.i, w.j, "+"))^2
  varkappa <- (sum(var.matrix) - chanceP^2)/(ns * (1 - chanceP)^2)
  SEkappa <- sqrt(varkappa)
  upper_ci <- value + z*SEkappa
  lower_ci <- value - z*SEkappa

  } else {
    # for clustered data
        set.seed(338)
        k2 <-
          c(irr::kappa2(cbind(r1,r2),
                       weight = weight)$value,
                boot::boot.ci(boot::boot(unique(clu.id),
                             kappa2.clu.boot,
                             n.boot,
                             datafull=cbind(r1,r2,clu.id),
                             weight2 = weight),
                             conf = ci.level,
                        type="bca")$bca[,4:5])

        value <- k2[1]
        lower_ci <- k2[2]
        upper_ci <- k2[3]
  }
  if (!table){
  rval <-
    structure(list(method = paste0("Cohen's Kappa for 2 Raters (Weights: ",
                                   paste(weight, collapse = ","), ")"),
                   subjects = ns,  raters = nr,
                   meas.name = "kappa",
                   cluster = cluster,
                   value = value,
                   lower_ci = lower_ci, upper_ci = upper_ci),
              class = "measfunlist")
  } else {
    rater <- names(ratings)
    data.frame(ratings) %>%
      gtsummary::tbl_summary(
        by = rater[2],
        statistic = list(all_categorial()~ "{n}")) %>%
      gtsummary::add_stat(
        fns = everything() ~ gtkappa2
      )
  }

  return(rval)
}

