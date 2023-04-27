#'
#' Calculate Cohen's kappa with the confidence interval
#'
#' The function estimates Cohen's kappa for 2 ratings with bootstrapped confidence interval.
#'
#' @param v1,v2 vectors, categorical ratings
#' @param weight one of c("equal","squared", "unweighted"), default is "equal"
#' @param cluster TRUE/FALSE, whether the data is clustered
#' @param clu.id subject id, when data is clustered
#'
#' @return a vector for the estimated kappa with the 95%CI
#'
#' @details
#' kappa2r calculates Cohen's kappa with 95% bootstrapped confidence intervals.
#' The confidence interval is the bias corrected ("bca") confidence interval. For clustered data, use kappa2r.clu
#' Check the categories in 2 ratings before running the function.
#' You should have them categoried in the same order.
#'
#' @examples
#' set.seed(234); x1 <- rbinom(n=50, size=2, prob=0.3)
#' set.seed(393); x2 <- rbinom(n=50, size=2, prob=0.4)
#' id <- c(1:25, 1:25)
#' kappa2r(x1,x2, weight="squared", nboot=500, cluster=TRUE, clu.id = id)
#'
#' @export
#'


# -----------------------------
# kappa and bootstrapped 95%CI
# -----------------------------
# agreement on categorical vars

# function output kappa for 2 ratings and bootstrapped confidence intervals
kappa2r <-
  function(v1, v2, weight = "equal", nboot = 500, cluster = FALSE, clu.id=NULL){
    if (!cluster){
      set.seed(338)
      k2 <- round(
        c(irr::kappa2(cbind(v1,v2),
                 weight = weight)$value,
          boot::boot.ci(boot::boot(cbind(v1,v2),
                       kappa2.boot, nboot, weight2 = weight),
                  type="bca")$bca[,4:5]),
        3)
      namekp <- "kappa (95%CI)"
    }
    if (cluster){
      set.seed(338)
      k2 <- round(
        c(irr::kappa2(cbind(v1,v2),
                 weight = weight)$value,
          boot::boot.ci(boot::boot(unique(clu.id),
                       kappa2.clu.boot, 1000,
                       datafull=cbind(v1,v2,clu.id),
                       weight2 = weight),
                  type="bca")$bca[,4:5]),
        3)
      kp <- paste0(k2[1], " (95%CI: ",
                      k2[2],", ",
                      k2[3],")")
      namekp <- "Clustered kappa (95%CI)"
    }
    kp <- paste0(k2[1], "(",
                k2[2],", ",
                k2[3],")")
    names(kp) <- namekp
    return(noquote(kp))
  }

