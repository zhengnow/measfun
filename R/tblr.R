
#' To estimate the diagnostic likelihood ratios (DLR)
#'
#' @param var a vector, test variable
#' @param reference a binary reference variable
#' @param thres a numeric value, the threshold value when the test variable is continuous
#' @param ci TRUE/FALSE, whether to include the confidence interval in the output
#' @param ci.level the confidence level of the confidence interval, the default is 0.95
#'
#' @details
#' tblr outputs estimates of DLR+ and DLR- for independent subjects.
#'
#' @return a data frame for the estimates, with the optional 95% confidence intervals
#'
#' @examples
#' y <- rbinom(n=50, size=1, prob=0.3) # reference variable
#' x <- rnorm(n=50, mean=20, sd=5)
#' tblr(x,y, thres = 20)
#'
#' @export
#'


# --------------------------
# report SE, SP, LR+ and LR-
# --------------------------
tblr <- function(var, reference, thres=NULL, ci=TRUE){

  # if the reference has more than 2 levels, quit and print error message
  assert_stopifnot(length(table(reference))==2)

  # show the reference and test variable levels
  cat("Reference levels: ",
      names(table(reference))[1], " (negative), ",
      names(table(reference))[2], " (positive) ",
      "\n")

  # var can be numeric, character or factor values
  if(is.numeric(var)){
    assert_stopifnot(length(table(var))>0)
    if(length(table(var))<=2){var.thres <- var} else {
      var.thres <-
        ifelse(var>thres, paste0(">",thres), paste0("<=",thres))
    }
  } else if (!is.numeric(var)) {
    assert_stopifnot(length(table(var))==2)
    var.thres <- var
  }

  cat("Test levels: ",
      names(table(var.thres))[1], " (negative), ",
      names(table(var.thres))[2], " (positive)", "\n")

  # crosstab numbers in table()
  tab1 <- table(var.thres,reference)
  clvl <- qnorm((1+ci.level)/2)
  # LR+
  lr1 <- s1/(1-s2)
  # SE of log(LR+)
  se1 <- sqrt(1/tab1[2,2]-1/sum(tab1[,2])+1/tab1[2,1]-1/sum(tab1[,1]))
  # confidence intervals of LR+
  lrci1 <- c(exp(log(lr1)-clvl*se1),exp(log(lr1)+clvl*se1))

  # LR-
  lr2 <- (1-s1)/s2
  # SE of log(LR-)
  se2 <- sqrt(1/tab1[1,2]-1/sum(tab1[,2])+1/tab1[1,1]-1/sum(tab1[,1]))
  # confidence intervals of LR-
  lrci2 <- c(exp(log(lr2)-clvl*se2),exp(log(lr2)+clvl*se2))

  # the estimates
  # LR+
  lrp<-c(round(lr1,3),round(lrci1,3))
  # LR-
  lrn<-c(round(lr2,3),round(lrci2,3))

  # output values
  if (ci==TRUE){
    c1<-c(paste0(lrp[1]," (",lrp[2],", ",lrp[3],")"),
          paste0(lrn[1]," (",lrn[2],", ",lrn[3],")"))
  } else {
    c1<-c(lrp[1],lrn[1])
  }
  names(c1) <- c("DLR+","DLR-")

  return(data.frame(Estimate = c1))
}

