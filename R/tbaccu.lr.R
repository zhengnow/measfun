
#' Table the measures of accuracy
#'
#' The function estimates the measures of accuracy (without adjustment)
#' for a binary reference variable and a test variable.
#' The test variable, when applicable, is dichotomized at the threshold value.
#'
#' @param var a vector, test variable
#' @param reference a binary reference variable
#' @param thres a numeric value, the threshold value when the test variable is continuous
#' @param ci TRUE/FALSE, whether to include the confidence interval in the output
#'
#' @details
#' tbaccu.sslr outputs the table including estimates of sensitivity, specificity, LR+ and LR-.
#' Clustered data is not accepted in this function yet.
#'
#' @return a data frame for the estimates, with the 95%CI (optional), and the counts
#'
#' @examples
#' y <- rbinom(n=50, size=1, prob=0.3) # reference variable
#' x <- rnorm(n=50, mean=20, sd=5)
#' tbaccu.lr(x,y, thres = 20)
#'
#' @export
#'


# --------------------------
# report SE, SP, LR+ and LR-
# --------------------------
tbaccu.lr <- function(var, reference, thres=NULL, ci=TRUE){

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

  s1 <- tab1[2,2]/sum(tab1[,2]) #sen
  s2 <- tab1[1,1]/sum(tab1[,1]) #spe
  # LR+
  lr1 <- s1/(1-s2)
  # SE of log(LR+)
  se1 <- sqrt(1/tab1[2,2]-1/sum(tab1[,2])+1/tab1[2,1]-1/sum(tab1[,1]))
  # confidence intervals of LR+
  lrci1 <- c(exp(log(lr1)-1.96*se1),exp(log(lr1)+1.96*se1))
  # LR-
  lr2 <- (1-s1)/s2
  # SE of log(LR-)
  se2 <- sqrt(1/tab1[1,2]-1/sum(tab1[,2])+1/tab1[1,1]-1/sum(tab1[,1]))
  # confidence intervals of LR-
  lrci2 <- c(exp(log(lr2)-1.96*se2),exp(log(lr2)+1.96*se2))

  # output values
  # sensitivity
  se<-c(round(tab1[2,2]/sum(tab1[,2]),3),round(binom.test(tab1[2,2],sum(tab1[,2]))$conf.int,3))
  # specificity
  sp<-c(round(tab1[1,1]/sum(tab1[,1]),3),round(binom.test(tab1[1,1],sum(tab1[,1]))$conf.int,3))
  # LR+
  lrp<-c(round(lr1,3),round(lrci1,3))
  # LR-
  lrn<-c(round(lr2,3),round(lrci2,3))

  # add numbers (counts)
  if (ci==TRUE){
    c1<-c(paste0(se[1]," (",se[2],", ",se[3],"), ",tab1[2,2],"/",sum(tab1[,2])),
         paste0(sp[1]," (",sp[2],", ",sp[3],"), ",tab1[1,1],"/",sum(tab1[,1])),
         paste0(lrp[1]," (",lrp[2],", ",lrp[3],")"),
         paste0(lrn[1]," (",lrn[2],", ",lrn[3],")"))
  } else {
    c1<-c(paste0(se[1],", ",tab1[2,2],"/",sum(tab1[,2])),
         paste0(sp[1],", ",tab1[1,1],"/",sum(tab1[,1])),
         paste0(lrp[1]),
         paste0(lrn[1]))
  }
  names(c1) <- c("Sensitivity","Specificity","DLR+","DLR-")

  return(data.frame(Estimate = c1))
}

# LR estimates in clustered data are pending
