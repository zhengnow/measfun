
#' Table the AUC, threshold, and measures of accuracy at the threshold
#'
#' This function estimate the AUC for a binary reference variable and a continuous predictor,
#' reporting the threshold value maximizing the Youden's index,
#' and the corresponding measures of accuracy (without adjustment).
#'
#' @param var a continuous predictor
#' @param reference a binary reference variable
#' @return a vector
#' @export
#'
# -----------------------------------------
# ROC analysis
# Estimate se and sp max Youden's index
# -----------------------------------------

youdentb <- function(var,reference){
  # if the reference had more than 2 levels, quit and print error message
  assert(length(table(reference))==2)
  # print the
  # var and reference variables need to have matched levels (at least in the same order)
  print(paste("The AUC is for the level of", names(table(reference))[2]))

  reference <- as.numeric(as.factor(reference))

  r1 <- pROC::roc(reference,var)
  sum1 <- r1$sensitivities+r1$specificities
  r1$sensitivities[which(sum1==max(sum1))]
  r1$specificities[which(sum1==max(sum1))]
  # thresholds
  (t1 <- r1$thresholds[which(sum1==max(sum1))])
  tab1 <- table(var>t1, reference)

  # sensitivity
  se <- c(round(tab1[2,2]/sum(tab1[,2]),3),
          round(binom.test(tab1[2,2], sum(tab1[,2]))$conf.int,3))
  # specificity
  sp <- c(round(tab1[1,1]/sum(tab1[,1]),3),
          round(binom.test(tab1[1,1], sum(tab1[,1]))$conf.int,3))
  # PPV
  ppv <- c(round(tab1[2,2]/sum(tab1[2,]),3),
           round(binom.test(tab1[2,2], sum(tab1[2,]))$conf.int,3))
  # npv
  npv <- c(round(tab1[1,1]/sum(tab1[1,]),3),
           round(binom.test(tab1[1,1], sum(tab1[1,]))$conf.int,3))

  # add numbers
  c1 <- c(paste0(round(r1$auc,3)," (",
              round(ci.auc(r1)[1],3),", ",
              round(ci.auc(r1)[3],3),")"),
       t1,
       paste0(se[1]," (",se[2],", ",se[3],"), ",
              tab1[2,2],"/",sum(tab1[,2])),
       paste0(sp[1]," (",sp[2],", ",sp[3],"), ",
              tab1[1,1],"/",sum(tab1[,1])),
       paste0(ppv[1]," (",ppv[2],", ",ppv[3],"), ",
              tab1[2,2],"/",sum(tab1[2,])),
       paste0(npv[1]," (",npv[2],", ",npv[3],"), ",
              tab1[1,1],"/",sum(tab1[1,])))
  names(c1) <-
    c("AUC","cutoff","Sensitivity","Specificity","PPV","NPV")
  return(c1)
}
