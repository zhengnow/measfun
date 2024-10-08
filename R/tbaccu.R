
#' Table the measures of accuracy
#'
#' The function estimate the measures of accuracy (without adjustment)
#' for a test variable and a binary reference variable, independent observations or clustered data.
#'
#' @param var a vector, test variable
#' @param reference a binary reference variable
#' @param thres a numeric value, the threshold value when the test variable is continuous
#' @param ci TRUE/FALSE, whether to include the confidence interval in the output
#' @param ci.level the confidence level for the confidence interval calculation, default is 0.95
#' @param cluster TRUE/FALSE, an indicator of whether the data is clustered. The data was assume to adopt simple random sampling with equal weights
#' @param id subject cluster id, when data is clustered
#'
#'@details
#'tbaccu outputs the table including estimates of sensitivity, specificity, PPV and NPV.
#' The test variable, when applicable, are dichotomized at the threshold value.
#' When the data is clustered, the survey package is used to estimate the confidence intervals,
#' assuming simple random sampling.
#'
#' @return a data frame on measure of accuracy estimates, with the 95%CI (optional), and the counts.
#'
#' @examples
#' y <- rbinom(n=50, size=1, prob=0.3) # reference variable
#' x <- rnorm(n=50, mean=20, sd=5) # test variable
#' tbaccu(x, y, thres = 20)
#'
#' @export
#'




# --------------------------------------------
# report sensitivity, specificity, PPV and NPV
# --------------------------------------------
tbaccu <-
  function(var, reference, thres=NULL, ci=TRUE, ci.level = 0.95, cluster=FALSE, id = NULL){

  # if the reference has more than 2 levels, quit and print error message
  assert_stopifnot(length(table(reference))==2)

  # show the levels of reference and test variables
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

  if(!cluster){
    # sensitivity
    se<-c(round(tab1[2,2]/sum(tab1[,2]),3),
          round(binom.test(tab1[2,2],sum(tab1[,2]), conf.level = ci.level)$conf.int,3))
    # specificity
    sp<-c(round(tab1[1,1]/sum(tab1[,1]),3),
          round(binom.test(tab1[1,1],sum(tab1[,1]), conf.level = ci.level)$conf.int,3))
    # PPV
    ppv<-c(round(tab1[2,2]/sum(tab1[2,]),3),
           round(binom.test(tab1[2,2],sum(tab1[2,]), conf.level = ci.level)$conf.int,3))
    # npv
    npv<-c(round(tab1[1,1]/sum(tab1[1,]),3),
           round(binom.test(tab1[1,1],sum(tab1[1,]), conf.level = ci.level)$conf.int,3))
  }
  if(cluster) {

    binvar <- ifelse(var>ct, 1,0)
    svydata <- data.frame(binvar, reference, id)
    vclus1 <- survey::svydesign(id = ~id, data=svydata)

    # se and sp
    svy1 <-
      survey::svyby(~binvar, by = ~reference,
                    design=vclus1, survey::svymean)

    # sensitivity
    se<-c(round(tab1[2,2]/sum(tab1[,2]),3),
          round(capci(confint(svy1, level = ci.level)["1",]),3))
    # specificity
    sp<-c(round(tab1[1,1]/sum(tab1[,1]),3),
          round(capci(1-confint(svy1, level = ci.level)["0",2:1]),3))

    # ppv and npv
    svy2 <-
      survey::svyby(~reference, by = ~binvar,
                    design=vclus1, survey::svymean)
    # PPV
    ppv<-c(round(tab1[2,2]/sum(tab1[2,]),3),
           round(capci(confint(svy2, level = ci.level)["1",]),3))
    # npv
    npv<-c(round(tab1[1,1]/sum(tab1[1,]),3),
           round(capci(1-confint(svy2, level = ci.level)["0",2:1]),3))
  }

  # format the estimates (95%CI) and add numbers
  if (ci==TRUE){
    c1=c(paste0(se[1]," (",se[2],", ",se[3],"), ",
                tab1[2,2],"/",sum(tab1[,2])),
         paste0(sp[1]," (",sp[2],", ",sp[3],"), ",
                tab1[1,1],"/",sum(tab1[,1])),
         paste0(ppv[1]," (",ppv[2],", ",ppv[3],"), ",
                tab1[2,2],"/",sum(tab1[2,])),
         paste0(npv[1]," (",npv[2],", ",npv[3],"), ",
                tab1[1,1],"/",sum(tab1[1,])))
  } else {
    c1=c(paste0(se[1],", ",tab1[2,2],"/",sum(tab1[,2])),
         paste0(sp[1],", ",tab1[1,1],"/",sum(tab1[,1])),
         paste0(ppv[1],", ",tab1[2,2],"/",sum(tab1[2,])),
         paste0(npv[1],",",tab1[1,1],"/",sum(tab1[1,])))
  }
  names(c1) <-
    c("Sensitivity","Specificity","PPV","NPV")

  out <- data.frame(Estimate = c1)
  return()
}
