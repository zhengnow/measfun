
#' Table the measures of accuracy
#'
#' The functions estimate the measures of accuracy (without adjustment)
#' for a binary reference variable and a test variable.
#' The test variable, when applicable, are dichotomized at the threshold value.
#' When the data is clustered, the survey package is used to estimate the confidence interval.
#'
#' @param var a vector, test variable
#' @param reference a binary reference variable
#' @param thres a numeric value, the threshold value when the test variable is continuous
#' @param ci TRUE/FALSE, whether to include the confidence interval in the output
#' @param cluster TRUE/FALSE, an indicator of whether the data is clustered. The data was assume to adopt simple random sampling with equal weights
#' @param id subject ID variable, when data is clustered
#'
#'@details
#'tbaccu.sspv outputs the table including estimates of sensitivity, specificity, PPV and NPV.
#'tbaccu.sslr outputs the table including estimates of sensitivity, specificity, LR+ and LR-.
#'Clustered data is not accepted in this function yet.
#'
#' @return a data frame for the estimates, with the 95%CI (optional), and the counts
#'
#' @examples
#' y <- rbinom(n=50, size=1, prob=0.3) # reference variable
#' x <- rnorm(n=50, mean=20, sd=5)
#' tbaccu.sspv(x,y, thres = 20)
#' tbaccu.sslr(x,y, thres = 20)
#'
#' @export
#'




# --------------------------------------------
# report sensitivity, specificity, PPV and NPV
# --------------------------------------------
tbaccu.sspv <-
  function(var, reference, thres=NULL, ci=TRUE, cluster=FALSE, id = NULL){

  # if the reference has more than 2 levels, quit and print error message
  assert_stopifnot(length(table(reference))==2)

  # show the levels of reference and test variables
  cat("The reference levels: ",
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

  cat("The test levels: ",
        names(table(var.thres))[1], " (negative), ",
        names(table(var.thres))[2], " (positive)", "\n")

  # crosstab numbers in table()
  tab1 <- table(var.thres,reference)

  if(!cluster){
    # sensitivity
    se<-c(round(tab1[2,2]/sum(tab1[,2]),3),
          round(binom.test(tab1[2,2],sum(tab1[,2]))$conf.int,3))
    # specificity
    sp<-c(round(tab1[1,1]/sum(tab1[,1]),3),
          round(binom.test(tab1[1,1],sum(tab1[,1]))$conf.int,3))
    # PPV
    ppv<-c(round(tab1[2,2]/sum(tab1[2,]),3),
           round(binom.test(tab1[2,2],sum(tab1[2,]))$conf.int,3))
    # npv
    npv<-c(round(tab1[1,1]/sum(tab1[1,]),3),
           round(binom.test(tab1[1,1],sum(tab1[1,]))$conf.int,3))
  }
  if(cluster) {

    binvar <- ifelse(var>ct, 1,0)
    negbinvar <- ifelse(var<=ct, 1, 0)
    vdata <- data.frame(binvar, reference, negbinvar, id)
    vclus1 <- survey::svydesign(id = ~id, data=vdata)

    # se and sp
    svyse <-
      survey::svyby(~binvar, by = ~reference,
                    design=vclus1, svymean)
    #svysp <- survey::svyby(~negbinvar, by = ~reference, design=vclus1, svymean)
    # confidence interval of reference=0 group = 1- confint(svyse)

    # sensitivity
    se<-c(round(tab1[2,2]/sum(tab1[,2]),3),
          round(capci(confint(svyse)["1",]),3))
    # specificity
    sp<-c(round(tab1[1,1]/sum(tab1[,1]),3),
          round(capci(1-confint(svyse)["0",2:1]),3))

    # ppv and npv svy
    svyppv <-
      survey::svyby(~reference, by = ~binvar,
                    design=vclus1, svymean)
    # PPV
    ppv<-c(round(tab1[2,2]/sum(tab1[2,]),3),
           round(capci(confint(svyppv)["1",]),3))
    # npv
    npv<-c(round(tab1[1,1]/sum(tab1[1,]),3),
           round(capci(1-confint(svyppv)["0",2:1]),3))
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
  return(data.frame(Estimate = c1))
}

# --------------------------
# report SE, SP, LR+ and LR-
# --------------------------
tbaccu.sslr <- function(var, reference, thres=NULL, ci=TRUE){

  # if the reference has more than 2 levels, quit and print error message
  assert_stopifnot(length(table(reference))==2)

  # show the reference and test variable levels
  cat("The reference levels: ",
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

  cat("The test levels: ",
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

  # add numbers
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
