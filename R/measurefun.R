





# =========================================
# A function to output measures of accuracy
# =========================================
# var = test result
# reference = reference standard
# ct = cutoff value. For binary, when 1 is positive, set as 0.
# ci = output confidence interval, set as T/F.

# for clustered set up, use survey package to calculate confidence intervals

# a function to cap the 95%CI of measures of accuracy (0-1)
capci <- function(x){
  # lower bound
  x[1] = max(x[1],0)
  # upper bound
  x[2] = min(x[2],1)
  return(x)
}

assert <-
  function(exprs) eval.parent(substitute(stopifnot(exprs = exprs)))

tb.sesp <- function(var,reference,ct=0,ci=TRUE, cluster=F, id = NULL){

  # if the reference had more than 2 levels, quit and print error message
  assert(length(table(reference))==2, length(table(var))==2)

  # print the
  # var and reference variables need to have matched levels (at least in the same order)
  print(paste("The positive level in the reference is ", names(table(reference))[2]))
  print(paste("The test levels are"), names(table(var))[1], "and", names(table(var))[2])


  tab1 <- table(var>ct,reference)
  # dim.tab1 <- dim(tab1)
  # if(sum(dim.tab1) != 4){
  #   if (dim.tab1[1] == 1){
  #     present.grp <- rownames(tab1)
  #
  #   } else if (dim.tab1[2]==1)
  # }
  if(cluster == FALSE){
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
  if(cluster == TRUE) {

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

  # add numbers
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
  names(c1)=c("Sensitivity","Specificity","PPV","NPV")
  return(c1)
}

# --------------------------
# report SE, SP, LR+ and LR-
LR.f = function(var,reference,ct=0,ci=TRUE){
  tab1=table(var>ct,reference)
  s1 = tab1[2,2]/sum(tab1[,2]) #sen
  s2 = tab1[1,1]/sum(tab1[,1]) #spe
  # LR+
  lr1 = s1/(1-s2)
  # SE of log(LR+)
  se1 = sqrt(1/tab1[2,2]-1/sum(tab1[,2])+1/tab1[2,1]-1/sum(tab1[,1]))
  # confidence intervals of LR+
  lrci1 = c(exp(log(lr1)-1.96*se1),exp(log(lr1)+1.96*se1))
  # LR-
  lr2 = (1-s1)/s2
  # SE of log(LR-)
  se2 = sqrt(1/tab1[1,2]-1/sum(tab1[,2])+1/tab1[1,1]-1/sum(tab1[,1]))
  # confidence intervals of LR-
  lrci2 = c(exp(log(lr2)-1.96*se2),exp(log(lr2)+1.96*se2))

  # output values
  # sensitivity
  se=c(round(tab1[2,2]/sum(tab1[,2]),3),round(binom.test(tab1[2,2],sum(tab1[,2]))$conf.int,3))
  # specificity
  sp=c(round(tab1[1,1]/sum(tab1[,1]),3),round(binom.test(tab1[1,1],sum(tab1[,1]))$conf.int,3))
  # LR+
  lrp=c(round(lr1,3),round(lrci1,3))
  # LR-
  lrn=c(round(lr2,3),round(lrci2,3))

  # add numbers
  if (ci==TRUE){
    c1=c(paste(se[1]," (",se[2],", ",se[3],"), ",tab1[2,2],"/",sum(tab1[,2]),sep=""),
         paste(sp[1]," (",sp[2],", ",sp[3],"), ",tab1[1,1],"/",sum(tab1[,1]),sep=""),
         paste(lrp[1]," (",lrp[2],", ",lrp[3],")",sep=""),
         paste(lrn[1]," (",lrn[2],", ",lrn[3],")",sep=""))
  } else {
    c1=c(paste(se[1],", ",tab1[2,2],"/",sum(tab1[,2]),sep=""),
         paste(sp[1],", ",tab1[1,1],"/",sum(tab1[,1]),sep=""),
         paste(lrp[1],sep=""),
         paste(lrn[1],sep=""))
  }
  names(c1)=c("Sensitivity","Specificity","DLR+","DLR-")
  return(c1)
}
# -----------------------------
# kappa and 95%CI
# -----------------------------
# agreement on categorical vars
# function to bootstrap for confidence interval
kappa2.boot <-
  function(data, x, weight2 = NULL){
    kappa2(data[x,],
           weight = weight2)$value
  }

# function output kappa and bootstrapped confidence intervals
kappa.f <-
  function(v1, v2, weight = "equal"){
    set.seed(338)
    k2 <- round(
      c(kappa2(cbind(v1,v2),
               weight = weight)$value,
        boot.ci(boot(cbind(v1,v2),
                     kappa2.boot, 1000, weight2 = weight),
                type="bca")$bca[,4:5]),
      3)
    kp <- paste(k2[1], "(95%CI: ",
                k2[2],", ",
                k2[3],")",
                sep="")
    names(kp) <- "kappa"
    return(kp)
  }

# -----------------------------
# kappa and 95%CI
# boot at cluster level
# -----------------------------
# agreement on categorical vars
# function to bootstrap for confidence interval
kappa2.clu.boot <-
  function(data, x, datafull, weight2 = NULL){
    # id indicator is listed in the 3rd column of the data
    # When a cluster id is sampled, include all observations of the cluster id.
    # data1 <- data[data[,3] %in% data[x,3], ]
    data1 <- datafull[datafull[,3] %in% x, ]
    kappa2(data1[,1:2],
           weight = weight2)$value
  }

# function output kappa and bootstrapped confidence intervals
kappa.clu.f <-
  function(v1, v2, clu.id = id, weight = "equal"){
    set.seed(338)
    k2 <- round(
      c(kappa2(cbind(v1,v2),
               weight = weight)$value,
        boot.ci(boot(unique(clu.id),
                     kappa2.clu.boot, 1000,
                     datafull=cbind(v1,v2,clu.id),
                     weight2 = weight),
                type="bca")$bca[,4:5]),
      3)
    ckappa <- paste(k2[1], " (95%CI: ",
                    k2[2],", ",
                    k2[3],")",
                    sep="")
    names(ckappa) <- "Clustered kappa"
    return(ckappa)
  }

# --------------------------
# Light's kappa, add weights
# --------------------------
# under construction
kappam2.light <-
  function(ratings, weight=c("unweighted", "equal", "squared")) {
    ratings <- as.matrix(na.omit(ratings))
    if (is.character(weight))
      weight = match.arg(weight)
    ns <- nrow(ratings)
    nr <- ncol(ratings)
    for (i in 1:(nr - 1)) for (j in (i + 1):nr) {
      if ((i == 1) & (j == (i + 1)))
        kappas <- kappa2(ratings[, c(i, j)], weight = weight)$value
      else kappas <- c(kappas, kappa2(ratings[, c(i, j)], weight = "u")$value)
    }
    value <- mean(kappas)
    lev <- levels(as.factor(ratings))
    levlen <- length(levels(as.factor(ratings)))
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
    u <- value/SEkappa
    p.value <- 2 * (1 - pnorm(abs(u)))
    rval <- structure(list(method = "Light's Kappa for m Raters",
                           subjects = ns, raters = nr, irr.name = "Kappa", value = value,
                           stat.name = "z", statistic = u, p.value = p.value), class = "irrlist")
    return(rval)
  }

# -----------------------------------
# concordance correlation coefficient
# -----------------------------------
# concordence correlation coefficient
# 09-09 from Chaya
ccc <-
  function(var1, var2)  {
    mat <- cbind(var1, var2)
    cov <- var(mat)[1, 2]
    s1 <- var(mat[,1])
    s2 <- var(mat[,2])
    m1 <- mean(mat[,1])
    m2 <- mean(mat[,2])
    corr <- cor(mat)[1, 2]
    n <- dim(na.omit(mat))[1]
    ccc <- (2 * cov)/(s1 + s2 + (m1 - m2)^2)
    omcccsq <- 1 - ccc^2
    omcorrsq <- 1 - corr^2
    u <- (m1 - m2)/sqrt(sqrt(s1) * sqrt(s2))
    sigmasq <- (((omcorrsq * (ccc^2))/(omcccsq * (corr^2))
    ) + ((4 * (ccc^3) * (1 - ccc) * (u^2))/(corr * (
      omcccsq^2))) - ((2 * (ccc^4) * (u^4))/((corr^2
      ) * (omcccsq^2))))/(n - 2)
    std <- sqrt((omcccsq^2) * sigmasq)
    cilow <- max(ccc - 1.96 * std,-1)
    cihigh <- min(ccc + 1.96 * std,1)
    p <- round((1 - pnorm(abs(ccc/std))) * 2, 4)
    cccout <- round(c(ccc, cilow, cihigh, p),3)
    names(cccout) <- c("CCC","ci95lower","ci95upper","p")
    return(cccout)
  }

# example
# set.seed(567); v1 <- rnorm(50);
# set.seed(448); v2 <- rnorm(50)
# ccc(v1, v2)

# overall concordence correlation coefficient
# for multiple readers (>2)
# 08-23-2010
occc <-
  function(mat){
    ncmat=ncol(mat)
    summ=sumc=sums=0
    for (j in 1:(ncmat-1)) {
      for (k in (j+1):ncmat) {
        mjk = (mean(mat[,j])-mean(mat[,k]))^2;
        summ=summ+mjk
        cov.jk = var(mat[,c(j,k)])[1, 2]
        sumc=sumc+cov.jk
      }}
    for (j in 1:ncmat) {
      sj= var(mat[,j]) ; sums=sums+sj
    }
    ccc <- (2 * sumc)/((ncmat-1)*sums + summ)
    ccc <- round(ccc,3)
    names(ccc) = "Overall CCC"
    return(ccc)
  }

# ----------------------------------------
# Compare AUCs in clustered data
# ----------------------------------------
# https://www.lerner.ccf.org/qhs/software/lib/funcs_clusteredROC.R
# The R script funcs_clusteredROC.R contains functions to perform the statistical methods in:
#   Obuchowski NA. Nonparametric analysis of clustered ROC curve data. Biometrics. 1997: 567-578.
# =================================
#  kern
# =================================
kern = function(x1, x0)
{
  if(x1  > x0)      return(1.0)
  else if(x1 == x0) return(0.5)
  else if(x1  < x0) return(0.0)
}
# =================================
#  V10
# =================================
V10 = function(Xi, Ys) { sum(sapply(Ys, FUN=kern, x1=Xi)) / length(Ys) }
# =================================
#  V01
# =================================
V01 = function(Xs, Yi) { sum(sapply(Xs, FUN=kern, x0=Yi)) / length(Xs) }
# =================================
#  getAUC
# =================================
getAUC = function(Xs, Ys) { sum(sapply(Xs, function(Xi) { sapply(Ys, FUN=kern, x1=Xi) })) / (length(Xs)*length(Ys)) }

# =================================
#  clusteredROC
# =================================
clusteredROC = function(predictor1, predictor2=NULL, response, clusterID, alpha=0.05, level=NULL, print.all=F)
{
  if(is.null(level))
  {
    absent  = levels(as.factor(response))[1]
    present = levels(as.factor(response))[2]
  } else
  {
    absent  = levels(as.factor(response))[levels(as.factor(response)) != level]
    present = level
  }

  # *** One ROC curve ***
  if(is.null(predictor2)==T)
  {
    missings = which(is.na(predictor1)==T)
    if(length(missings)>0)
    {
      predictor1 = predictor1[-missings]
      response   = response[-missings]
      clusterID  = clusterID[-missings]
    }

    I   = length(unique(clusterID))
    I01 = length(unique(clusterID[response==absent]))
    I10 = length(unique(clusterID[response==present]))
    m   = sapply(1:I, function(i) sum(clusterID==unique(clusterID)[i] & response==present))
    n   = sapply(1:I, function(i) sum(clusterID==unique(clusterID)[i] & response==absent))
    s   = m + n
    N   = sum(n)
    M   = sum(m)

    AUC1 = getAUC(predictor1[response==present], predictor1[response==absent])
    Xcomps.1 = rep(NA, I)
    Ycomps.1 = rep(NA, I)

    for(i in 1:I)
    {
      if(m[i]==0) { Xcomps.1[i] = 0 }  else
      {
        Xcomps.1[i] = sum(sapply(1:m[i], function(j) V10(Xi=predictor1[clusterID==unique(clusterID)[i] & response==present][j], Ys=predictor1[response==absent])))
      }
      if(n[i]==0) { Ycomps.1[i] = 0  }  else
      {
        Ycomps.1[i] = sum(sapply(1:n[i], function(j) V01(Yi=predictor1[clusterID==unique(clusterID)[i] & response==absent][j], Xs=predictor1[response==present])))
      }
    }

    S10_1 = (I10/((I10-1)*M)) * sum((Xcomps.1 - m*AUC1) * (Xcomps.1 - m*AUC1))
    S01_1 = (I01/((I01-1)*N)) * sum((Ycomps.1 - n*AUC1) * (Ycomps.1 - n*AUC1))
    S11_1 = (I/(I-1))         * sum((Xcomps.1 - m*AUC1) * (Ycomps.1 - n*AUC1))
    var_1 = S10_1/M + S01_1/N + (2*S11_1)/(M*N)
    AUC1.SE = sqrt(var_1)

    AUC1.CIlo = AUC1 - qnorm(1-alpha/2)*AUC1.SE
    AUC1.CIhi = AUC1 + qnorm(1-alpha/2)*AUC1.SE

    cat ("\n")
    cat ("Total # of clusters: ", I, "\n", sep='')
    cat ("Total # of observations: ", length(clusterID), "\n", sep='')
    cat ("Min # of observations per cluster: ", min(s), "\n", sep='')
    cat ("Max # of observations per cluster: ", max(s), "\n", sep='')
    cat ("AUC (SE) for ROC curve: ", round(AUC1,4), " (", round(AUC1.SE,4), ")\n", sep='')
    cat ((100*(1-alpha)),"% CI for AUC: ", "(", round(AUC1.CIlo,4), ", ", round(AUC1.CIhi,4), ")\n\n", sep='')

    if(print.all)
    {
      name = c("I", "I10", "I01", "M", "N", "S10", "S01", "S11")
      value = c(I, I10, I01, M, N, S10_1, S01_1, S11_1)
      print(data.frame(name, value))
    }
  }

  # *** Two correlated ROC curves ***
  if(is.null(predictor2)==F)
  {
    missings = which(is.na(predictor1)==T | is.na(predictor2)==T)
    if(length(missings)>0)
    {
      predictor1 = predictor1[-missings]
      predictor2 = predictor2[-missings]
      response   = response[-missings]
      clusterID  = clusterID[-missings]
    }

    I   = length(unique(clusterID))
    I01 = length(unique(clusterID[response==absent]))
    I10 = length(unique(clusterID[response==present]))
    m   = sapply(1:I, function(i) sum(clusterID==unique(clusterID)[i] & response==present))
    n   = sapply(1:I, function(i) sum(clusterID==unique(clusterID)[i] & response==absent))
    s   = m + n
    N   = sum(n)
    M   = sum(m)

    AUC1 = getAUC(predictor1[response==present], predictor1[response==absent])
    AUC2 = getAUC(predictor2[response==present], predictor2[response==absent])
    Xcomps.1 = rep(NA, I)
    Xcomps.2 = rep(NA, I)
    Ycomps.1 = rep(NA, I)
    Ycomps.2 = rep(NA, I)

    for(i in 1:I)
    {
      if(m[i]==0) { Xcomps.1[i] = 0;  Xcomps.2[i] = 0 }  else
      {
        Xcomps.1[i] = sum(sapply(1:m[i], function(j) V10(Xi=predictor1[clusterID==unique(clusterID)[i] & response==present][j], Ys=predictor1[response==absent])))
        Xcomps.2[i] = sum(sapply(1:m[i], function(j) V10(Xi=predictor2[clusterID==unique(clusterID)[i] & response==present][j], Ys=predictor2[response==absent])))
      }
      if(n[i]==0) { Ycomps.1[i] = 0;  Ycomps.2[i] = 0  }  else
      {
        Ycomps.1[i] = sum(sapply(1:n[i], function(j) V01(Yi=predictor1[clusterID==unique(clusterID)[i] & response==absent][j], Xs=predictor1[response==present])))
        Ycomps.2[i] = sum(sapply(1:n[i], function(j) V01(Yi=predictor2[clusterID==unique(clusterID)[i] & response==absent][j], Xs=predictor2[response==present])))
      }
    }

    S10_1 = (I10/((I10-1)*M)) * sum((Xcomps.1 - m*AUC1) * (Xcomps.1 - m*AUC1))
    S01_1 = (I01/((I01-1)*N)) * sum((Ycomps.1 - n*AUC1) * (Ycomps.1 - n*AUC1))
    S11_1 = (I/(I-1))         * sum((Xcomps.1 - m*AUC1) * (Ycomps.1 - n*AUC1))
    var_1 = S10_1/M + S01_1/N + (2*S11_1)/(M*N)

    S10_2  = (I10/((I10-1)*M)) * sum((Xcomps.2 - m*AUC2) * (Xcomps.2 - m*AUC2))
    S01_2  = (I01/((I01-1)*N)) * sum((Ycomps.2 - n*AUC2) * (Ycomps.2 - n*AUC2))
    S11_2  = (I/(I-1))         * sum((Xcomps.2 - m*AUC2) * (Ycomps.2 - n*AUC2))
    var_2  = S10_2/M + S01_2/N + (2*S11_2)/(M*N)

    S10_12 = (I10/((I10-1)*M)) * sum((Xcomps.1 - m*AUC1) * (Xcomps.2 - m*AUC2))
    S01_12 = (I01/((I01-1)*N)) * sum((Ycomps.1 - n*AUC1) * (Ycomps.2 - n*AUC2))
    S11_12 = (I/(I-1))         * sum((Xcomps.1 - m*AUC1) * (Ycomps.2 - n*AUC2))
    S11_21 = (I/(I-1))         * sum((Xcomps.2 - m*AUC2) * (Ycomps.1 - n*AUC1))
    cov_12 = S10_12/M + S01_12/N + S11_12/(M*N) + S11_21/(M*N)

    AUC1.SE = sqrt(var_1)
    AUC2.SE = sqrt(var_2)

    DIFF      = abs(AUC1 - AUC2)
    DIFF.SE   = sqrt(var_1 + var_2 - 2*cov_12)
    DIFF.CIlo = DIFF - qnorm(1-alpha/2)*DIFF.SE
    DIFF.CIhi = DIFF + qnorm(1-alpha/2)*DIFF.SE
    p = 2*(1-pnorm(DIFF/DIFF.SE))

    cat ("\n")
    cat ("Total # of clusters: ", I, "\n", sep='')
    cat ("Total # of observations: ", length(clusterID), "\n", sep='')
    cat ("Min # of observations per cluster: ", min(s), "\n", sep='')
    cat ("Max # of observations per cluster: ", max(s), "\n", sep='')
    cat ("AUC (SE) for ROC curve 1: ", round(AUC1,4), " (", round(AUC1.SE,4), ")\n", sep='')
    cat ("AUC (SE) for ROC curve 2: ", round(AUC2,4), " (", round(AUC2.SE,4), ")\n", sep='')
    cat ("Difference (SE): ", round(DIFF,4), " (", round(DIFF.SE,4), ")\n", sep='')
    cat ((100*(1-alpha)),"% CI for difference: ", "(", round(DIFF.CIlo,4), ", ", round(DIFF.CIhi,4), ")\n", sep='')
    cat ("Associated p-value: ", format(p,digits=4), "\n\n", sep='')

    if(print.all)
    {
      name = c("I", "I10", "I01", "M", "N", "reader 1 S10", "reader 1 S01", "reader 1 S11",
               "reader 2 S10", "reader 2 S01", "reader 2 S11", "S10_12", "S01_12", "S11_12", "S11_21")
      value = c(I, I10, I01, M, N, S10_1, S01_1, S11_1, S10_2, S01_2, S11_2, S10_12, S01_12, S11_12, S11_21)
      print(data.frame(name, value))
    }
  }
  if(is.null(predictor2)==T)
  {
    invisible(list("auc"=AUC1, "auc.se"=AUC1.SE, "ci.for.auc"=c(AUC1.CIlo, AUC1.CIhi)))
  } else
  {
    invisible(list("auc"=c(AUC1, AUC2), "auc.se"=c(AUC1.SE, AUC2.SE), "diff"=DIFF, "diff.se"=DIFF.SE,
                   "ci.for.diff"=c(DIFF.CIlo, DIFF.CIhi), "p.value"=p))
  }
}
