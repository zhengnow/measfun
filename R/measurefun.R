




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
