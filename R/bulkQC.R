##Dependencies
#install.packages('isotree')
#library(isotree)
#library(stddiff)

##Count number of unique values in a vector (to focus on cont variables)
nVal = function(x) {
 length(unique(x))
}

##Obtain IQR-based outliers in a vector
getOutliers = function(x, mult=1.5) {
 x = as.numeric(x)
 ss = summary(x)
 q1 = ss["1st Qu."]
 q3 = ss["3rd Qu."]
 iqr = q3-q1
 l = q1-mult*iqr
 u = q3+mult*iqr
 x[x<l | x>u]
}

##Function to identify individual univariate outliers in a data frame
ind_uni = function(d0, exclude=c("pid", "site"), n_uniq=10, m=1.5) {
 d1 = d0[, !names(d0) %in% exclude]
 if (length(dim(d1))==0) {
  d1 = as.data.frame(as.matrix(d1, 1, length(d1)))
  names(d1) = names(d0)[!names(d0) %in% exclude]
  }

 tt = apply(d1, 2, nVal) #only keep those with sufficient variation
 keep = names(tt)[tt>=n_uniq]
 d2 = d1[, names(d1) %in% keep]
 if (length(dim(d2))==0) {
  d2 = as.data.frame(as.matrix(d2, 1, length(d2)))
  names(d2) = keep
 }
 out = apply(d2, 2, function(x) getOutliers(x, mult=m)) #look for outliers
 if (typeof(out)=="double") {out=data.frame(out)}
 n_out = sum(unlist(lapply(out, function(x) sum(!is.na(x)))))

 if (n_out==0) {
  data = "No outliers found"
  } else if (n_out>0) {
   d3 = d0
   d3[, !names(d3) %in% exclude] = NA
   for (i in 1:length(names(out))) {
    vv = d0[, names(d0)==names(out)[i]]
    vv[ !vv %in% out[[i]]] = NA
    d3[, names(d3)==names(out)[i]] = vv
    }
   nvar = length(names(d0)[!names(d0) %in% exclude])
   d4 = d3[apply(d3, 1, function(x) sum(is.na(x)))!=nvar, ]
   d5 = d4[, apply(d4, 2, function(x) sum(is.na(x)))!=dim(d4)[1]]
   data = d5[, match(c(names(d5)[names(d5) %in% exclude], names(d5)[!names(d5) %in% exclude]), names(d5))]
   }
 nID = dim(d2)[1]
 nVar = dim(d2)[2]
 output.all = list(nID=nID, nVar=nVar, data=data)
 return(output.all)
}

##Test it out
#iris2 = iris
#iris2$pid = 1:dim(iris2)[1]
#ind_uni(iris2, exclude=c("pid", "Species"), m=1.5)
#ind_uni(iris2, exclude=c("pid", "Species"), m=3)


##Function to identify individual multivariate outliers
ind_multi = function(d0, exclude=c("pid", "site"), thresh=0.7, n_uniq=10) {
 d1 = d0[, !names(d0) %in% exclude]
 tt = apply(d1, 2, nVal) #only keep those with sufficient variation
 keep = names(tt)[tt>=n_uniq]
 d2 = d1[, names(d1) %in% keep]

 model = suppressWarnings(isolation.forest(d2))
 d0$outlier_score = predict(model, d2)
 n_out = length(d0$outlier_score[d0$outlier_score>thresh])
 if (n_out==0) {
  data = "No outliers found"
  } else if (n_out>0) {
    d3 = d0[d0$outlier_score>thresh, ]
    data = d3[, match(c(names(d3)[names(d3) %in% exclude], names(d3)[!names(d3) %in% exclude]), names(d3))]  
   }
 nID = dim(d2)[1]
 nVar = dim(d2)[2]
 output.all = list(nID=nID, nVar=nVar, data=data)
 return(output.all)
}

#ind_multi(iris2, exclude=c("pid", "Species"), thresh=0.7, n_uniq=10)
#ind_multi(iris2, exclude=c("pid", "Species"), thresh=0.6, n_uniq=10)


##Function to identify differences between sites
site_outliers = function(d0, exclude=c("pid"), siteID="site", covs=c("age"), threshG=0.001, thresh2=0.05, threshS=0.5, n_uniq=10, n_dec=4, n_decS=2) {
 d1 = d0[, !names(d0) %in% exclude]
 vars = d1[, names(d1)[!names(d1) %in% c(siteID, covs)]] #only keep those variables with sufficient variation
 tt = apply(vars, 2, nVal)
 keep = names(tt)[tt>=n_uniq]
 d2 = d0[, match(c(keep, siteID, covs), names(d0))]
 d2[, names(d2)==siteID] = as.factor(d2[, names(d2)==siteID])
 names(d2)[names(d2)==siteID] = "site" #rename if needed

 pp = data.frame(var=keep, p_unadj=NA, p_adj=NA)
 for (i in 1:length(keep)) {
  elim = keep[-i]
  d3 = d2[, !names(d2) %in% elim]
  names(d3)[1] = "outc"
  unadj = lm(outc ~ site, data=d3)
  adj = lm(outc ~ ., data=d3)
  pp$p_unadj[pp$var==keep[i]] = anova(unadj)["site", "Pr(>F)"]
  pp$p_adj[pp$var==keep[i]] = anova(adj)["site", "Pr(>F)"]
 }

 test2 = pp$var[pp$p_unadj<threshG | pp$p_adj<threshG]
 d2 = d2[, match(c(test2, "site", covs), names(d2))]
 if (length(test2)==0) {output.all="No outliers found"}
 else if (length(test2)>0) {
  ##For those variables with unadj and adj p-values<XX, then look at site-wise deviations
  un_Pval = matrix(NA, nrow=length(unique(d2$site)), ncol=length(test2))
  rownames(un_Pval) = unique(d2$site)
  colnames(un_Pval) = test2
  un_StDf = matrix(NA, nrow=length(unique(d2$site)), ncol=length(test2))
  rownames(un_StDf) = unique(d2$site)
  colnames(un_StDf) = test2
  a_Pval = matrix(NA, nrow=length(unique(d2$site)), ncol=length(test2))
  rownames(a_Pval) = unique(d2$site)
  colnames(a_Pval) = test2
  a_StDf = matrix(NA, nrow=length(unique(d2$site)), ncol=length(test2))
  rownames(a_StDf) = unique(d2$site)
  colnames(a_StDf) = test2
  for (i in 1:length(unique(d2$site))) {
   for (j in 1:length(test2)) {
    elim = test2[-j]
    d3 = d2[, !names(d2) %in% elim]
    names(d3)[1] = "outc"
    d3$site_ind = ifelse(d3$site==unique(d3$site)[i], 1, 0)
    t_test = t.test(outc ~ site_ind, data=d3)
    un_Pval[i, j] = round(t_test$p.value, n_dec)
    stdf = stddiff.numeric(d3, gcol=which(names(d3)=="site_ind"), vcol=which(names(d3)=="outc"))
    un_StDf[i, j] = round(stdf[colnames(stdf)=="stddiff"], n_decS)
    d4 = d3[, names(d3)!="site"]
    mod = lm(outc ~ ., data=d4)
    a_Pval[i, j] = round(summary(mod)$coef["site_ind", "Pr(>|t|)"], n_dec)
    a_StDf[i, j] = round(abs(summary(mod)$coef["site_ind", "Estimate"]/sqrt(anova(mod)["Residuals", "Mean Sq"])), n_decS) #abs(coeff/RMSE)
   }
  }

 ##Only print those with overall p<threshG
 pp2 = pp[pp$p_unadj<threshG | pp$p_adj<threshG, ]
 pp2$p_unadj = round(pp2$p_unadj, n_dec)
 pp2$p_adj = round(pp2$p_adj, n_dec)

 ##Only print those with p<thresh2
 un_Pval[un_Pval>=thresh2] = NA
 a_Pval[a_Pval>=thresh2] = NA

 ##Only print those with stddiff>=threshS
 un_StDf[un_StDf<=threshS] = NA
 a_StDf[a_StDf<=threshS] = NA

 ##Further reduce the matrices to remove missing data and preserve row names for those with only one observation
 i_p = apply(un_Pval, 1, function(x) sum(is.na(x)))!=length(test2)
 sitewise_P=un_Pval[i_p, ]
 if (length(dim(sitewise_P))==0) {
  nms = names(sitewise_P)
  sitewise_P = matrix(sitewise_P, 1, length(sitewise_P))
  rownames(sitewise_P) = rownames(un_Pval)[i_p]
  colnames(sitewise_P) = nms
 }

 i_pa = apply(a_Pval, 1, function(x) sum(is.na(x)))!=length(test2)
 sitewise_P_adj=a_Pval[i_pa, ]
 if (length(dim(sitewise_P_adj))==0) {
  nms = names(sitewise_P_adj)
  sitewise_P_adj = matrix(sitewise_P_adj, 1, length(sitewise_P_adj))
  rownames(sitewise_P_adj) = rownames(a_Pval)[i_pa]
  colnames(sitewise_P_adj) = nms
 }

 i_sd = apply(un_StDf, 1, function(x) sum(is.na(x)))!=length(test2)
 sitewise_StDf=un_StDf[i_sd, ]
 if (length(dim(sitewise_StDf))==0) {
  nms = names(sitewise_StDf)
  sitewise_StDf = matrix(sitewise_StDf, 1, length(sitewise_StDf))
  rownames(sitewise_StDf) = rownames(un_StDf)[i_sd]
  colnames(sitewise_StDf) = nms
 }

 i_sda = apply(a_StDf, 1, function(x) sum(is.na(x)))!=length(test2)
 sitewise_StDf_adj=a_StDf[i_sda, ]
 if (length(dim(sitewise_StDf_adj))==0) {
  nms = names(sitewise_StDf_adj)
  sitewise_StDf_adj = matrix(sitewise_StDf_adj, 1, length(sitewise_StDf_adj))
  rownames(sitewise_StDf_adj) = rownames(un_StDf)[i_sda]
  colnames(sitewise_StDf_adj) = nms
 }

 output.all = list(overall=pp2, sitewise_P=sitewise_P, sitewise_P_adj=sitewise_P_adj, sitewise_StDf=sitewise_StDf, sitewise_StDf_adj=sitewise_StDf_adj)
 }
 return(output.all)
}


#iris2$temp = rnorm(dim(iris2)[1]) #for covariate adjustment
#site_outliers(iris2, site="Species", covs=c("temp"))

##Create package
#package.skeleton(name="bulkQC", code_files="bulkQC.R")
