library(randomForest)

# function
rfcv1 <- function(trainx, trainy, cv.fold = 5, scale = "log", step = 0.5,output_dir=".", mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE,
  ipt = NULL, ...) {
  classRF <- is.factor(trainy)
  n <- nrow(trainx) # sample size
  p <- ncol(trainx) # variables
  if (scale == "log") {
    k <- floor(log(p, base = 1/step))
    n.var <- round(p * step^(0:(k - 1)))
    same <- diff(n.var) == 0
    if (any(same))  # only keep the distinct ones
      n.var <- n.var[-which(same)]
    if (!1 %in% n.var)
      n.var <- c(n.var, 1)
  } else {
    n.var <- seq(from = 1, to = p, by = step)
  }
  # Note: n.var is the size of sequential marker sets, from p to 1 in decreasing order
  k <- length(n.var)
  cv.pred <- vector(k, mode = "list")

  for (i in 1:k) cv.pred[[i]] <- trainy
  if (classRF) {
    f <- trainy
    if (is.null(ipt))
      ipt <- nlevels(trainy) + 1
  } else {
    f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])  #???
    if (is.null(ipt))
      ipt <- 1
  }
  nlvl <- table(f)
  idx <- numeric(n) # index of samples in each Y group
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, length = nlvl[i]))
  }
  res = list()
  for (i in 1:cv.fold) {
    all.rf <- randomForest(trainx[idx != i, , drop = FALSE], trainy[idx != i], trainx[idx == i, , drop = FALSE], trainy[idx ==
      i], mtry = mtry(p), importance = TRUE, ...)
    cv.pred[[1]][idx == i] <- all.rf$test$predicted

    impvar <- (1:p)[order(all.rf$importance[, ipt], decreasing = TRUE)]
    # write.table(impvar, paste0(output_dir,"/feature_importances_",i,".txt"),sep = "\t", quote = F)
    res[[i]] <- impvar

    for (j in 2:k) {
      cat("Fold:", i,"/",cv.fold, "-->> Gene set:", j,"/", k, "\n")
      imp.idx <- impvar[1:n.var[j]]
      sub.rf <- randomForest(trainx[idx != i, imp.idx, drop = FALSE], trainy[idx != i], trainx[idx == i, imp.idx, drop = FALSE],
        trainy[idx == i], mtry = mtry(n.var[j]), importance = recursive, ...)
      cv.pred[[j]][idx == i] <- sub.rf$test$predicted
      if (recursive) {
        impvar <- (1:length(imp.idx))[order(sub.rf$importance[, ipt], decreasing = TRUE)]
      }
      NULL
    }
    NULL
  }

  #print(cv.pred)
  if (classRF) {
    error.cv <- sapply(cv.pred, function(x) mean(trainy != x))
  } else {
    error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
  }
  names(error.cv) <- names(cv.pred) <- n.var
  list(n.var = n.var, error.cv = error.cv, predicted = cv.pred, res = res)
}

library(mRMRe)
mRMR_features <- function(train, output_dir, cv.fold = 5,cv.time=5, scale = "log", step = 2,marker.num=0){
set.thread.count(4)
  dd <- mRMR.data(data =  data.frame(train))
  filter <- mRMR.ensemble(data = dd, target_indices = ncol(train), levels=c(1:ncol(train)-1))
  return()
}


mRMR_features_old <- function(train, output_dir, cv.fold = 5,cv.time=5, scale = "log", step = 2,marker.num=0){
  set.seed(0)
  
  train.x <- train[, !names(train) %in% c("Y")]
  train.y <- as.factor(train$Y)
  
  train.cv <- replicate(cv.time, rfcv1(train.x, train.y, cv.fold = cv.fold, scale = scale, step = step,output_dir=output_dir), simplify = F)

  error.cv <- sapply(train.cv, "[[", "error.cv")
  error.cv.rm <- rowMeans(error.cv)
  # id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
  # Note by Xianjun Dong: instead of using overall sd, this should be changed to use the corresponding se (or sd) at the position with the minimal mean
  # Ref: see slide 13 of https://web.stanford.edu/class/stats202/content/lec11-cond.pdf
  # The method is actually from Maximum Relevance Minimum Redundancy (mRMR):
  # Peng et al., 2005, “Feature selection based on mutual information: criteria of max-dependency, max-relevance, and min-redundancy,” IEEE Trans Pattern Anal Mach Intell., 27, 1226-1238, doi:10.1109/TPAMI.2005.159
  # https://patentswarm.com/patents/US10526659B2

  error.cv.se <- apply(error.cv, 1, function(x) sd(x)/sqrt(length(x)))
  id <- error.cv.rm < min(error.cv.rm) + error.cv.se[which.min(error.cv.rm)]

  if (marker.num == 0) {
    marker.num <- min(as.numeric(names(error.cv.rm)[id]))
  }

  matplot(train.cv[[1]]$n.var, error.cv, type = "l", log = "x", col = rep(1, cv.time),
          main = paste("select", marker.num, "Vars"), xlab = "Number of vars",
    ylab = "CV Error", lty = 1)
  lines(train.cv[[1]]$n.var, error.cv.rm, lwd = 2)
  abline(v = marker.num, col = "pink", lwd = 2)

  # pick marker by coross-validation
  marker.t <- table(unlist(lapply(train.cv, function(x) {
    lapply(x$res, "[", 1:marker.num)
  })))
  
  marker.t <- sort(marker.t, d = T)
  names(marker.t) <- colnames(train.x)[as.numeric(names(marker.t))]
  marker.dir <- paste0(output_dir, "/marker.txt")
  write.table(marker.t, marker.dir, col.names = F, sep = "\t", quote = F)
  marker.p <- names(marker.t)[1:marker.num]

  return(marker.p)
}

library(glmnet)
LASSO_features <- function(train, output_dir){
  set.seed(0)
  
  train.x <- train[, !names(train) %in% c("Y")]
  train.y <- train$Y

  x <- as.matrix(train.x) # all X vars
  y <-  as.double(as.matrix(train.y)) # Only Class

  # Fit the LASSO model (Lasso: Alpha = 1)
  cv.lasso <- cv.glmnet(x, y, family='binomial', alpha=1, parallel=TRUE,nfolds=5,type.measure="mse")
  plot(cv.lasso)

  df_coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 4)
  # See all contributing variables
  df_coef.dir <- paste0(output_dir, "/lasso_coef.txt")
  write.table(df_coef, df_coef.dir, col.names = F, sep = "\t", quote = F)

  df_coef = df_coef[df_coef[, 1] != 0, ,drop=FALSE]
  df_coef = df_coef[ !rownames(df_coef)=='(Intercept)', ,drop=FALSE]
  # df_coef = df_coef[order(df_coef[, 1]), ,drop=FALSE]

  marker.p = variables<-row.names(df_coef)
  return(marker.p)
}


library(pROC)
# function
plot_roc <- function(response, predictor, direction = "auto",main="plot") {
  print(paste0("warning: levels of respone is ", paste(levels(as.factor(response)), collapse = ", "),
               " and should corresponding to controal and case, the default direction is auto"))
  roc.obj <- roc(response, predictor, percent = T, ci = T, plot = T, direction = direction)
  ci.se.obj <- ci.se(roc.obj, specificities = seq(0, 100, 5))
  plot(ci.se.obj, type = "shape", col = rgb(0, 1, 0, alpha = 0.2))
  plot(ci.se.obj, type = "bars")
  plot(roc.obj, col = 2, add = T,main=main)
  txt <- c(paste("AUC=", round(roc.obj$ci[2], 2), "%"), paste("95% CI:", round(roc.obj$ci[1], 2), "%-", round(roc.obj$ci[3], 2),
    "%"))
  legend("bottomright", txt)
}

