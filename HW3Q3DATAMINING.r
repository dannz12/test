library(glmnet)

prostate_cancer_data_raw = read.csv('C:/prostate.csv')
prostate_cancer_data_cleaned = prostate_cancer_data_raw[complete.cases(prostate_cancer_data_raw),]


sample_size <- floor(.33 * nrow(prostate_cancer_data_cleaned))
set.seed(42)
test_idx <- sample(seq_len(nrow(prostate_cancer_data_cleaned)), size = sample_size)

prostate_cancer_train <- prostate_cancer_data_cleaned[-test_idx,]
prostate_cancer_test <- prostate_cancer_data_cleaned[test_idx,]

dep <- c("lcavol")
prostate_cancer_tdep <- subset(prostate_cancer_train, select = dep)
prostate_cancer_tindep <-subset(prostate_cancer_train, select = names(prostate_cancer_train) != dep)

lasso_alphas <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
ridge_alphas <- c(0.00001, 0.0001,0.001, 0.005, 0.01, 0.05, 0.1, 1, 5, 10, 100)

#part 1
mprostate_cancer_tindep <- as.matrix(prostate_cancer_tindep)
mprostate_cancer_tdep <- as.matrix(prostate_cancer_tdep)

lassocv <- cv.glmnet(mprostate_cancer_tindep, mprostate_cancer_tdep, alpha=1, nfold=5, lambda=lasso_alphas)
ridgecv <- cv.glmnet(mprostate_cancer_tindep, mprostate_cancer_tdep, alpha=0, nfold=5, lambda=ridge_alphas)
lasso_best_alpha <- as.numeric(lassocv["lambda.min"])
ridge_best_alpha <- as.numeric(ridgecv["lambda.min"])
                                  
cat("Lasso best alpha Value: ", lasso_best_alpha)
cat("Ridge best alpha value: ", ridge_best_alpha)

#part 2

symbols = c(1:11)


lassomdl <- glmnet(mprostate_cancer_tindep, mprostate_cancer_tdep, alpha=1, lambda=lasso_alphas)

lasso_coef <- lassomdl[["beta"]]
row_names = rownames(lasso_coef)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(c(1:8), lasso_coef[,1], axes=FALSE, pch=symbols[1], ylim=c(min(lasso_coef), max(lasso_coef)))
axis(2)
axis(1, at=1:8, labels=row_names, las=2)
for(i in 2:nrow(lasso_coef)){
  points(c(1:8), lasso_coef[,i], pch=symbols[i])
}

legend(8.5,1, legend=lasso_alphas, pch=symbols)

print("As the regularization parameter is changed, the coefficients associated with each category 
      vary, where some coefficients increase then decrease as alpha increases or vice versa. Some coefficients
      are not influenced as shown in the plot")

#part 3
prostate_cancer_ttdep <- subset(prostate_cancer_test, select = dep)
prostate_cancer_ttindep <-subset(prostate_cancer_test, select = names(prostate_cancer_test) != dep)

mprostate_cancer_ttindep <- as.matrix(prostate_cancer_ttindep)
mprostate_cancer_ttdep <- as.matrix(prostate_cancer_ttdep)


#using best alpha values obtained from part 1
lasso_fit <- glmnet(mprostate_cancer_tindep, mprostate_cancer_tdep, alpha=1, lambda=0.01)
ridge_fit <- glmnet(mprostate_cancer_tindep, mprostate_cancer_tdep, alpha=1, lambda=0.1)
linear_fit <- lm(mprostate_cancer_tdep~., data=prostate_cancer_tindep)
lasso_actual = predict(lasso_fit, mprostate_cancer_ttindep)
ridge_actual = predict(ridge_fit, mprostate_cancer_ttindep)
linear_actual = predict(linear_fit, newdata=prostate_cancer_ttindep)
linear_actual= t(unname(linear_actual))

findMSE <- function(actuallist, expectedlist){
  sum = 0
  for(i in 1:nrow(actuallist)){
    t = (actuallist[i] - expectedlist[i])^2
    #print(t)
    sum = sum + t
  }
  res = sum/nrow(actuallist)
  return(res)
}


lasso_mse = findMSE(lasso_actual, mprostate_cancer_ttdep)
ridge_mse = findMSE(ridge_actual, mprostate_cancer_ttdep)
linear_mse = findMSE(linear_actual, mprostate_cancer_ttdep)

cat("Lasso prediction error:", lasso_mse)
cat("Ridge prediction error:", ridge_mse)
cat("Linear prediction error:", linear_mse)

#part 4

linear_fit <- lm(mprostate_cancer_tdep~., data=prostate_cancer_tindep)
linear_fit_best <- lm(mprostate_cancer_tdep~., data=prostate_cancer_tindep[,-c(2,3,7)])



linear_actual_best = predict(linear_fit_best, newdata=prostate_cancer_ttindep[,-c(2,3,7)])
linear_actual_best = t(unname(linear_actual_best))
linear_mse_best = findMSE(linear_actual_best, mprostate_cancer_ttdep)
cat("Linear prediction error with selective variables:", linear_mse_best)
