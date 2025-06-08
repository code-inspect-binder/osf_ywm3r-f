# function pcor2beta gives you data from a partial correlation/network
#-- input --
# pcor = a partial correlation matrix / network
# -- output --
# a matrix of betas, each column corresponds to a dependent variable
# so that you can get predicted values by a matrix multiplication
# in the form betas %*% data

pcor2beta <- function(pcor)
{
  require(psych)
  require(corpcor)
  diag(pcor) <- 1
  p <- ncol(pcor)
  betas <- matrix(0, ncol = p, nrow = p)
  
  for(i in 1:p)
    betas[-i,i] <- matReg(y = i, x = seq(p)[-i],
                          C = pcor2cor(pcor))$beta
  
  betas[abs(betas) < 1e-13] <- 0
  betas
}


# function R2 gives you two different types of R2 and of predicted values
#-- input --
# pcor = a partial correlation matrix / network
# -- output --
# a matrix of betas, each column corresponds to a dependent variable
# so that you can get predicted values by a matrix multiplication
# in the form betas %*% data

# this function gives you R2 and predicted values from betas + data
# two types of R2 and predicted values are considered:
# - R2_orig and predicted_orig use beta weights directly implied by the
#   graphical lasso regularization
# - R2_refit and predicted_refit use only the sparsity pattern of the network
#   and then refit linear regression using the network only for prediction
#   but not for shrinkage (this one should do better in terms of prediction)
# - refit: logical, regulates whether R2_refit and predicted_refit are computed.
#   set it to FALSE to speed up computations

R2 <- function(betas, dt, refit = TRUE)
{
  dt <- data.frame(scale(dt))
  out <- list()
  p <- ncol(betas)
  # predicted values
  predicted <- as.matrix(dt) %*% betas
  # R squared using the formula in Haslbeck & Waldorp (2017, BRM)
  R2 <- 1-apply(predicted - dt, 2,var)
  out$predicted_orig <- predicted
  out$R2_orig <- R2
  
  if(refit)
  {
    # refit the model considering the sparsity pattern inicated by the network
    # first fit regression using the sparsity indciated in the matrix of betas
    betas_refit <- matrix(0, ncol = p, nrow = p)
    for(i in 1:p)
    {
      if(any(betas[,i] != 0))
      {
        fit <- lm(dt[,i] ~ as.matrix(dt[,betas[,i] != 0]))
        betas_refit[betas[,i] != 0, i] <- fit$coefficients[-1]
      } else betas_refit[,i] <- 0
    }
    
    predicted <- as.matrix(dt) %*% betas_refit
    # R squared
    R2<- 1-apply(predicted - dt, 2,var)
    out$predicted_refit <- predicted
    out$R2_refit <- R2
  }
  out
}

# a wrapper that gives you directly the predictability
predictability <- function(net, dt, refit = FALSE)
{
  betas <- pcor2beta(net)
  predict <- R2(betas, dt = dt, refit = refit)
  predict$R2_orig
}
