MultStepwiseReg <- function(datas, params){
## Simple function used for mult-mult stepwise regression.
##
## Input:
##    data: a list contain the Y data and the X data
##        data$Y: Y data
##        data$X: X data
##
##    params: a list contains all the parameters used in the function
##          params$standardize: 0 means not standardize
##                              1 means standardize using the standard error
##          params$verbose: 0 means silence
##                          1 means output information at the first and the last iteration
##                          2 means output information at every iteration
##          params$ain: criterion for entrance of the variable
##          params$aout: criterion for clearance of the variable
##          params$max_iter: maximum iteration times for the inner loop
##
## Output: 
##    result: a list contains all the debug information 
##          result$iteration: total iteration times in the loop
##          result$chosenvar: variables chosen during the iteration
##          result$p_val: contains each ps calculated in the function
##          result$F_n1: each F parameter 2 used in the function
##          result$F_n2: each F paramater 1 used in the function
##          result$F_num: each F num used in the function
##          result$beta_est: estimate of the coefficient matriex beta, also given while none of the variable is included
##          result$sigma_est: estimate of the covariance matrix, also given while none of the variable is included
##          result$Rsquare:  sum of square residual, also given while none of the variable is included
##          result$multiple_cor: multiple correlation of y, also given while none of the variable is included
##
## by Liyifan Gan  April 17, 2022 copyright 2018-2022 LoveTianyi
  
  # extract and get the basic information of the data
  Y <- as.matrix(datas$Y)
  X <- as.matrix(datas$X)
  dim_y <- dim(Y)
  dim_x <- dim(X)
  if(dim_x[1] != dim_y[1])
  {
    stop("The length of the data is not the same.") 
  }
  
  # whether the data need to be standardize or not
  if(params$standardize)
  {
    for(i in dim_x[2])
    {
      X[,i] <- (X[,i] - mean(X[,i])) / sd(X[,i])
    }
    
    for(i in dim_y[2])
    {
      Y[,i] <- (Y[,i] - mean(Y[,i])) / sd(Y[,i])
    }
  }
  
  # preparation for the regression, generation of the L_0 matrix - centralized form
  X_tit <- X
  Y_tit <- Y
  for(i in seq(1,dim_x[2]))
  {
    X_tit[,i] <- X[,i] - mean(X[,i])
  }
  
  for(i in seq(1,dim_y[2]))
  {
    Y_tit[,i] <- Y[,i] - mean(Y[,i])
  }
  
  Lxx <- matrix(t(X_tit) %*% X_tit, dim_x[2], dim_x[2])
  Lxy <- matrix(t(X_tit) %*% Y_tit, dim_x[2], dim_y[2])
  Lyy <- matrix(t(Y_tit) %*% Y_tit, dim_y[2], dim_y[2])
  L_0 <- rbind(cbind(Lxx, Lxy), cbind(t(Lxy), Lyy))
  
  flag <- rep(0, dim_x[2])
  
  if(params$verbose >= 1)
  {
    cat("Stepwise regression of multiple dependent variables!\n")
  }
  
  # the first step
  steptime = 1
  
  # calculate the contribution of Xs
  V <- rep(0, length(flag))
  
  for(i in seq(1,dim_x[2]))
  {
    V[i] <- Lxy[i,] %*% solve(Lyy) %*% Lxy[i,] / L_0[i,i]
  }
  
  # get the maximum number and the position in the Vs
  V_max <- max(V)
  pos <- which.max(V)
  
  # hypothesis test for the chosen variable with F distribution
  p_val <- NaN
  F_n1 <- NaN
  F_n2 <- NaN
  F_num <- NaN
  F_n1[steptime] <- dim_y[1] - dim_y[2] - 1 
  F_n2[steptime] <- dim_y[2]
  F_num[steptime] <- F_n1 / F_n2 * V_max / (1 -  V_max)
  p_val[steptime] <- 1 - pf(F_num[steptime], F_n2[steptime], F_n1[steptime])
  
  # if fail for entrance at the very first beginning
  if(p_val[steptime] >= params$ain) 
  {
    warning("Iteration fails at the beginning, no variable is chosen!")
    
    L_r <- L_0
    beta_est <- matrix(0, 2, 1)
    sigma_est <- 1 / (dim_y[1] - sum(flag) - 1) * Lyy
    multiple_cor <- NaN
    for(i in 1:dim_y[2])
    {
      idx <- dim_x[2] + i
      beta_est[i,1] <- mean(Y[,1]) - L_r[chosenvar, idx] %*% apply(X, 2, mean)[chosenvar]
      Rsquare[i] <-  L_r[idx, idx]
      multiple_cor[i] <- sqrt(1 - L_r[idx, idx] / L_0[idx, idx])
    }
    beta_est <- cbind(beta_est, Lxy[chosenvar,])
    
    result <- list(
      steptime = steptime,
      chosenvar = flag,
      p_val = p_val,
      F_n1 = F_n1,
      F_n2 = F_n2,
      F_num = F_num,
      beta_est = beta_est,
      sigma_est = sigma_est,
      Rsquare = Rsquare,
      multiple_cor = multiple_cor
    )
    return(result)
  }
  
  # include the variable and preparation for the next iteration
  flag[pos] <- 1
  L_r <- sweepTransform(L_0, pos, pos)
  Lxx <- matrix(L_r[1:dim_x[2], 1:dim_x[2]],
                length(1:dim_x[2]), length(1:dim_x[2]))
  Lxy <- matrix(L_r[1:dim_x[2], (dim_x[2]+1):dim(L_r)[1]],
                   length(1:dim_x[2]), length((dim_x[2]+1):dim(L_r)[1]))
  Lyx <- matrix(L_r[(dim_x[2]+1):dim(L_r)[1], 1:dim_x[2]], 
                   length((dim_x[2]+1):dim(L_r)[1]), length(1:dim_x[2]))
  Lyy <- matrix(L_r[(dim_x[2]+1):dim(L_r)[1], (dim_x[2]+1):dim(L_r)[1]],
                   length((dim_x[2]+1):dim(L_r)[1]), length((dim_x[2]+1):dim(L_r)[1]))
  
  # main loop
  stop = 0
  while(!stop && steptime < params$max_iter)
  {
    steptime = steptime + 1
    
    # clearance of the variable
    V <- rep(Inf, length(flag))
    for(i in which(flag == 1))
    {
      V[i] <- -Lxy[i,] %*% solve(Lyy) %*% Lyx[,i] / L_r[i,i]
    }
    
    # get the minimum of V
    V_min <- min(V)
    pos <- which.min(V)
    
    # hypothesis test for clearance of the variable with F distribution
    F_n1[steptime] <- dim_y[1] - dim_y[2] - sum(flag)
    F_n2[steptime] <- dim_y[2]
    F_num[steptime] <- F_n1[steptime] / F_n2[steptime] * V_min
    p_val[steptime] <- 1 - pf(F_num[steptime], F_n2[steptime], F_n1[steptime])
    
    if(p_val[steptime] >= params$aout)
    {
      if(params$verbose >= 2)
      {
        cat("Variable",pos,"is excluded in the",steptime,"step\n")
      }
      
      # renew the L_r matrix
      L_r <- sweepTransform(L_r, pos, pos)
      Lxx <- matrix(L_r[1:dim_x[2], 1:dim_x[2]],
                    length(1:dim_x[2]), length(1:dim_x[2]))
      Lxy <- matrix(L_r[1:dim_x[2], (dim_x[2]+1):dim(L_r)[1]],
                    length(1:dim_x[2]), length((dim_x[2]+1):dim(L_r)[1]))
      Lyx <- matrix(L_r[(dim_x[2]+1):dim(L_r)[1], 1:dim_x[2]], 
                    length((dim_x[2]+1):dim(L_r)[1]), length(1:dim_x[2]))
      Lyy <- matrix(L_r[(dim_x[2]+1):dim(L_r)[1], (dim_x[2]+1):dim(L_r)[1]],
                    length((dim_x[2]+1):dim(L_r)[1]), length((dim_x[2]+1):dim(L_r)[1]))
      
      # check to see if there's no variable left, if so, give the warning and exit the program
      if(sum(flag == 0))
      {
        warning("All the variables are excluded!")
        
        sigma_est <- 1 / (dim_y[1] - sum(flag) - 1) * Lyy
        multiple_cor <- NaN
        for(i in 1:dim_y[2])
        {
          idx <- dim_x[2] + i
          beta_est[i,1] <- mean(Y[,1]) - L_r[chosenvar, idx] %*% apply(X, 2, mean)[chosenvar]
          Rsquare[i] <-  L_r[idx, idx]
          multiple_cor[i] <- sqrt(1 - L_r[idx, idx] / L_0[idx, idx])
        }
        beta_est <- cbind(beta_est, Lxy[chosenvar,])
        
        result <- list(
          steptime = steptime,
          chosenvar = flag,
          p_val = p_val,
          F_n1 = F_n1,
          F_n2 = F_n2,
          F_num = F_num,
          beta_est = beta_est,
          sigma_est = sigma_est,
          Rsquare = Rsquare,
          multiple_cor = multiple_cor
        )
        return(result)
      }
      next
    }
    
    steptime = steptime + 1
    
    # entrance of the variable
    V <- rep(0, length(flag))
    for(i in which(flag == 0))
    {
      V[i] <- Lxy[i,] %*% solve(Lyy) %*% Lyx[,i] / L_r[i,i]
    }
    
    # get the maximum of V
    V_max <- max(V)
    pos <- which.max(V)
    
    # hypothesis test for entrance of the variable with F distribution
    F_n1[steptime] <- dim_y[1] - dim_y[2] - sum(flag) - 1
    F_n2[steptime] <- dim_y[2]
    F_num[steptime] <- F_n1[steptime] / F_n2[steptime] * V_max / (1 - V_max)
    p_val[steptime] <- 1 - pf(F_num[steptime], F_n2[steptime], F_n1[steptime])
    
    if(p_val[steptime] < params$ain)
    {
      if(params$verbose >= 2)
      {
        cat("Variable",pos,"is included in the",steptime,"step\n")
      }
      
      # include the variable and prepare for the next iteration
      flag[pos] <- 1
      L_r <- sweepTransform(L_r, pos, pos)
      Lxx <- matrix(L_r[1:dim_x[2], 1:dim_x[2]],
                    length(1:dim_x[2]), length(1:dim_x[2]))
      Lxy <- matrix(L_r[1:dim_x[2], (dim_x[2]+1):dim(L_r)[1]],
                    length(1:dim_x[2]), length((dim_x[2]+1):dim(L_r)[1]))
      Lyx <- matrix(L_r[(dim_x[2]+1):dim(L_r)[1], 1:dim_x[2]], 
                    length((dim_x[2]+1):dim(L_r)[1]), length(1:dim_x[2]))
      Lyy <- matrix(L_r[(dim_x[2]+1):dim(L_r)[1], (dim_x[2]+1):dim(L_r)[1]],
                    length((dim_x[2]+1):dim(L_r)[1]), length((dim_x[2]+1):dim(L_r)[1]))
    }
    else
    {
      cat("The Stepwise Regression Process ended!\n")
      stop = 1
    }
  }
  
  # output the result
  if(params$verbose >= 1)
  {
    cat("The final choice of the variable is Variable", which(flag == 1),"\n")
  }
  
  # give the other necessary result
  beta_est <- matrix(0, 2, 1)
  sigma_est <- 1 / (dim_y[1] - sum(flag) - 1) * Lyy
  multiple_cor <- NaN
  Rsquare <- NaN
  chosenvar <- which(flag == 1)
  for(i in 1:dim_y[2])
  {
    idx <- dim_x[2] + i
    beta_est[i,1] <- mean(Y[,1]) - L_r[chosenvar, idx] %*% apply(X, 2, mean)[chosenvar]
    Rsquare[i] <-  L_r[idx, idx]
    multiple_cor[i] <- sqrt(1 - L_r[idx, idx] / L_0[idx, idx])
  }
  beta_est <- cbind(beta_est, Lxy[chosenvar,])
  
  # collect the debug information
  result <- list(
    steptime = steptime,
    chosenvar = flag,
    p_val = p_val,
    F_n1 = F_n1,
    F_n2 = F_n2,
    F_num = F_num,
    beta_est = beta_est,
    sigma_est = sigma_est,
    Rsquare = Rsquare,
    multiple_cor = multiple_cor
  )
  
  return(result)
}

sweepTransform <- function(A, i, j){
## Sweep transform for matrix A
##
## Input:
##      A: a matrix for sweep transform
##      i: row index for sweep transform, no more than the dim of A
##      j: column index for sweep transform, no more than the dim of A
##
## Output:
##      B: a matrix, transformed result
##
## by Liyifan Gan  April 18, 2022 copyright 2018-2022 LoveTianyi
  
  dim_A <- dim(A)
  B <- matrix(0, dim_A[1], dim_A[2])
  # generate the index
  if(i == 1)
  {
    id_i <- 2:dim_A[1]
  }
  else if(i == dim_A[1])
  {
    id_i <- 1:(dim_A[1] - 1)
  }
  else
  {
    id_i <- c(1:(i - 1), (i + 1):dim_A[1])
  }
  
  if(j == 1)
  {
    id_j <- 2:dim_A[2]
  }
  else if(j == dim_A[2])
  {
    id_j <- 1:(dim_A[2] - 1)
  }
  else
  {
    id_j <- c(1:(i - 1), (i + 1):dim_A[1])
  }
  
  # case 1
  B[i, j] <- 1 / A[i, j]
  
  # case 2
  for(k in id_i)
  {
    B[k, j] <- -A[k, j] / A[i, j] 
  }
  
  # case 3
  for(k in id_i)
  {
    B[i, k] <- A[i, k] / A[i, j]
  }
  
  # case 4
  for(k in id_i)
  {
    for(s in id_j)
    {
      B[k, s] <- A[k, s] - A[k, j] * A[i, s] / A[i, j]
    } 
  }
    
  return(B)
}