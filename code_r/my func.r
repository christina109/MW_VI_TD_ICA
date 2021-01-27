which.localMax <- function(mat, xVec = NaN, yVec = NaN){
  
  maxPositionCol <- array(0,ncol(mat))
  maxValCol      <- maxPositionCol
  for (ci in 1:ncol(mat)){
    maxPositionCol[ci] <- which.max(mat[,ci])
    maxValCol[ci] <- mat[maxPositionCol[ci], ci]
  }
  colId <- which.max(maxValCol)
  rowId <- maxPositionCol[colId]
  
  pos <- c(rowId, colId)
  if (!is.nan(xVec)[1]){
    pos <- c(yVec[rowId], xVec[colId])
  }
  
  return(pos)
  
}



mean_ci <- function(x, confidence=0.95) {
  if (length(x)==0){
    warning("'x' was empty", call.=FALSE)
    interval<-c(-Inf, Inf)
  }else{
    se <- sd(x) / sqrt(length(x)-1)
    alpha <- 1 - confidence
    interval<-mean(x) + se * qnorm(c(alpha / 2, 1 - alpha / 2))
  }
  return(interval)
}



smd <- function(x1, x2, s1, s2, n1, n2) {
  sp <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2)/(n1 + n2 -2))
  smd <- (x1 - x2)/sp
  return(smd)
}



t2smd <- function(t, n1, n2) {
  smd <- t*sqrt(1/n1 + 1/n2)
  return(smd)
}



se4smd <- function(smd, n1, n2) {
  se <- sqrt((n1+n2)/n1/n2 + smd^2/2/(n1+n2))
  return(se)
}


ci4smd <- function(smd, se, confidence = 0.95) {
  alpha <- 1 - confidence
  interval <- smd + se * qnorm(c(alpha/2, 1-alpha/2))
  return(interval)
}


