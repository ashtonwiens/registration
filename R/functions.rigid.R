

#' Translate the x and y-coordinates in the first three columns
#'
#' @param df the dataframe with coordinates in first three columns
#' @param tr a vector of length 3 with the x, y, and z translations
#'
#' @return the dataframe with translated coordinates
#' @export
#'
#' @examples
translate <- function(df, tr){
  df %>% mutate( X = .[[1]]+tr[1], Y = .[[2]]+tr[2], Z = .[[3]]+tr[3] )
}

#' Title
#'
#' @param df the dataframe with coordinates in first three columns
#' @param tr a rotation angle in radians
#'
#' @return the dataframe with rotated coordinates
#' @export
#'
#' @examples
rotate2d <- function(df, rt){
  phi <- rt
  rotat <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), nr=2, nc=2 )
  g <- t( rotat %*% t( data.matrix(df %>% dplyr::select(X,Y)) ) )
  df[,1:2] <- g
  return(df)
}

#' Log likelihood function with embedded registration
#'
#' Registration includes translations in the x, y, and z coordinates and a rotate in the xy plane
#'
#' @param p parameters defining Gaussian process covariance and registration
#' @param nu smoothness parameter of GP
#' @param grd three column matrix of point set 1
#' @param grd2 three column matrix of point set 2
#'
#' @return the value of the negative log likelihood
#' @export
#'
#' @examples
logLik.translate.rotate2d.Matern.allp <- function(p, nu, grd, grd2){
  # p = c(theta, sigma2, tau2, ksi)
  # ksi = c( x, y, z, rotation_phi)
  #cat(ksi[1], ', ', ksi[2], ', ', ksi[3], ', ', ksi[4])
  cat(p[1], p[2], p[3], p[4], p[5], p[6], p[7] )
  cat(", \n")
  y <- grd[,3]
  g <- data.matrix( grd[,1:2])
  names(g) <- NULL

  y2 <- grd2[,3] + p[6]
  phi <- p[7]
  rotat <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), nr=2, nc=2 )
  g2 <- t( rotat %*% t(cbind( grd2[,1] + p[4], grd2[,2] + p[5] )) )
  d <- rdist(rbind(g,g2))
  SS <- exp(p[2]) * Matern(d, range=exp(p[1]), smoothness=nu)
  W <- exp(p[3]) * diag(dim(SS)[1])
  C <- SS + W
  #Q <- solve(C)
  yy <- c(y,y2)
  C.c <- t(chol(C))
  out <- forwardsolve(C.c, yy)
  quad.form <- sum(out^2)
  det.part <- 2*sum(log(diag(C.c)))

  # up to the normalizing constants
  NLL <- 0.5*det.part + 0.5*quad.form
  #NLL <- dim(SS)[1]/2*log(2*pi) + 0.5*sum( log(eigen(C)$values)) +
  #    0.5*t(yy) %*% Q %*% yy
  cat(NLL)
  cat(", \n")
  return(NLL)
}
