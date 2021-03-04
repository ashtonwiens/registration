#' Registers two data sets using Gaussian process nonrigid registration
#'
#' This function registers using two data sets using the nonrigid registration method developed in
#' \cite{}. The two data sets are registered as follows: first, locally estimate transformation parameters using a
#' Gaussian process rigid registration model developed in \cite{} and implemented in the function
#' \code{gp.rigid.registration}; then fit spatial model(s) to the local estimates using \code{spatialProcess}
#' or \code{Tps}; finally predict the transformation parameters using these models at the two
#' dimensional observation locations in data2, the moving data set.
#'
#' Details: est.pars is a list containing n.subsamples (the number of subsamples to be taken from each data set
#' in each window), smoothness (of the Matern covariance function in the local estimation procedure),
#' lambd (value of tuning parameter in the penalty in the objective function, can be a vector of length four
#' specifying a different multiplicative tuning parameter for each transformation parameter),
#' penalty.type (one of 'L1', 'L2', 'elastic', or 'none'), upper.bounds and lower.bounds (box constraints
#' in the local likelihood parameter estimation), min.num.samples (the minimum number of samples needed
#' in both data sets to perform local estimation on a window of data), and est.type (a character vector
#' 'serial' or 'Rmpi' indicating how the estimation should be performed).
#'
#' The setup parameters list \code{setup.pars} includes the number of centroids in the grid in the x and y
#' dimensions (\code{nxx} and \code{nyy}), an ovrlp parameter indicating how much each window overlaps
#' (<0.5 indicates the windows are mutually exclusive, 0.5 indicates overlap only along the boundaries,
#' and >0.5 indicates overlap among adjacent windows), and a sqz (or inset) parameter which insets the
#' grid of centroids by (domainSize in one dimension)/sqz on each edge of the grid of centroids.
#'
#' @param pd1 either a character vector path pointing the a delimited file whose first three columns are 3d
#' locations of data set 1, a data frame containing data set 1 with the former specifications, or a
#' gpnl.obj
#' @param pd2 character vector path pointing the a delimited file whose first three columns are 3d
#' locations of data set 2, or a data frame containing data set 2
#' @param setup.pars a list containing components \code{nxx}, \code{nyy}, \code{ovrlp}, and \code{sqz}
#' @param est.pars a list containing n.subsamples, smoothness, lambd, upper.bounds and lower.bounds,
#' min.num.samples, and est.type
#' @return a gpnl.obj with data.list, data.subsets, est.pars, setup.pars, est.results, and fit.results
#' components
#' @examples
#' library(dplyr)
#' library(fields)
#' library(data.table)
#' data(largeCliffs)
#' d1 <- largeCliffs$data1
#' d2 <- largeCliffs$data2
#' setup.pars <- list(nxx = 20, nyy = 20, ovrlp = 1.9, sqz = 25)
#' trns.bwd <- 1
#' angl.bwd <- pi / 4
#' est.pars <- list(
#'   n.subsamples = 100,
#'   smoothness = 1,
#'   lambd = 20, # c(50,50,50,25)
#'   penalty.type = "L2",
#'   lower.bounds = c(-Inf, -Inf, -Inf, -trns.bwd, -trns.bwd, -trns.bwd, -angl.bwd + 0.1),
#'   upper.bounds = c(4, Inf, Inf, trns.bwd, trns.bwd, trns.bwd, angl.bwd),
#'   est.type = "serial",
#'   tps.lambda = 0.1,
#'   krig.theta = 1
#' ) # 'serial', 'Rmpi', 'parallel'
#' gpnl.obj <- gpnlreg.setup(d1, d2, setup.pars)
#' gpnl.obj <- gpnlreg.estim.def(gpnl.obj, est.pars)
#' gpnl.obj <- gpnlreg(gpnl.obj, setup.pars = setup.pars, est.pars = est.pars)
#' ## or all in one step:
#' gpnl.obj <- gpnlreg(d1, d2, setup.pars, est.pars)
#' ## now plot for diagnostics
#' library(rgl)
#' plot.data.subset(gpnl.obj, 40)
#' plot.windows.i(gpnl.obj, c(40, 41, 42))
#' plot.param(gpnl.obj, "x")
#' plot.param(gpnl.obj, "y")
#' plot.param(gpnl.obj, "z")
#' plot.param(gpnl.obj, "phi")
#' plot.data(gpnl.obj$data.list$data1, 100, 100)
#' plot.data(gpnl.obj$data.list$data2, 100, 100, add = TRUE)
#' plot.data(gpnl.obj$data.list$data2.trf, 100, 100, add = TRUE)
#' plot.data.3d(gpnl.obj$data.list$data1)
#' plot.data.3d(gpnl.obj$data.list$data2, add = TRUE, col = "blue")
#' plot.data.3d(gpnl.obj$data.list$data2.trf, add = TRUE, col = "green")
#' fields::quilt.plot(total.change[, 1:2], abs(gpnl.obj$fit.results$trnsf.That[, 1]))
#' fields::quilt.plot(total.change[, 1:2], abs(gpnl.obj$fit.results$trnsf.That[, 2]))
#' fields::quilt.plot(total.change[, 1:2], abs(gpnl.obj$fit.results$trnsf.That[, 3]))
#' fields::quilt.plot(total.change[, 1:2], abs(gpnl.obj$fit.results$trnsf.That[, 4]))
#' rs <- apply(gpnl.obj$fit.results$trnsf.That, 1, function(x) {
#'   sum(abs(x))
#' })
#' total.change <- cbind(data2[, 1:2], rs)
#' fields::quilt.plot(total.change[, 1:2], total.change[, 3])
#' @export
gpnlreg <- gp.nonrigid.registration <- nonrigid.registration <- function(pd1, pd2 = NULL, setup.pars, est.pars) {
  if (!is.data.frame(pd1) & !is.character(pd1)) {
    gpnl.obj <- pd1
  } else {
    gpnl.obj <- list()
  }
  if (is.character(pd1) & is.character(pd2)) {
    gpnl.obj <- gpnlreg.load.ssv(pd1, pd2)
  }
  if (is.data.frame(pd1) & is.data.frame(pd2)) { # need to check column names
    data.list <- list(data1 = pd1, data2 = pd2)
    gpnl.obj <- list()
    gpnl.obj$data.list <- data.list
  }
  if (is.null(gpnl.obj$data.list$sDomain)) {
    sx <- range(rbind(range(gpnl.obj$data.list$data1[, 1]), range(gpnl.obj$data.list$data2[, 1])))
    sy <- range(rbind(range(gpnl.obj$data.list$data1[, 2]), range(gpnl.obj$data.list$data2[, 2])))
    sDomain <- list(sx = sx, sy = sy)
    gpnl.obj$data.list$sDomain <- sDomain
  }
  colnames(gpnl.obj$data.list$data1)[1:3] <- c("x", "y", "z")
  colnames(gpnl.obj$data.list$data2)[1:3] <- c("x", "y", "z")

  if (is.null(gpnl.obj$data.subsets)) {
    gpnl.obj <- gpnlreg.setup(gpnl.obj, d2 = NULL, setup.pars)
    gpnl.obj$setup.pars <- setup.pars
  }
  if (is.null(gpnl.obj$est.results)) {
    gpnl.obj <- gpnlreg.estim.def(gpnl.obj, est.pars)
    gpnl.obj$est.pars <- est.pars
  }
  if (is.null(gpnl.obj$fit.results)) {
    print("fitting")
    gpnl.obj <- gpnlreg.krig.local.pars(gpnl.obj)
  }
  gpnl.obj <- gpnlreg.apply.transformation(gpnl.obj)
  return(gpnl.obj)
}

#' Compute the log likelihood function, which may include a penalty term on the transformation parameters
#'
#' Used as the objective function in a call to \code{optim}
#'
#' @param p covariance parameters (range, sill and nugget variances) and transformation parameters
#' (x, y, z translations and 2d rotation) to be optimized over
#' @param nu fixed smoothness parameter of the Matern covariance function
#' @param grd three column data set 1
#' @param grd2 three column data set 2
#' @param lambda multiplicative tuning parameter in the penalty term
#' @param pen.type type of penalty term to apply ('L1', 'L2', 'elastic', or 'none')
#' @return the value of the penalized likelihood function
#' @export
logLik.deformation.penalized <- function(p, nu, grd, grd2, lambda, int = NULL, kappa = NULL, pen.type = "L2") {
  # p = c(theta, sigma2, tau2, ksi)
  # ksi = c( x, y, z, rotation_phi)
  # cat(ksi[1], ', ', ksi[2], ', ', ksi[3], ', ', ksi[4])
  # cat(p[1], ', ', p[2], ', ', p[3], p[4], p[5], p[6], p[7] )
  # cat(", \n")
  y <- grd[, 3]
  g <- data.matrix(grd[, 1:2])
  names(g) <- NULL

  y2 <- grd2[, 3] + p[6]
  phi <- p[7]
  rotat <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), nrow = 2, ncol = 2)
  g2 <- t(rotat %*% t(cbind(grd2[, 1] + p[4], grd2[, 2] + p[5])))
  d <- fields::rdist(rbind(g, g2))
  SS <- exp(p[2]) * fields::Matern(d, range = exp(p[1]), smoothness = nu)
  W <- exp(p[3]) * diag(dim(SS)[1])
  C <- SS + W
  # Q <- solve(C)
  yy <- c(y, y2)
  C.c <- t(chol(C))
  out <- forwardsolve(C.c, yy)
  quad.form <- sum(out^2)
  det.part <- 2 * sum(log(diag(C.c)))

  if (!is.null(int)) {
    intk <- int[4]
    int <- int[1:3]
    # print(intk)
  } else {
    int <- c(0, 0, 0)
    inkt <- c(0)
  }

  if (pen.type == "L2") {
    pen <- 0.5 * lambda * sum((p[4:6] - int)^2) # lambda in sum?
  }
  if (pen.type == "L1") {
    pen <- lambda / 2 * sum(abs(p[4:6] - int))
  }
  if (pen.type == "elastic") {
    pen <- 0.5 * (sum(lambda * abs(p[4:6] - int)) + sum(lambda * (p[4:6] - int)^2)) # ?
  }
  if (pen.type == "none") {
    pen <- 0
  }
  if (!is.null(kappa)) {
    pen <- pen - kappa * cos(p[7]) + log(2 * pi * besselI(abs(kappa - intk), nu = 0))
  }
  # cat(p[4], p[5], p[6], p[7] )
  # up to the normalizing constants
  NLL <- 0.5 * det.part + 0.5 * quad.form + 0.5 * dim(SS)[1] * log(2 * pi) + pen
  # NLL <- dim(SS)[1]/2*log(2*pi) + 0.5*sum( log(eigen(C)$values)) +
  #    0.5*t(yy) %*% Q %*% yy
  # cat('LL',0.5*det.part + 0.5*quad.form + 0.5*dim(SS)[1]*log(2*pi))
  # cat(0.5*det.part)
  # cat(", \n")
  # cat(0.5*quad.form)
  # cat('    penalty:', pen )
  # cat(", \n")
  return(NLL)
}

#' Compute the log likelihood function, which may include a penalty term on the transformation parameters
#'
#' Used as the objective function in a call to \code{optim}
#'
#' @param p transformation parameters (x, y, z translations and 2d rotation) to be optimized over
#' @param nu fixed smoothness parameter of the Matern covariance function
#' @param grd three column data set 1
#' @param grd2 three column data set 2
#' @param lambda multiplicative tuning parameter in the penalty term
#' @param pen.type type of penalty term to apply ('L1', 'L2', 'elastic', or 'none')
#' @param fixed.p covariance parameters (range, sill and nugget variances), fixed not estimated
#' @return the value of the penalized likelihood function
#' @export
logLik.deformation.penalized.fixed <- function(p, nu, grd, grd2,
                                               lambda, pen.type = "L2",
                                               fixed.p) { # range sill nugget
  # p = c(theta, sigma2, tau2, ksi)
  # ksi = c( x, y, z, rotation_phi)
  # cat(ksi[1], ', ', ksi[2], ', ', ksi[3], ', ', ksi[4])
  # cat(p[1], ', ', p[2], ', ', p[3], p[4], p[5], p[6], p[7] )
  # cat(", \n")
  y <- grd[, 3]
  g <- data.matrix(grd[, 1:2])
  names(g) <- NULL

  y2 <- grd2[, 3] #+ p[6]
  phi <- 0 # p[7]
  rotat <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), nrow = 2, ncol = 2)
  g2 <- t(rotat %*% t(cbind(grd2[, 1] + p[1], grd2[, 2])))
  d <- fields::rdist(rbind(g, g2))
  SS <- fixed.p[2] * fields::Matern(d, range = fixed.p[1], smoothness = nu)
  W <- fixed.p[3] * diag(dim(SS)[1])
  C <- SS + W
  # Q <- solve(C)
  yy <- c(y, y2)
  C.c <- t(chol(C))
  out <- forwardsolve(C.c, yy)
  quad.form <- sum(out^2)
  det.part <- 2 * sum(log(diag(C.c)))

  if (pen.type == "L2") {
    pen <- sum(lambda * p^2)
  }
  if (pen.type == "L1") {
    pen <- sum(lambda * abs(p))
  }
  if (pen.type == "elastic") {
    pen <- 0.5 * (sum(lambda * abs(p)) + sum(lambda * p^2))
  }
  if (pen.type == "none") {
    pen <- 0
  }
  # cat(p[4], p[5], p[6], p[7] )
  # up to the normalizing constants
  NLL <- 0.5 * det.part + 0.5 * quad.form + 0.5 * dim(SS)[1] * log(2 * pi) + pen
  # NLL <- dim(SS)[1]/2*log(2*pi) + 0.5*sum( log(eigen(C)$values)) +
  #    0.5*t(yy) %*% Q %*% yy
  # cat('LL',0.5*det.part + 0.5*quad.form + 0.5*dim(SS)[1]*log(2*pi))
  # cat(0.5*det.part)
  # cat(", \n")
  # cat(0.5*quad.form)
  # cat('    penalty:', pen )
  # cat(", \n")
  return(NLL)
}

#' Compute the log likelihood function for a Gaussian process with Matern covariance
#'
#' Used as the objective function in a call to \code{optim}
#'
#' @param p covariance parameters (range, sill and nugget variances) to be optimized over
#' @param nu fixed smoothness parameter of the Matern covariance function
#' @param grd three column data set
#' @return the value of the negative log likelihood function
#' @export
logLikMatern <- function(p, nu, grd) {
  # p = c(theta, sigma2, tau2, ksi)
  # ksi = c( x, y, z, rotation_phi)
  # cat(ksi[1], ', ', ksi[2], ', ', ksi[3], ', ', ksi[4])
  # cat(p[1], ', ', p[2], ', ', p[3], p[4], p[5], p[6], p[7] )
  # cat(", \n")
  y <- grd[, 3]
  g <- data.matrix(grd[, 1:2])
  names(g) <- NULL
  d <- fields::rdist(g)
  SS <- exp(p[2]) * fields::Matern(d, range = exp(p[1]), smoothness = nu)
  W <- exp(p[3]) * diag(dim(SS)[1])
  C <- SS + W
  # Q <- solve(C)
  yy <- c(y)
  C.c <- t(chol(C))
  out <- forwardsolve(C.c, yy)
  quad.form <- sum(out^2)
  det.part <- 2 * sum(log(diag(C.c)))

  # cat(p[4], p[5], p[6], p[7] )
  # up to the normalizing constants
  NLL <- 0.5 * det.part + 0.5 * quad.form + 0.5 * dim(SS)[1] * log(2 * pi)
  # NLL <- dim(SS)[1]/2*log(2*pi) + 0.5*sum( log(eigen(C)$values)) +
  #    0.5*t(yy) %*% Q %*% yy
  # cat('LL',0.5*det.part + 0.5*quad.form + 0.5*dim(SS)[1]*log(2*pi))
  # cat(0.5*det.part)
  # cat(", \n")
  # cat(0.5*quad.form)
  # cat('    penalty:', pen )
  # cat(", \n")
  return(NLL)
}

#' Split the data into subsets based on inclusion in a grid box centered around each centroid
#'
#' This function creates a grid of centroids of size nx*ny, then creates square boundaries centered on
#' each centroid, and finally creates subsets of data based on inclusion in each window.
#'
#' @param data1 three column data set 1
#' @param data2 three column data set 2
#' @param nx number of centroids in the x dimension
#' @param ny number of centroids in the y dimension
#' @param overlap a positive real number indicating how much overlap among adjacent windows
#' @param squeeze a positve real number which sets the inset of the outermost centroids
#' @export
split.data.overlap <- function(data1, data2, nx, ny, overlap = 0.6, squeeze = 10) {
  # library(data.table)
  sx <- range(rbind(range(data1[, 1]), range(data2[, 1])))
  sy <- range(rbind(range(data1[, 2]), range(data2[, 2])))
  squeeze <- nx + 1
  rx <- (sx[2] - sx[1]) / squeeze
  squeeze <- ny + 1
  ry <- (sy[2] - sy[1]) / squeeze
  num.centroids.x <- nx
  num.centroids.y <- ny
  seq.x <- seq(sx[1] + rx, sx[2] - rx, length.out = num.centroids.x)
  seq.y <- seq(sy[1] + ry, sy[2] - ry, length.out = num.centroids.y)
  dx <- seq.x[2] - seq.x[1]
  dy <- seq.y[2] - seq.y[1]
  centroids <- expand.grid(seq.x, seq.y)
  centroids <- dplyr::mutate(centroids, grid = 1:nrow(centroids), x = Var1, y = Var2, Var1 = NULL, Var2 = NULL)
  centroids <- data.table::setDT(centroids)
  overlap.factor <- overlap # 0.625 # 1/2 for no overlap, must be > 1/2 or else deficient (unused data)
  bounds <- centroids[, .(
    xl = x - (dx * overlap.factor), xu = x + (dx * overlap.factor),
    yl = y - (dy * overlap.factor), yu = y + (dy * overlap.factor),
    grid
  )]

  D1_i <- D2_i <- list()
  data.size <- matrix(0, nrow = nrow(centroids), ncol = 2)
  colnames(data1) <- colnames(data2) <- c("x", "y", "z")
  for (i in 1:nrow(bounds)) {
    D1_i[[i]] <- dplyr::filter(
      data1, data.table::between(data1$x, bounds$xl[i], bounds$xu[i]),
      data.table::between(data1$y, bounds$yl[i], bounds$yu[i])
    ) # can also be bounds$grid[i]
    D2_i[[i]] <- dplyr::filter(
      data2, data.table::between(data2$x, bounds$xl[i], bounds$xu[i]),
      data.table::between(data2$y, bounds$yl[i], bounds$yu[i])
    )
    # D1_i[[i]] <- dplyr::filter(data1, x >= bounds$xl[i], x <= bounds$xu[i],  y <= bounds$yu[i], y >= bounds$yl[i] ) # i can also be bounds$grid[i]
    # D2_i[[i]] <- dplyr::filter(data2, x >= bounds$xl[i], x <= bounds$xu[i],  y <= bounds$yu[i], y >= bounds$yl[i] )
    data.size[i, 1] <- nrow(D1_i[[i]])
    data.size[i, 2] <- nrow(D2_i[[i]])
  }
  # plot(D_i[[1]][,2:3])
  return(list(D1_i = D1_i, D2_i = D2_i, centroids = centroids, data.size = data.size, bounds = bounds))
}


#' Wrapper for \code{\link{load.data}}
#'
#' \code{gpnlreg.load} takes two paths to delimited files and loads them as data frames using \code{load.data}.
#'
#' @param pth1 a character vector with a path to a delimited file containing the first data set in the first three columns
#' @param pth2 a character vector with a path to a delimited file containing the second data set in the first three columns
#' @param delim a single vector indicating how the data are delimited
#' @return a gpnl.obj object with component gpnl.obj$data.list
#' @seealso \code{\link{load.data}}
#' @export
gpnlreg.load <- gpnlreg.load.ssv <- function(pth1, pth2, delim = " ") {
  gpnl.obj <- list()
  data.list <- load.data(pth1, pth2, delim = delim)
  gpnl.obj$data.list <- data.list
  return(gpnl.obj)
}

#' Splits data into subsets using setup.pars
#'
#' \code{gpnlreg.setup} takes a gpnl.obj and a list of setup parameters, splits the data to be registered
#' into a list of subsets within each window, and returns a gpnl.obj with a data.subsets component.
#'
#' @param gpnl.obj a gpnl.obj containing a data.list
#' @param setup.pars a list containing components \code{nxx}, \code{nyy}, \code{ovrlp}, and \code{sqz}
#' @return a gpnl.obj containing a data.subsets component
#' @seealso \code{\link{split.data.overlap}}
#' @export
gpnlreg.setup <- function(gpnl.obj, d2 = NULL, setup.pars) {
  # setup.pars <- list(nxx=nxx, nyy=nyy, ovrlp=ovrlp, sqz=sqz)
  if (!is.null(d2)) {
    data1 <- gpnl.obj
    data2 <- d2
    data.list <- list(data1 = data1, data2 = data2)
    gpnl.obj <- list()
    gpnl.obj$data.list <- data.list
    sx <- range(rbind(range(gpnl.obj$data.list$data1[, 1]), range(gpnl.obj$data.list$data2[, 1])))
    sy <- range(rbind(range(gpnl.obj$data.list$data1[, 2]), range(gpnl.obj$data.list$data2[, 2])))
    sDomain <- list(sx = sx, sy = sy)
    gpnl.obj$data.list$sDomain <- sDomain
    colnames(gpnl.obj$data.list$data1)[1:3] <- c("x", "y", "z")
    colnames(gpnl.obj$data.list$data2)[1:3] <- c("x", "y", "z")
  } else {
    data1 <- gpnl.obj$data.list$data1
    data2 <- gpnl.obj$data.list$data2
    colnames(data1)[1:3] <- colnames(data2)[1:3] <- c("x", "y", "z")
  }
  start.time.split <- Sys.time()
  data.subsets <- split.data.overlap(data1, data2, # split data
    nx = setup.pars$nxx, ny = setup.pars$nyy,
    overlap = setup.pars$ovrlp, squeeze = setup.pars$sqz
  )
  end.time.split <- Sys.time()
  split.time <- end.time.split - start.time.split
  cat("Time to subset data: ", split.time, " \n")
  N <- setup.pars$nxx * setup.pars$nyy
  setup.pars$N <- N
  gpnl.obj$setup.pars <- setup.pars
  gpnl.obj$data.subsets <- data.subsets
  return(gpnl.obj)
}

#' Locally estimates Matern and rigid transformation parameters in serial or parallel
#'
#' \code{gpnlreg.estim.def} takes a gpnl.obj and a list of estimation and fitting parameters,
#' and using either a serial or an \code{Rmpi} implementation, calls \code{parallel.fcn} on every
#' window of data in gpnl.obj$data.subsets, returning a gpnl.obj with an est.results component
#' including all of the results of the estimation and model fitting.
#'
#' @param gpnl.obj a gpnl.obj containing a data.subsets component
#' @param est.pars a list containing n.subsamples (the number of subsamples to be taken from each data set
#' in each window), smoothness (of the Matern covariance function in the local estimation procedure),
#' lambd (value of tuning parameter in the penalty in the objective function), penalty.type (one of 'L1',
#' 'L2', 'elastic', or 'none'), upper.bounds and lower.bounds (box constraints in the local likelihood
#' parameter estimation), min.num.samples (the minimum number of samples needed in both data sets to
#' perform local estimation on a window of data), and est.type (a character vector 'serial' or 'Rmpi'
#' indicating how the estimation should be performed).
#' @return a gpnl.obj containing an est.results component
#' @seealso \code{\link{serial.implementation}} \code{\link{rmpi.implementation}} \code{\link{parallel.fcn}}
#' @export
gpnlreg.estim.def <- function(gpnl.obj, est.pars) {
  # est.pars <- list(n.subsamples=200, smoothness=1,
  #                  lambd=5, #c(50,50,50,25)
  #                 penalty.type='L2',
  #                  lower.bounds=c(-Inf,-Inf,-Inf, -trns.bwd, -trns.bwd, -trns.bwd, -angl.bwd+0.1),
  #                   upper.bounds=c(4, Inf, Inf, trns.bwd, trns.bwd, trns.bwd, angl.bwd),
  #                  min.num.samples=30,
  #                   est.type='Rmpi')
  # est.type=c('serial','Rmpi', 'parallel')
  est.type <- est.pars$est.type
  gpnl.obj$est.pars <- est.pars
  if (is.null(gpnl.obj$est.pars)) {
    gpnl.obj$est.pars <- est.pars
  }
  if (est.type == "serial") {
    est.results <- serial.implementation(gpnl.obj, est.pars)
  }
  if (est.type == "parallel") {

  }
  if (est.type == "Rmpi") {
    est.results <- rmpi.implementation(gpnl.obj, est.pars)
  }
  gpnl.obj$est.results <- est.results
  return(gpnl.obj)
}

#' Takes locally estimated transformation parameters and fits surface models independently
#'
#' \code{gpnlreg.krig.local.pars} takes a gpnl.obj with est.results, and fits a spatial model to each
#' transformation parameter independently using a Gaussian process or an approximate thin plate spline.
#' Then, using these surfaces models, the transformation parameters are predicted at the observation locations
#' of the moving data set (\code{data2}). \code{gpnl.obj$fit.results} contains the fitted models and the
#' predicted transformation parameters.
#'
#' @param gpnl.obj a gpnl.obj containing data.subsets, est.pars, and est.results components
#' @return a gpnl.obj containing an fit.results component
#' @seealso \code{\link{spatialProcess}} \code{\link{Tps}}
#' @export
gpnlreg.krig.local.pars <- function(gpnl.obj) {
  D1_i <- gpnl.obj$data.subsets$D1_i
  D2_i <- gpnl.obj$data.subsets$D2_i
  centroids <- gpnl.obj$data.subsets$centroids[, -1]
  ds <- gpnl.obj$data.subsets$data.size # calculate size of subsets of data
  local.pars <- gpnl.obj$est.results$local.pars
  est.params <- matrix(0, nrow = nrow(centroids), ncol = 7)

  for (i in 1:length(local.pars)) { # extract parameter estimates from optimization results
    if (!is.null(local.pars[[i]]$par)) {
      if (length(local.pars[[i]]$par) == 1) {
        est.params[i, ] <- rep(0, 7)
        est.params[i, 4] <- local.pars[[i]]$par
      } else {
        est.params[i, ] <- local.pars[[i]]$par
      }
    } else {
      est.params[i, ] <- rep(0, 7)
    }
  }
  est.params[, 1:3] <- exp(est.params[, 1:3]) # Matern parameters > 0 so estimated in log scale; transforms back to original scale
  # wghts <- apply(ds, 1, mean)
  # condition <-  (ds[,1] > min.num.samples) & (ds[,2] > min.num.samples)
  # est.params.rf <- est.params[condition,]  #filter data sets by whether parameters were estimated
  est.params.rf <- est.params[!dplyr::near(est.params[, 4], 0), ] # filter data sets by whether parameters were estimated
  if (!is.matrix(est.params.rf)) {
    est.params.rf <- matrix(est.params.rf, ncol = 7)
  }
  # centroids.rf <- centroids[condition,1:2]
  centroids.rf <- centroids[!dplyr::near(est.params[, 4], 0), 1:2]
  # wghts.rf <- wghts[condition]
  data2.prd.loc <- gpnl.obj$data.list$data2[, 1:2]
  # utils::str(data2.prd.loc)
  # fit a Gaussian process to each spatial field independently
  GP.fit <- list()
  trnsf.That <- gpnl.obj$data.list$data2
  trnsf.SE <- gpnl.obj$data.list$data2
  # utils::str(trnsf.That)
  prd.seq <- seq(1, nrow(data2.prd.loc), by = 10000)
  if (length(prd.seq) == 1) {
    prd.seq <- c(1, nrow(data2.prd.loc))
  }
  prd.seq[length(prd.seq)] <- nrow(data2.prd.loc)
  fit.time <- system.time(
    for (i in 4:ncol(est.params.rf)) {
      param.response <- est.params.rf[, i]
      # print(dput(param.response))
      # print(dput(data.matrix(centroids.rf)))
      # GP.fit[[i-3]] <- mod <- spatialProcess(x = data.matrix(centroids.rf), y = param.response, #, weights = wghts,
      #                                        mKrig.args = list(m=2), theta = gpnl.obj$est.pars$krig.theta ,
      #                                        cov.args=list(Covariance="Matern", smoothness=1) )
      dm <- data.matrix(centroids.rf)
      mdl <- gpnl.obj$est.pars$mdl
      mdl.name <- gpnl.obj$est.pars$mdl.name
      mdl.pars <- gpnl.obj$est.pars$mdl.pars
      # if( length(mdl.pars)>0 ){
      if (mdl.name == "Tps") {
        print("Tps")
        # utils::str(data.matrix(centroids.rf))
        # utils::str(param.response)
        GP.fit[[i - 3]] <- mod <- fields::Tps(x = data.matrix(centroids.rf), Y = param.response)
        # GP.fit[[i-3]] <- mod <- do.call(mdl.name, c( x = data.matrix(centroids.rf), Y= param.response, mdl.pars))
        #          mKrig.args = list(m = 2))#, theta=10)
        #                                  theta=tht, lambda=lbd)
      }
      if (mdl.name == "LatticeKrig") {
        print("LK")
        GP.fit[[i - 3]] <- mod <- LatticeKrig::LatticeKrig(x = data.matrix(centroids.rf), y = param.response)
        #          mKrig.args = list(m = 2))#, theta=10)
        #                                  theta=tht, lambda=lbd)
      }
      # }else{
      #   GP.fit[[i-3]] <- mod <- mdl(x = data.matrix(centroids.rf),  param.response)#,
      #   #          mKrig.args = list(m = 2))#, theta=10)
      #   #                                  theta=tht, lambda=lbd)
      # }
      # GP.fit[[i-3]] <- mod <- spatialProcess( x=dm, y=param.response,
      #                                         mKrig.args = list(m = 2))
      # GP.fit[[i-3]] <- mod <- LatticeKrig::LatticeKrig(x=dm, y=param.response)

      #                                      theta=gpnl.obj$est.pars$krig.theta, lambda=gpnl.obj$est.pars$tps.lambda)
      mod <- GP.fit[[i - 3]]
      # print(length(prd.seq))
      for (j in 1:(length(prd.seq) - 1)) {
        print(j)
        k1 <- prd.seq[j]
        k2 <- prd.seq[j + 1] - 1
        if (j == (length(prd.seq) - 1)) {
          k2 <- k2 + 1
        }
        # print(k1)
        # print(k2)
        # utils::str(data2.prd.loc[k1:k2,])
        # utils::str(predict(mod, xnew = data.matrix(data2.prd.loc[k1:k2,]) ))
        # utils::str(trnsf.That[k1:k2,i-3])
        if (mdl.name == "LatticeKrig") {
          trnsf.That[k1:k2, i - 3] <- LatticeKrig::predict.LKrig(mod, data2.prd.loc[k1:k2, ])
          trnsf.SE[k1:k2, i - 3] <- LatticeKrig::predictSE.LKrig(mod, data2.prd.loc[k1:k2, ])
        }
        if (mdl.name == "Tps") {
          trnsf.That[k1:k2, i - 3] <- fields::predict.Tps(mod, data2.prd.loc[k1:k2, ])
          trnsf.SE[k1:k2, i - 3] <- fields::predictSE(mod, data2.prd.loc[k1:k2, ])
        }
      }
      # trnsf.That[[i-3]] <- predict(mod, xnew = data2.prd.loc )
      # print(i)
    }
  )
  cat("Time to fit GP and predict transformation: ", fit.time[3], " \n")
  fit.results <- list(
    GP.fit = GP.fit, trnsf.That = trnsf.That, trnsf.SE = trnsf.SE,
    fit.time = fit.time,
    est.params.rf = est.params.rf, est.params = est.params,
    centroids.rf = centroids.rf, centroids = centroids
  )
  gpnl.obj$fit.results <- fit.results
  return(gpnl.obj)
}


#' Takes locally estimated covariance parameters and fits surface models independently
#'
#' \code{gpnlreg.krig.cov.pars} takes a gpnl.obj with est.results, and fits a spatial model to each
#' covariance parameter independently using a Gaussian process or an approximate thin plate spline.
#' Then, using these surfaces models, the transformation parameters are predicted at the observation locations
#' of the moving data set (\code{data2}). \code{gpnl.obj$fit.results} contains the fitted models and the
#' predicted covariance parameters.
#'
#' @param gpnl.obj a gpnl.obj containing data.subsets, est.pars, and est.results components
#' @return a gpnl.obj containing an fit.results component
#' @seealso \code{\link{spatialProcess}} \code{\link{Tps}}
#' @export
gpnlreg.krig.cov.pars <- function(gpnl.obj) {
  D1_i <- gpnl.obj$data.subsets$D1_i
  D2_i <- gpnl.obj$data.subsets$D2_i
  centroids <- gpnl.obj$data.subsets$centroids[, -1]
  ds <- gpnl.obj$data.subsets$data.size # calculate size of subsets of data
  local.pars <- gpnl.obj$est.results$local.pars
  est.params <- matrix(0, nrow = nrow(centroids), ncol = 7)

  for (i in 1:length(local.pars)) { # extract parameter estimates from optimization results
    if (!is.null(local.pars[[i]]$par)) {
      est.params[i, ] <- local.pars[[i]]$par
    } else {
      est.params[i, ] <- rep(0, 7)
    }
  }
  est.params[, 1:3] <- est.params[, 1:3] # Matern parameters > 0 so estimated in log scale; transforms back to original scale
  # wghts <- apply(ds, 1, mean)
  # condition <-  (ds[,1] > min.num.samples) & (ds[,2] > min.num.samples)
  # est.params.rf <- est.params[condition,]  #filter data sets by whether parameters were estimated
  est.params.rf <- est.params[!dplyr::near(est.params[, 4], 0), ] # &
  # !dplyr::near(est.params[,1], 4)
  # & est.params[,1]<=log(20),]  #filter data sets by whether parameters were estimated
  if (!is.matrix(est.params.rf)) {
    est.params.rf <- matrix(est.params.rf, ncol = 7)
  }
  # centroids.rf <- centroids[condition,1:2]
  centroids.rf <- centroids[!dplyr::near(est.params[, 4], 0), 1:2] # & !dplyr::near(est.params[,1],4)
  #  & est.params[,1]<=log(20)
  # wghts.rf <- wghts[condition]
  data2.prd.loc <- gpnl.obj$data.list$d1cv[, 1:2]
  # fit a Gaussian process to each spatial field independently
  GP.fit <- list()
  trnsf.That <- gpnl.obj$data.list$d1cv
  trnsf.SE <- gpnl.obj$data.list$d1cv

  prd.seq <- seq(1, nrow(data2.prd.loc), by = 10000)
  if (length(prd.seq) == 1) {
    prd.seq <- c(1, nrow(data2.prd.loc))
  }
  prd.seq[length(prd.seq)] <- nrow(data2.prd.loc)
  # print(prd.seq)
  fit.time <- system.time(
    for (i in 1:3) {
      param.response <- est.params.rf[, i]
      # print(utils::str(param.response))
      # print(utils::str(data.matrix(centroids.rf)))

      # GP.fit[[i]] <- mod <- spatialProcess(x = data.matrix(centroids.rf), y = param.response, #, weights = wghts,
      #                                        mKrig.args = list(m=2), theta = gpnl.obj$est.pars$krig.theta ,
      #                                        cov.args=list(Covariance="Matern", smoothness=1) )
      mdl <- gpnl.obj$est.pars$mdl
      mdl.name <- gpnl.obj$est.pars$mdl.name
      mdl.pars <- gpnl.obj$est.pars$mdl.pars
      # if( length(mdl.pars)>0 ){
      if (mdl.name == "Tps") {
        print("Tps")
        GP.fit[[i]] <- mod <- fields::Tps(x = data.matrix(centroids.rf), Y = param.response)
        # GP.fit[[i-3]] <- mod <- do.call(mdl.name, c(x = data.matrix(centroids.rf), Y= param.response, mdl.pars))
        #          mKrig.args = list(m = 2))#, theta=10)
        #                                  theta=tht, lambda=lbd)
      }
      if (mdl.name == "LatticeKrig") {
        print("LK")
        GP.fit[[i]] <- mod <- LatticeKrig::LatticeKrig(x = data.matrix(centroids.rf), y = param.response)
        #          mKrig.args = list(m = 2))#, theta=10)
        #                                  theta=tht, lambda=lbd)
      }
      # }else{
      #   GP.fit[[i]] <- mod <- mdl(x = data.matrix(centroids.rf),  param.response)#,
      #   #          mKrig.args = list(m = 2))#, theta=10)
      #   #                                  theta=tht, lambda=lbd)
      # }

      mod <- GP.fit[[i]]
      for (j in 1:(length(prd.seq) - 1)) {
        k1 <- prd.seq[j]
        k2 <- prd.seq[j + 1] - 1
        if (j == (length(prd.seq) - 1)) {
          k2 <- k2 + 1
        }
        # utils::str(trnsf.SE)
        # utils::str(trnsf.That)
        utils::str(predict(mod, data2.prd.loc[k1:k2, ]))
        # utils::str( data2.prd.loc[k1:k2,] )
        # print(j)
        trnsf.That[k1:k2, i] <- predict(mod, data2.prd.loc[k1:k2, ]) # xnew=
        trnsf.SE[k1:k2, i] <- predictSE(mod, data2.prd.loc[k1:k2, ])
        # print(j)
      }
      # trnsf.That[[i-3]] <- predict(mod, xnew = data2.prd.loc )
      # print(i)
    }
  )
  cat("Time to fit GP and predict transformation: ", fit.time[3], " \n")
  cov.fit.results <- list(
    GP.fit = GP.fit, trnsf.That = trnsf.That, fit.time = fit.time,
    est.params.rf = est.params.rf, est.params = est.params,
    centroids.rf = centroids.rf, centroids = centroids
  )
  gpnl.obj$cov.fit.results <- cov.fit.results
  return(gpnl.obj)
}

#' Applies transformation parameters to the moving data set
#'
#' \code{gpnlreg.apply.transformation} applies \code{gpnl.obj$fit.results$trnsf.That} to
#' \code{gpnl.obj$data.list$data2} (the moving data set).
#'
#' @param gpnl.obj a gpnl.obj with a data.list and fit.results component
#' @return a gpnl.obj with a data2.trf data set in the data.list component
#' @export
gpnlreg.apply.transformation <- function(gpnl.obj) {
  data2.trf <- data2 <- gpnl.obj$data.list$data2
  trnsf.That <- gpnl.obj$fit.results$trnsf.That
  # trf.time <- system.time(
  #   for(i in 1:nrow(data2)){
  #     data2.trf[i,4] <- data2[i,4] + trnsf.That[[3]][i] # z translation
  #     phi <- trnsf.That[[4]][i]
  #     rotat <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), nrow=2, ncol=2 )
  #     data2.trf[i,2:3] <- t( rotat %*% t(cbind( data2[i,2] + trnsf.That[[1]][i], data2[i,3] + trnsf.That[[2]][i] )) )
  #     print(i)
  #   }
  # )
  s.t <- Sys.time()
  data2.trf[, 3] <- data2[, 3] + trnsf.That[[3]] # z translation
  # tmp.x <- data2[,1] + trnsf.That[[1]]
  # tmp.y <- data2[,2] + trnsf.That[[2]]
  cs <- cos(trnsf.That[[4]])
  sn <- sin(trnsf.That[[4]])
  # data2.trf[,1] <- cs * tmp.x - sn * tmp.y
  # data2.trf[,2] <- sn * tmp.x + cs * tmp.y
  data2.trf[, 1] <- (data2[, 1] * cs + data2[, 2] * sn) + trnsf.That[[1]]
  data2.trf[, 2] <- (data2[, 1] * -sn + data2[, 2] * cs) + trnsf.That[[2]]
  s.t2 <- Sys.time()
  trf.time <- s.t2 - s.t
  cat("Time to apply transformation: ", trf.time, " \n")
  gpnl.obj$data.list$data2.trf <- data2.trf
  gpnl.obj$fit.results$prd.time <- trf.time
  return(gpnl.obj)
}

#' Use fitted models in \code{fit.results} to predict at locations in \code{pred}.
#'
#' The models which were fitted to the locally estimated transformation parameters in fit.results are used
#' to predict the transformation parameters at the two dimensional locations in the first two columns of
#' \code{pred}.
#'
#' @param gpnl.obj a gpnl.obj containing fit.results
#' @param pred a three column matrix containing xyz locations of data at which the
#' transformation parameters will be predicted at and to which the predicted transformatin parameters
#' will be applied to
#' @return a three column matrix containing the transformed coordinates of \code{pred}
#' @seealso \code{\link{gpnlreg.apply.transformation}} \code{\link{gpnlreg.krig.local.pars}}
#' @export
gpnlreg.predict <- function(gpnl.obj, pred) { # pred is a 3 column matrix with x, y, and z coordinates
  trnsf.That <- list()
  GP.fit <- gpnl.obj$fit.results$GP.fit
  for (i in 1:length(GP.fit)) { # order of lists GP.fit and trnsf.That: x, y, z, phi
    trnsf.That[[i]] <- predict(GP.fit[[i]], xnew = pred[, 1:2])
  }
  pred.trf <- pred
  s.t <- Sys.time()
  data2.trf[, 3] <- pred[, 3] + trnsf.That[[3]] # z translation
  tmp.x <- pred[, 1] + trnsf.That[[1]]
  tmp.y <- pred[, 2] + trnsf.That[[2]]
  cs <- cos(trnsf.That[[4]])
  sn <- sin(trnsf.That[[4]])
  data2.trf[, 1] <- cs * tmp.x - sn * tmp.y
  data2.trf[, 2] <- sn * tmp.x + cs * tmp.y
  s.t2 <- Sys.time()
  trf.time <- s.t2 - s.t
  cat("Time to predict transformation: ", trf.time[3], " \n")
  return(pred.trf)
}


#' Load two delimited files from two character vectors paths and coerce to data frames
#'
#' \code{load.data} takes two character vectors which are paths to two delimited files, reads the data,
#' takes the first three columns from each data set and coerces them into two data frames, and
#' finally creates an sDomain list with the spatial extent of the first two dimensions of both data sets.
#'
#' @param pth1 a character vector with a path to a delimited file containing the first data set in the first three columns
#' @param pth2 a character vector with a path to a delimited file containing the second data set in the first three columns
#' @param delim a single vector indicating how the data are delimited
#' @return a list which is added to the gpnl.obj as gpnl.obj$data.list
#' @seealso \code{\link{gpnlreg.load}}
#' @export
load.data <- function(pth1, pth2, delim) {
  ## reads in the data, coerces to data.frame
  data1 <- readr::read_delim(pth1, delim = delim, col_names = FALSE)
  data2 <- readr::read_delim(pth2, delim = delim, col_names = FALSE)
  data1 <- data.frame(data1)
  data2 <- data.frame(data2)
  ## Sets column names to 'x', 'y', and 'z',
  data1 <- data1[, 1:3]
  data2 <- data2[, 1:3]
  colnames(data1) <- colnames(data2) <- c("x", "y", "z")
  ## creates a list containing the 2d spatial domain
  sx <- range(rbind(range(data1$x), range(data2$x)))
  sy <- range(rbind(range(data1$y), range(data2$y)))
  sDomain <- list(sx = sx, sy = sy)
  return(list(data1 = data1, data2 = data2, sDomain = sDomain, path1 = pth1, path2 = pth2))
}

#' Apply a rigid registration tot data2 to align with the coordinate system of data1
#'
#' This function applied a rigid transformation, which is estimated as part of the Gaussian process
#' model into which it is embedded.
#'
#' @param data1 the fixed data set
#' @param data2 the moving data set
#' @param est.pars a list of parameters used in the likelihood estimation, including the fixed
#' \code{smoothness} of the Matern covariance function, the multiplicative penalty term tuning parameter
#' \code{lambd}, the penalty type \code{penalty.type}, and box constraints \code{lower.bounds} and
#' \code{upper.bounds}.
#' @return The optimal Gaussian process (and transformation) parameters
#' @export
rigid.registration <- gp.rigid.registration <- function(data1, data2, est.pars, init = NULL) { ### est.pars <- list(smoothness, lambda, penalty.type, lower.bounds, upper.bounds)
  D1_i <- data1
  D2_i <- data2
  d1i <- data.frame(D1_i)
  d2i <- data.frame(D2_i)
  x1 <- mean(d1i$x)
  y1 <- mean(d1i$y)
  z1 <- mean(d1i$z)
  cxyz <- c(x1, y1, z1)
  d1i <- dplyr::mutate(d1i, x = x - x1, y = y - y1, z = z - z1)
  d2i <- dplyr::mutate(d2i, x = x - x1, y = y - y1, z = z - z1)
  if (is.null(init)) {
    init <- c(log(3), log(1), log(0.01), 0, 0, 0, 0)
  }
  opt.res <- tryCatch(stats::optim(
    par = init,
    fn = logLik.deformation.penalized,
    grd = d1i[, 1:3], grd2 = d2i[, 1:3],
    nu = est.pars$smoothness, int = init,
    kappa = est.pars$kappa,
    lambda = est.pars$lambd,
    pen.type = est.pars$penalty.type,
    method = "L-BFGS-B",
    lower = est.pars$lower.bounds,
    upper = est.pars$upper.bounds,
    control = list(trace = 1, maxit = 500)
  ),
  error = function(e) {
    return(e)
  }
  )
  return(list(opt.res = opt.res, data1 = data1, data2 = data2, cxyz = cxyz, est.pars = est.pars))
}

#' Function for paralleled estimation of a Gaussian process rigid registration
#'
#' Using window i of data, subsamples the data and then performs estimation of a penalized Gaussian process
#' with Matern covariance function for two independent sets of data, where a rigid transformation is
#' embedded, applied to the second, moving data set. Each data set is coordinate-wise de-meaned by the
#' means of the first data set's coordinates.
#'
#' @param i an integer indicating which window/subset of data to apply the Gaussian process
#' rigid registration model to
#' @param gpnl.obj a gpnl.obj containing data.subsets and est.pars components
#' @return a list containing the optimization results, the timing results, and the means
#' of each coordinate used in de-meaning before optimization is performed
#' @seealso \code{\link{gp.rigid.registrataion}} \code{\link{gp.nonrigid.registration}}
#' @export
parallel.fcn <- function(i, gpnl.obj) { ### doTask, for parallel and Rmpi packages
  D1_i <- gpnl.obj$data.subsets$D1_i
  D2_i <- gpnl.obj$data.subsets$D2_i
  d1i <- data.frame(D1_i[[i]])
  d2i <- data.frame(D2_i[[i]])
  if (nrow(d1i) > gpnl.obj$est.pars$n.subsamples & nrow(d2i) > gpnl.obj$est.pars$n.subsamples) {
    # if(nrow(d1i)>30 & nrow(d2i)>30){
    x1 <- mean(d1i$x)
    y1 <- mean(d1i$y)
    z1 <- mean(d1i$z)
    cxyz <- c(x1, y1, z1)
    d1i <- dplyr::mutate(d1i, x = x - x1, y = y - y1, z = z - z1)
    d2i <- dplyr::mutate(d2i, x = x - x1, y = y - y1, z = z - z1)
    d1i <- dplyr::sample_n(d1i, gpnl.obj$est.pars$n.subsamples)
    d2i <- dplyr::sample_n(d2i, gpnl.obj$est.pars$n.subsamples)
    if (gpnl.obj$est.pars$fixed == TRUE) {
      print("fixed")
      init <- c(0)
      time.opt <- system.time(opt.res <- tryCatch(stats::optim(
        par = init,
        fn = logLik.deformation.penalized.fixed,
        fixed.p = data.gend$param.list$gen.p[2:4],
        grd = d1i[, 1:3], grd2 = d2i[, 1:3],
        nu = gpnl.obj$est.pars$smoothness,
        lambda = gpnl.obj$est.pars$lambd,
        pen.type = gpnl.obj$est.pars$penalty.type,
        method = "L-BFGS-B",
        lower = gpnl.obj$est.pars$lower.bounds,
        upper = gpnl.obj$est.pars$upper.bounds,
        hessian = TRUE,
        control = list(trace = 1, maxit = 1000)
      ),
      error = function(e) {
        return(NULL)
      }
      ))
    } else {
      print("full")
      init <- c(log(3), log(1), log(0.01), 0, 0, 0, 0)
      # init <- c(2.45049191710834, 1.49277821122548, -4.12095940249598,
      #            -1.01084936, -0.01962812, 0.22183046,0)
      #-0.875902411669938,  -0.0504028913352409, 0.205027609071423, -0.00227546051767845)
      time.opt <- system.time(opt.res <- tryCatch(stats::optim(
        par = init,
        fn = logLik.deformation.penalized,
        grd = d1i[, 1:3], grd2 = d2i[, 1:3],
        nu = gpnl.obj$est.pars$smoothness,
        lambda = gpnl.obj$est.pars$lambd,
        pen.type = gpnl.obj$est.pars$penalty.type,
        method = "L-BFGS-B",
        lower = gpnl.obj$est.pars$lower.bounds,
        upper = gpnl.obj$est.pars$upper.bounds,
        hessian = TRUE,
        control = list(trace = 1, maxit = 1000)
      ),
      error = function(e) {
        return(NULL)
      }
      ))
    }
  } else {
    time.opt <- system.time(opt.res <- NULL)
    cxyz <- rep(0, 3)
  }
  print(opt.res$par)
  return(list(opt.res = opt.res, time.opt = time.opt, cxyz = cxyz))
}

#' Local estimation implemented with the \code{Rmpi} package
#'
#' Set the number of workers to spawn
#'
#' @param gpnlobj a gpnl.obj with data.subsets and est.pars components
#' @param est.pars a list of estimation parameters (see \code{gp.nonrigid.registration})
#' @return a list with the optimization and timing results, and data set 1 means used in de-meaning
#' @seealso \code{\link{serial.implementation}}
#' @export
rmpi.implementation <- function(gpnl.obj, est.pars) {
  # options(echo=FALSE)
  start.time <- Sys.time()
  D1_i <- gpnl.obj$data.subsets$D1_i
  # D2_i <- gpnl.obj$data.subsets$D2_i
  # ns <- parallel::detectCores()
  ns <- 1 # (mpi.universe.size()-3)
  mpi.spawn.Rslaves(nslaves = ns)
  mpi.bcast.Robj2slave(logLik.deformation.penalized)
  # mpi.bcast.Robj2slave(smoothness)
  # mpi.bcast.Rfun2slave()
  namesLibraries <- c(
    "fields", "dplyr"
  )
  for (objName in namesLibraries) {
    mpi.bcast.cmd(
      cmd = "library",
      package = objName,
      # lib.loc = namesLibraryLocations,
      character.only = TRUE
    )
  }

  optim.results <- mpi.iapplyLB(1:length(D1_i), parallel.fcn, gpnl.obj)
  #    gpnl.obj=gpnl.obj) # mpi.iapplyLB
  mpi.close.Rslaves(dellog = FALSE)
  mpi.finalize()
  local.pars <- list()
  optim.times <- list()
  cxyz <- matrix(0, nrow = length(D1_i), ncol = 3)
  for (i in 1:length(optim.results)) { # extract optimization results, result also contains timing
    if (methods::is(optim.results[[i]][[1]], "list")) {
      local.pars[[i]] <- optim.results[[i]][[1]]
      optim.times[[i]] <- optim.results[[i]][[2]]
      cxyz[i, ] <- optim.results[[i]][[3]]
    } else {
      local.pars[[i]] <- NULL
    }
  }
  end.time <- Sys.time()
  est.time <- end.time - start.time
  cat("Total time of estimation: ", est.time, " \n")
  return(list(local.pars = local.pars, optim.times = optim.times, est.time = est.time, cxyz = cxyz))
}

#' Local estimation implemented with a naive for-loop
#'
#' Serial implementation for small problems, testing code, etc.
#'
#' @param gpnlobj a gpnl.obj with data.subsets and est.pars components
#' @param est.pars a list of estimation parameters (see \code{gp.nonrigid.registration})
#' @return a list with the optimization and timing results, and data set 1 means used in de-meaning
#' @seealso \code{\link{rmpi.implementation}}
#' @export
serial.implementation <- function(gpnl.obj, est.pars) {
  start.time <- Sys.time()
  D1_i <- gpnl.obj$data.subsets$D1_i
  D2_i <- gpnl.obj$data.subsets$D2_i
  assertthat::are_equal(length(D1_i), length(D2_i))
  optim.results <- list()
  init <- c(log(3), log(1), log(0.01), 0, 0, 0, 0)
  for (i in 1:length(D1_i)) {
    print(i / length(D1_i))
    optim.results[[i]] <- parallel.fcn(i, gpnl.obj)
  }
  local.pars <- list()
  optim.times <- list()
  cxyz <- matrix(0, nrow = length(D1_i), ncol = 3)
  for (i in 1:length(optim.results)) { # extract optimization results, result also contains timing
    if (methods::is(optim.results[[i]][[1]], "list")) {
      local.pars[[i]] <- optim.results[[i]][[1]]
      optim.times[[i]] <- optim.results[[i]][[2]]
      cxyz[i, ] <- optim.results[[i]][[3]]
    } else {
      local.pars[[i]] <- NULL
    }
  }
  end.time <- Sys.time()
  est.time <- end.time - start.time
  cat("Total time of estimation: ", est.time, " \n")
  return(list(local.pars = local.pars, optim.times = optim.times, est.time = est.time, cxyz = cxyz, init = init))
}


#' Quilt plot three column data frame
#'
#' @param datalist.data a three column data frame
#' @param n.x number of breaks in the x dimension
#' @param n.y number of breaks in the y dimension
#' @param add logical, whether to add to current plot
#' @export
plot.data <- function(datalist.data, n.x = 64, n.y = 64,
                      add = FALSE, zlim = NULL) {
  if (!is.null(zlim)) {
    fields::quilt.plot(datalist.data[, 1:2], datalist.data[, 3],
      nx = n.x, ny = n.y, add = add, zlim = zlim, xaxt = "n", yaxt = "n",
    )
  } else {
    fields::quilt.plot(datalist.data[, 1:2], datalist.data[, 3],
      nx = n.x, ny = n.y, add = add, xaxt = "n", yaxt = "n",
    )
  }
}

#' Plot three column data frame in 3D with the \code{rgl} package
#'
#' @param datalist.data a three column data frame
#' @param add logical, whether to add to current plot
#' @param col color of points
#' @export
plot.data.3d <- function(datalist.data, add = FALSE, col = 1, cex = 2) {
  rgl::plot3d(datalist.data[, 1], datalist.data[, 2], datalist.data[, 3],
    add = add, col = col,
    size = cex, xlab = "", ylab = "", zlab = ""
  )
}

#' 3D plot window and subsamples
#'
#' Plots window data and subsamples in \code{gpnl.obj$data.subsets[[window.i]]}
#'
#' @param gpnl.obj a gpnl.obj with data.subsets component
#' @param window.i an integer indicating which window to plot
#' @export
plot.data.subset <- function(gpnl.obj, window.i) {
  i <- window.i
  d1i <- data.frame(gpnl.obj$data.subsets$D1_i[[i]])
  d2i <- data.frame(gpnl.obj$data.subsets$D2_i[[i]])
  d1is <- dplyr::sample_n(d1i, gpnl.obj$est.pars$n.subsamples)
  d2is <- dplyr::sample_n(d2i, gpnl.obj$est.pars$n.subsamples)
  plot3d(d1i[, 1], d1i[, 2], d1i[, 3], col = "black")
  plot3d(d2i[, 1], d2i[, 2], d2i[, 3], col = "grey", add = TRUE)
  plot3d(d1is[, 1], d1is[, 2], d1is[, 3], col = "blue", add = TRUE, size = 10)
  plot3d(d2is[, 1], d2is[, 2], d2is[, 3], col = "red", add = TRUE, size = 10)
}

#' Plots the locally estimated transformation parameters and fitted model surface
#'
#' Plots the local estimates of any of the translation and rotation parameters,
#' as well as the surface predicted by the spatial model fitted to these data
#'
#' @param gpnl.obj a gpnl.obj with fit.results and setup.pars components
#' @param param a character vector, either 'x', 'y', 'z', or 'phi'
#' @export
plot.param <- function(gpnl.obj, param = "x", zlim = NULL) { #' y', 'z', 'phi'
  print(param)
  centroids.rf <- gpnl.obj$fit.results$centroids.rf
  est.params.rf <- gpnl.obj$fit.results$est.params.rf
  GP.fit <- gpnl.obj$fit.results$GP.fit
  if (param == "x") {
    param <- 1
  }
  if (param == "y") {
    param <- 2
  }
  if (param == "z") {
    param <- 3
  }
  if (param == "phi") {
    param <- 4
  }
  n.x <- gpnl.obj$setup.pars$nxx
  n.y <- gpnl.obj$setup.pars$nyy
  # graphics::par(mfrow=c(4,3), mar=c(2,2,1,1))
  if (!is.null(zlim)) {
    fields::quilt.plot(centroids.rf[, 1:2],
      est.params.rf[, param + 3],
      nx = n.x, ny = n.y, zlim = zlim,
      # xaxt='n',
      add.legend = FALSE
    )

    # surface(GP.fit[[param]], xaxt='n', yaxt='n',
    #         xlab='', ylab='', zlim=zlim, type='I')
    pr <- predictSurface(GP.fit[[param]])
    fields::quilt.plot(expand.grid(pr$x, pr$y), pr$z,
      # xaxt='n',
      yaxt = "n",
      xlab = "", ylab = "", zlim = zlim, # tck=seq(1,6,1),
      add.legend = FALSE
    )
  } else {
    fields::quilt.plot(centroids.rf[, 1:2],
      est.params.rf[, param + 3],
      nx = n.x, ny = n.y
    )
    # surface(GP.fit[[param]], xaxt='n', yaxt='n',
    #         xlab='', ylab='', type='I')
    pr <- predictSurface(GP.fit[[param]])
    fields::quilt.plot(expand.grid(pr$x, pr$y), pr$z,
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "",
      add.legend = FALSE
    )
  }
}

#' Plots data and boundaries of window.i
#'
#' Plots data1 and data2, and then for each integer in the vector window.i plots
#' the boundaries of window.i
#'
#' @param gpnl.obj a gpnl.obj with data.list and data.subsets components
#' @param windows.i a vector of integers
#' @export
plot.windows.i <- function(gpnl.obj, windows.i) {
  ## plot grid boxes to see size and overlap
  ii <- windows.i
  data1 <- gpnl.obj$data.list$data1
  data2 <- gpnl.obj$data.list$data2
  D1_i <- gpnl.obj$data.subsets$D1_i
  graphics::par(mfrow = c(1, 2), mar = c(2, 2, 1, 1))
  fields::quilt.plot(data1[, 1:2], data1[, 3])
  for (j in 1:length(ii)) {
    i <- ii[j]
    xl <- min(D1_i[[i]][, 1])
    xr <- max(D1_i[[i]][, 1])
    yb <- min(D1_i[[i]][, 2])
    yt <- max(D1_i[[i]][, 2])
    graphics::rect(xleft = xl, xright = xr, ytop = yb, ybottom = yt)
  }
  fields::quilt.plot(data2[, 1:2], data2[, 3])
  for (j in 1:length(ii)) {
    i <- ii[j]
    xl <- min(D1_i[[i]][, 1])
    xr <- max(D1_i[[i]][, 1])
    yb <- min(D1_i[[i]][, 2])
    yt <- max(D1_i[[i]][, 2])
    graphics::rect(xleft = xl, xright = xr, ytop = yb, ybottom = yt)
  }
}

#' Translates the 3D coordinates of a data frame
#'
#' Translates the first columns of \code{df} by \code{tr[1]}, and similarly for
#' the second and third columns.
#'
#' @param df a data frame
#' @param tr a vector with three real numbers to be applied to the x, y, and z
#' coordinates (columns) of \code{df}, respectively
#' @return a translated data frame
#' @export
translate <- function(df, tr) {
  # df %>% dplyr::mutate( X = .[[1]]+tr[1], Y = .[[2]]+tr[2], Z = .[[3]]+tr[3] )
  df[, 1] <- df[, 1] + tr[1]
  df[, 2] <- df[, 2] + tr[2]
  df[, 3] <- df[, 3] + tr[3]
  return(df)
}

#' Rotates the xy coordinates of a data frame
#'
#' Rotates the first and second columns of \code{df} by \code{rt}. The rotation in the
#' xy-plane is a two-dimensional linear transformation.
#'
#' @param df a data frame
#' @param rt a real number (angle in radians) to be applied to the x, y
#' coordinates of \code{df}
#' @return a rotated data frame
#' @export
rotate2d <- function(df, rt) {
  phi <- rt
  rotat <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), nrow = 2, ncol = 2)
  g <- t(rotat %*% t(data.matrix(df[, 1:2])))
  df[, 1:2] <- g
  return(df)
}

#' Rotates the 3D coordinates of a data frame
#'
#' Rotates the first, second, and third columns of \code{df} by \code{rt}. The 3D
#' rotation is a composition of three rotations applied in the xy, xz, and yz planes,
#' specified by the angles in \code{rt[3]}, \code{rt[2]}, and \code{rt[1]}, respectively.
#'
#' @param df a data frame
#' @param rt a vector of three real numbers (angles in radians) to be applied to the
#' x, y, and z coordinates of \code{df}
#' @return a rotated data frame
#' @export
rotate3d <- function(df, rt) {
  phi_x <- rt[1]
  phi_y <- rt[2]
  phi_z <- rt[3]
  rotat_x <- matrix(c(
    1, 0, 0,
    0, cos(phi_x), -sin(phi_x),
    0, sin(phi_x), cos(phi_x)
  ),
  nrow = 3, ncol = 3, byrow = TRUE
  )
  rotat_y <- matrix(c(
    cos(phi_y), 0, sin(phi_y),
    0, 1, 0,
    -sin(phi_y), 0, cos(phi_y)
  ),
  nrow = 3, ncol = 3, byrow = TRUE
  )
  rotat_z <- matrix(c(
    cos(phi_z), -sin(phi_z), 0,
    sin(phi_z), cos(phi_z), 0,
    0, 0, 1
  ),
  nrow = 3, ncol = 3, byrow = TRUE
  )
  gg <- data.frame(
    t(rotat_x %*% rotat_y %*% rotat_z %*%
      t(data.matrix(df[, 1:3])))
  )
  names(gg) <- c("X", "Y", "Z")
  df[, 1:3] <- gg
  return(df)
}

#' Generates two data sets to be registered
#'
#' \code{generate.data} makes two copies of a point set, and then de-registers one of
#' them using parameters specified in argument \code{dereg.list}.
#'
#' @param N integer specifying the number of points in both data sets combined
#' @param param.list a list of parameters specifying the deregistration applied
#' to the second, moving point set. \code{gen.p} contains the covariance parameters used
#' to generate the z-coordinates of the surfaces, including the smoothness, divisive range,
#' sill variance, and nugget variance. \code{dereg.type} can be either \code{'rigid'}
#' or \code{'nonrigid'}, indicating the type of deregistration to be applie to the
#' moving point set.
#' @return a list with two data frames
#' @examples
#' library(dplyr)
#' library(fields)
#' library(data.table)
#' library(rgl)
#' ### generate two rigidly deformed points sets
#' param.list <- list(
#'   gen.p = c(1, 2, 1, 0.1),
#'   dereg.type = "rigid",
#'   dereg.p = NULL, # c(1, 0.5, 0.1, 0.001),
#'   theta = c(0, 0.3, 0.3, 0.5)
#' )
#' data.gend <- generate.data(250, param.list)
#' d1 <- data.gend$d1
#' d2 <- data.gend$d2
#' plot.data.3d(d1)
#' plot.data.3d(d2, add = TRUE, col = "red")
#' model1 <- gp
#' ### generate two rigidly deformed points sets
#' param.list <- list(
#'   gen.p = c(1, 2, 1, 0.1),
#'   dereg.type = "rigid",
#'   dereg.p = NULL, # c(1, 0.5, 0.1, 0.001),
#'   theta = c(0, 0.3, 0.3, 0.5)
#' )
#' data.gend <- generate.data(250, param.list)
#' d1 <- data.gend$d1
#' d2 <- data.gend$d2
#' plot.data.3d(d1)
#' plot.data.3d(d2, add = TRUE, col = "red")
#' @export
generate.data <- function(N = 200, param.list = list(
                            gen.p = c(1, 2, 1, 0.1), gen.p2 = NULL, # c(15,0.002, 0.001, 0.0001)
                            gen.p3 = NULL, #
                            dereg.type = "rigid", #' nonrigid'
                            dereg.p = NULL, # c(1, 0.5, 0.1, 0.001),
                            theta = NULL, ## dereg.p makes deformation from Matern, theta makes constant (rigid) deformation, beta and beta2 make regression deformation, and thresh makes piecewise deformation
                            beta = NULL, beta2 = NULL, thresh = NULL,
                            param = c("x", "y", "z", "phi")
                          )) { # c(0,0.3,0.3, 0.5)
  p <- param.list$gen.p
  dereg.type <- param.list$dereg.type
  p2 <- param.list$dereg.p
  theta <- param.list$theta
  beta <- param.list$beta
  beta2 <- param.list$beta2
  thresh <- param.list$thresh
  param <- param.list$param
  if (is.null(p2) & is.null(theta) & is.null(beta) & is.null(thresh)) {
    warning("Need either dereg.p or theta to be a vector (not both NULL)")
  }
  xy <- cbind(stats::runif(N, 0, 6), stats::runif(N, 0, 6))
  d <- fields::rdist(xy)
  G <- p[3] * fields::Matern(d, range = p[2], smoothness = p[1])
  C <- G + p[4] * diag(nrow(G))
  z <- (C.c <- t(chol(C))) %*% (e <- stats::rnorm(N))
  if (!is.null(param.list$gen.p2)) {
    ## extremely smooth but short scale and low variance process to add "rocks"
    gp2 <- param.list$gen.p2
    G2 <- gp2[3] * fields::Matern(d, range = gp2[2], smoothness = gp2[1])
    C2 <- G2 + gp2[4] * diag(nrow(G2))
    z2 <- (C.c2 <- t(chol(C2))) %*% (e <- stats::rnorm(N))
    z <- z + z2
  }
  if (!is.null(param.list$gen.p3)) {
    ## extremely smooth but short scale and low variance process to add "rocks"
    gp3 <- param.list$gen.p3
    G3 <- gp3[3] * fields::Matern(d, range = gp3[2], smoothness = gp3[1])
    C3 <- G3 + gp3[4] * diag(nrow(G3))
    z3 <- (C.c3 <- t(chol(C3))) %*% (e <- stats::rnorm(N))
    z <- z + z3
  }
  # fields::quilt.plot(xy[,1], xy[,2], z)
  # source("~/Downloads/myColorRamp.R")
  # cols <- myColorRamp(tim.colors(alpha=0.75), values = z)
  # rgl::plot3d(cbind(xy,z), col=cols, box=FALSE, xlab='', ylab='', zlab='')
  # d <- data.frame(cbind(xy, z)); names(d) <- c('X','Y','Z')
  # p0 <- ggplot(d) + geom_point(aes(X,Y, color=Z)) + scale_color_viridis_c()
  # rgl::plot3d(d[,1:3] )
  d1 <- cbind(xy[1:floor(N / 2), ], z[1:floor(N / 2)])
  d2.b4 <- cbind(xy[(floor(N / 2) + 1):N, ], z[(floor(N / 2) + 1):N])
  # d1 <- data.frame(d1); d2.b4 <- data.frame(d2.b4)
  # names(d1) <- c('X','Y','Z'); names(d2.b4) <- c('X','Y','Z')

  if (dereg.type == "rigid" & !is.null(theta)) {
    # (theta <- -stats::runif(3, 0, 1) ) ### + c(0,0.3,0.3)
    # (ang <- -stats::runif(1, 0, pi/4)) ### + 0.5
    d2 <- translate(rotate2d(d2.b4, theta[4]), theta[1:3])
    return(list(d1 = d1, d2 = d2, d2.b4 = d2.b4, param.list = param.list))
  }
  # d1$f <- 1; d2$f <- 2
  # d1 <- cbind(stats::runif(100), stats::runif(100), stats::runif(100))
  # d2 <- d1 + matrix(rep((theta <- stats::runif(3)), each=100), nrow=100)
  # d3 <- rbind(d1, d2); d3$f <- as.factor(d3$f)
  # p1 <- ggplot(d3) + geom_point(aes(X, Y, color=f))
  # p2 <- ggplot(d3) + geom_point(aes(X, Z, color=f))
  # p3 <- ggplot(d3) + geom_point(aes(Y, Z, color=f))
  # grid.arrange(p0, p1, p2, p3, nrow=2)
  if ("x" %in% param) {
    print("x")
  }
  if ("y" %in% param) {
    print("y")
  }
  if ("z" %in% param) {
    print("z")
  }
  if ("phi" %in% param) {
    print("phi")
  }
  if (dereg.type == "nonrigid") {
    d2 <- d2.b4
    if (!is.null(p2)) {
      d <- fields::rdist(d2.b4[, 1:2])
      G <- p2[3] * fields::Matern(d, range = p2[2], smoothness = p2[1])
      C <- G + p2[4] * diag((N2 <- nrow(G)))
      C.c <- t(chol(C))
      if ("x" %in% param) {
        tx <- C.c %*% (e <- stats::rnorm(N2))
      }
      if ("y" %in% param) {
        ty <- C.c %*% stats::rnorm(N2)
      }
      if ("z" %in% param) {
        tz <- C.c %*% stats::rnorm(N2)
      }
      if ("phi" %in% param) {
        tphi <- C.c %*% stats::rnorm(N2)
      }
    }
    if (!is.null(beta)) {
      if (!is.list(beta)) {
        betal <- list()
        betal[[1]] <- betal[[2]] <- betal[[3]] <- betal[[4]] <- matrix(beta, ncol = 2)
        beta <- betal
      }
      # for linear
      if ("x" %in% param) {
        tx <- beta[[1]] %*% t(d2.b4[, 1:2])
      }
      if ("y" %in% param) {
        ty <- beta[[2]] %*% t(d2.b4[, 1:2])
      }
      if ("z" %in% param) {
        tz <- beta[[3]] %*% t(d2.b4[, 1:2])
      }
      if ("phi" %in% param) {
        tphi <- beta[[4]] %*% t(d2.b4[, 1:2])
      }
      param.list$beta <- beta
      if (!is.null(beta2)) {
        if (!is.list(beta2)) {
          betal <- list()
          betal[[1]] <- betal[[2]] <- betal[[3]] <- betal[[4]] <- matrix(beta2, ncol = 2)
          beta2 <- betal
        }
        # for quadratic
        if ("x" %in% param) {
          tx <- tx + beta2[[1]] %*% t(d2.b4[, 1:2])^2
        }
        if ("y" %in% param) {
          ty <- ty + beta2[[2]] %*% t(d2.b4[, 1:2])^2
        }
        if ("z" %in% param) {
          tz <- tz + beta2[[3]] %*% t(d2.b4[, 1:2])^2
        }
        if ("phi" %in% param) {
          tphi <- tphi + beta2[[4]] %*% t(d2.b4[, 1:2])^2
        }
        param.list$beta2 <- beta2
      }
    }
    if (!is.null(thresh)) {
      # for bump function
      if ("x" %in% param) {
        tx <- ifelse(d2.b4[, 1] > thresh, 0.1, 0)
      }
      if ("y" %in% param) {
        ty <- ifelse(d2.b4[, 1] > thresh, 0.1, 0)
      }
      if ("z" %in% param) {
        tz <- ifelse(d2.b4[, 1] > thresh, 0.1, 0)
      }
      if ("phi" %in% param) {
        tphi <- ifelse(d2.b4[, 1] > thresh, 0.1, 0)
      }
    }

    # tphi <- matrix(0, nrow=nrow(tphi), ncol=ncol(tphi))
    graphics::par(mfrow = c(2, 2))
    if (!exists("tx")) {
      tx <- matrix(0, nrow = nrow(d2.b4), ncol = 1)
    }
    if (!exists("ty")) {
      ty <- matrix(0, nrow = nrow(d2.b4), ncol = 1)
    }
    if (!exists("tz")) {
      tz <- matrix(0, nrow = nrow(d2.b4), ncol = 1)
    }
    if (!exists("tphi")) {
      tphi <- matrix(0, nrow = nrow(d2.b4), ncol = 1)
    }
    # fields::quilt.plot(d2.b4[,1:2], tx )
    # fields::quilt.plot(d2.b4[,1:2], ty )
    # fields::quilt.plot(d2.b4[,1:2], tz )
    # fields::quilt.plot(d2.b4[,1:2], tphi )
    trnsf.That <- cbind(c(tx), c(ty), c(tz), c(tphi))
    d2[, 3] <- d2.b4[, 3] + trnsf.That[, 3] # z translation
    tmp.x <- d2.b4[, 1] + trnsf.That[, 1]
    tmp.y <- d2.b4[, 2] + trnsf.That[, 2]
    cs <- cos(trnsf.That[, 4])
    sn <- sin(trnsf.That[, 4])
    d2[, 1] <- cs * tmp.x - sn * tmp.y
    d2[, 2] <- sn * tmp.x + cs * tmp.y
    return(list(
      d1 = d1, d2 = d2, d2.b4 = d2.b4, param.list = param.list,
      txyzphi = trnsf.That, beide = cbind(xy[, 1], xy[, 2], z)
    ))
  }
}

#' Process residuals from registration
#'
#' \code{process.residuals} calculates mean squared and mean absolute error summary
#' statistics, and the plotting functionality can be used as model diagnostics
#'
#' @param d1 data set 1
#' @param d2.obs the observed (deregistered) data set 2, which need(ed) to be registered
#' @param d2.true the true values of data set 2
#' @param d2.hat the registration model fitted values for data set 2 (i.e. the registered
#' data set 2)
#' @return a list of summary statistics and residual matrices
#' @export
#'
#'
process.residuals <- function(d1, d2.obs, d2.true, d2.hat, data.gend,
                              plot = FALSE,
                              gpnl.obj,
                              zlim = NULL) {
  # d1=datagendlist[[i]]$d1
  # d2.obs= datagendlist[[i]]$d2
  # d2.true=datagendlist[[i]]$d2.b4
  # d2.hat = gpobjlist[[i]]$data.list$data2.trf
  # data.gend=datagendlist[[i]]
  # gpnl.obj=gpobjlist[[i]]
  # zlim=c(-0.3,0.3)

  residuals <- d2.hat - d2.true

  abs.resid.sum <- apply(abs(residuals), 1, sum)
  resid.abs <- cbind(d2.true[, 1:2], abs.resid.sum)
  # plot.data(resid.abs, zlim=range(resid.abs[,3]))
  # plot.data.3d(resid.abs)
  mae <- mean(abs.resid.sum)
  print("MAE: ")
  print(mae)

  sqrd.resid.sum <- apply(residuals^2, 1, sum)
  resid.sqrd <- cbind(d2.true[, 1:2], sqrd.resid.sum)
  # plot.data(resid.sqrd, zlim=range(resid.sqrd[,3]))
  # plot.data.3d(resid.sqrd)
  rmse <- sqrt(mean(sqrd.resid.sum))
  print("RMSE: ")
  print(rmse)

  resid.x <- cbind(d2.true[, 1:2], residuals[, 1])
  resid.y <- cbind(d2.true[, 1:2], residuals[, 2])
  resid.z <- cbind(d2.true[, 1:2], residuals[, 3])

  target <- data.gend$txyzphi # d2.true - d2.obs
  target.x <- cbind(d2.true[, 1:2], -target[, 1])
  target.y <- cbind(d2.true[, 1:2], -target[, 2])
  target.z <- cbind(d2.true[, 1:2], -target[, 3])
  target.phi <- cbind(d2.true[, 1:2], -target[, 4])
  zlimx <- range(c(
    target.x[, 3],
    gpnl.obj$fit.results$est.params.rf[, 4]
  ))
  zlimy <- range(c(
    target.y[, 3],
    gpnl.obj$fit.results$est.params.rf[, 5]
  ))
  zlimz <- range(c(
    target.z[, 3],
    gpnl.obj$fit.results$est.params.rf[, 6]
  ))
  zlimphi <- range(c(
    target.phi[, 3],
    gpnl.obj$fit.results$est.params.rf[, 7]
  ))
  if (!is.null(zlim)) {
    zlimx <- zlimy <- zlimz <- zlimphi <- zlim
  }
  if (plot == TRUE) {
    graphics::par(mfrow = c(4, 3), mar = c(1, 1, 1, 2), oma = c(1, 1, 1, 1))
    plot.param(gpnl.obj, param = "x", zlim = zlimx)
    # mtext('Estimates and TPS fit', adj=7)
    plot.data(target.x, zlim = zlimx)
    # mtext('Estimation target values')
    # plot.data(resid.x, zlim=range(resid.x[,3]))
    # mtext('x residuals')
    plot.param(gpnl.obj, param = "y", zlim = zlimy)
    plot.data(target.y, zlim = zlimy)
    # plot.data(resid.y, zlim=range(resid.y[,3]))
    plot.param(gpnl.obj, param = "z", zlim = zlimz)
    plot.data(target.z, zlim = zlimz)
    # plot.data(resid.z, zlim=range(resid.z[,3]))
    plot.param(gpnl.obj, param = "phi", zlim = zlimphi)
    plot.data(target.phi, zlim = zlimphi)
  }
  return(list(target = target, residuals = residuals, rmse = rmse, mae = mae))
}

#' Plots fitted vs true registration
#'
#' @param d1 x, y, z coordinates of data set 1
#' @param d2.obs x, y, z coordinates of data set 2, pre-registration
#' @param d2.true x, y, z coordinates of data set 2, the true target alignment
#' @param d2.hat x, y, z coordinates of data set 2, post-registration
#'
#' @return plots
#' @export
#'
#' @examples
plot.fitted.vs.truth <- function(d1, d2.obs, d2.true, d2.hat) {
  plot.data.3d(d2.true, cex = 3)
  plot.data.3d(d2.hat, add = TRUE, col = "blue", cex = 4)
  plot.data.3d(d1, add = TRUE, col = "red")
  plot.data.3d(d2.obs, add = TRUE, col = "red", cex = 0.1)
}



#' Perform local kriging using 1000 nearest neighbors using covariance parameters in pars
#'
#' @param data available data
#' @param x locations at which to predict
#' @param pars list of data frames with covariance parameters
#'
#' @return the predicted surface and prediction variances
#' @export
#'
#' @examples
local.krig.rigid <- function(data, x, pars) {
  pp <- pars
  pred <- c()
  for (i in 1:nrow(x)) {
    print(i)
    # indx <- apply( fields::rdist(x[i,], data[,1:2]), 1, function(x) order(x, decreasing=F) )[1:500,]
    # datatouse <- data[indx, ]
    indx <- FNN::get.knnx(data = data.matrix(data[, 1:2]), query = data.matrix(x[i, 1:2]), k = 100)
    indx <- indx$nn.index
    # indx <- sample(indx, 100)
    dind <- duplicated(data[indx, 1:2])
    datatouse <- data[indx, ]
    datatouse <- datatouse[!dind, ]
    d <- fields::rdist(datatouse[, 1:2])
    G <- pp[2] * fields::Matern(d, range = pp[1], smoothness = 1) + pp[3] * diag(nrow(d))
    d0 <- fields::rdist(x[i, ], datatouse[, 1:2])
    G0 <- pp[2] * fields::Matern(d0, range = pp[1], smoothness = 1)
    pred[i] <- G0 %*% chol2inv(chol(G)) %*% c(datatouse[, 3])
  }
  return(pred)
}

#' Perform local kriging using 1000 nearest neighbors using covariance parameters in pars
#'
#' @param data available data
#' @param x locations at which to predict
#' @param pars list of data frames with covariance parameters
#'
#' @return the predicted surface and prediction variances
#' @export
#'
#' @examples
local.krig.nonrigid <- function(data, x, pars) {
  pred <- c()
  for (i in 1:nrow(x)) {
    print(i)
    pp <- c(data.matrix(pars[i, ]))
    # indx <- apply( fields::rdist(x[i,], data[,1:2]), 1, function(x) order(x, decreasing=F) )[1:500,]
    indx <- FNN::get.knnx(data = data.matrix(data[, 1:2]), query = data.matrix(x[i, 1:2]), k = 1000)
    print("Found NN")
    indx <- indx$nn.index
    # indx <- sample(indx, 10)
    dind <- duplicated(data[indx, 1:2])
    datatouse <- data[indx, ]
    datatouse <- datatouse[!dind, ]
    d <- fields::rdist(datatouse[, 1:2])
    G <- pp[2] * fields::Matern(d, range = pp[1], smoothness = 1) + pp[3] * diag(nrow(d))
    d0 <- fields::rdist(x[i, ], datatouse[, 1:2])
    G0 <- pp[2] * fields::Matern(d0, range = pp[1], smoothness = 1)
    pred[i] <- G0 %*% chol2inv(chol(G)) %*% c(datatouse[, 3])
  }
  return(pred)
}

# local.krig.nonrigid2 <- function(data, x, centroids, pars){
#   pred <- c()
#   for(i in 1:nrow(x)){
#     print(i)
#     indx1 <- apply( fields::rdist(x[i,], centroids), 1, which.min)
#     pp <- pars[indx1, ]
#     indx <- apply( fields::rdist(x[i,], data[,1:2]), 1, function(x) order(x, decreasing=F) )[1:500,]
#     datatouse <- distinct(data[indx, ])
#     d <- fields::rdist(datatouse[,1:2])
#     G <- pp[2] * fields::Matern(d, range=pp[1], smoothness = 1) + pp[3]*diag(nrow(d))
#     d0 <- fields::rdist(x[i,], datatouse[,1:2])
#     G0 <- pp[2] * fields::Matern(d0, range=pp[1], smoothness = 1)
#     pred[i] <- G0 %*% chol2inv(chol(G)) %*% c(datatouse[,3])
#   }
#   return(pred)
# }


#' Perform local kriging using 1000 nearest neighbors using covariance parameters in pars
#'
#' @param data available data
#' @param x locations at which to predict
#' @param pars list of data frames with covariance parameters
#'
#' @return the predicted surface and prediction variances
#' @export
#'
#' @examples
local.krig <- function(data, x, pars) { # pars is a list of data frames
  pred <- pred.var <- matrix(NA, nrow = nrow(x), ncol = length(pars))
  indxtt <- FNN::get.knnx(data = data.matrix(data[, 1:2]), query = data.matrix(x[, 1:2]), k = 1000)
  indxt <- indxtt$nn.index
  dists <- indxtt$nn.dist
  print("Found NN")
  for (i in 1:nrow(x)) {
    print(i)
    # indx <- apply( fields::rdist(x[i,], data[,1:2]), 1, function(x) order(x, decreasing=F) )[1:500,]
    # indx <- sample(indx, 10)
    indx <- indxt[i, ]
    # dind <- duplicated(data[indx, 1:2])
    datatouse <- data[indx, ]
    # datatouse <- datatouse[!dind,]
    d <- fields::rdist(datatouse[, 1:2])
    d0 <- dists[i, ] # fields::rdist(x[i,], datatouse[,1:2])
    for (j in 1:length(pars)) {
      pp <- c(data.matrix(pars[[j]][i, ]))
      G <- pp[2] * fields::Matern(d, range = pp[1], smoothness = 1) + pp[3] * diag(nrow(d))
      G0 <- pp[2] * fields::Matern(d0, range = pp[1], smoothness = 1)
      pred[i, j] <- G0 %*% (sig.inv <- chol2inv(chol(G))) %*% c(datatouse[, 3])
      pred.var[i, j] <- pp[2] - t(G0) %*% sig.inv %*% G0
    }
  }
  return(list(pred = pred, pred.var = pred.var))
}


#' Calculate log score for normal distribution
#'
#' @param mu mean of normal distribution
#' @param sigma2 variance of normal distribution
#' @param x realization from process
#'
#' @return the log score
#' @export
#'
#' @examples
log.score <- function(mu, sigma2, x) {
  # p = c(theta, sigma2, tau2)
  # cat(p[1], ', ', p[2], ', ', p[3], p[4], p[5], p[6], p[7] )
  # cat(", \n")
  y <- x
  NLL <- (y - mu)^2 / (2 * sigma2) + length(y) / 2 * (log(2 * pi) + log(sigma2))
  return(NLL)
}


#' Calculate crps score for normal distribution
#'
#' @param mu mean of normal distribution
#' @param sigma2 variance of normal distribution
#' @param x realization from process
#'
#' @return the crps score
#' @export
#'
#' @examples
crps.normal <- function(mu, sigma2, x) {
  ### mu, sigma2, and x are vectors of equal length
  ### mu and sigma2 are the mean and variance of the normal predictive distribution
  ### x is the actual realization that came to fruition
  s <- sqrt(sigma2)
  Z <- (x - mu) / s
  score <- -s * ((1 / sqrt(pi)) - 2 * stats::dnorm(Z) - Z * (2 * stats::pnorm(Z) - 1))
  return(score)
}

#' Calculate crps score
#'
#' @param X1 matrix with ncol independent realizations of the random vector X
#' @param X2 matrix with another ncol independent realizations of the random vector X
#' @param x the actual realization that came to fruition
#'
#' @return the crps score
#' @export
#'
#' @examples
crps <- function(X1, X2, x) {
  # X1 and X2 are matrices with ncol independent realizations of the random vector X
  # the number of rows of X1 and X2 is equal to the number of spatial locations
  # x is the actual realization that came to fruition
  term1 <- 0.5 * apply(abs(X1 - X2), 1, mean)
  xx <- replicate(ncol(X1), x) # repeats column vector x ncol(X1) times
  term2 <- apply(abs(X1 - xx), 1, mean)
  -(term1 - term2)
}
