#' Get Basis splines data to run JAGS country model
#'
#' Function to get basis splines data to run country model in JAGS.
#'
#' @param nyears.c vector with number of years of observations for each country
#' @param niso number of countries
#' @param order order of splines penalization; must be 1 or 2.
#' @param I interval length between two knots during observation period. Default is 2.5 years.
#' @export
#' @return A list of data for input to JAGS, which includes matrix of basis splines for each country.


getSplinesData <- function(nyears.c,
                           niso,
                           order,
                           I = 2.5){
  if(order!=1&order!=2){
    stop("Order of penalization must either be 1 or 2.")
  }
  max.nyears <- max(nyears.c)
  # get max dimensions
  x.t <- 1:max.nyears
  sp <- GetSplines(x.t, I = I)
  K <- length(sp$knots.k)
  B.tk <- sp$B.ik
  # this is the maximum dimensions of the basis splines.
  max.dim.B <- dim(B.tk)
  len.t <- max.dim.B[1]
  len.k <- max.dim.B[2]
  # for order =1 , H = k - 1
  # for order = 2 H = k-2
  # following from this, for order = 1, dim(Z.ih) will be t x (k-1)
  # for order = 2, dim(Z.ih) will be t x (k-2)

  K.c <- rep(NA, niso)
  H1.c <-rep(NA, niso)
  H2.c <-rep(NA, niso)
  B.tkc<- array(NA, c(len.t, len.k, niso))
  Z1.tkc<- array(0, c(len.t, len.k-1, niso)) #first diff
  Z2.tkc<- array(0, c(len.t, len.k-2, niso)) #second diff
  Delta1comb.khc <- array(NA, c(len.k, len.k-1, niso)) #first diff
  Delta2comb.khc <- array(NA, c(len.k, len.k -2, niso)) #second diff
  G.kdc <-  array(NA, c(len.k, 2, niso)) #second diff
  BG.tdc <-  array(NA, c(len.t, 2, niso)) #second diff

  for(i in 1:niso){
    x.t <- 1:nyears.c[i]
    sp <- GetSplines(x.t, I = I)
    K <- length(sp$knots.k)
    B.tk <- sp$B.ik
    B.tkc[1:nrow(B.tk),1:ncol(B.tk),i] <- B.tk
    K.c[i] <- K
    H1.c[i] <- K-1
    H2.c[i] <- K-2
    # stuff for reparameterization
    ## first difference
    Delta.hk <- diff(diag(K), diff = 1)
    Delta1comb.kh <- t(Delta.hk)%*%solve(Delta.hk%*%t(Delta.hk))
    Delta1comb.khc[1:K, 1:(K-1), i] <- Delta1comb.kh
    Z.ih <-B.tk%*%Delta1comb.kh
    Z1.tkc[1:nrow(B.tk),1:(ncol(B.tk)-1),i ] <- Z.ih
    ## second difference
    Delta.hk <- diff(diag(K), diff = 2) # difference matrix
    Delta2comb.kh <- t(Delta.hk)%*%solve(Delta.hk%*%t(Delta.hk))
    Delta2comb.khc[1:K, 1:(K-2), i] <- Delta2comb.kh
    Z.ih <-B.tk%*%Delta2comb.kh
    Z2.tkc[1:nrow(B.tk),1:(ncol(B.tk)-2),i ] <- Z.ih
    G.kd <- cbind(rep(1, K), seq(1, K)-K/2)
    G.kdc[1:K,,i] <-  G.kd
    BG.td <- B.tk%*%G.kd
    BG.tdc[1:nrow(B.tk), ,i] <- BG.td
  }
  if(order==1){
    return(list(K.c = K.c, H.c = H1.c, H = max(H1.c),B.tkc = B.tkc, Z.tkc = Z1.tkc))
  }
  if(order==2){
    return(list(K.c = K.c, H.c = H2.c, H = max(H2.c), D = 2, B.tkc = B.tkc, Z.tkc = Z2.tkc, BG.tdc = BG.tdc))
  }
}
