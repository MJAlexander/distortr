plotCovariance <- function(method, params, order = NULL, nyears = 10){
  for (i in 1:length(params)){
    assign(names(params)[i], params[[i]])
  }
  if(method=="ar"){
   cov.fn <- function(t){
     out <- sigma^2/(1-rho^2)*rho^(abs(t))
     return(out)
   }
  }

  if(method == "arma"){
    cov.fn <- function(t){
      if(t == 0){
        out <- (sigma.ar^2*(1 + 2*phi*theta + theta^2))/(1-phi^2)
      }
      else{
        out <- (sigma.ar^2*(1 + phi*theta)*(phi + theta))/(1-phi^2)*phi^(abs(t))
      }
      return(out)
    }
  }

  if(method == "gp"){
    cov.fn <- function(t){
      out <- sigma.g^2 * p^(abs(t)^2)
      return(out)
    }
  }

  if(method=="splines"){
    if(is.null(order)){
      stop("Order of penalization must be specified.")
    }
    if(order==0){
      cov.fn <- function(nyears){
        x.i <- seq(-nyears,nyears, by = 0.1)
        sp <- GetSplines(x.i)
        Ca<- sigma^2*solve(t(sp$B.ik)%*%sp$B.ik)
        out <- (sp$B.ik%*%Ca%*%t(sp$B.ik))[which(x.i==0),]
        return(out)
      }
    }
    else{
      cov.fn <- function(nyears){
        x.i <- seq(-nyears,nyears, by = 0.1)
        sp <- GetSplines(x.i)
        D <- diff(diag(length(sp$knots.k)), diff = order)
        Ca<- sigma^2*solve(t(sp$B.ik)%*%sp$B.ik + (sigma.alpha^-2)*t(D)%*%D)
        out <- (sp$B.ik%*%Ca%*%t(sp$B.ik))[which(x.i==0),]
        return(out)
      }
    }

  }

  if(method=="splines"){
    res <- cov.fn(nyears)
    y <- res/max(res)
    pt <- ggplot(data = NULL, aes(x = seq(-nyears,nyears, by = 0.1), y = y)) +
      geom_line() + theme_bw() + xlab("x") + ylab("covariance") + ggtitle(paste0("Covariance between points for ", method))
  }
  else{
    res <- sapply(seq(-nyears,nyears, by = 0.1), function(i) cov.fn(i))
    y <- res/max(res)
    pt <- ggplot(data = NULL, aes(x = seq(-nyears,nyears, by = 0.1), y = y)) +
      geom_line() + theme_bw() + xlab("x") + ylab("covariance") + ggtitle(paste0("Covariance between points for ", method))

  }
  return(pt)
}
