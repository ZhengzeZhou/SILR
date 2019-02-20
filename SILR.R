library(cladoRcpp)

SILR <- function(hnum = 1, alpha = 0.05, pnull = 0.2, ri = 0.01, printallps = 0, hitfreq = c(6, 10, 8, 5, 4, 2)) {
  # 1. Enter parameters and data, and find n

  # hnum <- 2 # 1 for an upper-tail test;2 for a lower-tail test; 3 for a some-tail test
  # usefft <- 0 # 1 to use FFT; 0 otherwise
  # alpha <- 0.05 # significance level, often 0.05 or 0.01
  # pnull <- 0.2 # the probability of a correct answer on each trial if the numm H is true
  # t <- 5 # the number of trials each participant undergoes
  # ri <- 0.01 # "rounding interval", usually set to 0.01
  # printallps <- 0 # whether to print all values of ps

  # # a vector of t+1 values showing the number of participants who got each number of hits from 0 to t
  # # hitfreq <- c(6, 10, 8, 5, 4, 2) 
  # hitfreq <- c(60, 109, 80, 50, 40, 20) 
  # # hitfreq <- c(200, 100, 80, 50, 40, 20)

  # if (alpha < 1e-12){
  #   alpha <- 1e-12
  # }

  # t1 <- hsec

  n <- sum(hitfreq) # the number of participants


  # 2. Construct p1
  t = length(hitfreq) - 1
  
  xup <- 0:t
  # calculate the vector showing the probability under the null that 
  # any one participant will get each possible number of hits from 0 to t
  probH0 <- factorial(t) / (factorial(xup) * factorial(t - xup)) * pnull ^ xup * (1 - pnull) ^ (t - xup) 
  # entry j is the probability of j-1 hits under the H which maximized that probability
  probH1 <- factorial(t) / (factorial(xup) * factorial(t - xup)) * (xup / t) ^ xup * (1 - xup / t) ^ (t - xup)

  gvec <- log(probH1 / probH0) # logarithm of the reciprocal of a likelihood ratio

  # set some values in gvec to 0 based on which test
  if (hnum == 1){
    nbelow <- sum(xup / t < pnull - 0.00001)
    gvec[1:nbelow] <- 0
  } else if (hnum == 2) {
    nabove <- sum(xup / t > pnull + 0.00001)
    gvec[(t - nabove + 2):(t + 1)] <- 0
  }

  gvec <- ri * round(gvec/ri)
  g <- hitfreq %*% gvec
  gcn <- round(g / ri + 1)

  order_index <- order(gvec)

  gcum <- rep(0, n)
  start <- 1
  for (index in order_index) {
    gcum[start:(start + hitfreq[index] - 1)] <- gvec[index]
    start <- start + hitfreq[index]
  }
  gcum <- cumsum(gcum)

  p1cellnum <- round(gvec/ri + 1)
  numcellsp1 <- max(p1cellnum)
  p1 <- rep(0, numcellsp1)
  for (i in 1:(t + 1)) {
    p1[p1cellnum[i]] <- p1[p1cellnum[i]] + probH0[i]
  }

  store <- list()
  store[["1"]] <- p1

  # 3. Doubling convolutions
  nt <- 1
  currentvec <- p1
  ntoolarge <- n + 1
  pstoolow <- -1

  while (2 * nt <= n) { 
    # convolve
    # newvec <- convolve(currentvec, rev(currentvec), type = "open")
    newvec <- rcpp_convolve(currentvec, currentvec)
    # trim
    if (length(newvec) > gcn) {
      finalp <- sum(newvec[gcn:length(newvec)])
      newvec[gcn] <- finalp
      newvec <- newvec[1:gcn]
    }
    
    # find and print ps
    nt <- 2 * nt
    currentcell <- round(gcum[nt]/ri) + 1
    newps <- sum(newvec[currentcell:length(newvec)])
    if (printallps == 1) {
      print(nt)
      print(newps)
    }
    if (newps < alpha) {
      ntoolarge <- nt
      pstoolow <- newps
    } else {
      pstoohigh <- newps
    }
    
    
    if (newps > alpha) {
      maxnsn <- nt
      maxnsvec <- newvec
      store[[toString(nt)]] <- newvec
      currentvec <- newvec
    }
  } 

  # Phase 2

  while (ntoolarge - maxnsn > 1) {
    maxntoadd <- ntoolarge - maxnsn - 1
    
    for (i in rev(names(store))) {
      if (as.integer(i) <= maxntoadd) {
        ntoadd <- as.integer(i)
        break
      } 
    }
    
    vectoadd <- store[[toString(ntoadd)]]
    # newvec <- convolve(maxnsvec, rev(vectoadd), type = "open")
    newvec <- rcpp_convolve(maxnsvec, vectoadd)
    nt <- maxnsn + ntoadd
    currentcell <- round(gcum[nt] / ri) + 1
    newps <- sum(newvec[currentcell:length(newvec)])
    
    # trim 
    if (length(newvec) > gcn) {
      finalp <- sum(newvec[gcn:length(newvec)])
      newvec[gcn] <- finalp
      newvec <- newvec[1:gcn]
    }
    
    if (newps < alpha) {
      ntoolarge <- nt
      pstoolow <- newps
    } else {
      pstoohigh <- newps
      maxnsn <- nt
      maxnsvec <- newvec
    }
  }

  my_list <- list("maxnsn" = maxnsn, "pstoohigh" = pstoohigh, "ntoolarge" = ntoolarge, "pstoolow" = pstoolow, "withtrait" = n - maxnsn)

  return(my_list)

}


findPmax <- function(hnum = 1, alpha = 0.05, pnull = 0.04, ri = 0.01, printallps = 0, hitfreq = c(1, 5, 10, 10, 6, 5), tol = 1e-2) {
  
  ps <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh
  
  if (ps > alpha) {
    print("Current pnull already yields p value larger than alpha!")
    return()
  } 
  
  pg <- pnull
  pnull <- (pg + 1) / 2 
  ps <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh
  
  while (ps < alpha) {
    pg <- pnull
    pnull <- (pg + 1) / 2 
    ps <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh
  }
  
  pg_e <- SILR(hnum, alpha = 0, pg, ri, printallps, hitfreq)$pstoohigh - alpha
  pnull_e <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh - alpha
  
  if (abs(pg_e) < tol)
    return(pg)
  
  if (abs(pnull_e) < tol)
    return(pull)
  
  xv <- c(pg, pnull)
  ev <- c(pg_e, pnull_e)
  
  xv_log <- c(pg, pnull)
  ev_log <- c(pg_e, pnull_e)
  
  order <- 1
  
  while (abs(ev[length(ev)]) >= tol) {
    
    if (order > 10) {
      index <- which.max(ev)
      ev <- ev[-index]
      xv <- xv[-index]
      order <- 10
    }
    fit <- lm(xv~poly(ev, order, raw = TRUE))
    pnull <- fit$coefficients[[1]]
    pnull_e <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh - alpha
    xv <- append(xv, pnull)
    ev <- append(ev, pnull_e)
    xv_log <- append(xv_log, pnull)
    ev_log <- append(ev_log, pnull_e)
    
    order <- order + 1
    
  }
  return(list(xv = xv_log, ev = ev_log))
}

findPmaxlog <- function(hnum = 1, alpha = 0.05, pnull = 0.04, ri = 0.01, printallps = 0, hitfreq = c(1, 5, 10, 10, 6, 5), tol = 1e-2) {
  
  ps <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh
  
  if (ps > alpha) {
    print("Current pnull already yields p value larger than alpha!")
    return()
  } 
  
  pg <- pnull
  pnull <- (pg + 1) / 2 
  ps <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh
  
  while (ps < alpha) {
    pg <- pnull
    pnull <- (pg + 1) / 2 
    ps <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh
  }
  
  # pg_e <- SILR(hnum, alpha = 0, pg, ri, printallps, hitfreq)$pstoohigh - alpha
  # pnull_e <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh - alpha
  
  pg_e <- log(SILR(hnum, alpha = 0, pg, ri, printallps, hitfreq)$pstoohigh / alpha)
  pnull_e <- log(SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh / alpha)
  
  if (abs(pg_e) < tol)
    return(pg)
  
  if (abs(pnull_e) < tol)
    return(pull)
  
  xv <- c(pg, pnull)
  ev <- c(pg_e, pnull_e)
  
  xv_log <- c(pg, pnull)
  ev_log <- c(pg_e, pnull_e)
  
  order <- 1
  
  while (abs(ev[length(ev)]) >= tol) {
    
    if (order > 10) {
      index <- which.max(ev)
      ev <- ev[-index]
      xv <- xv[-index]
      order <- 10
    }
    fit <- lm(xv~poly(ev, order, raw = TRUE))
    pnull <- fit$coefficients[[1]]
    # pnull_e <- SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh - alpha
    pnull_e <- log(SILR(hnum, alpha = 0, pnull, ri, printallps, hitfreq)$pstoohigh / alpha)
    xv <- append(xv, pnull)
    ev <- append(ev, pnull_e)
    xv_log <- append(xv_log, pnull)
    ev_log <- append(ev_log, pnull_e)
    
    order <- order + 1
    
  }
  return(list(xv = xv_log, ev = ev_log))
}


