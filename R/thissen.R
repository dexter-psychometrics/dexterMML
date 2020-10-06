



# Thissen stuff


# Bock and Aitkin IRT
# using the Stouffer-Toby (1951)
test_thissen = function()
{
x1 <- c(rep(T,4),F,rep(T,3),rep(F,3),T,rep(F,4))
x2 <- c(rep(T,3),F,T,T,F,F,T,T,F,F,T,rep(F,3))
x3 <- c(T,T,F,T,T,F,T,F,T,F,T,F,F,T,F,F)
x4 <- c(T,F,T,T,T,F,F,T,F,T,T,F,F,F,T,F)
x <- cbind(x1,x2,x3,x4)
r <- c(42,23,6,6,1,24,25,7,4,2,1,38,9,6,2,20)

stouffer <- data.frame(I(x),r)

dat = as.matrix(stouffer)[rep(1:nrow(stouffer),stouffer$r),1:4]

est(dat)$item_em_step
thissen_main(stouffer,seq(-4,4,0.5))

}

Gaussian.pts <-
  function(mu,sigma,theta) {
    curve <- exp(-0.5*((theta - mu)/sigma)^2)
    curve <- curve/sum(curve)
  }

trace.line.pts <-
  function(a,b,theta) {
    traceline <- 1/(1+exp(-a*(theta-b)))
  }

# note that this loglikelihood function is of the same sort as
# was used for the logistic analysis of Constant Method data
ll.2pl <-
  function(p2,r1,r0,theta) {
    a <- p2[1]
    b <- p2[2]

    itemtrace <- trace.line.pts(a,b,theta)

    l <- (-1)*(sum(r1*(log(itemtrace))) + (sum(r0*(log(1.0-itemtrace)))))
  }

Estep.2pl <-
  function(p,testdata,theta) {
    a=double(length(p)/2)
    b=double(length(p)/2)

    for (i in 1:ncol(testdata$x)) {
      a[i] <- p[2*(i-1) + 1]
      b[i] <- p[2*i]
    }

    # the following three blocks "make space" for the trace lines, and the E-tables
    itemtrace <- matrix(0,nrow=ncol(testdata$x),ncol=length(theta))
    r1 <- matrix(0,nrow=ncol(testdata$x),ncol=length(theta))
    r0 <- matrix(0,nrow=ncol(testdata$x),ncol=length(theta))

    # compute the trace lines
    for (i in 1:length(a)) {
      itemtrace[i,] <- trace.line.pts(a[i],b[i],theta)
    }
    # loop over response patterns and compute posteriors
    for (i in 1:length(testdata$r)) {
      posterior <- Gaussian.pts(0,1,theta)
      for (item in 1:length(a)) {
        x <- I(testdata$x[i,item])
        if (x)
          posterior <- posterior*itemtrace[item,]
        else
          posterior <- posterior*(1-itemtrace[item,])
      }
      # normalize posterior
      # to area equals number of persons with this response pattern
      expd <- sum(posterior)
      posterior <- (posterior/expd)*testdata$r[i]
      # put this response pattern's people in the r1 and r0 tables
      for (item in 1:length(a)) {
        x <- I(testdata$x[i,item])
        if (x)
          r1[item,] <- r1[item,] + posterior
        else
          r0[item,] <- r0[item,] + posterior
      }
    } # end loop over response patterns
    rlist <- list(r1,r0)
  } # end of E-step


thissen_main = function(stouffer,theta,a = c(1,1,1,1),b = c(0,0,0,0))
{
  # the "main" routine:
  # initialize the item parameters (starting values)


  p <- rep(0,2*length(a))
  for (i in 1:length(a)) {
    p[2*(i-1) + 1] <- a[i]
    p[2*i] <- b[i]
  }
  # lastp, previousp, rate, and Hessian are storage for
  # the rate-of-convergence approximations for the standard errors
  lastp <- p
  Hessian <- matrix(0, nrow=length(p), ncol=length(p))

  # cycles loop
  for (cycles in 1:25) {
    # E-step:
    rlist <- Estep.2pl(p,stouffer,theta)

    # M-step:
    # for each item (loop)
    for (item in 1:length(a)) {
      r1 <- rlist[[1]][item,]
      r0 <- rlist[[2]][item,]
      p2 <- c(a[item], b[item])
      # note that this call to nlm is a two-variable logit regression problem:
      Mout <- nlm(f=ll.2pl,p=p2,hessian=TRUE,r1=r1,r0=r0,theta=theta)
      # new a,b are in Mout$estimate[1] and [2]
      a[item] <- Mout$estimate[1]
      b[item] <- Mout$estimate[2]
    }
    # update parameter list
    lastp <- p
    for (i in 1:length(a)) {
      p[2*(i-1) + 1] <- a[i]
      p[2*i] <- b[i]
    }

  }
  list(a=p[seq(1,8,2)],b=p[seq(2,8,2)])
}



