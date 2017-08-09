## extension of ns() to include different boundary derivatives,
## centering and cure
nsx <- 
function (x, df = NULL, knots = NULL, intercept = FALSE,
          Boundary.knots = range(x),
          derivs = if (cure) c(2,1) else c(2,2),
          log=FALSE, # deprecated: only used in rstpm2:::stpm2Old
          centre = FALSE, cure = FALSE, stata.stpm2.compatible=FALSE) 
{
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
            Boundary.knots[2L])
    }
    else outside <- FALSE
    if (!missing(df) && missing(knots)) {
        nIknots <- df - 1 - intercept + 4 - sum(derivs) 
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used ", 1 + intercept)
        }
        knots <- if (nIknots > 0) {
          knots <- if (!cure)
            seq.int(0, 1, length.out = nIknots + 2L)[-c(1L, 
                            nIknots + 2L)]
          else c(seq.int(0, 1, length.out = nIknots + 1L)[-c(1L, 
                                 nIknots + 1L)], 0.95)
          if (!stata.stpm2.compatible)
            stats::quantile(x[!outside], knots)
          else stats::quantile(x[!outside], round(knots,2), type=2)
        }
    }
    else nIknots <- length(knots)
    Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
    if (any(outside)) {
        basis <- array(0, c(length(x), nIknots + 4L))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1L]
            xl <- cbind(1, x[ol] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
                1))$design
            basis[ol, ] <- xl %*% tt
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2L]
            xr <- cbind(1, x[or] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
                1))$design
            basis[or, ] <- xr %*% tt
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 
                4)$design
    }
    else basis <- spline.des(Aknots, x, 4)$design
    const <- splineDesign(Aknots, rep(Boundary.knots, 3-derivs), 4, c(derivs[1]:2, derivs[2]:2))
    if (!intercept) {
        const <- const[, -1, drop = FALSE]
        basis <- basis[, -1, drop = FALSE]
    }
    qr.const <- qr(t(const))
    q.const <- qr.Q(qr.const, complete=TRUE)[, -(1L:2L), drop = FALSE] # NEW
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:nrow(const)), drop = FALSE])
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    if (centre) {
      centreBasis <- nsx(centre,
                         knots=if (is.null(knots)) numeric(0) else knots,
                         Boundary.knots=Boundary.knots, 
                         intercept=intercept, derivs=derivs, centre=FALSE, log=log)
      oldAttributes <- attributes(basis)
      basis <- t(apply(basis,1,function(x) x-centreBasis))
      attributes(basis) <- oldAttributes
    }
    a <- list(degree = 3, knots = if (is.null(knots)) numeric(0) else knots, 
        Boundary.knots = Boundary.knots, intercept = intercept, derivs=derivs,
              centre=centre, log=log, q.const=q.const)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("nsx", "basis", "matrix")
    basis
}
makepredictcall.nsx <- 
function (var, call) 
{
    if (as.character(call)[1L] != "nsx") 
        return(call)
    at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
                            "derivs", "centre", "log")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}
predict.nsx <- 
function (object, newx, ...) 
{
    if (missing(newx)) 
        return(object)
    a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
        "intercept", "derivs", "centre", "log")])
    do.call("nsx", a)
}
Shat <- function(obj)
  {
    ## predicted survival for individuals (adjusted for covariates)
    newobj = survfit(obj,se.fit=FALSE)
    surv = newobj$surv
    rr = try(predict(obj,type="risk"),silent=TRUE)
    if (inherits(rr,"try-error"))
        rr <- 1
    surv2 = surv[match(obj$y[,ncol(obj$y)-1],newobj$time)]
    return(surv2^rr)
  }
replaceCall=function(obj,old,new) {
  if (is.atomic(obj) && length(obj)>1)
    return(as.call(c(quote(c),lapply(as.list(obj),replaceCall,old,new))))
  if (is.name(obj) || is.symbol(obj) || (is.atomic(obj) && length(obj)==1)) {
    if (obj==old) return(new)
    else return(obj)
  }
##   if (length(obj)==1 && length(obj[[1]])==1) {
##     if (obj==old) return(new)
##     else return(obj)
##   }
  as.call(lapply(obj,replaceCall,old,new))
}
replaceFormula=function(...) as.formula(replaceCall(...))
## replaceFormula(~f(a+b),quote(f),quote(g))
allCall=function(obj) {
  if (is.atomic(obj) && length(obj)==1) return(obj)
  if (is.atomic(obj) && length(obj)>1) return(as.call(c(quote(c),as.list(obj))))
  if (is.name(obj) || is.symbol(obj)) return(obj)
  as.call(lapply(obj,allCall))
}
## allCall(as.call(c(quote(ns),list(df=3,knots=c(1,2)))))[[2]]
vector2call=function(obj) {
  if (is.atomic(obj) && length(obj)==1) return(obj)
  if (is.atomic(obj) && length(obj)>1) return(as.call(c(quote(c),as.list(obj))))
  if (is.name(obj) || is.symbol(obj)) return(obj)
  lapply(obj,allCall) # is this correct?
}
## vector2call(list(df=3,knots=c(1,2)))
findSymbol <- function(obj,symbol) {
  if (is.symbol(obj) && obj==symbol) TRUE else
  if (is.symbol(obj)) FALSE else
  if (is.atomic(obj)) FALSE else
  Reduce(`|`,lapply(obj,findSymbol,symbol),FALSE)
}
rhs=function(formula) 
  if (length(formula)==3) formula[[3]] else formula[[2]]
lhs <- function(formula) 
  if (length(formula)==3) formula[[2]] else NULL
"rhs<-" = function(formula,value) {
  newformula <- formula
  newformula[[length(formula)]] <- value
  newformula
}
"lhs<-" <- function(formula,value) {
  if (length(formula)==2)
    as.formula(as.call(c(formula[[1]],value,formula[[2]])))
  else {
    newformula <- formula
    newformula[[2]] <- value
    newformula
  }
}

## numerically calculate the partial gradient \partial func_j \over \partial x_i
## (dim(grad(func,x)) == c(length(x),length(func(x)))
grad <- function(func,x,...) # would shadow numDeriv::grad()
  {
    h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
    temp <- x+h
    h.hi <- temp-x
    temp <- x-h
    h.lo <- x-temp
    twoeps <- h.hi+h.lo
    nx <- length(x)
    ny <- length(func(x,...))
    if (ny==0L) stop("Length of function equals 0")
    df <- if(ny==1L) rep(NA, nx) else matrix(NA, nrow=nx,ncol=ny)
    for (i in 1L:nx) {
      hi <- lo <- x
      hi[i] <- x[i] + h.hi[i]
      lo[i] <- x[i] - h.lo[i]
      if (ny==1L)
        df[i] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
      else df[i,] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
      }
    return(df)
  }
## numerically calculate the gradient \partial func_i \over \partial x_i
## length(grad(func,x)) == length(func(x)) == length(x)
grad1 <- function(func,x,...)
  {
    h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
    temp <- x+h
    h.hi <- temp-x
    temp <- x-h
    h.lo <- x-temp
    twoeps <- h.hi+h.lo
    ny <- length(func(x,...))
    if (ny==0L) stop("Length of function equals 0")
    (func(x+h, ...) - func(x-h, ...))/twoeps
  }
## predict lpmatrix for an lm object
lpmatrix.lm <- 
  function (object, newdata, na.action = na.pass) {
    tt <- terms(object)
    if (!inherits(object, "lm")) 
      warning("calling predict.lm(<fake-lm-object>) ...")
    if (missing(newdata) || is.null(newdata)) {
      X <- model.matrix(object)
    }
    else {
      Terms <- delete.response(tt)
      m <- model.frame(Terms, newdata, na.action = na.action, 
                       xlev = object$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses"))) 
        .checkMFClasses(cl, m)
      X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    }
    X
  }
## fun: takes coef as its first argument
## requires: coef() and vcov() on the object
numDeltaMethod <- function(object,fun,gd=NULL,...) {
  coef <- coef(object)
  est <- fun(coef,...)
  Sigma <- vcov(object)
  if (is.null(gd))
      gd <- grad(fun,coef,...)
  ## se.est <- as.vector(sqrt(diag(t(gd) %*% Sigma %*% gd)))
  se.est <- as.vector(sqrt(colSums(gd* (Sigma %*% gd))))
  data.frame(Estimate = est, SE = se.est)
}
"coef<-" <- function (x, value) 
  UseMethod("coef<-")
predictnl <- function (object, ...) 
  UseMethod("predictnl")
"coef<-.default" <- function(x,value) {
    x$coefficients <- value
    x
}
predictnl.default <- function(object,fun,newdata=NULL,gd=NULL,...)
  {
    ## link=c(I,log,sqrt),invlink=NULL
    ## link <- match.arg(link)
    ## if (is.null(invlink))
    ##       invlink <- switch(deparse(substitute(link)),I=I,log=exp,sqrt=function(x) x^2)
    if (is.null(newdata) && !is.null(object$data))
      newdata <- object$data
    localf <- function(coef,...)
      {
        if ("coefficients" %in% names(object)) {
            object$coefficients <- coef
        } else if ("coef" %in% names(object)) {
            object$coef <- coef
        } else coef(object) <- coef
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
  }
setMethod("predictnl", "mle2", function(object,fun,newdata=NULL,gd=NULL,...)
  {
    if (is.null(newdata) && !is.null(object@data))
      newdata <- object@data
    localf <- function(coef,...)
      {
        object@fullcoef <- coef # changed from predictnl.default()
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,gd=gd,...)
  })
## setMethod("predictnl", "mle", function(object,fun,gd=NULL,...)
##   {
##     localf <- function(coef,...)
##       {
##         object@fullcoef = coef # changed from predictnl.default()
##         fun(object,...)
##       }
##     numDeltaMethod(object,localf,gd=gd,...)
##   })
predict.formula <- function(formula,data,newdata,na.action,type="model.matrix") 
{
  mf <- match.call(expand.dots = FALSE)
  type <- match.arg(type)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  xlevels <-.getXlevels(mt, mf)
  mfnew <- model.frame(mt, newdata, na.action=na.action, xlev=xlevels)
  if (!is.null(cl <- attr(mt, "dataClasses"))) .checkMFClasses(cl, mfnew)
  model.matrix(mt, mfnew, contrasts=contrasts)
}
`%call+%` <- function(left,right) call("+",left,right)
##
bread.stpm2 <- function (x, ...) {
  rval <- vcov(x) * nrow(x@y)
  dimnames(rval) <- list(names(coef(x)), names(coef(x)))
  return(rval)
}
estfun.stpm2 <- function(obj, weighted=FALSE, ...) {
  rr <- t(grad(obj@logli,coef(obj)))
  colnames(rr) <- names(coef(obj))
  if (weighted)
    rr <- rr * obj@weights
  rr
}
applyTapplySum <- function(X,index) apply(X, 2, function(col) tapply(col, index, sum))
meat.stpm2 <- function (x, adjust = FALSE, cluster=NULL, ...) 
{
    psi <- estfun.stpm2(x, ...)
    k <- NCOL(psi)
    n <- NROW(psi)
    if (!is.null(cluster))
        psi <- applyTapplySum(as.matrix(psi),cluster)
    rval <- crossprod(as.matrix(psi))/n
    if (adjust) 
        rval <- n/(n - k) * rval
    rownames(rval) <- colnames(rval) <- colnames(psi)
    return(rval)
}
sandwich.stpm2 <- 
function (x, bread. = bread.stpm2, meat. = meat.stpm2, ...) 
{
    if (is.function(bread.)) 
        bread. <- bread.(x)
    if (is.function(meat.)) 
        meat. <- meat.(x, ...)
    n <- NROW(estfun.stpm2(x))
    return(1/n * (bread. %*% meat. %*% bread.))
}
incrVar <- function(var,increment=1) {
  ##var <- deparse(substitute(var))
  ##function(data) "$<-"(data,var,"$"(data,var)+increment) # FAILS
  n <- length(var)
  if (n>1 && length(increment)==1)
    increment <- rep(increment,n)
  function(data) {
    for (i in 1:n) {
      data[[var[i]]] <- data[[var[i]]] + increment[i]
    }
    data
  }
}
cloglog <- function(x) log(-log(x))
cexpexp <- function(x) exp(-exp(x))
setOldClass("terms")
setClassUnion("listOrNULL",c("list","NULL"))
setClassUnion("nameOrcall",c("name","call"))
setClassUnion("nameOrcallOrNULL",c("name","call","NULL"))
##setClassUnion("numericOrNULL",c("numeric","NULL"))
setOldClass("Surv")
setOldClass("lm")
expit <- function(x) {
    ifelse(x==-Inf, 0, ifelse(x==Inf, 1, 1/(1+exp(-x))))
}
logit <- function(p) {
    ifelse(p==0, -Inf, ifelse(p==1, Inf, log(p/(1-p))))
} # numerical safety for large values?
## check: weights
##
## adapted from ordinal::drop.coef
which.dim <- function (X, silent = TRUE) 
{
    stopifnot(is.matrix(X))
    silent <- as.logical(silent)[1]
    qr.X <- qr(X, tol = 1e-07, LAPACK = FALSE)
    if (qr.X$rank == ncol(X)) 
        return(TRUE)
    if (!silent) 
        message(gettextf("design is column rank deficient so dropping %d coef", 
            ncol(X) - qr.X$rank))
    return(qr.X$pivot[1:qr.X$rank])
}

## link families
link.PH <- list(link=function(S) log(-log(as.vector(S))),
                ilink=function(eta) exp(-exp(as.vector(eta))),
                gradS=function(eta,X) -exp(as.vector(eta))*exp(-exp(as.vector(eta)))*X,
                h=function(eta,etaD) as.vector(etaD)*exp(as.vector(eta)),
                H=function(eta) exp(as.vector(eta)),
                gradh=function(eta,etaD,obj) obj$XD*exp(as.vector(eta))+obj$X*as.vector(etaD)*exp(as.vector(eta)),
                gradH=function(eta,obj) obj$X*exp(as.vector(eta)))
link.PO <- list(link=function(S) -logit(as.vector(S)),
                ilink=function(eta) expit(-as.vector(eta)),
                gradS=function(eta,X) -(exp(as.vector(eta))/(1+exp(as.vector(eta)))^2)*X,
                H=function(eta) log(1+exp(as.vector(eta))),
                h=function(eta,etaD) as.vector(etaD)*exp(as.vector(eta))*expit(-as.vector(eta)),
                gradh=function(eta,etaD,obj) {
                    as.vector(etaD)*exp(as.vector(eta))*obj$X*expit(-as.vector(eta)) -
                        exp(2*as.vector(eta))*obj$X*as.vector(etaD)*expit(-as.vector(eta))^2 +
                            exp(as.vector(eta))*obj$XD*expit(-as.vector(eta))
                    },
                gradH=function(eta,obj) obj$X*exp(as.vector(eta))*expit(-as.vector(eta)))
link.probit <-
    list(link=function(S) -qnorm(as.vector(S)),
         ilink=function(eta) pnorm(-as.vector(eta)),
         gradS=function(eta,X) -dnorm(-as.vector(eta))*X,
         H=function(eta) -log(pnorm(-as.vector(eta))),
         h=function(eta,etaD) dnorm(as.vector(eta))/pnorm(-as.vector(eta))*as.vector(etaD),
         gradh=function(eta,etaD,obj) {
             -as.vector(eta)*obj$X*dnorm(as.vector(eta))*as.vector(etaD)/pnorm(-as.vector(eta)) +
                 obj$X*dnorm(as.vector(eta))^2/pnorm(-as.vector(eta))^2*as.vector(etaD) +
                     dnorm(as.vector(eta))/pnorm(-as.vector(eta))*obj$XD
         },
         gradH=function(eta,obj) obj$X*dnorm(as.vector(eta))/pnorm(-as.vector(eta)))
link.AH <- list(link=function(S) -log(S),
                ilink=function(eta) exp(-as.vector(eta)),
                gradS=function(eta,X) -as.vector(exp(-as.vector(eta)))*X,
                h=function(eta,etaD) as.vector(etaD),
                H=function(eta) as.vector(eta),
                gradh=function(eta,etaD,obj) obj$XD,
                gradH=function(eta,obj) obj$X)
link.AO <- function(theta) { # Aranda-Ordaz
    if (theta==0) {
        return(link.PH)
        } else {
            list(link = function(S) log((S^(-theta)-1)/theta),
                 ilink = function(eta) exp(-log(theta*exp(as.vector(eta))+1)/theta),
                 gradS = function(eta,X) -as.vector(exp(as.vector(eta))*exp(-log(theta*exp(as.vector(eta))+1)/theta)/(1+theta*exp(as.vector(eta))))*X,
                 H = function(eta) log(theta*exp(as.vector(eta))+1)/theta,
                 h = function(eta,etaD) exp(as.vector(eta))*as.vector(etaD)/(theta*exp(as.vector(eta))+1),
                 gradH = function(eta,obj) exp(as.vector(eta))*obj$X/(1+theta*exp(as.vector(eta))),
                 gradh = function(eta,etaD,obj) {
                     eta <- as.vector(eta)
                     etaD <- as.vector(etaD)
                     ((theta*exp(2*eta)+exp(eta))*obj$XD+exp(eta)*etaD*obj$X) /
                         (theta*exp(eta)+1)^2
                 })
        }
    }
## fd <- function(f,x,eps=1e-5) (f(x+eps)-f(x-eps))/2/eps
fd <- function(f,x,eps=1e-5)
    t(sapply(1:length(x),
             function(i) {
                 upper <- lower <- x
                 upper[i]=x[i]+eps
                 lower[i]=x[i]-eps
                 (f(upper)-f(lower))/2/eps
             }))
## test code for the link functions
if (FALSE) {
    Xstar <- cbind(1,1:3) # beta[1] + beta[2]*t
    betastar <- c(-4,0.5)
    XDstar <- cbind(0,Xstar[,2])
    etastar <- as.vector(Xstar %*% betastar)
    etaDstar <- as.vector(XDstar %*% betastar)
    obj <- list(X=Xstar,XD=XDstar)
    for(link in list(rstpm2:::link.PH,rstpm2:::link.PO,rstpm2:::link.probit,rstpm2:::link.AH,rstpm2:::link.AO(.3))) {
        print(rstpm2:::fd(function(beta) link$ilink(Xstar%*%beta), betastar)-t(link$gradS(etastar,Xstar)))
        print(rstpm2:::fd(function(beta) link$h(Xstar%*%beta, XDstar%*%beta), betastar)-t(link$gradh(etastar,etaDstar,obj)))
        print(rstpm2:::fd(function(beta) link$H(Xstar%*%beta), betastar)-t(link$gradH(etastar,obj)))
    }
}

## copy of bbmle:::strwrapx
strwrapx <-
function (x, width = 0.9 * getOption("width"), indent = 0, exdent = 0, 
          prefix = "", simplify = TRUE, parsplit = "\n[ \t\n]*\n", 
          wordsplit = "[ \t\n]") 
{
  if (!is.character(x)) 
    x <- as.character(x)
  indentString <- paste(rep.int(" ", indent), collapse = "")
  exdentString <- paste(rep.int(" ", exdent), collapse = "")
  y <- list()
  plussplit = function(w) {
    lapply(w, function(z) {
      plusloc = which(strsplit(z, "")[[1]] == "+")
      plussplit = apply(cbind(c(1, plusloc + 1), c(plusloc, 
                                                   nchar(z, type = "width"))), 1, function(b) substr(z, 
                                                                                                     b[1], b[2]))
      plussplit
    })
  }
  z <- lapply(strsplit(x, parsplit), function(z) {
    lapply(strsplit(z, wordsplit), function(x) unlist(plussplit(x)))
  })
  for (i in seq_along(z)) {
    yi <- character(0)
    for (j in seq_along(z[[i]])) {
      words <- z[[i]][[j]]
      nc <- nchar(words, type = "w")
      if (any(is.na(nc))) {
        nc0 <- nchar(words)
        nc[is.na(nc)] <- nc0[is.na(nc)]
      }
      if (any(nc == 0)) {
        zLenInd <- which(nc == 0)
        zLenInd <- zLenInd[!(zLenInd %in% (grep("\\.$", 
                                                words) + 1))]
        if (length(zLenInd) > 0) {
          words <- words[-zLenInd]
          nc <- nc[-zLenInd]
        }
      }
      if (length(words) == 0) {
        yi <- c(yi, "", prefix)
        next
      }
      currentIndex <- 0
      lowerBlockIndex <- 1
      upperBlockIndex <- integer(0)
      lens <- cumsum(nc + 1)
      first <- TRUE
      maxLength <- width - nchar(prefix, type = "w") - 
        indent
      while (length(lens) > 0) {
        k <- max(sum(lens <= maxLength), 1)
        if (first) {
          first <- FALSE
          maxLength <- maxLength + indent - exdent
        }
        currentIndex <- currentIndex + k
        if (nc[currentIndex] == 0) 
          upperBlockIndex <- c(upperBlockIndex, currentIndex - 
                                 1)
        else upperBlockIndex <- c(upperBlockIndex, currentIndex)
        if (length(lens) > k) {
          if (nc[currentIndex + 1] == 0) {
            currentIndex <- currentIndex + 1
            k <- k + 1
          }
          lowerBlockIndex <- c(lowerBlockIndex, currentIndex + 
                                 1)
        }
        if (length(lens) > k) 
          lens <- lens[-(1:k)] - lens[k]
        else lens <- NULL
      }
      nBlocks <- length(upperBlockIndex)
      s <- paste(prefix, c(indentString, rep.int(exdentString, 
                                                 nBlocks - 1)), sep = "")
      for (k in (1:nBlocks)) {
        s[k] <- paste(s[k], paste(words[lowerBlockIndex[k]:upperBlockIndex[k]], 
                                  collapse = " "), sep = "")
      }
      s = gsub("\\+ ", "+", s)
      yi <- c(yi, s, prefix)
    }
    y <- if (length(yi)) 
      c(y, list(yi[-length(yi)]))
    else c(y, "")
  }
  if (simplify) 
    y <- unlist(y)
  y
}
