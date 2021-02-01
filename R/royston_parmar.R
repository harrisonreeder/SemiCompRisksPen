ns_d <- function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x))
{

  #inspired by rstpm2 package nsxD function, extending the ns function from splines
  #generate derivative basis of natural cubic spline function
  #I purposefully write this function as the derivative of ns, for use with Royston parmar model

  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)

  #temporarilty set aside all NA values
  if (nas <- any(nax))
    x <- x[!nax]

  #create a boolean flag for whether any of the given x values are outside of the boundaries
  outside <- if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  } else {
    if (length(x) == 1L)
      Boundary.knots <- x * c(7, 9)/8
    FALSE
  }

  #if df is given instead of actual knot locations, fill in knots based on quantiles
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - 1L - intercept
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d",
                       1L + intercept), domain = NA)
    }
    #notice this funny nesting, 'knots' inside is sequence of percentiles, and then reassigned as the quantiles of internal data
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots +
                         2L)[-c(1L, nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    }
  } else{ nIknots <- length(knots) }
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots)) #this code essentially puts two copies of each boundary on the ends, forming intervals
  if (any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      tt <- splines::splineDesign(Aknots, rep(k.pivot, sum(ol)), 4,
                         1)
      basis[ol, ] <- tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      tt <- splines::splineDesign(Aknots, rep(k.pivot, sum(or)), 4,
                         1)
      basis[or, ] <- tt
    }
    if (any(inside <- !outside))
      basis[inside, ] <- splines::splineDesign(Aknots, x[inside],
                                      4,derivs = 1)
  } else{ basis <- splines::splineDesign(Aknots, x, ord = 4L,derivs = 1) }
  const <- splines::splineDesign(Aknots, Boundary.knots, ord = 4L, derivs = c(2L,
                                                                     2L))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L),
                                                     drop = FALSE])#uses QR decomp to orthogonalize the basis
  n.col <- ncol(basis)

  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3L, knots = if (is.null(knots)) numeric() else knots,
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns_d", "basis", "matrix")
  basis
}
