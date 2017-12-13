###############################################################################
##
## toIdx - create integer index for Image object from two points
##
###############################################################################

	toIdx <- function(p) {
		p <- lapply(p, round)
		return(lapply(p, function(v) min(v[1:2]):max(v[1:2])))
	}

###############################################################################
##
## extract - extract and mark rectangular selections of an image
## return single image or list of images
## 
###############################################################################

extract <- function(x, plot = TRUE, col = "white") {
# process and check arguments
	if (!interactive())	# need an operator
		return(NULL)
	if(!is(x, "Image"))	# only Image objects
		stop ("'x' must be of class 'Image'")
	if (plot)			# plot image
		plot(x)
	if (is.null(col))	# accept either NULL or NA for no color
		col <- NA
# alert and repeat with operator in loop
	alarm()
	message("Define selections, right-click to stop")
	flush.console()
	ans <- list()	# hold extracted images in a list 
	repeat {
	# get and mark 2 corners, return when NULL
		p <- locator(2, type = "p", pch = 3, col = col)
		if (is.null(p))
			break
	# round values and draw rectangle
		p <- lapply(p, round)
		pp <- do.call(rbind, p)
		do.call(rect, c(as.list(pp), list(border = col)))
		idx <- toIdx(p)
	# assemble argument for `[` function
		arg <- list(x, idx$x, idx$y)
		other.dims <- dim(x)[-(1:2)]	# for variable dimensions
		for (n in other.dims)
			arg <- c(arg, list(seq_len(n)))
		sel <- do.call(`[`, arg)		# extracted image
		ans <- c(ans, list(sel))		# collect extracted images
	}
	if (length(ans) == 0)
		return(NULL)
	else if (length(ans) == 1)
		return(ans[[1]])
	else
		return(ans)
}

###############################################################################
##
## Improved perimeter calculator and related functions
## Based on code posted by Vincent Zoonekynd on Stack Overflow, April 9, 2013
##
## Functions:
##   circularity & perimeter (meant to be used)
##   S1, S2, S3, S4, aRea, eDge (support functions)
##
###############################################################################
# Necessary improvements on EBImage computeFeatures.shape
# considerably slower than existing code or ImageJ functions

# Functions to shift the image in z in one direction
	S1 <- function(z) cbind(rep(0,nrow(z)), z[,-ncol(z)] )
	S2 <- function(z) cbind(z[,-1], rep(0,nrow(z)) )
	S3 <- function(z) rbind(rep(0,ncol(z)), z[-nrow(z),] )
	S4 <- function(z) rbind(z[-1,], rep(0,ncol(z)) )

# EBimage functions to shift the image in z in one direction (not much faster)
#	S1 <- function(z) translate(z, c(1, 0))
#	S2 <- function(z) translate(z, c(-1, 0))
#	S3 <- function(z) translate(z, c(0, 1))
#	S4 <- function(z) translate(z, c(0, -1))

# aRea, eDge helper functions
	aRea <- function(z) sum(z)
	eDge <- function(z) z & !(S1(z) & S2(z) & S3(z) & S4(z))

# Improved perimeter function
	perimeter <- function(z) {
		e <- eDge(z)
	# add horizontal and vertical segments
		segs1 <- sum(e & S1(e)) + sum(e & S2(e) ) + sum(e & S3(e)) + sum(e & S4(e))
	# add diagonal segments
		segs2 <- sqrt(2)*(sum(e & S1(S3(e))) + sum(e & S1(S4(e))) +
		sum(e & S2(S3(e))) + sum(e & S2(S4(e))))
	# each segmented was counted twice, divide by 2
		return((segs1 + segs2)/2)
	}

# Improved circularity (eccentricity) function
	circularity <- function(z) {
		4*pi*aRea(z) / perimeter(z)^2
	}

###############################################################################
##
## crop - crop grayscale or binary image to non-zero values plus border
##
###############################################################################

crop <- function(img, border = 1, fill = 0) {
	img <- cbind(0, img, 0)   # ensure an edge exists!
	img <- rbind(0, img, 0)
	mask <- img > 0

	xb <- apply(mask, 1, Negate(any))	# blank rows
	xr <- rle(xb)
	nx <- xr$lengths * xr$values
	ix <- seq_len(dim(mask)[1])
	ix <- intersect(tail(ix, -head(nx, 1)), head(ix, -tail(nx, 1)))

	yb <- apply(mask, 2, Negate(any))	# blank columns
	yr <- rle(yb)
	ny <- yr$lengths * yr$values
	iy <- seq_len(dim(mask)[2])
	iy <- intersect(tail(iy, -head(ny, 1)), head(iy, -tail(ny, 1)))

	dm <- c(length(ix), length(iy))
	z <- Image(fill, dim = dm + 2 * border) 
	z[border + seq_along(ix), border + seq_along(iy)] <- img[ix, iy]
	return(z)
}

###############################################################################
##
## pnpoly - point inclusion in polygon test
##
## return logical vector indicating which points are within the
## (convex) polygon defined by vertices
##
## Arguments:
##   points - points to be tested in a form appropriate for xy.coords 
##   vertices - vertices of the polygon, also processed by xy.coords
##
## Based on original FORTRAN code by W. Randolph Franklin
##
###############################################################################

pnpoly <- function(points, vertices)
{
	pp <- xy.coords(points)
	vv <- xy.coords(vertices)
	nvert <- length(vv$x)

  # working function
	.fun <- function(x, y, vx, vy, nvert) {
		inside <- FALSE
		j <- nvert
		for (i in seq_len(nvert)) {
			if (((vy[i] > y) != (vy[j] > y)) &&
					(x < (vx[j] - vx[i]) * (y - vy[i]) / (vy[j] - vy[i]) + vx[i]))
				inside <- !inside
			j <- i
		}
		return(inside)
	}

  # apply to each pair of coordinates in points
	sapply(seq_len(lengths(points)[1]),
		function(k) .fun(pp$x[k], pp$y[k], vv$x, vv$y, nvert))
}