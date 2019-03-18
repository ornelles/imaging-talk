###############################################################################
##
## Set 'home' to path with 'imaging-talk'
## The example here is specific for my office computer...
##
###############################################################################
  home <- "~/Documents/github/imaging-talk"

# library and helper functions
  if (!require(EBImage))
    stop("EBImage must be installed to proceed")
  setwd(home)
  source("helpers.R")

###############################################################################
##
## Gel labeling exercise with straightening, cropping
##
## This loads a skewed image and asks the interacts with the user to choose
## two points that define the horizontal. Next, the user chooses two points
## to crop the image. Next the leading edge of the 1st and 4th lane are 
## selected in order to calculate labels.  
##
###############################################################################

# read RNA gel image
  setwd(file.path(home, "rna"))
  img <- readImage("gel.tif")
  dev.new(width = 7.5, height = 4, xpos = 5, ypos = 5)
  plot(img)

# straighten image
  p <- locator(2, type = "l")	# choose 2 points along wells

  p <- lapply(p, diff)        # compute delta-y, delta-x
  a <- 180*atan2(p$y, p$x)/pi # compute angle in degrees
  img <- rotate(img, -a)      # rotate image
  plot(img)

# interactively crop (1st four lanes) and replot with a margin
  p <- locator(2, type = "p", pch = 3)  # top left, bottom right corner

  p <- toIdx(p)        # round, sort, convert to integer sequence
  img <- img[p$x, p$y]
  dev.off()
  dev.new(width = 6.5, height = 5, xpos = 5, ypos = 5)
  plot(img, margin = 80)

# select leading edge of 1st and 4th lanes to locate 4 labels
  w <- diff(locator(2)$x) / 3 # leading edge of 1st and 4th lanes

# create labels, compute coordinates, place labels and lines
  labs <- c("MW", "blank", "A", "B")
  xx <- 0:3 * w + w/2  # center of each lane 
  text(xx, -20, labs, cex = 1.5, xpd = TRUE)
  lines(w*c(1,4) + 20*c(1,-1), c(-40, -40), lwd = 2, xpd = TRUE)
  text(xx[3], -60, "RNA Samples", cex = 1.5, xpd = TRUE)

# draw lines around image based on dimensions of image
  dm <- dim(img) + 1
  rect(0, dm[2], dm[1], 0, border = 1, lwd = 4, lend = "square", xpd = TRUE)

# save as JPEG file
  dev.print(jpeg, file="rna figure.jpg", res=300, unit="in",
    width=dev.size()[1])
  dev.off()

###############################################################################
##
## Flatfield correction for microscope images
##
###############################################################################

# collect names of files in /flatfield
  setwd(file.path(home, "flatfield"))
  ff <- list.files()
  ff  # first image is image with no specimen in place (brightfield)

# read images from vector of filenames and note dimensions (1567 x 1567 x 3)
  img <- readImage(ff)
  img

# background remains constant
  dev.new(width = 6.5, height = 3.2, xpos=5, ypos=5)
  plot(img[,,2:3], all = T)  # two views of yeast

  plot(img[,,1], all = TRUE) # brightfield image with no specimen
  dev.off()

# flatfield adjustment by simply dividing image with "blank" brightfield
  adj <- img[,,2:3] / img[,,c(1,1)] # need to replicate blank with c(1,1)
  apply(adj, 3, range)    # note values above 1 are preserved
  adj <- normalize(adj)   # rescale values in Image data
  apply(adj, 3, range)
  z <- combine(img[,,2:3], adj)
  dev.new(width = 6.5, height = 6.5, xpos = 5, ypos = 5)
  plot(z, all = TRUE)
  dev.off()

###############################################################################
##
## Sea urchin embryo segmentation and analysis
##
###############################################################################

# read image from ImageJ website at NIH
  setwd(file.path(home, "embryos"))
  img.url <- "https://imagej.nih.gov/ij/images/embryos.jpg"
  img <- readImage(img.url)
# NOTE: sample is JPEG. JPEG files are not generally appropriate for analysis
  dev.new(width = 7.5, height = 6, xpos = 5, ypos = 5)
  plot(img)

# create grayscale image to extract features
  x <- channel(img, "gray")
  plot(x)

# the objective is to isolate individual embryos
# a simple threshold is often appropriate...but not always
  hist(x, log = "y")        # two populations evident 
  abline(v = 0.45, col = 2)  # approximate break point at 0.45
  dev.off()

# note blurred spot that is not an embryo but identified by global thresholding
  dev.new(width = 9, height = 4, xpos = 5, ypos = 5)
  z <- combine(x, x < 0.45)  # identify pixels darker than 0.45
  plot(z, all = TRUE)
  mtext("original", line=3, adj=0.2)
  mtext("global threshold", line=3, adj=0.8)
  mark <- c(x0 = 575, y0 = 938, x1 = 670, y1 = 849) # arrow coordinates
  do.call(arrows, c(as.list(mark), len = 0.1, col = "yellow", lwd = 2))
  mark <- mark + c(dim(x)[1], 0, dim(x)[1], 0)
  do.call(arrows, c(as.list(mark), len = 0.1, col = "yellow", lwd = 2))
  
##
## solution is local threshold along with fundamental image processing
##
# median filter to smooth noise
  xm <- medianFilter(x, 5)
  plot(combine(x, xm), all = TRUE)
  mtext("original", line=3, adj=0.2)
  mtext("median filter", line=3, adj=0.8)

# local threshold with default 5 x 5 rectangle and threshold of 0.01
# here, a local threshold applied to dark objects creates a white border
# around the object, creating a mask larger than the object
  xt <- thresh(xm)
  plot(combine(xm, xt), all = TRUE) # bright values above local background
  mtext("median filter", line=3, adj=0.2)
  mtext("local threshold", line=3, adj=0.8)

# fill holes in binary objects produced by thresh()
  xf <- fillHull(xt)
  plot(combine(xt, xf), all = TRUE)
  mtext("local threshold", line=3, adj=0.2)
  mtext("fill hull", line=3, adj=0.8)

# opening = erode pixels on edges followed by dilating pixels on edges 
  xo <- opening(xf, makeBrush(11, "disc"))  # erode then dilate
  plot(combine(xf, xo), all = TRUE)
  mtext("fill hull", line=3, adj=0.2)
  mtext("opening", line=3, adj=0.8)

# distance map - replace each pixel with the distance to nearest background value
  xd <- distmap(xo)         
  plot(combine(xo, normalize(xd)), all = TRUE)
  mtext("opening", line=3, adj=0.2)
  mtext("distance map", line=3, adj=0.8)

# separate joined objects by 'watershed' algorithm
  mask <- watershed(xd)
  plot(combine(toRGB(normalize(xd)), colorLabels(mask)), all = TRUE)
  mtext("distance map", line=3, adj=0.2)
  mtext("watershed", line=3, adj=0.8)

# if we want to further expand the mask, one could do so with the following code
  xo <- dilate(mask, makeBrush(21, "disc"))
  xd <- distmap(xo)
  big.mask <- watershed(xd)
  plot(combine(colorLabels(mask),colorLabels(big.mask)), all = TRUE)
  mtext("original mask", line=3, adj=0.2)
  mtext("expanded", line=3, adj=0.8)

# examine mask and a couple objects defined by mask
  print(mask, short = TRUE) # image object
  range(mask)       # range is no longer 0 - 1
  table(mask)       # number of pixels for each is the area
  plot(mask == 10)  # show specific objects (10 and 12)
  plot(mask == 12)
	dev.off()

# EBImage feature: stack objects then and display objects with default background
  dev.new(width = 8, height = 5.5, xpos = 5, ypos = 5)
  plot(img)
  stk <- stackObjects(mask, img)  # apply mask in 'mask' to color image in 'img'
  plot(stk, all = TRUE, nx = 7)

# better background color from selected region
# surprisingly, bg.col accepts a vector of 3 values as RGB values
  plot(img)
  p <- locator(2, type = "p", pch = 3)  # select corners for background

  idx <- toIdx(p)
  bg <- apply(img[idx$x, idx$y, ], 3, median)
  bg  # returns three values for red, green and blue
  stk <- stackObjects(mask, img, bg.col = bg)
  plot(stk, all = TRUE, nx = 7)

# even better background from median value of entire image
  plot(img)
  bg <- apply(img, 3, median) # most of the image is background
  bg  # nearly same as selected
  stk <- stackObjects(mask, img, bg.col = bg) # hidden 'feature' in EBImage
  plot(stk, all = TRUE, nx = 7)
  dev.off()

# extract features with computeFeatures family of functions in EBImage
# use the features of the objects to classify
  FS <- computeFeatures.shape(mask)
  FB <- computeFeatures.basic(mask, x, basic.quantile = 0.95)
  FM <- computeFeatures.moment(mask, x)
  FH <- computeFeatures.haralick(mask, x)

# examine matrices returned by computeFeatures...
  FS[1:3,]  # quantifies object (mask) shape properties
  FB[1:3,]  # quantifies basic intensity values in target
  FM[1:3,]  # quantifies image moments with or without reference to intensity
  FH[1:3,]  # Haralick features characterize pixel "texture" 26 parameters

# calculate correlation for Haralick features and display with heatmap 
  dev.new(width = 5, height = 4, xpos = 5, ypos = 5)
  heatmap(cor(FH))

# some correlated, some not...combine all 35 properties linked to intensity
  F <- cbind(FB, FM, FH)

# group "related" objects by clustering principle components
  pc <- prcomp(F, scale. = TRUE)
  summary(pc) # says 7 components explain >95% of variance
  plot(pc)    # shows that 6 may be enough but we'll go with 7
  centers <- 7
  k <- kmeans(pc$x, centers, nstart = 101)

# examine clusters
  dev.off()
  if (require(rgl)) {  # 3D interactive plotting may not work in R Studio
    open3d(windowRect=c(25, 25, 25+512, 25+512))
    plot3d(pc$x[,1:3], col = seq_len(centers)[k$cl], type = "s")
  }

  rgl.close()  # put away 3D window

# establish order based on clusters 
  k$cluster
  (ord <- order(k$cluster))

# for asthetics. create blank place holder from background value
  (bgCol <- rgb(bg[1], bg[2], bg[3])) # create RGB color string
  blank <- Image(bgCol, dim = dim(stk)[1:2], colormode = "Color")

# show original and then plot final form with labels
  dev.new(width = 8, height = 5.5, xpos = 5, ypos = 5)
  plot(img)
  plot(combine(stk[,,,ord], blank), all = TRUE, nx = 7)
  xx <- 10 + 0:6 * dim(stk)[1]        # positions along left columns
  yy <- 10 + 0:3 * dim(stk)[2]        # positions along top rows
  pos <- expand.grid(x = xx, y = yy)  # all top left corners
  pos <- head(pos, -1)
  text(pos$x, pos$y, k$cluster[ord], col = (1:6)[k$cluster[ord]])
  mtext("Sea Urchin Embryo Classification", 3, line = 2.5)
  mtext(paste("Image from:", img.url), side=1, line=3, cex=0.8, adj=0.95)

# save
  dev.print(jpeg, "embryo figure.jpg", width = 7, unit = "in", res = 200)
  dev.off()

###############################################################################
##
## Create contour map from raw hand written data
##
###############################################################################

# set working directory and source helper functions
  setwd(file.path(home, "contour map"))
  source(file.path(home, "helpers.R"))

# read image, convert to grayscale, crop
  img0 <- readImage("map assignment.jpg")
  dev.new(width = 5.5, height = 7, xpos = 5, ypos = 5)
  plot(img0)

  img <- channel(img0, "gray")
  img <- img[350:2154, 434:2278]  # predetermined crop coordinates
  dev.off()
  dev.new(width = 5.5, height = 5.5, xpos = 5, ypos = 5)
  plot(img)

# perform image segmentation with EBImage tools
  x <- normalize(1 - img)   # invert image to 'typical' photo
  xf <- medianFilter(x, 3)  # median
  xb <- gblur(xf, 1.5)      # Gaussian blur
  xt <- thresh(xb, 24, 24, 0.05)  # Local threshold to binary
  xo <- opening(xt)   # Clean up pixels at edges
  xm <- bwlabel(xo)   # Label segmented objects

# examine binary/thresholded image
  plot(xo)
  plot(xo[730:900, 960:1050], interp = FALSE) # want to find "circles"

################################################################################
##
## Switch to discussion of perimeter problem in EBImage
##
  squares <- readImage("squares.tif")
  plot(squares)
  text(c(256, 769), 256, labels = 2:1, cex = 3)
  sm <- bwlabel(squares)
  computeFeatures.shape(sm) # perimeter is goofy!
##
## compare with ImageJ (area and shape descriptors)
##

# look more closely at rotated square
  plot(squares[581:596, 248:262], interp = FALSE)
  726 * sqrt(2)   # shows failure to measure Euclidean distance
  dev.off()
################################################################################

##
## relevant to the problem at hand, which is to find little circles
##
  FM <- computeFeatures.moment(xm)  # for x, y position and (bad) eccentricity
  FS <- computeFeatures.shape(xm)   # for area
  df <- data.frame(x = FM[,"m.cx"], y = FM[,"m.cy"])
  df[1:2] <- lapply(df, round)  # keep only integer positions
  df$area <- FS[,"s.area"]      # add area measurement

# Explicit code to calculate circularity will be used because of 
# shortcomings in perimeter (and eccentricity) calculations with EBImage.
# The circularity of each object will be determined after cropping.
# Circularity will be used to select small, circular points in the image.
# NOTE: The native R code is very slow...minutes slow
  
  df$circ <- sapply(1:max(xm), function(i) circularity(crop(xm==i)))
  plot(circ ~ area, df)         # show circularity and area
  v <- locator(type = "l")      # draw polygon around points, right-click stops 

# identify the points in polygon and check against plot by looking for
# a tight cluster of small but not too small points close to circ = 1
  sel <- pnpoly(df[c("area", "circ")], v)  # identify points within polygon
  plot(1 - xt)      # show black on white image
  points(df[sel,], col = "red", pch = 16, cex = 1.2) # add red dots to check
  res <- df[sel,]   # keep only selected in 'res'

#
# interact with plot to highlight point and enter value associated with point
# for misidentified points, enter the value assigned to 'BAD' 
#
  BAD <- 0        # use this for misidentified points
  res$z <- NA     # create variable in result data frame
  plot(1 - xt)    # work interactively with the thresholded image
  for (i in 1:nrow(res)) { # run through each selected point
    points(res[i,], col = "red", cex = 2, lwd = 2)
    res$z[i] <- scan(n = 1, quiet = TRUE)
    points(res[i,], col = "white", cex = 2, lwd = 2) 
  }
  res <- res[!res$z == BAD, ]	# exclude bad points (if any)

# (To save time, use previously saved values if 52 were found...)
  if (sum(sel) == 52) {	# have to link "}" and "else" in scripts
    message("Found exactly 52 points!")
    res$z <- readRDS("z values.Rda")
  } else {
      warning("Something is not right. Try again...")
  }
  head(res)

# invert y coordinates in res to correspond to image coordinates
  last <- max(res$x, res$y) # images read top-down, plots bottom-up
  res$y <- last - res$y + min(res$y)
  opar <- par(mar = c(1, 1, 3, 1))
  plot(y ~ x, res, pch = 16, frame = FALSE, axes = FALSE, xlab = "", ylab = "")
  text(res$x, res$y, res$z, pos = 3, offset = 0.5, cex=2/3, xpd = TRUE)

##
## Two of many methods to interpolate spatial values
##
# method 1: local polynomial regression with loess()
  fm1 <- loess(z ~ x * y, data = res, span = 1/5)

# create prediction data frame with every 20th position
  xp <- yp <- seq(1, last, 20)  
  new.data <- expand.grid(x = xp, y = yp)

# predict and plot contours from this direct model
  zp <- predict(fm1, new.data)  # predict yields matrix
  surf1 <- list(x = xp, y = yp, z = zp) # object for contour()

# make pretty features and plot
  plevels <- seq(230, 320, 5)   # contour lines at every 5 meters
  opar <- par(mar = c(1, 1, 3, 1))
  contour(surf1, levels = plevels, lty = c(1, 3), col = "red", asp = 1)
  title(main = "Raw loess fit")
  points(res$x, res$y, col = "steelblue")
  text(res$x, res$y, res$z, adj = c(0.5, -1), cex = 2/3, col = "steelblue")

# method 2: probably the best method with Akima's bivariate interpolation and
# smooth surface fitting algorithm from 1974 (need to install akima package)
  if(!require(akima))
  	warning("'akima' package must be installed to proceed")

# create prediction values at every 20 "pixels" in xp, yp
  lim <- range(res$x, res$y, na.rm = TRUE)
  xp <- yp <- seq(lim[1], lim[2], 20)

# Akima's function (provides rather few choices for "best" fit)
# using non-linear fit with extrapolation here
  surf2 <- interp(res$x, res$y, res$z, xp, yp, linear = FALSE, extrap = TRUE)

# plot as before with alternating line styles
  opar <- par(mar=c(1, 1, 3, 1))
  contour(surf2, levels = plevels, lty = c(1, 3), col = 2, asp = 1)
  points(res$x, res$y, col = "steelblue")
  text(res$x, res$y, res$z, cex = 0.6, col = "steelblue", adj = c(0.5, -1))
  title(main = "Akima fit (nonlinear) with extrapolation")

  dev.off()

###############################################################################
##
## Virus titer - two color immunofluorescence staining
##
###############################################################################

# This requires installation of the virustiter package
  if (!require("virustiter"))
    stop("the 'virustiter' package must be installed to proceed")

# set working directory
  setwd(home)

# load selected image pairs of DNA and viral protein
  f <- system.file("extdata", "by_stack/file005.tif", package = "virustiter")
  img <- readImage(f)  # two pairs of images

# Two pairs of images showing infected cells stained by DAPI for the nucleus and
# and by a fluorescent antibody for the e1a viral protein. The objective is to
# determine how many of the cells (identified by DAPI) also express the viral
# protein and are therefore infected.
  z <- normalize(img)  # 12-bit images will appear dim
  dev.new(width = 7, height = 5.3, xpos = 5, ypos = 5)
  plot(z, all = TRUE) 

# separate nuc and e1a
  nuc <- img[,,seq(1,4,2)]
  e1a <- img[,,seq(2,4,2)]

# reproduce color micrographs
  x <- normalize(nuc)
  y <- normalize(e1a, separate = FALSE) 
  N <- rgbImage(blue = x, green = 0.3 * x)  # mixture of blue and green
  E <- rgbImage(red = 0.2 * y, green = y)   # mixture of green and red
  z <- combine(N, E)    # combine single color and merged
  z <- z[,,,c(1,3,2,4)] # rorder to place on same row
  plot(z, all = T)

# one function to get nuclear masks (median, blur, threshold, fill, etc.)
  nm0 <- nucMask(nuc)
  z <- combine(N, colorLabels(nm0))[,,,c(1,3,2,4)]
  plot(z, all = T)  # not obvious, but small fragments captured

# calculate area of nuclei to identify excessively small fragments
  area <- lapply(1:2, function(i) computeFeatures.shape(nm0[,,i])[,1])
  plot(density(unlist(area)))
  rug(unlist(area))
  abline(v = 75, col = 2) # reasonable break point
  small <- lapply(area, function(v) which(v < 75))
  large <- lapply(area, function(v) which(v >= 75))
  lengths(small) # how many small nuclei?
  lengths(large) # how many large nuclei?

# remove small objects and save in new mask 
  nm <- rmObjects(nm0, small)

# examine the nuclei that were excluded for being too small
  nms <- rmObjects(nm0, large)
  obj <- paintObjects(nms, N, col = c(2, 2))
  plot(obj, all = T, nx = 1) # all fragments around the edges

# a more effective way to trim the nuclear mask of small and large
# objects and remove those touching the edge is with trimMask() function:
  dev.new(width = 7, height = 5.3, xpos = 5, ypos = 5)
	nmx <- trimMask(nm0, cutoff = c(75, 500), border = 1)
	plot(colorLabels(nmx), all = TRUE, nx = 1)
	dev.off()

# nuclear masks used to determine the staining intensity in e1a image
# faint fluorescence is enhanced by gamma and brightness adjustments
  obj <- paintObjects(nm, E^0.4 - 0.1, col = "magenta")
  plot(obj, all = TRUE)
  mtext("Nuclei on E1A (B/C adjusted)", 3)

# loop through each image pair and collect properties in a data frame
  df <- data.frame()
  for (i in seq_len(dim(nm)[3])) {
    area <- computeFeatures.shape(nm[,,i])[,1]
    XY <- computeFeatures.moment(nm[,,i])[,1:2]
    val <- computeFeatures.basic(nm[,,i], e1a[,,i])[,1]
    df <- rbind(df, data.frame(x = XY[,1], y = XY[,2], area, val))
  }

# plot the intensity as a densityplot and then interact with the plot
# to click on the point that separates the negative from positive cells
  plot(density(log(df$val)))
  rug(log(df$val))
  cutoff <- locator(1, type = "p", pch = 3)$x # waits for the operator
  (cutoff <- exp(cutoff))

# score positive samples in data frame
  df$pos <- df$val > cutoff
  xtabs(~ pos, df)

# use this value to create a mask for positive cells
  val <- lapply(1:2, function(i) computeFeatures.basic(nm[,,i], e1a[,,i])[,1])
  neg <- lapply(val, function(v) which(v < cutoff))
  nm.pos <- rmObjects(nm, neg)

# use gamma adjustment of 0.4 with brightness adjustment again
  obj <- paintObjects(nm.pos, y^0.4 - 0.1, col = "gray")
  plot(obj, all = TRUE)
  mtext("Positive Cells (B/C adjusted)", 3)

##
## All assembled in collection of code to do this rather automatically
##
  fpd <- system.file("extdata", "by_stack/phenoData.csv", package = "virustiter")
	path <- dirname(f)          # parental directory
  list.images(path)           # shows 8 total multi-layered TIF files
  img <- getImages(path)      # read paired images from each file 
  df <- parseImages(img)      # read all images in enclosing directory
  (pd <- read.csv(fpd))       # read phenotype data (about infection)
  df <- mergePdata(pd, df)    # join phenotype and image data
  bg <- getBgnd(df)           # programmatically determine cutoff
  df <- score(df, bg)         # add positives to data
  (res <- tally(df))          # assembly results data frame
  fm <- getFit(res)           # perform Poisson fit with glm()
  fm                          # glm fit
  plotFit(fm)                 # show fit with results: but problem with high moi

##
## Done!
## David Ornelles, Winston-Salem, December 11, 2017
##