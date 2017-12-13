###############################################################################
##
## Winston-Salem R Users Group Imaging talk
## David Ornelles
## December 11, 2017
##
## EBImage and ImageJ installation
##
# Download ImageJ from: imagej.nih.gov/ij/download.html
#
###############################################################################
#
# Install EBImage from Bioconductor with biocLite installer:
#  source("https://bioconductor.org/biocLite.R")
#  biocLite()
#  biocLite("EBImage")
#
###############################################################################
##
## Set 'home' to path with 'imaging-talk'
## Here, it runs from my local github directories 
##
###############################################################################
  home <- gsub("\\", "/", path.expand("~"), fixed = TRUE)
  home <- sub("/Documents", "", home)   # adjust for Windows
  home <- file.path(home, "Documents/github/imaging-talk")

# library and helper functions
  if (!require(EBImage))
    stop("EBImage must be installed to proceed")
  setwd(home)

###############################################################################
##
## Image processing basics
##
## Plot windows are opened at the upper left with a fixed size to accommodate
## demonstrating this on a little laptop computer. The position specification
## (xpos, ypos) may only work on Windows platforms. 
##
###############################################################################
#
# EBImage reads jpeg, tiff and png files
  img <- readImage("images/rabbit island.jpg")

# two methods of displaying images
  display(img) # JavaScript viewer in web browser
  dev.new(width = 7.5, height = 6, xpos = 5, ypos = 5)
  plot(img)    # same as display(img, method = "raster")

# images are matrix objects with special features from S4 system
  str(img)
  img                      # default for print(img)
  print(img, short = TRUE) # more compact representation
  apply(img, 3, range)     # apply to 3rd dimension
  img0 <- img
  colorMode(img0) <- Grayscale  # or "grayscale" or 0
  print(img0, short = TRUE)

# plot makes images accessible with base graphics 
  plot(img0, all = TRUE)
  text(x = c(101, 3301, 101), y = c(101, 101, 2468), adj = c(0, 1),
    label = c("Red channel", "Green channel", "Blue channel"),
    col = 2:4, cex = 1.2)

# images can be treated as matrix with few special features
  hist(img0[,,1])   # typical histogram
  hist(img)         # hist() has method for "Image"
  dev.off()         # windows sometimes benefits from clean device

###############################################################################
##
## Standard image manipulation
##
## mtext() is used throughout assuming that the default margins are in place
## i.e., par("mar") = c(5.1, 4.1, 4.1, 2.1)
##
###############################################################################
#
# Brightness - add or subtract
  dev.new(width = 7.5, height = 6, xpos = 5, ypos = 5); plot(img)
  plot(img + 0.25); mtext("Brighten by addition (+0.25)", side = 3, line = 3)
  plot(img - 0.25); mtext("Darken by subtraction (-0.25)", 3, 3)

# Contrast - multiply or divide
  plot(img); mtext("Original image", side = 3, line = 3)
  plot(img / 2); mtext("Divide (by 2) to decrease contrast", 3, 3)
  plot(img * 2); mtext("Multiply (by 2) to increase contrast", 3, 3)

# functions and values to help illustrate gamma adjustment
  x <- seq(1, dim(img)[1], len = 51) # points to illustrate effect 
  fy <- function(x) {                # map y values to screen
    x <- dim(img)[1] * x / max(x)
    dim(img)[2] - x * dim(img)[2]/dim(img)[1] + 1
  }
  dash <- paste(rep("-", 20), collapse = "")  # labels
  ylab <- paste0(dash, "  Displayed Value  ", dash, ">")
  xlab <- paste0(dash, "  Actual Value  ", dash, ">")

# Gamma adjustment (with illustration of transfer function)
  # gamma = 1 (no change)
  gamma = 1
  plot(img); mtext(paste("Gamma =", gamma), 3, 3)
  points(x, fy(x), col = "white", pch = 16, xpd = TRUE)
  mtext(xlab, 1, 2.5, col = "white"); mtext(ylab, 2, 3, col = "white")

  # gamma < 1 selectively lightens darker pixels  
  gamma <- 0.5
  plot(img^gamma); mtext(paste("Gamma =", gamma), 3, line = 3) 
  points(x, fy(x^gamma), col = "white", pch = 16, xpd = TRUE)
  mtext(xlab, 1, 2.5, col = "white"); mtext(ylab, 2, 3, col = "white")

  # gamma > 1 selectively darkens darker pixels  
  gamma <- 2.2
  plot(img^gamma); mtext(paste("Gamma =", gamma), 3) 
  points(x, fy(x^gamma), col = "white", pch = 16, xpd = TRUE)
  mtext(xlab, 1, 2.5, col = "white"); mtext(ylab, 2, 3, col = "white")

# Cropping - same methods used to subset an array
  dev.off()
  dev.new(width = 7.5, height = 6, xpos = 5, ypos = 5); plot(img)
  dim(img)
  island <- img[510:1900, 930:1490, ]
  plot(island)

# Rotate with adjustment to background color (bg.col)
  x <- rotate(island, 25, bg.col = "white")
  plot(x); mtext("Rotate +25 degrees", 3)

# Scale (resize)
  x <- resize(island, w = 512, h = 512)
  plot(x); mtext("Resized with identical X and Y dimensions", 2, 1)

# Flip about horizon or "flop" about vertical
  plot(flip(island)); mtext("Flip (horizontal flip)", 3)
  plot(flop(island)); mtext("Flop (vertical flop)", 3)
  x <- island + flop(island) # some values will be greater than 1
  x <- normalize(x)          # remaps to [0,1]
  plot(x); mtext("Fake new(s) island", 3)

# Filters (gaussian, median, all others)
  plot(gblur(island, 5)); mtext("Low Pass or Gaussian Blur, sigma = 5", 3)
  plot(medianFilter(island, 5)); mtext("Median Filter, radius = 5", 3)
  plot(island - gblur(island, 15)); mtext("High Pass filter effect", 3)

# Thresholding - matrix of logical values
  x <- channel(img, "luminance")  # grayscale based on luminance
  plot(x); mtext("Grayscale image based on luminance (perception)", 3, 3)
  xt <- x > 0.6 # logical Image
  xt
  plot(xt); mtext("Thresholded or Binary Image, pixel value > 0.6", 3, 3)
  dev.off()

###############################################################################
##
## 	Interactive Photo processing with EBImage
##
###############################################################################
#
# straighten with user interaction and locator()
  dev.new(width = 7.5, height = 6, xpos = 5, ypos = 5); plot(img)
  p <- locator(2, type = "l") # two points along horizon

  delta <- lapply(p, diff)    # need delta-y over delta-x
  angle <- 180 * atan2(delta$y, delta$x) / pi # convert radians to degrees
  img <- rotate(img, angle = - angle)
  plot(img)
  abline(h = locator(1)$y)    # test rotation with horizontal line

  (p <- locator(4))           # identify 4 corners of image
  (p <- lapply(p, round))     # round to integers
  (p <- lapply(p, sort))      # sort x, y values
  (p <- lapply(p, "[", 2:3))  # extract middle values
  idx <- lapply(p, function(v) seq(v[1], v[2])) # create sequence
  img2 <- img[idx$x, idx$y, ] # crop
  plot(img2)

# select top left, bottom right of little bay to create inset
  p <- locator(2)  # top left, bottom right

  p <- lapply(p, round) # convert to integers
  idx <- lapply(p, function(v) seq(v[1], v[2])) # index of pixels
  inset <- img2[idx$x, idx$y, ] # extracted inset
  plot(inset)

# enlarge inset by factor of 3, other dimension scales proportionately
  inset <- resize(inset, w = dim(inset)[1] * 3)

# choose point for lower right corner to place inset in original image
  plot(img2)
  p2 <- locator(1)  # choose lower right corner

  p2 <- lapply(p2, round)
  p2 <- unlist(p2)  # convert list to vector to use with dim()
  p2

# calculate upper left corner for inset
  dim(inset)
  p1 <- p2 - dim(inset)[1:2] + 1
  p1  # upper left corner
  p2  # lower right corner

# create index and check that target on original image is acceptable
  idx <- mapply(seq, from = p1, to = p2, SIMPLIFY = FALSE)
  plot(img2[idx$x, idx$y, ])

##
## discuss quirks of mapply()
##
  mapply(seq, c(1, 11, 21), c(3, 13, 23))
## vs
  mapply(seq, c(1, 11, 21), c(3, 13, 25))
## fix inconsistent (unexpected) output with SIMPLIFY Option 
  mapply(seq, c(1, 11, 21), c(3, 13, 23), SIMPLIFY = FALSE)
##
## return to talk
##

# copy image to new object and replace selected area
  img3 <- img2
  img3[idx$x, idx$y, ] <- inset
  plot(img3)

# dress up image with vector elements
  p1 <- as.list(p1) # for the convenience of using x, y names
  p2 <- as.list(p2)
  rect(p1$x - 1, p2$y - 1, p2$x, p1$y, border = "white",
    lwd = 3, lend = "square")

# use locator() twice with options to add lines
  for(i in 1:2) locator(2, type = "l", col = "white", lty = 2)

# save image as jpeg from device (not Image)
  dev.print(jpeg, file = "./images/rabbit island figure.jpg",
    width = dev.size()[1], unit = "in", res = 300)

  dev.off()