## Synopsis
This repository contains the brief PowerPoint file, scripts, support code, and data files used for a talk on image analysis at the Winston-Salem R Users Group MeetUp held on December 11, 2017. The talk demonstrated image processing with the `EBImage` package from the bioconductor. This repository has been revised on March 2019 to make use of my `virustiter` package, version 0.0.2.0.

First, be sure that the latest version of R is installed before installing the `EBImage` package from the Bioconductor. This is accomplished with the latest version of `biocLite.R` through the following code. If you have used `biocLite.R` in the distant past, it may be necessary to uninstall the old version and re-install the latest version. 
```
source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")
```

To run the final example in `case studies script.R`, the `virustiter` package need to be downloaded and installed. This is most easily accomplished with `install` from the `devtools` package. The `latticeExtra` package is also installed for use by `virstiter`. 
```
install.packages("devtools")
install.packages("latticeExtra")
library(devtools)
install_github("ornelles/virustiter")
```
If the `install_github` function fails. Download and unzip the `virustiter` package. Set the working director to resulting folder and execute the `install` command from `devtools`:
```
setwd("<virustiter>") # directory containing the unzipped virustiter repository.
library("devtools")
install()
```

After downloading and unzipping the repository from GitHub [imaging-talk] (https://github.com/ornelles/imaging-talk), the two R scripts (`imaging_demo script.R` and `case studies script.R`) can be used to recreate the talk. In both of these scripts, the variable named `home` needs to changed to the path of the local copy of this repository. Note that these scripts are meant to be executed line-by-line and are not meant to be sourced. They require frequent input from the operator. A few additional packages are required. All of these are wrapped in a call to `require()` with a message indicating which package needs to be installed. 

Many thanks are due to the creators and maintainers of `EBImage` (Gregoire Pau, Florian Fuchs, Oleg Sklyar, Michael Boutros and Wolfgang Huber) as well as Wayne Rasband, creator of ImageJ. And special thanks to Mason DeCamillis for looking after the Winston-Salem R Users group (www.meetup.com/Winston-Salem-R-Users-Group)!

revised March 3, 2019

## Contact
David Ornelles  
Wake Forest School of Medicine  
Department of Microbiology and Immunology  
575 Patterson Ave  
Winston-Salem, NC 27101

## License
GPL-3
