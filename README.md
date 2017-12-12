## Synopsis
This repository contains the brief PowerPoint file, scripts, support code, and data files used for a talk on image analysis at the Winston-Salem R Users Group MeetUp held on December 11, 2017. The talk demonstrated image processing with the `EBImage` package from the bioconductor.

After downloading the repository, the two R scripts (`imaging_demo script.R` and `case studies script.R`) can be used to recreate the talk. It will be necessary to adjust the path to point to the local directory representing this repository. Note that these scripts are meant to be executed line-by-line and are not meant to be sourced. They require frequent input from the operator. A few additional packages are required. These are wrapped in a call to `require()` with a message indicating which package needs to be installed. 

Also included here is a directory representing a package that is still under development, `virustiter`. The package will be loaded locally with `load_all()` of the `devtools` package. This is necessary only to run the final example in `case studies script.R`.

Many thanks are due to the creators and maintainers of `EBImage` (Gregoire Pau, Florian Fuchs, Oleg Sklyar, Michael Boutros and Wolfgang Huber) as well as Wayne Rasband, creator of ImageJ. And special thanks to Mason DeCamillis for looking after the Winston-Salem R Users group (www.meetup.com/Winston-Salem-R-Users-Group)!

## Contact
David Ornelles  
Wake Forest School of Medicine  
Department of Microbiology and Immunology  
575 Patterson Ave  
Winston-Salem, NC 27101

## License
GPL-3
