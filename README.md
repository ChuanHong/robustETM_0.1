#Welcome to version 0.1 of robustETM package!

Installation
=============

Run following commands in R to install the source package directly from github: 

install.packages("devtools")

library(devtools)

install_github("robustETM_0.1", username="ChuanHong")

Note: The package contains C code that requires compilation. The compilers are installed by default on Linux systems. Windows users need to have "Rtools" installed (see https://www.biostat.wisc.edu/~kbroman/Rintro/Rwinpack.html for details). Mac users, depends on the system settings, likely need to install Xcode. 


Sample code
=============
library("robustETM")

data(pseudo_dat) 

myresult=PLEMT(dat.case[c(1:2),], dat.ctrl[c(1:2),], cc=2, niter=3, distn="beta") 

myresult

