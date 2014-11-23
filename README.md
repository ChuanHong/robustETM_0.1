robustETM_0.1
=============
Welcome to version 0.1 of robustETM package!

Installation
-------------

1. Automated Installation (get the current development version from github using devtools)
----------------------------------
install.packages("devtools")

library(devtools)

dev_mode(on=T)

install_github("robustETM_0.1", username="ChuanHong")

2. Alternative Installation 
-----------------------------------
Visit https://github.com/ChuanHong/robustETM.gz

Manually dowload robustETM_0.1.tar.gz from https://github.com/ChuanHong/robustETM.gz

Install the package: R CMD INSTALL robustETM_0.1.tar.gz

Sample code
--------------
library("robustETM")
data(methyl) 
myresult=PLEMT(dat.case[c(1:2),], dat.ctrl[c(1:2),], cc=2, niter=3, distn="beta") 
