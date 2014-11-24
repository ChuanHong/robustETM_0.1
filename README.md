
Welcome to version 0.1 of robustETM package!

Installation
=============

#Windows users need to install "Rtools"; MAC users need to install "Xcode". 
install.packages("devtools")

library(devtools)

dev_mode(on=T)

install_github("robustETM_0.1", username="ChuanHong")



Sample code
=============
library("robustETM")

data(pseudo_dat) 

myresult=PLEMT(dat.case[c(1:2),], dat.ctrl[c(1:2),], cc=2, niter=3, distn="beta") 
