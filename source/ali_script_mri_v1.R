# |----------------------------------------------------------------------------------|
# | Project: MRI processing in xenograft mice breast canser tumor model              |
# | Script: DICOM mage processing                                                    |
# | Coordinator: Ali                                                                 |
# | Author: Davit Sargsyan                                                           |
# | Created: 10/17/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# Sources: 
# 1. VivoQuant for drawing contours on DICOM files
# http://www.vivoquant.com/

# sink(file = "tmp/log_ali_dcemri_v1.txt")
date()

# Load packages----
library(help = oro.dicom)
require(oro.dicom)
require(oro.nifti)
?create3D

require(dcemriS4)
library(help = "dcemriS4")


# dcmList <- dicomSeparate(system.file("hk-40", package="oro.dicom"))

dcmList <- dicomSeparate(path = "C:\\ali_images")
dcmImage <- create3D(dcmList)

?oro.dicom::create4D

# bladder
image(dcmImage[110:160,60:120,10],
      col = grey(0:64/64), 
      axes = FALSE,
      xlab = "",
      ylab = "")

tiff(filename = "img1.tiff",
     height = 10,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
image(dcmImage[,,10],
      col = grey(0:64/64), 
      axes = FALSE,
      xlab = "",
      ylab = "")
graphics.off()

# imagePositionPatient <- attributes(dcmImage)$ipp
# 
# dSL <- abs(diff(imagePositionPatient[,3]))
# plot(dSL, ylim=range(range(dSL) * 1.5, 0, 10), xlab="Image", ylab="mm",
#      main="Difference in Slice Location")
# 
# dcmList <- dicomSeparate(path = "C:\\ali_images",
#                          pixelData=FALSE)
# dcmImage <- create3D(dcmList, pixelData=FALSE)
# tiff(filename = "img2.tiff",
#      height = 10,
#      width = 8,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# image(dcmImage[,,10],
#       col = grey(0:64/64), 
#       axes = FALSE,
#       xlab = "",
#       ylab = "")
# graphics.off()

fname <- paste("C:\\ali_images\\",
               dir("C:\\ali_images"),
                   sep = "")
fname
abdo <- dicomInfo(fname[1])
names(abdo)
head(abdo$hdr)

# require(DATforDCEMRI)
# library(help = "DATforDCEMRI")
# ?DATforDCEMRI::DAT
# tmp <- DAT()
# 
# DAT.checkData()

class(dcmList)
class(dcmImage)
?oro.dicom::dicom2nifti
niftiImage <- oro.dicom::dicom2nifti(dcm = dcmList,
                                     datatype = 4,
                                     mode = "integer")
dim(niftiImage)
tiff(filename = "img3.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
image(niftiImage)
graphics.off()
?dcemri.lm
mask <- array(data = FALSE,
              dim = c(256,
                      256,
                      24))

image(niftiImage[110:160, 
                 140:200,
                 13],
      col = grey(1:64/64))

mask[110:160,
     140:200,
     ] <- TRUE
sum(mask)/(24*256^2)

dim(niftiImage)
fit.lm <- dcemri.lm(conc = niftiImage,
                    time = 1:24,
                    img.mask = mask)
fit.lm
summary(fit.lm)

tiff(filename = "img4.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
overlay(niftiImage,
        ifelse(mask, 
               fit.lm$ktrans,
               NA))
graphics.off()

gc()

# NOTE: Ali will send correct masks that he created manually.

library(reshape2)
library(rgl)
dim(dcmImage)
class(dcmImage)
M=melt(dcmImage)
head(M)

summary(M$value)
M <- M[M$value > 300, ]
M$color <- 1 - M$value/max(M$value)


points3d(M$Var1,M$Var2,M$Var3,
         col = grey(M$color),
         alpha = 0.5)
aspect3d(1, 1, 1)
         