# |----------------------------------------------------------------------------------|
# | Project: MRI processing in xenograft mice breast canser tumor model              |
# | Script: DICOM mage processing                                                    |
# | Coordinator: Ali                                                                 |
# | Author: Davit Sargsyan                                                           |
# | Created: 10/17/2018                                                              |
# | Modified (DS, 10/17/2018): imported and processed 4D data (previously only 3D)   |
# |----------------------------------------------------------------------------------|
# Sources: 
# 1. VivoQuant for drawing contours on DICOM files
# http://www.vivoquant.com/
# 2. https://github.com/cran/DATforDCEMRI/blob/master/DESCRIPTION
# 3. Working with the DICOM and NIfTI Data Standards in R
# file:///C:/Users/dsargsy/Downloads/v44i06.pdf

# Can we use ImageJ (http://imagej.net/Welcome) as a free alternative for VivoQuant?

# sink(file = "tmp/log_ali_dcemri_v1.txt")
date()

# Load packages----
require(oro.dicom)
require(oro.nifti)
require(dcemriS4)

# Load data----
dcmList <- dicomSeparate(path = "data/dicom_combined")
t1 <- dicomTable(dcmList$hdr)

dcmImage <- create4D(dcmList)
dim(dcmImage)
# 256 256  24   6
# 24 slides at 6 timepoints, 256X256 pixels each

tiff(filename = "tmp/img1.tiff",
     height = 10,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
image(dcmImage[,,10, 1],
      col = grey(0:64/64), 
      axes = FALSE,
      xlab = "",
      ylab = "")
graphics.off()

fname <- paste(getwd(),
               "data/dicom_combined",
               dir("data/dicom_combined"),
               sep = "/")
fname
abdo <- dicomInfo(fname[1])
names(abdo)
head(abdo$hdr)

# Convert DICOM to Nifti format
niftiImage <- oro.dicom::dicom2nifti(dcm = dcmList,
                                     datatype = 4,
                                     mode = "integer",
                                     DIM = 4)
dim(niftiImage)
# 256 256  24   6
class(niftiImage)

# By time----
for (i in 1:dim(niftiImage)[4]) {
  tiff(filename = paste("tmp/nifti_by_time",
                        i,
                        ".tiff",
                        sep = ""),
       height = 10,
       width = 10,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  image(niftiImage,
        w = i)
  graphics.off()
}

Modeling
mask <- array(data = FALSE,
              dim = c(256,
                      256,
                      24,
                      6))

image(niftiImage[110:160, 
                 140:200,
                 13,
                 1],
      col = grey(0:64/64))

mask[110:160,
     140:200,,] <- TRUE
sum(mask)/(6*24*256^2)

dim(niftiImage)
fit.lm <- dcemri.lm(conc = niftiImage,
                    img.mask = mask)
gc()
fit.lm
summary(fit.lm)

tiff(filename = "tmp/overlay_time1.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
overlay(niftiImage,
        ifelse(mask[, , , 1], 
               fit.lm$ktrans,
               NA),
        w = 1)
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
         