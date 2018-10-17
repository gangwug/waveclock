This repository contains a copy of waveclock package maintained by Tom Price. 

## Introduction
waveclock is an R function designed to assess the period and amplitude of cycling cell luminescence data. The function reconstructs the modal frequencies from a continuous wavelet decomposition of the luminescence data using the 'crazy climbers'algorithm described in "Practical Time-Frequency Analysis: Gabor and Wavelet Transforms with an Implementation in S", by Rene Carmona, Wen L. Hwang and Bruno Torresani, Academic Press, 1998.

## Usage
```r
# install 'Rwave' package (if 'Rwave' is not installed yet)
install.packages("Rwave")

# install 'devtools' in R(>3.0.2) (if 'devtools' is not installed yet) 
install.packages("devtools")

# install waveclock
devtools::install_github('gangwug/waveclock')

```

## For more information
Price T.S., Baggs J.E., Curtis A.M., Fitzgerald G.A., Hogenesch J.B., WAVECLOCK: wavelet analysis of circadian oscillation.Bioinformatics. 2008, 24(23):2794-5. 
