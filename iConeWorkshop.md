---
title: "iCone_Workshop"
author: "Esteban del Castillo"
date: 'August, 2026'
output: html_document
---

## Demo data
> This R markdown script show the fundamentals of iCone data processing. The demo MSI data used here can be obtained in the link 
[Raw_data](https://doi.org/10.34810/data1744)

> After the download, the following files should appear:

For a tolerance of 33.33 ppm:
```
    231211_Au_P_MBr_cblm_30k.imzML
    231211_Au_P_MBr_cblm_30k.ibd 
```
For a tolerance of 16.66 ppm:
```
    231211_Au_P_MBr_cblm_60k.imzML
    231211_Au_P_MBr_cblm_60k.ibd 
```
For a tolerance of 8.33 ppm:
```
    231211_Au_P_MBr_cblm_120k.imzML 
    231211_Au_P_MBr_cblm_120k.ibd
```
> From now on we will assume that they have been copied to the /home/MSI/ folder.

## Obtaining the peak matrix

> The first step is to create a list with the desired parameters. 
> For example, to obtain the peak matrix corresponding to the 30k sample, with a signal-to-noise ratio of 3, using the 'estnoise_mad' algorithm for noise estimation, 
with a tolerance of 33.33 ppm and a minimum number of supporting pixels of 10%.
```
> params30 <-list("SNR"=3, "tolerance"=33.33, "minPixelsSupport"=10, "noiseMethod"="estnoise_mad")
```

> The second step is to call the getPeakMatrix() function.
```
> pk30k <-getPeakMatrix(dataFiles="/MSI/231211_Au_P_MBr_cblm_30k.imzML", params=params30);
> #It returns a peak matrix assuming the maximum mass range, for all pixels in the sample and making use of all available CPU cores (minus one).
>
> #Alternatively, we can limit the mass range. For example, between 500 and 1000 Da:
> params30 <-list("SNR"=3, "tolerance"=33.33, "minPixelsSupport"=10, "noiseMethod"="estnoise_mad", initMass=500, finalMass=1000)
> pk30k <-getPeakMatrix(dataFiles="/MSI/231211_Au_P_MBr_cblm_30k.imzML", params=params30);
>
> #It is also possible to pass a list with the pixels of interest:
> params30 <-list("SNR"=3, "tolerance"=33.33, "minPixelsSupport"=10, "noiseMethod"="estnoise_mad", pxList=c(1:100, 200, 300, 400:500))
> pk30k <-getPeakMatrix(dataFiles="/MSI/231211_Au_P_MBr_cblm_30k.imzML", params=params30);
```


## Precision of the centroids
> To determine the precision of the m/z measurement reported by the getPeakMatrix() function, the results are compared with the same tissue samples analyzed at higher resolutions. Specifically, samples with 30k and 60k resolutions are compared with samples analyzed at 120k resolution.

### obtaining the peak matrices
```
> params30 <-list("SNR"=3, "tolerance"=33.33, "minPixelsSupport"=10, "noiseMethod"="estnoise_mad")
> pk30k <-getPeakMatrix(dataFiles="/MSI/231211_Au_P_MBr_cblm_30k.imzML", params=params30);

> params60 <-list("SNR"=3, "tolerance"=16.66, "minPixelsSupport"=10, "noiseMethod"="estnoise_mad")
> pk60k <-getPeakMatrix(dataFiles="/MSI/231211_Au_P_MBr_cblm_60k.imzML", params=params60);

> params120 <-list("SNR"=3, "tolerance"=8.33, "minPixelsSupport"=1, "noiseMethod"="estnoise_mad")
> pk120k <-getPeakMatrix(dataFiles="/MSI/231211_Au_P_MBr_cblm_120k.imzML", params=params120);
```
### Precision

> To compare the results, we use the `statisticalQuality()` function. It reports statistical values regarding the deviation between the two m/z vectors being compared. This deviation is determined from the nearest centroids in the high-resolution sample to those in the low-resolution sample.
```
> statisticalQuality(mzRef, mzTest, tolerance) 
>arguments: 
>  mzRef      is the high-resolution centroid vector.
>  mzTest     is the low-resolution centroid vector.
>  tolerance  indicates the minimum distance for the result to be considered a false positive (expressed in ppm).

>report: 
>  refSize    size of the high-resolution centroid vector.
>  testSize   size of the low-resolution centroid vector.
>  mean       mean of the deviations between centroids.
>  sigma      standard deviation of the deviations between centroids.
>  median     median of the deviations between centroids.
>  FP         False Positive: number of low-resolution centroids with a deviation exceeding 1.5 times the tolerance relative to the nearest high-resolution centroid.
```
> Use
```
> statisticalQuality(pk120k$mass, pk30k$mass, 33.33)
  refSize=7842  testSize=1013     mean=3.4347  sigma=4.1666  median=1.7812  FP:1 (0.1%)

> statisticalQuality(pk120k$mass, pk30k$mass, 16.66)
  refSize=7842  testSize=1496     mean=1.6092  sigma=2.3784  median=0.6295  FP:0 (0.0%)
``` 
### Summary
To compare, for example, the 30k sample with the 120k sample, we will do the following:
```
> params30 <-list("SNR"=3, "tolerance"=33.33, "minPixelsSupport"=10, "noiseMethod"="estnoise_mad")
> pk30k <-getPeakMatrix(dataFiles="/MSI/231211_Au_P_MBr_cblm_30k.imzML", params=params30);
> params120 <-list("SNR"=3, "tolerance"=8.33, "minPixelsSupport"=1, "noiseMethod"="estnoise_mad")
> pk120k <-getPeakMatrix(dataFiles="/MSI/231211_Au_P_MBr_cblm_120k.imzML", params=params120);
> statisticalQuality(pk120k$mass, pk30k$mass, 33.33)
  refSize=7842  testSize=1013     mean=3.4347  sigma=4.1666  median=1.7812  FP:1 (0.1%)
```
