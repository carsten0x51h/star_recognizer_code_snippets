# star_recognizer
The C++ "star-recognizer" automatically determines stars in an image and calculates centoid, HFD and FWHM.


## Introduction
For some time now I am looking into the field of night sky image processing. It is a huge topic and there are already a lot of different approaches and solutions to many of the existing problems.

The overall result of this project in part 7 is a program which takes an astro image as input and detects all the stars in the image and calculates the centroid, the HFD and the FWHM of each star. All algorithms and routines developed in the parts 1-6 are used in part 7. Still, all parts are completely independent and can be compiled separately. The very simple Makefile shows the required libraries for each part.

A more detailed description of this project can be found here: https://www.lost-infinity.com/night-sky-image-processing-part-7-automatic-star-recognizer/


### Part 1
Part 1 starts with a problem every night sky photographer and astronomer is aware of: noise. I did some research and came across “anisotropic diffusion”. In short, anisotropic diffusion is a non-linear and space-variant transformation i.e. the transformation depends on the image content. The effect is that the resulting images preserve linear structures while at the same time smoothing is made along these structures. Luckily the CImg library already contains an implementation for anisotropic diffusion. For reading and writing the FITS file the CCFits library is used.


### Part 2
Part 2 is about image binarization using the Otsu method. In this step the algorithm decides which pixels belong to the background and which pixels belong to a star. This process is called thresholding i.e. converting a grey-scale image into a binary image. There are many different thresholding algorithms out there. All have there advantages and drawbacks. I found ImageJ very helpful to test different thresholding algorithms.


### Part 3
Part 3 is about finding all the pixels which belong to a star – also known as star clustering. The binary image from part 2 is now further examined in a process called clustering. In simple words the pixels which belong together (i.e. the neighbours) are grouped together. The result of this clustering process is a list of pixel groups which belong together – i.e. a list of stars and the pixels belonging to each star. The implemented algorithm is based on the paper “Improving night sky star image processing algorithm for star sensors“ from Mohammad Vali Arbabmir et. al. However, in detail it is slightly different.


### Part 4
The goal of part 4 is to find the center of each star – the star centroid. o keep it simple, in this example I apply the algorithm to calculate the star centroid only for one single star which I load from a FITS file. Mohammad Vali Arbabmir et. al.. proposed this algorithm in “Improving night sky star image processing algorithm for star sensors”. The paper mentioned above will probably help a lot to get a better understanding of the implementation.


### Part 5
In part 5 the FWHM determination of the star takes place. To do so the centroid position from part 4 is required.  The FWHM (Full Width Half Maximum) value is a measurement for the star width. In this part I write about the so called curve-fitting which is helpful to determine the FWHM value from such an image. The C++ implementation below uses the GSL (GNU scientific library) to make use of the Levenberg-Marquart algorithm. The implementation makes use of the Traits concept. The idea is to supply the curve details (derivations, parameters, …) as template parameters to the CurveFitTmplT template. I made this design decision since the curve type (e.g. a Gaussian curve) usually can be made already at compile-time.


### Part 6
Part 6 looks into another measure for the star focus - the Half Flux Diameter (HFD). It was invented by Larry Weber and Steve Brady. The main two arguments for using the HFD is robustness and less computational effort compared to the FWHM approach. The HFD is defined as the diameter of a circle that is centered on the unfocused star image in which half of the total star flux is inside the circle and half is outside.


### Part 7
The previous parts show the different steps which are required to finally create a "star recognizer" application. In this part I put all the pieces together and develop an “automatic star recognizer” which takes an astro-image as input and outputs a list of the recognized stars with their HFD and FWHM values.



## Required libraries

In order build the binaries the following libraries are required

 * CImg-devel.x86_64
 * CCfits-devel.x86_64
 * cfitsio-devel.x86_64
 * gsl-devel.x86_64


## Build
To build all mini-projects, just run

```make all```

to build only a particular part simply run

```make part<N>```


## Program execution

 * Part 1: ```./part1_anisotropic_diffusion test_data/part1_part2_anisotropic_diffusion.fits```
 * Part 2: ```./part2_otsu_thresholding test_data/part1_part2_anisotropic_diffusion.fits```
 * Part 3: ```./part3_star_clustering test_data/part3_bin_img.fits```
 * Part 4: ```./part4_star_centroid test_data/part4_star.fits```
 * Part 5: ```./part5_fwhm_levenmberg_marquart_curve_fitting```
 * Part 6: ```./part6_hfd test_data/part6_focus_star1.fits```
 * Part 7: ```./part7_star_recognizer test_data/part7_star_recognizer.fits```


## Further documentation

The specific parts are documented inthe following blog posts:

 * Part 1: http://www.lost-infinity.com/night-sky-image-processing-part-1-noise-reduction-using-anisotropic-diffusion-with-c-using-the-cimg-library/
 * Part 2: http://www.lost-infinity.com/night-sky-image-processing-part-2-image-binarization-using-otsus-thresholding-algorithm-a-simple-c-implementation-using-the-cimg-library/
 * Part 3: http://www.lost-infinity.com/night-sky-image-processing-part-3-star-clustering-a-simple-c-implementation/
 * Part 4: http://www.lost-infinity.com/night-sky-image-processing-part-4-calculate-the-star-centroid-with-sub-pixel-accuracy/
 * Part 5: http://www.lost-infinity.com/night-sky-image-processing-part-5-measuring-fwhm-of-a-star-using-curve-fitting-a-simple-c-implementation/
 * Part 6: http://www.lost-infinity.com/night-sky-image-processing-part-6-measuring-the-half-flux-diameter-hfd-of-a-star-a-simple-c-implementation/
 * Part 7: https://www.lost-infinity.com/night-sky-image-processing-part-7-automatic-star-recognizer/
