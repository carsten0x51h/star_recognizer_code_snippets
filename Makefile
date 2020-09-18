all: part1 part2 part3 part4 part5 part6 part7

part1:
	g++ part1_anisotropic_diffusion.cpp -o part1_anisotropic_diffusion -lCCfits -lX11 -lpthread -lcfitsio

part2:
	g++ part2_otsu_thresholding.cpp -o part2_otsu_thresholding -std=c++11 -lcfitsio -lCCfits -lX11 -lpthread

part3:
	g++ part3_star_clustering.cpp -o part3_star_clustering -std=c++0x -lX11 -lpthread -lCCfits -lcfitsio

part4:
	g++ part4_star_centroid.cpp -o part4_star_centroid -std=c++0x -lCCfits -lcfitsio -lX11 -lpthread

part5:
	g++ part5_fwhm_levenmberg_marquart_curve_fitting.cpp -o part5_fwhm_levenmberg_marquart_curve_fitting -std=c++0x -lgsl -lgslcblas

part6:
	g++ part6_hfd.cpp -o part6_hfd -lCCfits -lcfitsio -lX11 -lpthread

part7:
	g++ part7_star_recognizer.cpp -o part7_star_recognizer -std=c++0x -lX11 -lpthread -lCCfits -lcfitsio -lgsl -lgslcblas

clean:
	-rm part1_anisotropic_diffusion
	-rm part2_otsu_thresholding
	-rm part3_star_clustering
	-rm part4_star_centroid
	-rm part5_fwhm_levenmberg_marquart_curve_fitting
	-rm part6_hfd
	-rm part7_star_recognizer
