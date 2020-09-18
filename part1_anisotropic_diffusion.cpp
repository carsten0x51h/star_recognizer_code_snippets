#include <iostream>
#include <CImg.h>
#include <CCfits/CCfits>
 
using namespace cimg_library;
using namespace CCfits;
using namespace std;
 
void readFile(CImg<float> & inImg, const string & inFilename) {
  std::unique_ptr<FITS> pInfile(new FITS(inFilename, Read, true));
  PHDU & image = pInfile->pHDU(); 
  inImg.resize(image.axis(0) /*x*/, image.axis(1) /*y*/, 1/*z*/, 1 /*1 color*/);
  
  // NOTE: At this point we assume that there is only 1 layer.
  std::valarray<unsigned long> imgData;
  image.read(imgData);
  cimg_forXY(inImg, x, y) { inImg(x, inImg.height() - y - 1) = imgData[inImg.offset(x, y)]; }  
}
 
void writeFile(const CImg<float> & inImg, const string & inFilename) {
  long naxis = 2;
  long naxes[2] = { inImg.width(), inImg.height() };
  std::unique_ptr<FITS> pFits;
  pFits.reset(new FITS(string("!") + inFilename , USHORT_IMG , naxis , naxes) );
  
  // NOTE: At this point we assume that there is only 1 layer.
  long nelements = std::accumulate(& naxes[0], & naxes[naxis], 1, std::multiplies<long>());
  std::valarray<int> array(nelements);
  cimg_forXY(inImg, x, y) { array[inImg.offset(x, y)] = inImg(x, inImg.height() - y -1); }
 
  long fpixel(1);
  pFits->pHDU().write(fpixel, nelements, array);
}
 
int main(int argc, char *argv[]) {
  CImg<float> img;
 
  try {
    readFile(img, argv[1]);
  } catch (FitsException & exc) {
    cerr << "Read FITS failed." << endl;
    return 1;
  }
 
  CImgDisplay imgDisp(img, "Original image");
  while (! imgDisp.is_closed()) {
    imgDisp.wait();
  }
  
  // http://cimg.sourceforge.net/reference/structcimg__library_1_1CImg.html
  CImg<float> & resImg = img.blur_anisotropic(30.0f, /*amplitude*/
                                               0.7f, /*sharpness*/
                                               0.3f, /*anisotropy*/
                                               0.6f, /*alpha*/
                                               1.1f, /*sigma*/
                                               0.8f, /*dl*/
                                               30,   /*da*/
                                               2,    /*gauss_prec*/
                                               0,    /*interpolation_type*/
                                               false /*fast_approx*/
                                               );
 
  CImgDisplay imgDisp2(resImg, "Result image");
  while (! imgDisp2.is_closed()) {
    imgDisp2.wait();
  }
  
  try {
    writeFile(resImg, "test_res.fits");
  } catch (FitsException &) {
    cerr << "Writing FITS failed." << endl;
    return 1;
  }
  
  return 0;
}
