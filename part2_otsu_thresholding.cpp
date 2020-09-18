#include <iostream>
#include <algorithm>
#include <assert.h> 
#include <CImg.h>
#include <CCfits/CCfits>
 
using namespace cimg_library;
using namespace CCfits;
using namespace std;
 
void readFile(CImg<float> & inImg, const string & inFilename, long * outBitPix = 0) {
  std::unique_ptr<FITS> pInfile(new FITS(inFilename, Read, true));
  PHDU & image = pInfile->pHDU(); 
 
  if (outBitPix) {
    *outBitPix = image.bitpix();
  }
 
  inImg.resize(image.axis(0) /*x*/, image.axis(1) /*y*/, 1/*z*/, 1 /*1 color*/);
  
  // NOTE: At this point we assume that there is only 1 layer.
  std::valarray<unsigned long> imgData;
  image.read(imgData);
  cimg_forXY(inImg, x, y) { inImg(x, inImg.height() - y - 1) = imgData[inImg.offset(x, y)]; }  
}
 
int main(int argc, char *argv[]) {
  CImg<float> img;
  long bitPix = 0;
  readFile(img, argv[1], & bitPix);
 
  // Display initial image
  CImgDisplay imgDisp(img, "Click a point");
  while (! imgDisp.is_closed()) {
    imgDisp.wait();
  }
  
  CImg<> hist = img.get_histogram(pow(2.0, bitPix));
 
  float sum = 0;
  cimg_forX(hist, pos) { sum += pos * hist[pos]; }
 
  float numPixels = img.width() * img.height();
  float sumB = 0, wB = 0, max = 0.0;
  float threshold1 = 0.0, threshold2 = 0.0;
  
  cimg_forX(hist, i) { 
    wB += hist[i];
 
    if (! wB) { continue; }    
 
    float wF = numPixels - wB;
    
    if (! wF) { break; }
    
    sumB += i * hist[i];
 
    float mF = (sum - sumB) / wF;
    float mB = sumB / wB;
    float diff = mB - mF;
    float bw = wB * wF * pow(diff, 2.0);
    
    if (bw >= max) {
      threshold1 = i;
      if (bw > max) {
         threshold2 = i;
      }
      max = bw;            
    }
  } // end loop
  
  float th = (threshold1 + threshold2) / 2.0;
  CImg<> & thImg = img.threshold(th); 
 
  // Display result
  CImgDisplay imgDisp2(thImg, "Result");
  while (! imgDisp2.is_closed()) {
    imgDisp2.wait();
  }
  return 0;
}
