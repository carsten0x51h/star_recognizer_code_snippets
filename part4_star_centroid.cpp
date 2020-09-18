#include <iostream>
#include <assert.h>
#include <CImg.h>
#include <CCfits/CCfits>
 
using namespace cimg_library;
using namespace CCfits;
using namespace std;
 
// See http://heasarc.gsfc.nasa.gov/fitsio/ccfits/html/cookbook.html
void readFile(CImg<float> & inImg, const string & inFilename) {
 
  std::unique_ptr<FITS> pInfile(new FITS(inFilename, Read, true));
  PHDU & image = pInfile->pHDU(); 
  inImg.resize(image.axis(0) /*x*/, image.axis(1) /*y*/, 1/*z*/, 1 /*1 color*/);
  
  // NOTE: At this point we assume that there is only 1 layer.
  std::valarray<unsigned long> imgData;
  image.read(imgData);
  cimg_forXY(inImg, x, y) { inImg(x, inImg.height() - y - 1) = imgData[inImg.offset(x, y)]; }  
}
 
float calcIx2(const CImg<float> & img, int x) {
  float Ix = 0;
  cimg_forY(img, y) { Ix += pow(img(x, y), 2.0) * (float) x; }
  return Ix;
}
 
float calcJy2(const CImg<float> & img, int y) {
  float Iy = 0;
  cimg_forX(img, x) { Iy += pow(img(x, y), 2.0) * (float) y; }
  return Iy;
}
 
// Calculate Intensity Weighted Center (IWC) 
void calcIntensityWeightedCenter(const CImg<float> & inImg, float * outX, float * outY) {
  assert(outX && outY);
  assert(inImg.width() == inImg.height());
  const size_t L = inImg.width();
  
  // Determine weighted centroid - See http://cdn.intechopen.com/pdfs-wm/26716.pdf
  float Imean2 = 0, Jmean2 = 0, Ixy2 = 0;
  
  for(size_t i = 0; i < L; ++i) {
    Imean2 += calcIx2(inImg, i);
    Jmean2 += calcJy2(inImg, i);
    cimg_forY(inImg, y) { Ixy2 += pow(inImg(i, y), 2.0); }
  }
 
  *outX = Imean2 / Ixy2;
  *outY = Jmean2 / Ixy2;
}
 
void calcSubPixelCenter(const CImg<float> & inImg, float * outX, float * outY,
                         size_t inNumIter = 10 /*num iterations*/) {
  // Sub pixel interpolation
  float c, a1, a2, a3, a4, b1, b2, b3, b4;
  float a1n, a2n, a3n, a4n, b1n, b2n, b3n, b4n;
 
  assert(inImg.width() == 3 && inImg.height() == 3);
 
  b1 = inImg(0, 0); a2 = inImg(1, 0); b2 = inImg(2, 0);
  a1 = inImg(0, 1);  c = inImg(1, 1); a3 = inImg(2, 1);
  b4 = inImg(0, 2); a4 = inImg(1, 2); b3 = inImg(2, 2);
 
  for (size_t i = 0; i < inNumIter; ++i) {
    float c2 = 2 * c;
    float sp1 = (a1 + a2 + c2) / 4;
    float sp2 = (a2 + a3 + c2) / 4;
    float sp3 = (a3 + a4 + c2) / 4;
    float sp4 = (a4 + a1 + c2) / 4;
    
    // New maximum is center
    float newC = std::max({ sp1, sp2, sp3, sp4 });
    
    // Calc position of new center
    float ad = pow(2.0, -((float) i + 1));
 
    if (newC == sp1) {
      *outX = *outX - ad; // to the left
      *outY = *outY - ad; // to the top
 
      // Calculate new sub pixel values
      b1n = (a1 + a2 + 2 * b1) / 4;
      b2n = (c + b2 + 2 * a2) / 4;
      b3n = sp3;
      b4n = (b4 + c + 2 * a1) / 4;
      a1n = (b1n + c + 2 * a1) / 4;
      a2n = (b1n + c + 2 * a2) / 4;
      a3n = sp2;
      a4n = sp4;
 
    } else if (newC == sp2) {
      *outX = *outX + ad; // to the right
      *outY = *outY - ad; // to the top
 
      // Calculate new sub pixel values
      b1n = (2 * a2 + b1 + c) / 4;
      b2n = (2 * b2 + a3 + a2) / 4;
      b3n = (2 * a3 + b3 + c) / 4;
      b4n = sp4;
      a1n = sp1;
      a2n = (b2n + c + 2 * a2) / 4;
      a3n = (b2n + c + 2 * a3) / 4;
      a4n = sp3;
    } else if (newC == sp3) {
      *outX = *outX + ad; // to the right
      *outY = *outY + ad; // to the bottom
 
      // Calculate new sub pixel values
      b1n = sp1;
      b2n = (b2 + 2 * a3 + c) / 4;
      b3n = (2 * b3 + a3 + a4) / 4;
      b4n = (2 * a4 + b4 + c) / 4;
      a1n = sp4;
      a2n = sp2;
      a3n = (b3n + 2 * a3 + c) / 4;
      a4n = (b3n + 2 * a4 + c) / 4;
    } else {
      *outX = *outX - ad; // to the left
      *outY = *outY + ad; // to the bottom   
 
      // Calculate new sub pixel values
      b1n = (2 * a1 + b1 + c) / 4;
      b2n = sp2;
      b3n = (c + b3 + 2 * a4) / 4;
      b4n = (2 * b4 + a1 + a4) / 4;
      a1n = (b4n + 2 * a1 + c) / 4;
      a2n = sp1;
      a3n = sp3;
      a4n = (b4n + 2 * a4 + c) / 4;
    }
 
    c = newC; // Oi = Oi+1
 
    a1 = a1n;
    a2 = a2n;
    a3 = a3n;
    a4 = a4n;
 
    b1 = b1n;
    b2 = b2n;
    b3 = b3n;
    b4 = b4n;
  }
}
 
int main(int argc, char *argv[]) {
  CImg<float> img;
  float xc, yc;
 
  try {
    readFile(img, argv[1]);
  } catch (FitsException &) {
    cerr << "Read FITS failed." << endl;
    return 1;
  }
 
  // 1. Calculate the IWC
  calcIntensityWeightedCenter(img, & xc, & yc);
 
  // 2. Round xc, yc to nearest integer and then iteratively improve.
  int xi = floor(xc + 0.5);
  int yi = floor(yc + 0.5);
  
  CImg<float> img3x3 = img.get_crop(xi - 1 /*x0*/, yi - 1 /*y0*/, xi + 1 /*x1*/, yi + 1 /*y1*/);
 
  // 3. Interpolate using sub-pixel algorithm
  float xsc = xi, ysc = yi;
  calcSubPixelCenter(img3x3, & xsc, & ysc, 10 /*num iterations*/);
 
  cerr << "xc: " << xc << " --> xi: " << xi << ", yc: " << yc << " --> yi: " << yi << endl;
  cerr << "xsc: " << xsc << ", ysc: " << ysc << endl;
    
  return 0;
}
