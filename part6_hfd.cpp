#include <iostream>
#include <assert.h>
#include <CImg.h>
#include <CCfits/CCfits>
 
using namespace std;
using namespace cimg_library;
using namespace CCfits;
 
void readFile(CImg<float> & inImg, const string & inFilename) {
  std::unique_ptr<FITS> pInfile(new FITS(inFilename, Read, true));
  PHDU & image = pInfile->pHDU(); 
  inImg.resize(image.axis(0) /*x*/, image.axis(1) /*y*/, 1 /*z*/, 1 /*1 color*/);
  
  // NOTE: At this point we assume that there is only 1 layer.
  std::valarray<unsigned long> imgData;
  image.read(imgData);
  cimg_forXY(inImg, x, y) { inImg(x, inImg.height() - y - 1) = imgData[inImg.offset(x, y)]; }
}
 
/**
 * Get all pixels inside a radius: http://stackoverflow.com/questions/14487322/get-all-pixel-array-inside-circle
 * Algorithm: http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
 */
bool insideCircle(float inX /*pos of x*/, float inY /*pos of y*/, float inCenterX, float inCenterY, float inRadius) {
  return (pow(inX - inCenterX, 2.0) + pow(inY - inCenterY, 2.0) <= pow(inRadius, 2.0));
}
 
/**
 * Expects star centered in the middle of the image (in x and y) and mean background subtracted from image.
 *
 * HDF calculation: http://www005.upp.so-net.ne.jp/k_miyash/occ02/halffluxdiameter/halffluxdiameter_en.html
 *                  http://www.cyanogen.com/help/maximdl/Half-Flux.htm
 *
 * NOTE: Currently the accuracy is limited by the insideCircle function (-> sub-pixel accuracy).
 * NOTE: The HFD is estimated in case there is no flux (HFD ~ sqrt(2) * inOuterDiameter / 2).
 * NOTE: The outer diameter is usually a value which depends on the properties of the optical
 *       system and also on the seeing conditions. The HFD value calculated depends on this
 *       outer diameter value.
 */
float calcHfd(const CImg<float> & inImage, unsigned int inOuterDiameter) {
  // Sum up all pixel values in whole circle
  float outerRadius = inOuterDiameter / 2;
  float sum = 0, sumDist = 0;
  int centerX = ceil(inImage.width() / 2.0);
  int centerY = ceil(inImage.height() / 2.0);
 
  cimg_forXY(inImage, x, y) {
    if (insideCircle(x, y, centerX, centerY, outerRadius)) {
      sum += inImage(x, y);
      sumDist += inImage(x, y) * sqrt(pow((float) x - (float) centerX, 2.0f) + pow((float) y - (float) centerY, 2.0f));
    }
  }
  // NOTE: Multiplying with 2 is required since actually just the HFR is calculated above
  return (sum ? 2.0 * sumDist / sum : sqrt(2.0) * outerRadius);
}
 
 
int main(int argc, char *argv[]) {
  // For all filenames specified
  const size_t numCols = 4;
  const size_t numRows = ceil((float) argc / (float) numCols);
  const size_t imgWidth = 65, imgHeight = 85; // Hardcoded image size
  CImg<unsigned char> containerImg(numCols * imgWidth, numRows * imgHeight, 1/*depth*/, 3 /*3 channels - RGB*/);
 
  for (size_t i = 1; i < argc; ++i) {
    CImg<float> img;
  
    try {
      readFile(img, argv[i]);
    } catch (FitsException &) {
      cerr << "Read FITS failed." << endl;
      return 1;
    }
 
    // Subtract mean value from image which is required for HFD calculation
    CImg<float> img2(img);
    double mean = img.mean();
    
    cimg_forXY(img, x, y) {
      if (img(x, y) < mean) {
        img2(x, y) = 0;
      } else {
        img2(x, y) = img(x, y) - mean;
      }
    }
 
    // Calc the HFD
    const unsigned int outerDiameter = 60;
    float hfd = calcHfd(img2, outerDiameter /*outer diameter in px*/);
 
    cerr << "File " << argv[i] << " -> hfd: " << hfd << endl;
 
    
    // Create RGB image from fits file to paint circle (just for visualization)
    CImg<unsigned char> rgbImg(img2.width(), img2.height(), 1 /*depth*/, 3 /*3 channels - RGB*/);
    float min = img2.min(), mm = img2.max() - min;
    
    cimg_forXY(img2, x, y) {
      int value = 255.0 * (img2(x,y) - min) / mm;
      rgbImg(x, y, 0 /*red*/) = value;
      rgbImg(x, y, 1 /*green*/) = value;
      rgbImg(x, y, 2 /*blue*/) = value;
    }
    
    // Draw circles and text
    ostringstream oss;
    oss.precision(4);
    oss << "HFD" << endl << "~" << hfd;
    
    const unsigned char red[3] = { 255, 0, 0 }, green[3] = { 0, 255, 0 }, yellow[3] = { 255, 255, 0 }, black[3] = { 0, 0, 0 };
    rgbImg.draw_circle(img2.width() / 2, img2.height() / 2, outerDiameter / 2, red, 1 /*pattern*/, 1 /*opacity*/);
    rgbImg.draw_circle(img2.width() / 2, img2.height() / 2, hfd / 2, green, 1 /*pattern*/, 1 /*opacity*/);
    rgbImg.draw_text(0 /*x0*/, 0 /*y0*/, oss.str().c_str(), yellow /*fg color*/, black /*bg color*/, 0.7 /*opacity*/, 14 /*font-size*/);
 
    // Add rgb image to "container" image
    size_t rowIdx = floor((float)(i - 1) / 4.0f);
    size_t colIdx = (i - 1) % 4;
    containerImg.draw_image(imgWidth * colIdx, imgHeight * rowIdx, rgbImg);
  }
 
  // Display the container image
  CImgDisplay imgDisp(containerImg, "Image");
  while (! imgDisp.is_closed()) {
    imgDisp.wait();
  }
  
  return 0;
}
