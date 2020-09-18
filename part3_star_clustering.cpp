#include <iostream>
#include <CImg.h>
#include <CCfits/CCfits>
 
#include <list>
#include <set>
#include <algorithm>
#include <tuple>
 
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
 
 
typedef tuple<int /*x*/,int /*y*/> PixelPosT;
typedef set<PixelPosT> PixelPosSetT;
typedef list<PixelPosT> PixelPosListT;
 
/**
 * Removes all white neighbours arond pixel from whitePixels
 * if they exist and adds them to pixelsToBeProcessed and to
 * pixelsinCluster.
 */
void
getAndRemoveNeighbours(PixelPosT inCurPixelPos, PixelPosSetT * inoutWhitePixels, 
                       PixelPosListT * inoutPixelsToBeProcessed,
                       PixelPosListT * outPixelCluster)
{
  const size_t _numPixels = 8, _x = 0, _y = 1;
  const int offsets[_numPixels][2] = { { -1, -1 }, { 0, -1 }, { 1, -1 },
                                       { -1, 0 },              { 1, 0 },
                                       { -1, 1 }, { 0, 1 }, { 1, 1 } };
  
  for (size_t p = 0; p < _numPixels; ++p) {
    PixelPosT curPixPos(
                        std::get<0>(inCurPixelPos) + offsets[p][_x],
                        std::get<1>(inCurPixelPos) + offsets[p][_y]
    );
    PixelPosSetT::iterator itPixPos = inoutWhitePixels->find(curPixPos);
 
    if (itPixPos != inoutWhitePixels->end()) {
      const PixelPosT & curPixPos = *itPixPos;
      inoutPixelsToBeProcessed->push_back(curPixPos);
      outPixelCluster->push_back(curPixPos);
      // Remove white pixel from "white set" since it has
      // been now processed
      inoutWhitePixels->erase(itPixPos);
    }
  }
  return;
}
 
template<typename T>
void cluster_stars(const CImg<T> & inImg, vector<PixelPosListT> * outRecognizedClusters) {
  PixelPosSetT whitePixels;
 
  cimg_forXY(inImg, x, y) {
    if (inImg(x, y)) {
      whitePixels.insert(whitePixels.end(), PixelPosT(x, y));
    }
  }
 
  // Iterate over white pixels as long as set is not empty
  while (whitePixels.size()) {
    PixelPosListT pixelCluster;
    PixelPosListT pixelsToBeProcessed;
 
    PixelPosSetT::iterator itWhitePixPos = whitePixels.begin();
    pixelsToBeProcessed.push_back(*itWhitePixPos);
    whitePixels.erase(itWhitePixPos);
 
    while(! pixelsToBeProcessed.empty()) {
      PixelPosT curPixelPos = pixelsToBeProcessed.front();
      getAndRemoveNeighbours(curPixelPos, & whitePixels, & pixelsToBeProcessed, & pixelCluster);
      pixelsToBeProcessed.pop_front();
    }
 
    // Finally, append the cluster
    outRecognizedClusters->push_back(pixelCluster);
  }
}
 
int main(int argc, char *argv[]) {
  CImg<float> binImg;
  vector<PixelPosListT> recognizedPixelClusters;
 
  try {
    readFile(binImg, argv[1]);
    cluster_stars(binImg, & recognizedPixelClusters);
 
    cout << "Recognized " << recognizedPixelClusters.size() << " stars..." << endl;
  } catch(...) {
    cerr << "Problem occured!" << endl;
  }
  return 0;
}
