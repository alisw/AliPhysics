#ifndef ALIL3_HOUGHTEST
#define ALIL3_HOUGHTEST

#include "AliL3RootTypes.h"

struct SimData {
  Int_t npads;
  Int_t pads[10];//maximum 10 pads width
  Int_t minpad;
};

class AliL3Histogram;
class TH2;

class AliL3HoughTest {
  
 private:
  SimData *fData;
  Int_t fCurrentPatch;

 public:
  
  AliL3HoughTest(); 
  virtual ~AliL3HoughTest();
  
  Bool_t GenerateTrackData(Double_t pt,Double_t psi,Int_t sign,Int_t patch,Int_t minhits);
  void FillImage(TH2 *hist);
  void Transform2Circle(AliL3Histogram *hist);
  void Transform2CircleC(AliL3Histogram *hist);
  void Transform2Line(AliL3Histogram *hist,Int_t *rowrange);
  void Transform2LineC(AliL3Histogram *hist,Int_t *rowrange);
  
  ClassDef(AliL3HoughTest,1) //Hough transform base class
};

#endif







