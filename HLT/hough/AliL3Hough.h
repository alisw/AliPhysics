#ifndef ALIL3_HOUGH
#define ALIL3_HOUGH

#include "AliL3RootTypes.h"

class AliL3HoughMaxFinder;
class AliL3HoughTransformer;
class TH2F;

class AliL3Hough : public TObject {
  
 private:

  TH2F *fParamSpace;  //!
  Char_t fInputFile[100];
  
  AliL3HoughTransformer *fHoughTransformer;
  AliL3HoughMaxFinder *fPeakFinder;


 public:

  AliL3Hough(); 
  AliL3Hough(Char_t *rootfile,TH2F *hist);
  AliL3Hough(Char_t *rootfile,Int_t xbin,Double_t *xrange,Int_t ybin,Double_t *yrange);
  virtual ~AliL3Hough();
  
  void ProcessSlice(Int_t slice);
  void ProcessPatch(Int_t patch);
  void ProcessEtaSlice(Int_t patch,Double_t *eta);

  ClassDef(AliL3Hough,1)

};

#endif
