#ifndef ALIL3_HOUGH
#define ALIL3_HOUGH

#include "AliL3RootTypes.h"

class AliL3HoughMaxFinder;
class AliL3HoughTransformer;
class AliL3Histogram;
class AliL3FileHandler;
class AliL3HoughEval;
class AliL3Transform;
class AliL3TrackArray;

class AliL3Hough : public TObject {
  
 private:

  Char_t fPath[256];
  Int_t fNEtaSegments;
  AliL3Histogram **fHistos; //!
  AliL3FileHandler *fMemHandler; //!
  AliL3HoughMaxFinder *fMaxFinder; 
  AliL3HoughEval *fEval;
  AliL3HoughTransformer *fHoughTransformer;
  AliL3Transform *fTransform; //!
  Bool_t fUseBinary;
  Bool_t fDeleteTrack;
  AliL3TrackArray *fTracks; //!

  Int_t fNxbin;
  Int_t fNybin;
  Double_t fXmin;
  Double_t fXmax;
  Double_t fYmin;
  Double_t fYmax;

 public:

  AliL3Hough(); 
  AliL3Hough(Int_t n_eta_segments,Int_t xbin,Double_t *xrange,Int_t ybin,Double_t *yrange);
  virtual ~AliL3Hough();
  
  void SetInput(Char_t *input,Bool_t binary);
  void ProcessSlice(Int_t slice);
  void ProcessPatch(Int_t slice,Int_t patch);
  void SetDeleteTrack(Bool_t f) {fDeleteTrack = (Bool_t)f;}
  
  AliL3TrackArray *GetTracks() {return fTracks;}

  ClassDef(AliL3Hough,1)

};

#endif
