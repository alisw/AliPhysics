#ifndef ALIL3_HOUGH_Eval
#define ALIL3_HOUGH_Eval

#include "AliL3RootTypes.h"

class AliL3HoughTransformer;
class AliL3Transform;
class AliL3HoughTrack;
class AliL3DigitRowData;
class AliL3Histogram;

class AliL3HoughEval : public TObject {
  
 private:

  Int_t fSlice;
  Int_t fPatch;
  Int_t fNrows;
  Int_t fNEtaSegments;
  Double_t fEtaMin;
  Double_t fEtaMax;
  Int_t fNumOfPadsToLook;
  Int_t fNumOfRowsToMiss;
  
  //Flags
  Bool_t fRemoveFoundTracks;
  
  AliL3Transform *fTransform; //!
  AliL3HoughTransformer *fHoughTransformer; //!
  AliL3DigitRowData **fRowPointers; //!
  
 public:
  AliL3HoughEval(); 
  AliL3HoughEval(AliL3HoughTransformer *transform);
  virtual ~AliL3HoughEval();
  
  void GenerateLUT();
  void DisplayEtaSlice(Int_t eta_index,AliL3Histogram *hist);
  Bool_t LookInsideRoad(AliL3HoughTrack *track,Int_t eta_index,Bool_t remove=kFALSE);
  
  //Setters:
  void RemoveFoundTracks() {fRemoveFoundTracks = kTRUE;}
  
  ClassDef(AliL3HoughEval,1)

};

#endif
