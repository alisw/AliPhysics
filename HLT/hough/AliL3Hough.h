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
class TFile;

class AliL3Hough : public TObject {
  
 private:
  Char_t fPath[256];
  Bool_t fBinary;
  Int_t fNEtaSegments;
  AliL3FileHandler **fMemHandler; //!
  AliL3HoughTransformer **fHoughTransformer; //!
  TFile *fRootFile; //!

  void DeleteTransformers();
  void DeleteMemory();
  void Init();
  
 public:
  
  AliL3Hough(); 
  AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments=100);
  virtual ~AliL3Hough();
  
  void TransformSlice(Int_t slice);
  AliL3Histogram *AddHistograms();
  void Evaluate(AliL3Histogram *hist);

  //Setters
  void SetNEtaSegments(Int_t i) {fNEtaSegments = i;}
  
  //Getters
  AliL3HoughTransformer *GetTransformer(Int_t i) {if(!fHoughTransformer[i]) return 0; return fHoughTransformer[i];}

  ClassDef(AliL3Hough,1)

};

#endif
