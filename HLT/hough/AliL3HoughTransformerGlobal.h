// @(#) $Id$

#ifndef ALIL3_HOUGHTRANSFORMERGLOBAL
#define ALIL3_HOUGHTRANSFORMERGLOBAL

#include "AliL3RootTypes.h"
#include "AliL3HoughTransformer.h"

class AliL3Histogram;
class AliL3MemHandler;
class AliL3TrackArray;
class AliL3HoughMaxFinder;
class TH2F;

struct Seed {
  Float_t fPhi;
  Float_t fRadius;
  Float_t fEta;
  Int_t fIndex;
};

class AliL3HoughTransformerGlobal : public AliL3HoughTransformer {
  
 private:
  Int_t fEvent;
  Char_t fPath[1024]; 
  Int_t *fPadMin;  //!
  Int_t *fPadMax;  //!
  Int_t fNActiveSlice;
  Float_t fPsi;
  Float_t fPtMin;
  Int_t fSeedPadRow;
  AliL3MemHandler **fMemHandler; //!
  AliL3TrackArray *fTracks; //!
  AliL3HoughMaxFinder *fPeakFinder; //!
  
  void Rotate(Float_t *xyz,Int_t rel_slice);
  void FindPeaks(AliL3Histogram *hist,Float_t eta);
  Int_t LoadClusterSeeds(Seed *seeds);
  Float_t CalculateBorder(Float_t *xyz,Int_t charge);
  
 public:
  AliL3HoughTransformerGlobal(); 
  AliL3HoughTransformerGlobal(Char_t *path,Int_t event=0); 
  virtual ~AliL3HoughTransformerGlobal();
  
  void CreateHistograms(Float_t ptmin,Int_t nxbin,Int_t nybin);
  void VerifyTracks(AliL3TrackArray *tracks,Int_t &index);
  void TransformCircle();
  void TransformCircleC();
  void DefineRegion(Float_t ptmin,Float_t linephi,Int_t seedpadrow);
  void LoadActiveSlices(Bool_t binary=kTRUE);
  void UnloadActiveSlices();
  void DisplayActiveRegion(TH2F *hist,Int_t eta_index=0);
  
  void SetSeedPadRow(Int_t i) {fSeedPadRow=i;}
  
  AliL3TrackArray *GetTracks() {return fTracks;}
  
  ClassDef(AliL3HoughTransformerGlobal,1) 

};

#endif




