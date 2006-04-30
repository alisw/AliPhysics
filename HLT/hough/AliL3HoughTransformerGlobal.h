// @(#) $Id$

#ifndef ALIL3HOUGHTRANSFORMERGLOBAL_H
#define ALIL3HOUGHTRANSFORMERGLOBAL_H

#include "AliL3RootTypes.h"
#include "AliL3HoughTransformer.h"

class AliL3Histogram;
class AliL3MemHandler;
class AliL3TrackArray;
class AliL3HoughMaxFinder;
class TH2F;

struct AliL3Seed {
  Float_t fPhi;//Track emission angle
  Float_t fRadius;//Track radius
  Float_t fEta;//Track eta
  Int_t fIndex;//index?
};

class AliL3HoughTransformerGlobal : public AliL3HoughTransformer {

 public:
  AliL3HoughTransformerGlobal(); 
  AliL3HoughTransformerGlobal(Char_t *path,Int_t event=0); 
  virtual ~AliL3HoughTransformerGlobal();
  
  void CreateHistograms(Float_t ptmin,Int_t nxbin,Int_t nybin);
  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t pres,Int_t nybin,Float_t psi) {
    AliL3HoughTransformer::CreateHistograms(ptmin,ptmax,pres,nybin,psi);
  }
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax) {
    AliL3HoughTransformer::CreateHistograms(nxbin,ptmin,nybin,phimin,phimax);
  }
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax) {
    AliL3HoughTransformer::CreateHistograms(nxbin,xmin,xmax,nybin,ymin,ymax);
  }
  void VerifyTracks(AliL3TrackArray *tracks,Int_t &index);
  void TransformCircle();
  void TransformCircle(Int_t *row_range,Int_t every) {
    AliL3HoughTransformer::TransformCircle(row_range,every);
  }
  void TransformCircleC();
  void TransformCircleC(Int_t *rowrange,Int_t every){
    AliL3HoughTransformer::TransformCircleC(rowrange,every);
  }
  void DefineRegion(Float_t ptmin,Float_t linephi,Int_t seedpadrow);
  void LoadActiveSlices(Bool_t binary=kTRUE);
  void UnloadActiveSlices();
  void DisplayActiveRegion(TH2F *hist,Int_t etaindex=0);
  
  void SetSeedPadRow(Int_t i) {fSeedPadRow=i;}
  
  AliL3TrackArray *GetTracks() {return fTracks;}
  
 private:
  Int_t fEvent;//Event index
  Char_t fPath[1024];//path to the input files
  Int_t *fPadMin;  //!
  Int_t *fPadMax;  //!
  Int_t fNActiveSlice;//number of active slices
  Float_t fPsi;//psi range
  Float_t fPtMin;//minimum pt
  Int_t fSeedPadRow;//?
  AliL3MemHandler **fMemHandler; //!
  AliL3TrackArray *fTracks; //!
  AliL3HoughMaxFinder *fPeakFinder; //!
  
  void Rotate(Float_t *xyz,Int_t relslice);
  void FindPeaks(AliL3Histogram *hist,Float_t eta);
  Int_t LoadClusterSeeds(AliL3Seed *seeds);
  Float_t CalculateBorder(Float_t *xyz,Int_t charge);
  
  ClassDef(AliL3HoughTransformerGlobal,1) 

};

#endif




