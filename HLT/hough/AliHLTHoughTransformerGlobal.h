// @(#) $Id$

#ifndef ALIL3HOUGHTRANSFORMERGLOBAL_H
#define ALIL3HOUGHTRANSFORMERGLOBAL_H

#include "AliHLTRootTypes.h"
#include "AliHLTHoughTransformer.h"

class AliHLTHistogram;
class AliHLTMemHandler;
class AliHLTTrackArray;
class AliHLTHoughMaxFinder;
class TH2F;

struct AliHLTSeed {
  Float_t fPhi;//Track emission angle
  Float_t fRadius;//Track radius
  Float_t fEta;//Track eta
  Int_t fIndex;//index?
};

class AliHLTHoughTransformerGlobal : public AliHLTHoughTransformer {

 public:
  AliHLTHoughTransformerGlobal(); 
  AliHLTHoughTransformerGlobal(Char_t *path,Int_t event=0); 
  virtual ~AliHLTHoughTransformerGlobal();
  
  void CreateHistograms(Float_t ptmin,Int_t nxbin,Int_t nybin);
  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t pres,Int_t nybin,Float_t psi) {
    AliHLTHoughTransformer::CreateHistograms(ptmin,ptmax,pres,nybin,psi);
  }
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax) {
    AliHLTHoughTransformer::CreateHistograms(nxbin,ptmin,nybin,phimin,phimax);
  }
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax) {
    AliHLTHoughTransformer::CreateHistograms(nxbin,xmin,xmax,nybin,ymin,ymax);
  }
  void VerifyTracks(AliHLTTrackArray *tracks,Int_t &index);
  void TransformCircle();
  void TransformCircle(Int_t *row_range,Int_t every) {
    AliHLTHoughTransformer::TransformCircle(row_range,every);
  }
  void TransformCircleC();
  void TransformCircleC(Int_t *rowrange,Int_t every){
    AliHLTHoughTransformer::TransformCircleC(rowrange,every);
  }
  void DefineRegion(Float_t ptmin,Float_t linephi,Int_t seedpadrow);
  void LoadActiveSlices(Bool_t binary=kTRUE);
  void UnloadActiveSlices();
  void DisplayActiveRegion(TH2F *hist,Int_t etaindex=0);
  
  void SetSeedPadRow(Int_t i) {fSeedPadRow=i;}
  
  AliHLTTrackArray *GetTracks() {return fTracks;}
  
 private:
  Int_t fEvent;//Event index
  Char_t fPath[1024];//path to the input files
  Int_t *fPadMin;  //!
  Int_t *fPadMax;  //!
  Int_t fNActiveSlice;//number of active slices
  Float_t fPsi;//psi range
  Float_t fPtMin;//minimum pt
  Int_t fSeedPadRow;//?
  AliHLTMemHandler **fMemHandler; //!
  AliHLTTrackArray *fTracks; //!
  AliHLTHoughMaxFinder *fPeakFinder; //!
  
  void Rotate(Float_t *xyz,Int_t relslice);
  void FindPeaks(AliHLTHistogram *hist,Float_t eta);
  Int_t LoadClusterSeeds(AliHLTSeed *seeds);
  Float_t CalculateBorder(Float_t *xyz,Int_t charge);
  
  ClassDef(AliHLTHoughTransformerGlobal,1) 

};

typedef AliHLTHoughTransformerGlobal AliL3HoughTransformerGlobal; // for backward compatibility

#endif




