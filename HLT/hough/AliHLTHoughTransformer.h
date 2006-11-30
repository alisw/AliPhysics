// @(#) $Id$

#ifndef ALIL3HOUGHTRANSFORMER_H
#define ALIL3HOUGHTRANSFORMER_H

//-------------------------------------------------------------------------
//                Class AliHLTHoughTransformer
//   This is one of the possible implementations of the Hough Transform
//   TPC tracking algorithms for HLT. 
//-------------------------------------------------------------------------

#include "AliHLTRootTypes.h"
#include "AliHLTHoughBaseTransformer.h"

class AliHLTHistogram;

class AliHLTHoughTransformer : public AliHLTHoughBaseTransformer {

 public:
  AliHLTHoughTransformer(); 
  AliHLTHoughTransformer(Int_t slice,Int_t patch,Int_t netasegments,Bool_t DoEtaOverlap=kFALSE,Bool_t DoMC=kFALSE);
  virtual ~AliHLTHoughTransformer();
  
  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t ptres,Int_t nybin,Float_t psi);
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();
  void TransformCircle();
  void TransformCircle(Int_t *rowrange,Int_t every) {
    AliHLTHoughBaseTransformer:: TransformCircle(rowrange,every);
  }
  void TransformCircleC(Int_t *rowrange,Int_t every=1);
  void TransformLine(Int_t *rowrange=0,Float_t *phirange=0);
  void TransformLineC(Int_t *rowrange,Float_t *phirange);

  Int_t GetEtaIndex(Double_t eta) const;
  void GetEtaIndexes(Double_t eta,Int_t *indexes) const;
  AliHLTHistogram *GetHistogram(Int_t etaindex);
  Double_t GetEta(Int_t etaindex,Int_t slice) const;
  Int_t GetTrackID(Int_t etaindex,Double_t kappa,Double_t psi) const;
  
 private:
  
  AliHLTHistogram **fParamSpace; //!
#ifdef do_mc
  AliHLTTrackIndex **fTrackID; //!
#endif
  Bool_t fDoMC; // Calculate mc labels or not
  Bool_t fEtaOverlap; // Allow overlapping of eta slice or not
  
  void DeleteHistograms();

  ClassDef(AliHLTHoughTransformer,1) //Normal Hough transformation class

};

typedef AliHLTHoughTransformer AliL3HoughTransformer; // for backward compatibility

#endif




