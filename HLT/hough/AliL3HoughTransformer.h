// @(#) $Id$

#ifndef ALIL3HOUGHTRANSFORMER_H
#define ALIL3HOUGHTRANSFORMER_H

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

class AliL3Histogram;

class AliL3HoughTransformer : public AliL3HoughBaseTransformer {

 public:
  AliL3HoughTransformer(); 
  AliL3HoughTransformer(Int_t slice,Int_t patch,Int_t netasegments,Bool_t DoEtaOverlap=kFALSE,Bool_t DoMC=kFALSE);
  virtual ~AliL3HoughTransformer();
  
  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t ptres,Int_t nybin,Float_t psi);
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();
  void TransformCircle();
  void TransformCircleC(Int_t *rowrange,Int_t every=1);
  void TransformLine(Int_t *rowrange=0,Float_t *phirange=0);
  void TransformLineC(Int_t *rowrange,Float_t *phirange);

  Int_t GetEtaIndex(Double_t eta) const;
  void GetEtaIndexes(Double_t eta,Int_t *indexes) const;
  AliL3Histogram *GetHistogram(Int_t etaindex);
  Double_t GetEta(Int_t etaindex,Int_t slice) const;
  Int_t GetTrackID(Int_t etaindex,Double_t kappa,Double_t psi);
  
 private:
  
  AliL3Histogram **fParamSpace; //!
#ifdef do_mc
  AliL3TrackIndex **fTrackID; //!
#endif
  Bool_t fDoMC; // Calculate mc labels or not
  Bool_t fEtaOverlap; // Allow overlapping of eta slice or not
  
  void DeleteHistograms();

  ClassDef(AliL3HoughTransformer,1) //Normal Hough transformation class

};

#endif




