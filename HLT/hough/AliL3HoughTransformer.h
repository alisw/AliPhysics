// @(#) $Id$

#ifndef ALIL3_HOUGHTRANSFORMER
#define ALIL3_HOUGHTRANSFORMER

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

class AliL3Histogram;

class AliL3HoughTransformer : public AliL3HoughBaseTransformer {
  
 private:
  
  AliL3Histogram **fParamSpace; //!
#ifdef do_mc
  TrackIndex **fTrackID; //!
#endif
  Bool_t fDoMC;
  Bool_t fEtaOverlap;
  
  void DeleteHistograms();

 public:
  AliL3HoughTransformer(); 
  AliL3HoughTransformer(Int_t slice,Int_t patch,Int_t n_eta_segments,Bool_t DoEtaOverlap=kFALSE,Bool_t DoMC=kFALSE);
  virtual ~AliL3HoughTransformer();
  
  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t ptres,Int_t nybin,Float_t psi);
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();
  void TransformCircle();
  void TransformCircleC(Int_t *row_range,Int_t every=1);
  void TransformLine(Int_t *rowrange=0,Float_t *phirange=0);
  void TransformLineC(Int_t *rowrange,Float_t *phirange);

  Int_t GetEtaIndex(Double_t eta);
  void GetEtaIndexes(Double_t eta,Int_t *indexes);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  Double_t GetEta(Int_t eta_index,Int_t slice);
  Int_t GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi);
  
  ClassDef(AliL3HoughTransformer,1) //Normal Hough transformation class

};

#endif




