// @(#) $Id$

#ifndef ALIL3_HOUGHTRANSFORMERGAP
#define ALIL3_HOUGHTRANSFORMERGAP

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

class AliL3Histogram;

#define COUNTINGROWS
//#define SUBSLICES

class AliL3HoughTransformerGap : public AliL3HoughBaseTransformer {
  
 private:
  
  AliL3Histogram **fParamSpace; //!

#ifdef do_mc
  TrackIndex **fTrackID; //!
#endif
  Bool_t fDoMC;

  Float_t fZVertex;

  void DeleteHistograms();

#ifdef SUBSLICES
  Int_t **fEtaSpace; //!
  Int_t fEtaSpaceSize;
#endif
#ifdef COUNTINGROWS
  static UChar_t **fRowCount; //!
  static UChar_t **fGapCount; //!
  static UChar_t **fCurrentRowCount; //!

  Float_t *fLUT2sinphi0up; //!   
  Float_t *fLUT2cosphi0up; //!
  Float_t *fLUT2sinphi0low; //!   
  Float_t *fLUT2cosphi0low; //!
#endif

 public:
  AliL3HoughTransformerGap(); 
  AliL3HoughTransformerGap(Int_t slice,Int_t patch,Int_t n_eta_segments,Bool_t DoEtaOverlap=kFALSE,Bool_t DoMC=kFALSE,Float_t zvertex=0.0);
  virtual ~AliL3HoughTransformerGap();

  //void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t ptres,Int_t nybin,Float_t psi);
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();
  void TransformCircle();

  Int_t GetEtaIndex(Double_t eta);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  Double_t GetEta(Int_t eta_index,Int_t slice);

#ifdef SUBSLICES
  Int_t GetEtaSubIndex(Double_t eta,Int_t eta_index);
  Double_t GetMaxSubEta(Float_t kappa,Float_t phi0,Int_t eta_index,Int_t slice);
#endif

  //Int_t GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi);

  UChar_t *GetRowCount(Int_t eta_index);
  UChar_t *GetGapCount(Int_t eta_index);

#ifdef COUNTINGROWS
  //void CalculateNRows();
#endif

  ClassDef(AliL3HoughTransformerGap,1) //Hough transformation class

};

#endif




