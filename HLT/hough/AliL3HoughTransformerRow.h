// @(#) $Id$

#ifndef ALIL3_HOUGHTRANSFORMERROW
#define ALIL3_HOUGHTRANSFORMERROW

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

class AliL3Histogram;

class AliL3HoughTransformerRow : public AliL3HoughBaseTransformer {
  
 private:
  
  AliL3Histogram **fParamSpace; //!

#ifdef do_mc
  static TrackIndex **fTrackID; //!
#endif
  Bool_t fDoMC;

  void DeleteHistograms();

  static UChar_t **fRowCount; //!
  static UChar_t **fGapCount; //!
  static UChar_t **fCurrentRowCount; //!

  Float_t *fLUT2sinphi0up; //!   
  Float_t *fLUT2cosphi0up; //!
  Float_t *fLUT2sinphi0low; //!   
  Float_t *fLUT2cosphi0low; //!

  Float_t *fLUTforwardZ; //!
  Float_t *fLUTforwardZ2; //!
  Float_t *fLUTbackwardZ; //!
  Float_t *fLUTbackwardZ2; //!

 public:
  AliL3HoughTransformerRow(); 
  AliL3HoughTransformerRow(Int_t slice,Int_t patch,Int_t n_eta_segments,Bool_t DoMC=kFALSE,Float_t zvertex=0.0);
  virtual ~AliL3HoughTransformerRow();

  //void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t ptres,Int_t nybin,Float_t psi);
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();
  void TransformCircle();

  Int_t GetEtaIndex(Double_t eta);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  Double_t GetEta(Int_t eta_index,Int_t slice);
  Int_t GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi);
  UChar_t *GetRowCount(Int_t eta_index);
  UChar_t *GetGapCount(Int_t eta_index);

  ClassDef(AliL3HoughTransformerRow,1) //TPC Rows Hough transformation class

};

#endif




