// @(#) $Id$

#ifndef ALIL3_HOUGHTRANSFORMERROW
#define ALIL3_HOUGHTRANSFORMERROW

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

#define MAX_N_GAPS 4
#define MAX_GAP_SIZE 4
#define MIN_TRACK_LENGTH 70
#define MAX_MISS_ROWS 2

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

  static UChar_t *fTrackNRows; //!
  static UChar_t *fTrackFirstRow; //!
  static UChar_t *fTrackLastRow; //!

  Float_t *fLUTforwardZ; //!
  Float_t *fLUTforwardZ2; //!
  Float_t *fLUTbackwardZ; //!
  Float_t *fLUTbackwardZ2; //!

  static Float_t fBeta1,fBeta2;

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
  UChar_t *GetCurrentRowCount(Int_t eta_index);
  static UChar_t *GetTrackNRows();
  static UChar_t *GetTrackFirstRow();
  static UChar_t *GetTrackLastRow();
  static Float_t GetBeta1() {return fBeta1;}
  static Float_t GetBeta2() {return fBeta2;}

  ClassDef(AliL3HoughTransformerRow,1) //TPC Rows Hough transformation class

};

#endif




