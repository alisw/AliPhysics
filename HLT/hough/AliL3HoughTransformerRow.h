// @(#) $Id$

#ifndef ALIL3HOUGHTRANSFORMERROW_H
#define ALIL3HOUGHTRANSFORMERROW_H

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

#define MAX_N_GAPS 4
#define MAX_GAP_SIZE 4
#define MIN_TRACK_LENGTH 70
#define MAX_MISS_ROWS 2

class AliL3Histogram;

class AliL3HoughTransformerRow : public AliL3HoughBaseTransformer {

 public:
  AliL3HoughTransformerRow(); 
  AliL3HoughTransformerRow(Int_t slice,Int_t patch,Int_t netasegments,Bool_t DoMC=kFALSE,Float_t zvertex=0.0);
  virtual ~AliL3HoughTransformerRow();

  //void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t ptres,Int_t nybin,Float_t psi);
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();
  void TransformCircle();

  Int_t GetEtaIndex(Double_t eta) const;
  AliL3Histogram *GetHistogram(Int_t etaindex);
  Double_t GetEta(Int_t etaindex,Int_t slice) const;
  Int_t GetTrackID(Int_t etaindex,Double_t kappa,Double_t psi);
  UChar_t *GetRowCount(Int_t etaindex);
  UChar_t *GetGapCount(Int_t etaindex);
  UChar_t *GetCurrentRowCount(Int_t etaindex);
  static UChar_t *GetTrackNRows();
  static UChar_t *GetTrackFirstRow();
  static UChar_t *GetTrackLastRow();
  static Float_t GetBeta1() {return fgBeta1;}
  static Float_t GetBeta2() {return fgBeta2;}

 private:
  
  AliL3Histogram **fParamSpace; //!

#ifdef do_mc
  static AliL3TrackIndex **fgTrackID; //!
#endif
  Bool_t fDoMC; // Do MC labels or not

  void DeleteHistograms(); //Method to clean up the histograms containing Hough space

  static UChar_t **fgRowCount; //!
  static UChar_t **fgGapCount; //!
  static UChar_t **fgCurrentRowCount; //!

  static UChar_t *fgTrackNRows; //!
  static UChar_t *fgTrackFirstRow; //!
  static UChar_t *fgTrackLastRow; //!

  Float_t *fLUTforwardZ; //!
  Float_t *fLUTforwardZ2; //!
  Float_t *fLUTbackwardZ; //!
  Float_t *fLUTbackwardZ2; //!

  static Float_t fgBeta1,fgBeta2; // Two curves which define the Hough space

  ClassDef(AliL3HoughTransformerRow,1) //TPC Rows Hough transformation class

};

#endif




