// @(#) $Id$

#ifndef ALIL3_HOUGHTRANSFORMERLUT
#define ALIL3_HOUGHTRANSFORMERLUT

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

class AliL3Histogram;

class AliL3HoughTransformerLUT : public AliL3HoughBaseTransformer {
  
 protected:
  AliL3Histogram **fParamSpace; //!
#ifdef do_mc
  TrackIndex **fTrackID; //!
#endif
  Bool_t fDoMC;

  void DeleteHistograms();
  
  Int_t fMinRow;
  Int_t fMaxRow;
  Int_t fNRows;
  Int_t fNEtas;
  Int_t fNPhi0;
  Int_t fSlice;
  Int_t fSector;
  Int_t fSectorRow;
  Int_t fZSign;
  Float_t fZLengthPlusOff;
  Float_t fTimeWidth;
  Float_t fPadPitch;
  Float_t fEtaSlice;

  Float_t *fLUTX; //!
  Float_t *fLUTY; //!
  Float_t *fLUTEta; //!
  Float_t *fLUTEtaReal; //!
  Float_t *fLUTphi0; //!
  Float_t *fLUT2sinphi0; //!   
  Float_t *fLUT2cosphi0; //!
  Float_t *fLUTKappa; //!
  
  Int_t fLastPad;
  Int_t fLastIndex;
  Int_t fAccCharge;

  Float_t CalcRoverZ2(Float_t eta);
  Float_t CalcEta(Float_t roverz2);
  Float_t CalcX(Int_t row);
  Float_t CalcY(Int_t pad, Int_t row);
  Float_t CalcZ(Int_t time);  

  Int_t FindIndex(Float_t rz2, Int_t start=-100);
  void AddCurveToHistogram(Int_t new_eta_index=-1);

 public:

  AliL3HoughTransformerLUT(); 
  AliL3HoughTransformerLUT(Int_t slice,Int_t patch,Int_t n_eta_segments);
  virtual ~AliL3HoughTransformerLUT();
  
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();

  void TransformCircle();
  void TransformCircleC(Int_t *row_range,Int_t every){STDCERR<<"TransformCircleC is not defined for this transformer!"<<STDENDL;}
  void TransformLine(Int_t *rowrange=0,Float_t *phirange=0) {STDCERR<<"TransformLine is not defined for this transformer!"<<STDENDL;}

  Int_t GetEtaIndex(Double_t eta);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  Double_t GetEta(Int_t eta_index,Int_t slice);
  Int_t GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi);
  
  void Print();
  void Init(Int_t slice=0,Int_t patch=0,Int_t n_eta_segments=100,Int_t n_seqs=-1);

  ClassDef(AliL3HoughTransformerLUT,1) //Normal Hough transformation class

};
#endif
