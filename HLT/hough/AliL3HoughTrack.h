// @(#) $Id$

#ifndef ALIL3_HOUGH_Track
#define ALIL3_HOUGH_Track

#include "AliL3Track.h"

class AliL3HoughTrack : public AliL3Track {
  
 private:
  
  Double_t fMinDist;
  Int_t fWeight;
  Int_t fEtaIndex;
  Double_t fEta;
  Int_t fSlice; //The slice where this track was found

  Double_t fDLine;
  Double_t fPsiLine;
 
  Bool_t fIsHelix;

  Float_t fBinX,fBinY,fSizeX,fSizeY;

 public:
  AliL3HoughTrack(); 
  virtual ~AliL3HoughTrack();
  
  virtual void Set(AliL3Track *track);
  virtual Int_t Compare(const AliL3Track *track) const;
  
  Bool_t IsHelix() {return fIsHelix;}
  void UpdateToFirstRow();
  void SetTrackParameters(Double_t kappa,Double_t eangle,Int_t weight);  
  void SetLineParameters(Double_t psi,Double_t D,Int_t weight,Int_t *rowrange,Int_t ref_row);

  Int_t GetWeight()  const {return fWeight;}
  Double_t GetPsiLine() const {return fPsiLine;}
  Double_t GetDLine() const {return fDLine;}

  Int_t GetEtaIndex() const {return fEtaIndex;}
  Double_t GetEta() const {return fEta;}
  Int_t GetSlice()  const {return fSlice;}
  void GetLineCrossingPoint(Int_t padrow,Float_t *xy);

  Float_t    GetBinX()   const {return fBinX;}
  Float_t    GetBinY()   const {return fBinY;}
  Float_t    GetSizeX()  const {return fSizeX;}
  Float_t    GetSizeY()  const {return fSizeY;}
  
  void SetHelixTrue() {fIsHelix=kTRUE;}
  void SetSlice(Int_t slice) {fSlice=slice;}
  void SetEta(Double_t f);
  void SetWeight(Int_t i,Bool_t update=kFALSE) {if(update) fWeight+= i; else fWeight = i;}
  void SetEtaIndex(Int_t f) {fEtaIndex = f;}
  void SetBestMCid(Int_t f,Double_t min_dist);
  void SetDLine(Double_t f) {fDLine=f;}
  void SetPsiLine(Double_t f) {fPsiLine=f;}

  void SetBinXY(Float_t binx,Float_t biny,Float_t sizex,Float_t sizey) {fBinX = binx; fBinY = biny; fSizeX = sizex; fSizeY = sizey;}

  ClassDef(AliL3HoughTrack,1) //Track class for Hough tracklets

};

#endif
