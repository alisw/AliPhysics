#ifndef ALIL3_HOUGH_Track
#define ALIL3_HOUGH_Track

#include "AliL3Track.h"

class AliL3Transform;

class AliL3HoughTrack : public AliL3Track {
  
 private:

  AliL3Transform *fTransform; //!
  Double_t fMinDist;
  Int_t fWeight;
  Int_t fEtaIndex;
  Double_t fEta;

  Double_t fDLine;
  Double_t fPsiLine;
 
  Bool_t fIsHelix;

 public:
  AliL3HoughTrack(); 
  virtual ~AliL3HoughTrack();
  
  virtual void Set(AliL3Track *track);
  virtual Int_t Compare(const AliL3Track *track) const;

  void SetTrackParameters(Double_t kappa,Double_t phi,Int_t weight);  
  void SetLineParameters(Double_t psi,Double_t D,Int_t weight,Int_t *rowrange,Int_t ref_row);
  Int_t GetWeight()  const {return fWeight;}
  Double_t GetPsiLine() const {return fPsiLine;}
  Double_t GetDLine() const {return fDLine;}

  Int_t GetEtaIndex() const {return fEtaIndex;}
  Double_t GetEta() const {return fEta;}
  void GetLineCrossingPoint(Int_t padrow,Double_t *xy);
  
  void SetEta(Double_t f) {fEta = f;}
  void SetWeight(Int_t i,Bool_t update=kFALSE) {if(update) fWeight+= i; else fWeight = i;}
  void SetEtaIndex(Int_t f) {fEtaIndex = f;}
  void SetBestMCid(Int_t f,Double_t min_dist);
  void SetDLine(Double_t f) {fDLine=f;}
  void SetPsiLine(Double_t f) {fPsiLine=f;}

  ClassDef(AliL3HoughTrack,1)

};

#endif
