#ifndef ALIL3_HOUGH_Track
#define ALIL3_HOUGH_Track

#include "AliL3Track.h"

class AliL3Transform;

class AliL3HoughTrack : public AliL3Track {
  
 private:

  AliL3Transform *fTransform; //!
  Double_t fMinDist;
  Int_t fWeight;

  Double_t fDLine;
  Double_t fPsiLine;

  Bool_t fIsHelix;

 public:
  AliL3HoughTrack(); 
  virtual ~AliL3HoughTrack();
  
  void Set(AliL3Track *track);
  
  void SetTrackParameters(Double_t kappa,Double_t phi,Int_t weight);  
  void SetLineParameters(Double_t psi,Double_t D,Int_t weight,Int_t *rowrange,Int_t ref_row);
  Int_t GetWeight()  {return fWeight;}
  Double_t GetPsiLine() {return fPsiLine;}
  Double_t GetDLine() {return fDLine;}

  void GetLineCrossingPoint(Int_t padrow,Double_t *xy);
  Double_t GetCrossingAngle(Int_t padrow);
  Bool_t GetCrossingPoint(Int_t slice,Int_t padrow,Float_t *xyz);
  Bool_t GetCrossingPoint(Int_t padrow,Float_t *xyz);
  
  void SetBestMCid(Int_t f,Double_t min_dist);
  void SetDLine(Double_t f) {fDLine=f;}
  void SetPsiLine(Double_t f) {fPsiLine=f;}

  ClassDef(AliL3HoughTrack,1)

};

#endif
