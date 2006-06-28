#ifndef ALIMULTIPLICITY_H
#define ALIMULTIPLICITY_H

#include<TObject.h>

////////////////////////////////////////////////////////
////   Class containing multiplicity information      //
////   to stored in the ESD                           //
////////////////////////////////////////////////////////

class AliMultiplicity : public TObject {

 public:

  AliMultiplicity();               // default constructor
  // standard constructor
  AliMultiplicity(Int_t ntr,Float_t *t=NULL, Float_t *ph=NULL, Float_t *df=NULL); 
  AliMultiplicity(const AliMultiplicity& m);
  AliMultiplicity& operator=(const AliMultiplicity& m);
  virtual ~AliMultiplicity();
  Int_t GetNumberOfTracklets() const {return fNtracks;}
  Float_t GetTheta(Int_t i) const { if(i>=0 && i<fNtracks) {return fTh[i];}
  else {Error("GetTheta","Invalid track number %d",i); return -9999.;}}
  Float_t GetPhi(Int_t i) const { if(i>=0 && i<fNtracks) {return fPhi[i];}
  else {Error("GetTheta","Invalid track number %d",i); return -9999.;}}
  Float_t GetDeltaPhi(Int_t i) const {if(i>=0 && i<fNtracks) {return fDeltPhi[i];}
  else {Error("GetDeltaPhi","Invalid track number %d",i); return -9999.;}}

  protected:
  void Duplicate(const AliMultiplicity &m);  // used by copy ctr.
  Int_t fNtracks;         // Number of tracklets (=-1 when mult is not determined)
  Float_t *fTh;           //[fNtracks] array with theta values
  Float_t *fPhi;          //[fNtracks] array with phi values
  Float_t *fDeltPhi;      //[fNtracks] array with delta phi values

  ClassDef(AliMultiplicity,1);
};

#endif
