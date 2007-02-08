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
  AliMultiplicity(Int_t ntr,Float_t *t, Float_t *ph, Float_t *df,
                  Int_t ns, Float_t *ts, Float_t *ps); 
  AliMultiplicity(const AliMultiplicity& m);
  AliMultiplicity& operator=(const AliMultiplicity& m);
  virtual ~AliMultiplicity();
// methods to access tracklet information
  Int_t GetNumberOfTracklets() const {return fNtracks;}
  Float_t GetTheta(Int_t i) const { if(i>=0 && i<fNtracks) {return fTh[i];}
  else {Error("GetTheta","Invalid track number %d",i); return -9999.;}}
  Float_t GetPhi(Int_t i) const { if(i>=0 && i<fNtracks) {return fPhi[i];}
  else {Error("GetPhi","Invalid track number %d",i); return -9999.;}}
  Float_t GetDeltaPhi(Int_t i) const {if(i>=0 && i<fNtracks) {return fDeltPhi[i];}
  else {Error("GetDeltaPhi","Invalid track number %d",i); return -9999.;}}
// methods to access single cluster information
  Int_t GetNumberOfSingleClusters() const {return fNsingle;}
  Float_t GetThetaSingle(Int_t i) const { if(i>=0 && i<fNsingle) {return fThsingle[i];}
  else {Error("GetThetaSingle","Invalid cluster number %d",i); return -9999.;}}
  Float_t GetPhiSingle(Int_t i) const { if(i>=0 && i<fNsingle) {return fPhisingle[i];}
  else {Error("GetPhisingle","Invalid cluster number %d",i); return -9999.;}}

  protected:
  void Duplicate(const AliMultiplicity &m);  // used by copy ctr.
  Int_t fNtracks;         // Number of tracklets 
  Float_t *fTh;           //[fNtracks] array with theta values
  Float_t *fPhi;          //[fNtracks] array with phi values
  Float_t *fDeltPhi;      //[fNtracks] array with delta phi values
  Int_t fNsingle;         // Number of clusters on SPD layer 1, not associated
                          // with a tracklet on SPD layer 2
  Float_t *fThsingle;     //[fNsingle] array with theta values of L1 clusters
  Float_t *fPhisingle;    //[fNsingle] array with phi values of L2 clusters

  ClassDef(AliMultiplicity,2);
};

#endif
