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
  AliMultiplicity(Int_t ntr,Float_t *t, Float_t *ph, Float_t *df, Int_t *labels,
                  Int_t ns, Float_t *ts, Float_t *ps);
  AliMultiplicity(const AliMultiplicity& m);
  AliMultiplicity& operator=(const AliMultiplicity& m);
  virtual ~AliMultiplicity();
// methods to access tracklet information
  Int_t GetNumberOfTracklets() const {return fNtracks;}
  Double_t GetTheta(Int_t i) const { if(i>=0 && i<fNtracks) {return fTh[i];}
  else {Error("GetTheta","Invalid track number %d",i); return -9999.;}}
  Double_t GetPhi(Int_t i) const { if(i>=0 && i<fNtracks) {return fPhi[i];}
  else {Error("GetPhi","Invalid track number %d",i); return -9999.;}}
  Double_t GetDeltaPhi(Int_t i) const {if(i>=0 && i<fNtracks) {return fDeltPhi[i];}
  else {Error("GetDeltaPhi","Invalid track number %d",i); return -9999.;}}
  Int_t GetLabel(Int_t i) const {if(i>=0 && i<fNtracks) {return fLabels[i];}
  else {Error("GetLabel","Invalid track number %d",i); return -9999;}}
// methods to access single cluster information
  Int_t GetNumberOfSingleClusters() const {return fNsingle;}
  Double_t GetThetaSingle(Int_t i) const { if(i>=0 && i<fNsingle) {return fThsingle[i];}
  else {Error("GetThetaSingle","Invalid cluster number %d",i); return -9999.;}}
  Double_t GetPhiSingle(Int_t i) const { if(i>=0 && i<fNsingle) {return fPhisingle[i];}
  else {Error("GetPhisingle","Invalid cluster number %d",i); return -9999.;}}

  Short_t GetFiredChips(Int_t layer) { return fFiredChips[layer]; }
  void SetFiredChips(Int_t layer, Short_t firedChips) { fFiredChips[layer] = firedChips; }

  protected:
  void Duplicate(const AliMultiplicity &m);  // used by copy ctr.

  Int_t fNtracks;            // Number of tracklets
  Int_t fNsingle;            // Number of clusters on SPD layer 1, not associated
  Int_t *fLabels;            //[fNtracks] array with labels of tracklets
                             // with a tracklet on SPD layer 2
  Double32_t *fTh;           //[fNtracks] array with theta values
  Double32_t *fPhi;          //[fNtracks] array with phi values
  Double32_t *fDeltPhi;      //[fNtracks] array with delta phi values
  Double32_t *fThsingle;     //[fNsingle] array with theta values of L1 clusters
  Double32_t *fPhisingle;    //[fNsingle] array with phi values of L2 clusters
  Short_t fFiredChips[2]; // number of fired chips in the two SPD layers

  ClassDef(AliMultiplicity,6);
};

#endif
