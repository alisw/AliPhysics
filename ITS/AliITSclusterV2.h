#ifndef ALIITSCLUSTERV2_H
#define ALIITSCLUSTERV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                        ITS Cluster Class
//
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include "TMath.h"
#include "AliCluster.h"
#include "AliITSrecoV2.h"

//_____________________________________________________________________________
class AliITSclusterV2 : public AliCluster {
public:
  AliITSclusterV2() : AliCluster() {fQ=0; fLayer=0; fNz=1; fNy=1; fType=0;fDeltaProb=0;}
  AliITSclusterV2(Int_t *lab,Float_t *hit, Int_t *info) : AliCluster(lab,hit) {
    fIndex=lab[3];
    fQ=hit[4];
    fNy    = info[0];
    fNz    = info[1];
    fLayer = info[2];
    fChargeRatio = 0;
    fType=0;
    fDeltaProb=0.;
  }
  void Use() {fQ=-fQ;}
  void SetQ(Float_t q) {fQ=q;}
  void SetDetectorIndex(Int_t i) { fIndex=i; }
  void SetLayer(Int_t layer) {fLayer=layer;}
  void SetNz(Int_t nz) {fNz =nz;}
  void SetNy(Int_t ny){fNy=ny;}
  void SetChargeRatio(Float_t ratio) { fChargeRatio = ratio;}
  void SetType(Int_t type){ fType=type;}
  void SetDeltaProbability(Float_t prob){fDeltaProb = prob;}

  Int_t IsUsed() const {return (fQ<0) ? 1 : 0;}
  Float_t GetQ() const {return TMath::Abs(fQ);}
  Int_t GetDetectorIndex() const { return 0x3FF&fIndex; }
  Int_t GetLayer() const {return fLayer;}
  Int_t GetNz() const {return fNz;}
  Int_t GetNy() const {return fNy;}
  Float_t GetChargeRatio() const {return fChargeRatio;}
  Int_t GetPindex() const { return 0xFFF00000&fIndex; }  //SSD clusters only
  Int_t GetNindex() const { return 0xFFC00&fIndex; }  //SSD clusters only
  Int_t GetType() const {return fType;}  // type of the cluster
  Float_t GetDeltaProbability() const{return fDeltaProb;} //probability to belong to the delta ray
private:
  Int_t    fIndex;    // detector index
  Float_t  fQ ;       // Q of cluster (in ADC counts)
  Char_t   fLayer;    // layer number
  Short_t   fNz;       //number of digits in Z direction
  Short_t   fNy;       //number of digits in y direction 
  Float_t  fChargeRatio; //charge ratio
  Int_t    fType;         //quality factor of the cluster
  Float_t  fDeltaProb;    // probability to be deleta electron
  ClassDef(AliITSclusterV2,2)  // ITS clusters
};

#endif


