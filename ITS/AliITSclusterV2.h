#ifndef ALIITSCLUSTERV2_H
#define ALIITSCLUSTERV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                        ITS Cluster Class
//
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include "AliCluster.h"
#include "AliITSrecoV2.h"

//_____________________________________________________________________________
class AliITSclusterV2 : public AliCluster {
public:
  AliITSclusterV2() : AliCluster() {fQ=0;}
  AliITSclusterV2(Int_t *lab,Float_t *hit) : AliCluster(lab,hit) {
    fIndex=lab[3];
    fQ=hit[4];
  }
  void Use() {fSigmaY2=-fSigmaY2;}
  void SetQ(Float_t q) {fQ=q;}
  void SetDetectorIndex(Int_t i) { fIndex=i; }

  Int_t IsUsed() const {return (fSigmaY2<0) ? 1 : 0;}
  Float_t GetQ() const {return fQ;}
  Int_t GetDetectorIndex() const { return 0x3FF&fIndex; }

  Int_t GetPindex() const { return 0xFFF00000&fIndex; }  //SSD clusters only
  Int_t GetNindex() const { return 0xFFC00&fIndex; }  //SSD clusters only

private:
  Int_t    fIndex;    // detector index
  Float_t  fQ ;       // Q of cluster (in ADC counts)
  
  ClassDef(AliITSclusterV2,1)  // ITS clusters
};

#endif


