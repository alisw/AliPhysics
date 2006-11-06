#ifndef ALITPCCLUSTER_H
#define ALITPCCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                    TPC Cluster Class
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include "TMath.h"
#include "AliCluster.h"

//_____________________________________________________________________________
class AliTPCcluster : public AliCluster {
public:
  AliTPCcluster():AliCluster(),fQ(0){}
  AliTPCcluster(Int_t *lab, Float_t *hit);
  void Use(Int_t = 0) {fQ=-fQ;}
  void SetQ(Float_t q) {fQ=q;}

  Int_t IsUsed() const {return (fQ<0) ? 1 : 0;}
  Float_t GetQ() const {return TMath::Abs(fQ);}

  static Double_t SigmaY2(Double_t r, Double_t tgl, Double_t pt);
  static Double_t SigmaZ2(Double_t r, Double_t tgl);

private:
  Float_t   fQ ;       //Q of cluster (in ADC counts)
  
  ClassDef(AliTPCcluster,1)  // Time Projection Chamber clusters
};

#endif


