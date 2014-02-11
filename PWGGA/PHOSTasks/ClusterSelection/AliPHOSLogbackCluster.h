#ifndef AliAODCaloCluster_H
#define AliAODCaloCluster_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD calorimeter cluster class (for PHOS and EMCAL)
//     Author: Markus Oldenburg, CERN, 
//             Gustavo Conesa, INFN
//-------------------------------------------------------------------------

#include "TObject.h"
#include "TLorentzVector.h"

class AliPHOSLogbackCluster : public TObject {
 public:
  AliPHOSLogbackCluster(AliVCluster* cluster);
  
  virtual ~AliPHOSLogbackCluster();
  
  virtual Double_t E() const {return fE;}
  virtual Double_t CoreE() const {return fCoreE;}
  
  virtual TLorentzVector GetMomentum(TLorentzVector &momentum, Double_t * vertex);

 protected:
  AliPHOSLogbackCluster(const AliPHOSLogbackCluster& clus); 
  AliPHOSLogbackCluster& operator=(const AliPHOSLogbackCluster& clus);

  Double_t fE;
  Double_t fCoreE;
  Double_t fPosition[3];
  
  ClassDef(AliPHOSLogbackCluster,1);
};

#endif
