#ifndef ALIPMDRECOPARAM_H
#define ALIPMDRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with PMD reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliDetectorRecoParam.h"

class AliPMDRecoParam : public AliDetectorRecoParam
{
 public: 
  AliPMDRecoParam();
  AliPMDRecoParam(const AliPMDRecoParam &source); //Copy Ctor 
  AliPMDRecoParam& operator=(const AliPMDRecoParam &source); // ass. op.
  virtual ~AliPMDRecoParam();

  virtual void PrintParameters() const;

  Int_t GetClusteringParam() const { return fPmdClusteringParam;}
  void  SetClusteringParam(Int_t cluspar) {fPmdClusteringParam = cluspar;}

  static   AliPMDRecoParam *GetPbPbParam();   // reco param for PbPb.
  static   AliPMDRecoParam *GetPPParam();     // reco param for PP
  static   AliPMDRecoParam *GetCosmicParam(); // reco param for cosmic muons

  static AliPMDRecoParam *GetLowFluxParam();// make reco parameters for low flux env.
  static AliPMDRecoParam *GetHighFluxParam();// make reco parameters for high flux env. 
  
  

 private:

  Int_t fPmdClusteringParam;  // Clustering switch to decide crude or refined

  ClassDef(AliPMDRecoParam, 2)
};

#endif
