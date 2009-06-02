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


  Float_t GetNoiseCut(Int_t imod) const { return fNoiseCut[imod];}
  void    SetNoiseCut(Int_t imod, Float_t cut) {fNoiseCut[imod] = cut;}

  static   AliPMDRecoParam *GetPbPbParam();   // reco param for PbPb.
  static   AliPMDRecoParam *GetPPParam();     // reco param for PP
  static   AliPMDRecoParam *GetCosmicParam(); // reco param for cosmic muons
 private:

  Float_t fNoiseCut[48]; //Noise cut

  ClassDef(AliPMDRecoParam, 0)
};

#endif
