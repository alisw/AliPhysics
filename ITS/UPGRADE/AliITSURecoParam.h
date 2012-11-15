#ifndef ALIITSURECOPARAM_H
#define ALIITSURECOPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliITSURecoParam.h 57215 2012-06-17 14:47:08Z masera $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with ITS reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliDetectorRecoParam.h"

class AliITSURecoParam : public AliDetectorRecoParam
{
 public: 
  AliITSURecoParam();
  AliITSURecoParam(Int_t nLr);
  virtual ~AliITSURecoParam();

  static AliITSURecoParam *GetLowFluxParam();// make reco parameters for low flux env.
  static AliITSURecoParam *GetHighFluxParam();// make reco parameters for high flux env. 
  static AliITSURecoParam *GetCosmicTestParam();// special setting for cosmic  

  Double_t    GetTanLorentzAngle(Int_t lr)                 const;
  //
  void        SetNLayers(Int_t n);
  void        SetTanLorentzAngle(Int_t lr, Double_t v=0.);


 protected:
  Int_t          fNLayers;          // number of layers 
  Double_t*      fTanLorentzAngle;  //[fNLayers] Optional Lorentz angle for each layer
  //
 private:
  AliITSURecoParam(const AliITSURecoParam & param);
  AliITSURecoParam & operator=(const AliITSURecoParam &param);

  ClassDef(AliITSURecoParam,1) // ITS reco parameters
};

//_____________________________________________________________________________
inline Double_t AliITSURecoParam::GetTanLorentzAngle(Int_t lr) const 
{
  // get tg of Lorentz Angle for the layer
  return lr<fNLayers ? fTanLorentzAngle[lr]:0;
}

#endif


