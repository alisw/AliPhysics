#ifndef ALIMEVSIMCONFIG_H
#define ALIMEVSIMCONFIG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//__________________________________________________________________
////////////////////////////////////////////////////////////////////
//
// class AliMevSimConfig
//
// Class containing configuation inforamtion for MeVSim generator
// --------------------------------------------
// --------------------------------------------
// --------------------------------------------
// author: radomski@if.pw.edu.pl
//
////////////////////////////////////////////////////////////////////

#include "TObject.h"

class AliGenMevSim;

class AliMevSimConfig : public TObject {

 protected:

  static const Int_t fgkMAX_MODEL = 4; //Maximum number of available models
  static const Int_t fgkMAX_CTRL = 4;//Maximum number of available controls

  Int_t fModelType;  //current type of model

  Int_t fReacPlaneCntrl; //reaction plane simulation model
  Float_t fPsiRMean; //fPsiRMean mean psi
  Float_t fPsiRStDev; //fPsiRStDev  psi variance

  Float_t fMultFacMean;//fMultFacMean Mean multiplicity
  Float_t fMultFacStDev;//fMultFacStDev multiplicity variance

  Float_t fNStDevMult;//fNStDevMult
  Float_t fNStDevTemp;//fNStDevTemp
  Float_t fNStDevSigma;//fNStDevSigma
  Float_t fNStDevExpVel;//fNStDevExpVel
  Float_t fNStdDevPSIr;//fNStdDevPSIr
  Float_t fNStDevVn;//fNStDevVn
  Float_t fNStDevMultFac;//fNStDevMultFac

  Int_t fNIntegPts;//fNIntegPts
  Int_t fNScanPts;//fNScanPts

  void Init();

 public:

  AliMevSimConfig();
  AliMevSimConfig(Int_t modelType);

  ~AliMevSimConfig();

  void SetModelType(Int_t modelType);
  Int_t  GetModelType() const {return fModelType;}

  void SetRectPlane(Int_t ctrl, Float_t psiRMean = 0, Float_t psiRStDev = 0);
  void GetRectPlane(Int_t& ctrl, Float_t& psiRMean, Float_t& psiRStDev ) const
   {ctrl  = fReacPlaneCntrl; psiRMean = fPsiRMean; psiRStDev = fPsiRStDev;}
  
  void SetMultFac(Float_t mean, Float_t stDev);
  void GetMultFac(Float_t& mean, Float_t& stDev) const {mean = fMultFacMean ;stDev = fMultFacStDev;}

  void SetStDev(Float_t mult, Float_t temp, Float_t sigma,
                Float_t expVel, Float_t psiR, Float_t Vn, Float_t multFac);
  void GetStDev(Float_t& mult, Float_t& temp, Float_t& sigma,
                Float_t& expVel, Float_t& psiR, Float_t& Vn, Float_t& multFac) const;
  void SetGrid(Int_t integr, Int_t scan);
  void GetGrid(Int_t& integr, Int_t& scan) const {scan=fNScanPts;integr=fNIntegPts;}

  ClassDef(AliMevSimConfig,1)

};


#endif
