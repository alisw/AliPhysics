#ifndef ALIMEVSIMCONFIG_H
#define ALIMEVSIMCONFIG_H

#include "TObject.h"

class AliGenMevSim;

class AliMevSimConfig : public TObject {

  friend class AliGenMevSim;

 protected:

  static const Int_t kMAX_MODEL = 4;
  static const Int_t kMAX_CTRL = 4;

  Int_t fModelType;

  Int_t fReacPlaneCntrl;
  Float_t fPsiRMean, fPsiRStDev;

  Float_t fMultFacMean, fMultFacStDev;

  Float_t fNStDevMult, fNStDevTemp, fNStDevSigma, fNStDevExpVel, fNStdDevPSIr, fNStDevVn, fNStDevMultFac;

  Int_t fNIntegPts;
  Int_t fNScanPts;

  void Init();

 public:

  AliMevSimConfig();
  AliMevSimConfig(Int_t modelType);

  ~AliMevSimConfig();

  void SetModelType(Int_t modelType);

  void SetRectPlane(Int_t ctrl, Float_t psiRMean = 0, Float_t psiRStDev = 0);
  void SetMultFac(Float_t mean, Float_t stDev);

  void SetStDev(Float_t mult, Float_t temp, Float_t sigma,
		Float_t expVel, Float_t psiR, Float_t Vn, Float_t multFac);

  void SetGrid(Int_t integr, Int_t scan);


  ClassDef(AliMevSimConfig,1)

};


#endif
