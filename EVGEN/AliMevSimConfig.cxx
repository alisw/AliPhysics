

#include "AliMevSimConfig.h"

ClassImp(AliMevSimConfig)


//////////////////////////////////////////////////////////////////////////////////////////////////

AliMevSimConfig::AliMevSimConfig() {

  Init();
}

//////////////////////////////////////////////////////////////////////////////////////////////////

AliMevSimConfig::AliMevSimConfig(Int_t modelType) {

  Init();
  SetModelType(modelType);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

AliMevSimConfig::~AliMevSimConfig() {
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::Init() {

  // Default Values

  fModelType = 1;
  fReacPlaneCntrl = 4;
  fPsiRMean = fPsiRStDev = 0;

  fMultFacMean  = 1.0;
  fMultFacStDev = 0.0;

  fNStDevMult = fNStDevTemp = fNStDevSigma = 3.0;
  fNStDevExpVel = fNStdDevPSIr = fNStDevVn = fNStDevMultFac = 3.0;
  
  fNIntegPts = fNScanPts = 100;

}

//////////////////////////////////////////////////////////////////////////////////////////////////
 
void AliMevSimConfig::SetModelType(Int_t modelType) {

  if (modelType < 0 || modelType > kMAX_MODEL)
    Error("SetModelType","Wrog Model Type indentifier (%d)",modelType);

  fModelType = modelType;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::SetRectPlane(Int_t ctrl, Float_t psiRMean, Float_t psiRStDev) {

  if (ctrl < 0 || ctrl > kMAX_CTRL)
    Error("SetReactPlane","Wrong Control Parameter (%d)", ctrl);

  fReacPlaneCntrl = ctrl;
  fPsiRMean = psiRMean;
  fPsiRStDev = psiRStDev;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::SetMultFac(Float_t mean, Float_t stDev) {
  
  fMultFacMean = mean;
  fMultFacStDev = stDev;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::SetStDev(Float_t mult, Float_t temp, Float_t sigma,
  Float_t expVel, Float_t psiR, Float_t Vn, Float_t multFac) {

  fNStDevMult = mult;
  fNStDevTemp = temp;
  fNStDevSigma = sigma;
  fNStDevExpVel = expVel;
  fNStdDevPSIr = psiR;
  fNStDevVn = Vn;
  fNStDevMultFac =multFac;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void AliMevSimConfig::SetGrid(Int_t integr, Int_t scan) {

  fNIntegPts = integr;
  fNScanPts = scan;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
 
