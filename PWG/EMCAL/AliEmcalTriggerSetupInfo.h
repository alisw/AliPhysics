#ifndef ALIEMCALTRIGGERSETUPINFO_H
#define ALIEMCALTRIGGERSETUPINFO_H

// $Id$

#include "TNamed.h"

static const Double_t kEMCL1ADCtoGeV = 0.07874;
static const Double_t kEMCL1ADCtoADCSum = 3.40;

class AliEmcalTriggerSetupInfo: public TNamed {
 public:
  AliEmcalTriggerSetupInfo();
  AliEmcalTriggerSetupInfo(const AliEmcalTriggerSetupInfo &p); 
  AliEmcalTriggerSetupInfo &operator=(const AliEmcalTriggerSetupInfo &p);
  virtual ~AliEmcalTriggerSetupInfo();

  Int_t GetThresholdJetLow() const { return fThresholds[2]; }
  Int_t GetThresholdJetHigh() const { return fThresholds[0]; }
  Int_t GetThresholdJetLowSimple() const { return fThresholdsSimple[2]; }
  Int_t GetThresholdJetHighSimple() const { return fThresholdsSimple[0]; }
  
  Int_t GetThresholdGammaLow() const { return fThresholds[3]; }
  Int_t GetThresholdGammaHigh() const { return fThresholds[1]; }
  Int_t GetThresholdGammaLowSimple() const { return fThresholdsSimple[3]; }
  Int_t GetThresholdGammaHighSimple() const { return fThresholdsSimple[1]; }


   Double_t GetThresholdGeVRoughJetLow() const { return ((Double_t)fThresholds[2])*kEMCL1ADCtoGeV; }
   Double_t GetThresholdGeVRoughJetHigh() const { return ((Double_t)fThresholds[0])*kEMCL1ADCtoGeV; }
   Double_t GetThresholdGeVRoughJetLowSimple() const { return ((Double_t)fThresholdsSimple[2])*kEMCL1ADCtoGeV; }
   Double_t GetThresholdGeVRoughJetHighSimple() const { return ((Double_t)fThresholdsSimple[0])*kEMCL1ADCtoGeV; }
  
  void SetThresholds( Int_t i0, Int_t i1, Int_t i2, Int_t i3 ) {
            fThresholds[0] = i0; fThresholds[1] = i1; fThresholds[2] = i2; fThresholds[3] = i3;}
  void SetThresholdsSimple( Int_t i0, Int_t i1, Int_t i2, Int_t i3 ) {
            fThresholdsSimple[0] = i0; fThresholdsSimple[1] = i1; fThresholdsSimple[2] = i2; fThresholdsSimple[3] = i3;}
            
  void Clean();


 protected:
  Int_t             fThresholds[4];                 // per event L1 online thresholds in ADC counts
  Int_t             fThresholdsSimple[4];           // per event L1 simple offline thresholds

  ClassDef(AliEmcalTriggerSetupInfo, 2) // Emcal trigger setup class
};
#endif
