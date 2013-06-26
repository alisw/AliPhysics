#ifndef ALIEMCALTRIGGERSETUPINFO_H
#define ALIEMCALTRIGGERSETUPINFO_H

// $Id$

#include "TNamed.h"

static const Double_t kEMCL1ADCtoGeV = 0.07874;

class AliEmcalTriggerSetupInfo: public TNamed {
 public:
  AliEmcalTriggerSetupInfo();
  AliEmcalTriggerSetupInfo(const AliEmcalTriggerSetupInfo &p); 
  AliEmcalTriggerSetupInfo &operator=(const AliEmcalTriggerSetupInfo &p);
  virtual ~AliEmcalTriggerSetupInfo();

  Int_t GetThresholdJetLow() const { return fThresholds[2]; }
  Int_t GetThresholdJetHigh() const { return fThresholds[0]; }
  
   Double_t GetThresholdGeVRoughJetLow() const { return ((Double_t)fThresholds[2])*kEMCL1ADCtoGeV; }
   Double_t GetThresholdGeVRoughJetHigh() const { return ((Double_t)fThresholds[0])*kEMCL1ADCtoGeV; }
  
  void SetThresholds( Int_t i0, Int_t i1, Int_t i2, Int_t i3 ) {
            fThresholds[0] = i0; fThresholds[1] = i1; fThresholds[2] = i2; fThresholds[3] = i3;}
            
  void Clean();


 protected:
  Int_t             fThresholds[4];                 // per event L1 online thresholds in ADC counts

  ClassDef(AliEmcalTriggerSetupInfo, 1) // Emcal trigger setup class
};
#endif
