#ifndef ALIEMCALCORRECTIONCLUSTERNONLINEARITY_H
#define ALIEMCALCORRECTIONCLUSTERNONLINEARITY_H

#include "AliEmcalCorrectionComponent.h"

class AliEmcalCorrectionClusterNonLinearity : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionClusterNonLinearity();
  virtual ~AliEmcalCorrectionClusterNonLinearity();

  // Sets up and runs the task
  Bool_t Initialize();
  Bool_t Run();
  
protected:
  TH1F                  *fEnergyDistBefore;          //!<!energy distribution before
  TH2F                  *fEnergyTimeHistBefore;      //!<!energy/time distribution before
  TH1F                  *fEnergyDistAfter;           //!<!energy distribution after
  TH2F                  *fEnergyTimeHistAfter;       //!<!energy/time distribution after
  
 private:
  AliEmcalCorrectionClusterNonLinearity(const AliEmcalCorrectionClusterNonLinearity &);               // Not implemented
  AliEmcalCorrectionClusterNonLinearity &operator=(const AliEmcalCorrectionClusterNonLinearity &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterNonLinearity> reg;

  ClassDef(AliEmcalCorrectionClusterNonLinearity, 1) // EMCal cluster non-linearity correction component
};

#endif /* ALIEMCALCORRECTIONCLUSTERNONLINEARITY_H */
