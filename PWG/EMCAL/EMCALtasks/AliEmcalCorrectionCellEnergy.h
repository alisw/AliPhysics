#ifndef ALIEMCALCORRECTIONCELLENERGY_H
#define ALIEMCALCORRECTIONCELLENERGY_H

#include "AliEmcalCorrectionComponent.h"

class AliEmcalCorrectionCellEnergy : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellEnergy();
  virtual ~AliEmcalCorrectionCellEnergy();

  // Sets up and runs the task
  Bool_t Initialize();
  Bool_t Run();
  
protected:
  TH1F* fCellEnergyDistBefore;
  TH1F* fCellEnergyDistAfter;

private:
  Int_t                  InitRecalib();
  Int_t                  InitRunDepRecalib();
  
  // Change to false if experts
  Bool_t                 fUseAutomaticRecalib;       // On by default the check in the OADB of the energy recalibration
  Bool_t                 fUseAutomaticRunDepRecalib; // On by default the check in the OADB of the run dependent energy recalibration
  
  AliEmcalCorrectionCellEnergy(const AliEmcalCorrectionCellEnergy &);               // Not implemented
  AliEmcalCorrectionCellEnergy &operator=(const AliEmcalCorrectionCellEnergy &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellEnergy> reg;

  ClassDef(AliEmcalCorrectionCellEnergy, 1) // EMCal cell energy correction component
};

#endif /* ALIEMCALCORRECTIONCELLENERGY_H */
