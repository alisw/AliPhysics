#ifndef ALIEMCALCORRECTIONCELLBADCHANNEL_H
#define ALIEMCALCORRECTIONCELLBADCHANNEL_H

#include "AliEmcalCorrectionComponent.h"

class AliEmcalCorrectionCellBadChannel : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellBadChannel();
  virtual ~AliEmcalCorrectionCellBadChannel();

  // Sets up and runs the task
  Bool_t Initialize();
  Bool_t Run();
  
protected:
  TH1F* fCellEnergyDistBefore;
  TH1F* fCellEnergyDistAfter;
  
private:
  Int_t      InitBadChannels();

  AliEmcalCorrectionCellBadChannel(const AliEmcalCorrectionCellBadChannel &);             // Not implemented
  AliEmcalCorrectionCellBadChannel &operator=(const AliEmcalCorrectionCellBadChannel &);   // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellBadChannel> reg;

  ClassDef(AliEmcalCorrectionCellBadChannel, 1) // EMCal cell bad channel correction component
};

#endif /* ALIEMCALCORRECTIONCELLBADCHANNEL_H */
