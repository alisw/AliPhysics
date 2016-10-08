#ifndef ALIEMCALCORRECTIONCLUSTEREXOTICS_H
#define ALIEMCALCORRECTIONCLUSTEREXOTICS_H

#include "AliEmcalCorrectionComponent.h"

class TH2F;

class AliEmcalCorrectionClusterExotics : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionClusterExotics();
  virtual ~AliEmcalCorrectionClusterExotics();

  // Sets up and runs the task
  Bool_t Initialize();
  Bool_t Run();

protected:
  TH2F                  *fEtaPhiDistBefore;          //!<!eta/phi distribution before
  TH2F                  *fEtaPhiDistAfter;           //!<!eta/phi distribution after
  TH1F                  *fEnergyExoticClusters;      //!<!energy of exotic clusters
    
 private:
  AliEmcalCorrectionClusterExotics(const AliEmcalCorrectionClusterExotics &);               // Not implemented
  AliEmcalCorrectionClusterExotics &operator=(const AliEmcalCorrectionClusterExotics &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterExotics> reg;

  ClassDef(AliEmcalCorrectionClusterExotics, 1) // EMCal cluster exotics correction component
};

#endif /* ALIEMCALCORRECTIONCLUSTEREXOTICS_H */
