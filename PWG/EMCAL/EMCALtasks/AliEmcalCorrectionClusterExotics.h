#ifndef ALIEMCALCORRECTIONCLUSTEREXOTICS_H
#define ALIEMCALCORRECTIONCLUSTEREXOTICS_H

#include "AliEmcalCorrectionComponent.h"

class TH2F;

/**
 * @class AliEmcalCorrectionClusterExotics
 * @brief Exotic cluster removal in the EMCal correction framework
 * @ingroup EMCALCOREFW
 *
 * Based on code in AliEmcalClusterMaker.
 *
 * @author Constantin Loizides, LBNL, AliEmcalClusterMaker
 * @author Salvatore Aiola, LBNL, AliEmcalClusterMaker
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Jul 8, 2016
 */


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

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterExotics, 1); // EMCal cluster exotics correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTEREXOTICS_H */
