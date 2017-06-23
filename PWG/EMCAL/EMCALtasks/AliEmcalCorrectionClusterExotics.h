#ifndef ALIEMCALCORRECTIONCLUSTEREXOTICS_H
#define ALIEMCALCORRECTIONCLUSTEREXOTICS_H

#include "AliEmcalCorrectionComponent.h"

class TH2F;

/**
 * @class AliEmcalCorrectionClusterExotics
 * @ingroup EMCALCOREFW
 * @brief Exotic cluster removal in the EMCal correction framework.
 *
 * "Exotic" clusters are energetic clusters where most energy deposition is concentrated in one single cell. This clusters are not reproduced in MC simulations and are believed to arise from neutrons showering directly into the APD. These clusters need to be flagged, so that they can be easily rejected during the analysis.
 
 The "exotic" flag can be retrieved using `cluster->GetIsExotic()`. "Exotic" clusters can be easily rejected if clusters are accessed using an AliClusterContainer object. "Exotic" cluster removal is switched on by default in AliClusterContainer, however **it is necessary to run the ClusterExotics component (via AliEmcalCorrectionTask) to flag "exotic" cluster beforehand**.
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
  void UserCreateOutputObjects();
  Bool_t Run();

protected:
  TH2F                  *fEtaPhiDistBefore;          //!<!eta/phi distribution before
  TH2F                  *fEtaPhiDistAfter;           //!<!eta/phi distribution after
  TH1F                  *fEnergyExoticClusters;      //!<!energy of exotic clusters
  Float_t               fExoticMinCellAmplitude;     ///< Min energy of leading cell in order for exotic cut to be attempted
  Float_t               fMaxFcross;                  ///< Max value of Fcross = 1-Ecross/ecell allowed for clusters to pass exotic cut
  Float_t               fCellCrossMaxTimeDiff;       ///< Max time difference allowed between leading cell and cross cells (in ns)

 private:
  AliEmcalCorrectionClusterExotics(const AliEmcalCorrectionClusterExotics &);               // Not implemented
  AliEmcalCorrectionClusterExotics &operator=(const AliEmcalCorrectionClusterExotics &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterExotics> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterExotics, 2); // EMCal cluster exotics correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTEREXOTICS_H */
