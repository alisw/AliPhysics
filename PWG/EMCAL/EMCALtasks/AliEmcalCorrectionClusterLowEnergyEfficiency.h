#ifndef ALIEMCALCORRECTIONCLUSTERLOWENERGYEFFICIENCY_H
#define ALIEMCALCORRECTIONCLUSTERLOWENERGYEFFICIENCY_H

#include "AliEmcalCorrectionComponent.h"

#include "AliEMCALRecoUtils.h"

/**
 * @class AliEmcalCorrectionClusterLowEnergyEfficiency
 * @ingroup EMCALCORRECTIONFW
 * @brief Cluster number of cells correction component in the EMCal correction framework.
 *
 * This task flags a fraction of the 1 cell clusters as 2 cell clusters. This is necessary because the fraction of 1 cell clusters does not match between data and MC.
   If a cut on the number of cells is applied this correction is necessary to get correct results.
   The artificially widened clusters will have the flag chi2 set to 1, other 1 cell clusters are set to chi2 = 0

 *
 * Based on code in AliEmcalClusterMaker.
 *
 * @author Constantin Loizides, LBNL, AliEmcalClusterMaker
 * @author Salvatore Aiola, LBNL, AliEmcalClusterMaker
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Jul 8, 2016
 */


class AliEmcalCorrectionClusterLowEnergyEfficiency : public AliEmcalCorrectionComponent {
 public:
  /// Relates string to the number of cells efficiency function enumeration for %YAML configuration
  static const std::map <std::string, AliEMCALRecoUtils::NCellEfficiencyFunctions> fgkNCellEfficiencyFunctionMap; //!<!

  AliEmcalCorrectionClusterLowEnergyEfficiency();
  virtual ~AliEmcalCorrectionClusterLowEnergyEfficiency();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();

protected:
  TH1F                  *fNCellDistBefore;          //!<!number of cells distribution before
  TH1F                  *fNCellDistAfter;           //!<!number of cells distribution after

  Bool_t                fRejectNextToClus;          //!<!Switch waether to reject 1 cell clusters with direct neighbours

 private:
  AliEmcalCorrectionClusterLowEnergyEfficiency(const AliEmcalCorrectionClusterLowEnergyEfficiency &);               // Not implemented
  AliEmcalCorrectionClusterLowEnergyEfficiency &operator=(const AliEmcalCorrectionClusterLowEnergyEfficiency &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterLowEnergyEfficiency> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterLowEnergyEfficiency, 1); // EMCal cluster non-linearity correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERLOWENERGYEFFICIENCY_H */
