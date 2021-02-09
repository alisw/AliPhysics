#ifndef ALIEMCALCORRECTIONCLUSTERENERGYVARIATION_H
#define ALIEMCALCORRECTIONCLUSTERENERGYVARIATION_H

#include "AliEmcalCorrectionComponent.h"

#include "TRandom3.h"

/**
 * @class AliEmcalCorrectionClusterEnergyVariation
 * @ingroup EMCALCORRECTIONFW
 * @brief Cluster energy variation component in the EMCal correction framework, for EMCal systematics.
 *
 * This component allows the cluster energy to be scaled by a constant factor and smeared by a random value drawn
 * from a Gaussian with specified width. Note that you must configure the input cluster container (in the YAML file)
 * with the desired default cluster energy type (nonlincorr, hadcorr, etc) to be shifted/smeared.
 * The shifted/smeared energy is written into the user-defined cluster energy field AliVCluster::kUserDefEnergy1,
 * which can then be retrieved in later tasks through AliVCluster::GetUserDefEnergy() or through the usual cluster
 * container machinery (see http://alidoc.cern.ch/AliPhysics/master/READMEcontainers.html#emcalContainerClusterEnergyCorrections).
 *
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University
 * @date Oct 30 2018
 */

class AliEmcalCorrectionClusterEnergyVariation : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionClusterEnergyVariation();
  virtual ~AliEmcalCorrectionClusterEnergyVariation();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  void ExecOnce();
  Bool_t Run();
  
protected:
  Double_t               fEnergyScaleShift;               ///< Fraction of cluster energy to shift (positive means upward shift)
  Double_t               fSmearingWidth;                  ///< Width of Gaussian used to smear cluster energy
  
  TRandom3               fRandom;                         //!<! Random number generator for cluster energy smearing

 private:
  AliEmcalCorrectionClusterEnergyVariation(const AliEmcalCorrectionClusterEnergyVariation &);               // Not implemented
  AliEmcalCorrectionClusterEnergyVariation &operator=(const AliEmcalCorrectionClusterEnergyVariation &);    // Not implemented

  // Allows the registration of the class so that it is available to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterEnergyVariation> reg;
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterEnergyVariation, 1); // EMCal cluster energy variation component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERENERGYVARIATION_H */
