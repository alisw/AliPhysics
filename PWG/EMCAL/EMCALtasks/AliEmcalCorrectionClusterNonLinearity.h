#ifndef ALIEMCALCORRECTIONCLUSTERNONLINEARITY_H
#define ALIEMCALCORRECTIONCLUSTERNONLINEARITY_H

#include "AliEmcalCorrectionComponent.h"

#include "AliEMCALRecoUtils.h"

/**
 * @class AliEmcalCorrectionClusterNonLinearity
 * @ingroup EMCALCORRECTIONFW
 * @brief Cluster energy non-linearity correction component in the EMCal correction framework.
 *
 * Non-linearity correction to the cluster energy is necessary because the response of the calorimeter is not linear for very low momentum particles or very high momentum (shower leakage).
 
 The energy of the cluster **after** the non-linearity correction can be retrieved using the method `cluster->GetNonLinCorrEnergy()`.
 *
 * Based on code in AliEmcalClusterMaker.
 *
 * @author Constantin Loizides, LBNL, AliEmcalClusterMaker
 * @author Salvatore Aiola, LBNL, AliEmcalClusterMaker
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Jul 8, 2016
 */


class AliEmcalCorrectionClusterNonLinearity : public AliEmcalCorrectionComponent {
 public:
  /// Relates string to the non-linearity function enumeration for %YAML configuration
  static const std::map <std::string, AliEMCALRecoUtils::NonlinearityFunctions> fgkNonlinearityFunctionMap; //!<!

  AliEmcalCorrectionClusterNonLinearity();
  virtual ~AliEmcalCorrectionClusterNonLinearity();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();

protected:
  TH1F                  *fEnergyDistBefore;          //!<!energy distribution before
  TH2F                  *fEnergyTimeHistBefore;      //!<!energy/time distribution before
  TH1F                  *fEnergyDistAfter;           //!<!energy distribution after
  TH2F                  *fEnergyTimeHistAfter;       //!<!energy/time distribution after
  
  Bool_t                 fSetForceClusterE;          ///< Only for backwards compatibility, force cluster->E() to be set to the cluster non-linearity corrected energy. Off by default. For the standard methods, see: http://alidoc.cern.ch/AliPhysics/master/READMEcontainers.html#emcalContainerClusterEnergyCorrections
  
 private:
  AliEmcalCorrectionClusterNonLinearity(const AliEmcalCorrectionClusterNonLinearity &);               // Not implemented
  AliEmcalCorrectionClusterNonLinearity &operator=(const AliEmcalCorrectionClusterNonLinearity &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterNonLinearity> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterNonLinearity, 3); // EMCal cluster non-linearity correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERNONLINEARITY_H */
