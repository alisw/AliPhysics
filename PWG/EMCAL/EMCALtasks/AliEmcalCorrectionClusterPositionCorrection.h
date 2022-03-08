#ifndef ALIEMCALCORRECTIONCLUSTERPOSITIONCORRECTION_H
#define ALIEMCALCORRECTIONCLUSTERPOSITIONCORRECTION_H

#include "AliEmcalCorrectionComponent.h"

#include "AliEMCALRecoUtils.h"

#include <AliEMCALGeometry.h>

/**
 * @class AliEmcalCorrectionClusterPositionCorrection
 * @ingroup EMCALCORRECTIONFW
 * @brief Cluster position correction in eta and phi
 *
 * This task corrects a small missalignemt between the central detectors and the emcal in eta and phi for each SM.
 * The shifts were evaluated by calculating the mean difference between eta/phi of a matched track and the corresponding cluster

 *
 * Based on code in AliEmcalClusterMaker.
 *
 * @author Constantin Loizides, LBNL, AliEmcalClusterMaker
 * @author Salvatore Aiola, LBNL, AliEmcalClusterMaker
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Jul 8, 2016
 */


class AliEmcalCorrectionClusterPositionCorrection : public AliEmcalCorrectionComponent {
 public:

  AliEmcalCorrectionClusterPositionCorrection();
  virtual ~AliEmcalCorrectionClusterPositionCorrection();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();

protected:
  TH1F                  *fClusterPositionDiffEta;         //!<! shift of cluster position  in eta
  TH1F                  *fClusterPositionDiffPhi;         //!<! shift of cluster position  in phi

  std::vector<Double_t> fSMEtaValues;                     ///< vector containing position shifts in eta
  std::vector<Double_t> fSMPhiValues;                     ///< vector containing position shifts in phi

  Bool_t                fApplyToMC;                       ///< switch if correction should be applied on data or MC

 private:
  AliEmcalCorrectionClusterPositionCorrection(const AliEmcalCorrectionClusterPositionCorrection &);               // Not implemented
  AliEmcalCorrectionClusterPositionCorrection &operator=(const AliEmcalCorrectionClusterPositionCorrection &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterPositionCorrection> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterPositionCorrection, 2); // EMCal cluster low energy efficiency correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERPOSITIONCORRECTION_H */
