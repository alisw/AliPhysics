#ifndef ALIEMCALCORRECTIONCELLSINGLECHANNELCALIBRATION_H
#define ALIEMCALCORRECTIONCELLSINGLECHANNELCALIBRATION_H

#include "AliEmcalCorrectionComponent.h"

/**
 * @class AliEmcalCorrectionCellSingleChannelCalibration
 * @ingroup EMCALCOREFW
 * @brief Energy calibration correction using single channel calibration in the EMCal correction framework.
 *
 * Performs energy calibration of cells, using OADB calibration. The original cell information in the event **will be overwritten**.
 *
 * Based on code in AliEMCALTenderSupply.
 *
 * @author Deepa Thomas (Utrecht University), AliEMCALTenderSupply
 * @author Jiri Kral (University of Jyvaskyla), AliEMCALTenderSupply mods/rewrite
 * @author Salvatore Aiola, make AliEMCALTenderSupply work for AODs
 * @author C. Loizides, make AliEMCALTenderSupply work for AODs
 * @author Gustavo Conesa, LPSC-Grenoble, AliEMCALTenderSupply several mods.
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author Dhruv Dixit <dhruvdixit@berkeley.edu>, University of California, Berkeley, AliEmcalCorrectionCellSingleChannelCalibration
 * @date March 17, 2019
 */

class AliEmcalCorrectionCellSingleChannelCalibration : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellSingleChannelCalibration();
  virtual ~AliEmcalCorrectionCellSingleChannelCalibration();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  Bool_t CheckIfRunChanged();
  
 protected:
  TH1F* fCellSingleChannelEnergyDistBefore;        //!<! cell energy distribution, before energy calibration
  TH1F* fCellSingleChannelEnergyDistAfter;         //!<! cell energy distribution, after energy calibration

private:
  Int_t                  InitRecalib();
  
  // Change to false if experts
  Bool_t                 fUseAutomaticRecalib;       ///< On by default the check in the OADB of the energy recalibration
  TString                fCustomRecalibFilePath;     ///< Empty string by default the path to the OADB file of the custom energy recalibration
  
  AliEmcalCorrectionCellSingleChannelCalibration(const AliEmcalCorrectionCellSingleChannelCalibration &);               // Not implemented
  AliEmcalCorrectionCellSingleChannelCalibration &operator=(const AliEmcalCorrectionCellSingleChannelCalibration &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellSingleChannelCalibration> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellSingleChannelCalibration, 4); // EMCal cell energy correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCELLSINGLECHANNELCALIBRATION_H */
