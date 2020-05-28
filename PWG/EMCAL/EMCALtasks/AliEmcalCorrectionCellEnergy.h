#ifndef ALIEMCALCORRECTIONCELLENERGY_H
#define ALIEMCALCORRECTIONCELLENERGY_H

#include "AliEmcalCorrectionComponent.h"

/**
 * @class AliEmcalCorrectionCellEnergy
 * @ingroup EMCALCORRECTIONFW
 * @brief Energy calibration correction component in the EMCal correction framework.
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
 * @date Jul 8, 2016
 */

class AliEmcalCorrectionCellEnergy : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellEnergy();
  virtual ~AliEmcalCorrectionCellEnergy();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  Bool_t CheckIfRunChanged();
  
protected:
  TH1F* fCellEnergyDistBefore;        //!<! cell energy distribution, before energy calibration
  TH1F* fCellEnergyDistAfter;         //!<! cell energy distribution, after energy calibration

private:
  Int_t                  InitRecalib();
  Int_t                  InitRunDepRecalib();
  
  // Change to false if experts
  Bool_t                 fUseAutomaticRecalib;       ///< On by default the check in the OADB of the energy recalibration
  Bool_t                 fUseAutomaticRunDepRecalib; ///< On by default the check in the OADB of the run dependent energy recalibration
  Bool_t                 fUseNewRunDepTempCalib;     ///< Off by default the check in the OADB of the new run dependent temp calib Run1/Run2
  Bool_t                 fDisableTempCalib;          ///< Off by default, disables temp calibration totally
  Bool_t                 fUseShaperCorrection;       ///< Off by default the correction for the shaper nonlinearity
  TString                fCustomRecalibFilePath;     ///< Empty string by default the path to the OADB file of the custom energy recalibration
  Bool_t                 fLoad1DRecalibFactors;      ///< Flag to load 1D energy recalibration factors
  
  AliEmcalCorrectionCellEnergy(const AliEmcalCorrectionCellEnergy &);               // Not implemented
  AliEmcalCorrectionCellEnergy &operator=(const AliEmcalCorrectionCellEnergy &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellEnergy> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellEnergy, 7); // EMCal cell energy correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCELLENERGY_H */
