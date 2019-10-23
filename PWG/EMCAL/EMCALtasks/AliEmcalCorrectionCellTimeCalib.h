#ifndef AliEmcalCorrectionCellTimeCalib_H
#define AliEmcalCorrectionCellTimeCalib_H

#include "AliEmcalCorrectionComponent.h"

/**
 * @class AliEmcalCorrectionCellTimeCalib
 * @ingroup EMCALCORRECTIONFW
 * @brief Time calibration correction component in the EMCal correction framework.
 *
 * Performs time calibration of cells, using OADB calibration. The original cell information in the event **will be overwritten**.
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

class AliEmcalCorrectionCellTimeCalib : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellTimeCalib();
  virtual ~AliEmcalCorrectionCellTimeCalib();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  Bool_t CheckIfRunChanged();
  
protected:
  TH1F* fCellTimeDistBefore;            //!<! cell energy distribution, before time calibration
  TH1F* fCellTimeDistAfter;             //!<! cell energy distribution, after time calibration

private:
  Int_t      InitTimeCalibration();
  Int_t      InitTimeCalibrationL1Phase();
  
  Bool_t                 fCalibrateTime;          ///< flag cell time calibration
  Bool_t                 fCalibrateTimeL1Phase;   ///< flag cell time calibration with L1phase shift
  Bool_t                 fDoMergedBCs;            ///< flag to use one histogram for all BCs
  Bool_t                 fDoCalibrateLowGain;     ///< flag to calibrate the low gain cells
  Bool_t                 fDoCalibMergedLG;        ///< flag to calibrate the low gain cells using a period merged histogram for LG calibration
  
  // Change to false if experts
  Bool_t                 fUseAutomaticTimeCalib;     ///< On by default the check in the OADB of the time recalibration
  
  AliEmcalCorrectionCellTimeCalib(const AliEmcalCorrectionCellTimeCalib &);               // Not implemented
  AliEmcalCorrectionCellTimeCalib &operator=(const AliEmcalCorrectionCellTimeCalib &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellTimeCalib> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellTimeCalib, 4); // EMCal cell time calibration component
  /// \endcond
};

#endif /* AliEmcalCorrectionCellTimeCalib_H */
