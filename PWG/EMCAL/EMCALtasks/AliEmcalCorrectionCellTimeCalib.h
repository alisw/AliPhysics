#ifndef AliEmcalCorrectionCellTimeCalib_H
#define AliEmcalCorrectionCellTimeCalib_H

#include "AliEmcalCorrectionComponent.h"

class AliEmcalCorrectionCellTimeCalib : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellTimeCalib();
  virtual ~AliEmcalCorrectionCellTimeCalib();

  // Sets up and runs the task
  Bool_t Initialize();
  Bool_t Run();
  
protected:
  TH1F* fCellTimeDistBefore;
  TH1F* fCellTimeDistAfter;

private:
  Int_t      InitTimeCalibration();
  Int_t      InitTimeCalibrationL1Phase();
  
  Bool_t                 fCalibrateTime;          // flag cell time calibration
  Bool_t                 fCalibrateTimeL1Phase;   // flag cell time calibration with L1phase shift
  
  // Change to false if experts
  Bool_t                 fUseAutomaticTimeCalib;     // On by default the check in the OADB of the time recalibration
  
  AliEmcalCorrectionCellTimeCalib(const AliEmcalCorrectionCellTimeCalib &);               // Not implemented
  AliEmcalCorrectionCellTimeCalib &operator=(const AliEmcalCorrectionCellTimeCalib &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellTimeCalib> reg;

  ClassDef(AliEmcalCorrectionCellTimeCalib, 1) // EMCal cell energy correction component
};

#endif /* AliEmcalCorrectionCellTimeCalib_H */
