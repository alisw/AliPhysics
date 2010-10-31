#ifndef AliZDCTowerCalibTask_h
#define AliZDCTowerCalibTask_h

// analysis task for ZN towers intercalibration
// Author: Alessandro De Falco, INFN Cagliari

#include "TMatrixD.h"
#include "TVectorD.h"

class TH1F;
class AliESDEvent;

#include "AliAnalysisTask.h"

class AliZDCTowerCalibTask : public AliAnalysisTask {
 public:
  AliZDCTowerCalibTask();
  AliZDCTowerCalibTask(const char *name);
  virtual ~AliZDCTowerCalibTask() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   SetADCMin(Double_t adcmin) { fADCMin = adcmin; } 

  
 private:
  AliESDEvent *fESD;   //ESD object
  TMatrixD fAZNA;      // coefficient matrix for ZNA calorimeter
  TMatrixD fAZNC;      // coefficient matrix for ZNC calorimeter
  TVectorD fBZNA;      // vector of known terms for ZNA calorimeter
  TVectorD fBZNC;      // vector of known terms for ZNC calorimeter
  Double_t fADCMin; 
  AliZDCTowerCalibTask(const AliZDCTowerCalibTask&); // not implemented
  AliZDCTowerCalibTask& operator=(const AliZDCTowerCalibTask&); // not implemented
  
  ClassDef(AliZDCTowerCalibTask, 1); 
};

#endif
