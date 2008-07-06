#ifndef ALIEMCALCALIBTIMEDEP_H
#define ALIEMCALCALIBTIMEDEP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEMCALCalibTimeDep.h $ */

/////////////////////////////////////////////////
// class for EMCAL time-dependence calibration //
/////////////////////////////////////////////////

#include "TObject.h"

class AliEMCALSensorTempArray;

class AliEMCALCalibTimeDep : public TObject {

 public:
  AliEMCALCalibTimeDep(); //! ctor
  AliEMCALCalibTimeDep(const AliEMCALCalibTimeDep &calibt); //! copy ctor
  AliEMCALCalibTimeDep& operator= (const AliEMCALCalibTimeDep &calibt); //! 
  virtual ~AliEMCALCalibTimeDep(); //! dtor
  virtual void Print() const; 

  void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);//!

  // access pointer to the basic temperature object: AliEMCALSensorTempArray 
  AliEMCALSensorTempArray  * GetTempArray() const { return fTempArray; } //

  // simple getters
  int GetRunNumber() const { return fRun; } //
  int GetStartTime() const { return fStartTime; } // 
  int GetEndTime() const { return fEndTime; } // 
  double GetLengthOfRunInHours() const; // 

  // range of temperature readings/values during the run 
  double GetMinTemp() const { return fMinTemp; } // 
  double GetMaxTemp() const { return fMaxTemp; } // 
  int GetMinTime() const { return fMinTime; } // 
  int GetMaxTime() const { return fMaxTime; } // 
  double GetRangeOfTempMeasureInHours() const; // 
  double GetRangeOfTempMeasureInDegrees() const; // 

  // reference temperature, that we normalize to
  double GetRefTemp() const { return fRefTemp; } //
  void SetRefTemp(double refTemp) { fRefTemp = refTemp; } //

  // basic calibration info
  double GetTemperature(int secSinceStart) const; // for all sensors, all SuperModules
  double GetTemperatureSM(int imod, int secSinceStart) const; // for all sensors, in a SuperModule
  double GetTemperatureSMSensor(int imod, int isens, int secSinceStart) const; // for a sensor, in a SuperModule
  double GetCorrection(double temperature) const; //

 private:

  void Reset();
  void GetTemperatureInfo(); // pick up Preprocessor output
  //
  int fRun;
  int fStartTime;
  int fEndTime;
  double fMinTemp;
  double fMaxTemp;
  int fMinTime;
  int fMaxTime;
  double fRefTemp;
  //
  AliEMCALSensorTempArray  *fTempArray;     // CDB class for temperature sensors
  //
  ClassDef(AliEMCALCalibTimeDep,1)    // EMCAL Calibration data
};

#endif
