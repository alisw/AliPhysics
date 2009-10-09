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
class AliCaloCalibSignal;
class AliEMCALBiasAPD;
class AliEMCALCalibMapAPD;
class AliEMCALCalibReference; 
class AliEMCALCalibTimeDepCorrection; 

class AliEMCALCalibTimeDep : public TObject {

 public:
  AliEMCALCalibTimeDep(); //! ctor
  AliEMCALCalibTimeDep(const AliEMCALCalibTimeDep &calibt); //! copy ctor
  AliEMCALCalibTimeDep& operator= (const AliEMCALCalibTimeDep &calibt); //! 
  virtual ~AliEMCALCalibTimeDep(); //! dtor
  virtual void PrintInfo() const; 

  void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);//!

  // simple getters
  Int_t GetRunNumber() const { return fRun; } //
  Int_t GetStartTime() const { return fStartTime; } // 
  Int_t GetEndTime() const { return fEndTime; } // 
  Double_t GetLengthOfRunInHours() const; // 
  Double_t GetLengthOfRunInBins() const; // 

  // Temperature Section
  // access pointer to the basic temperature object: AliEMCALSensorTempArray 
  AliEMCALSensorTempArray  * GetTempArray() const { return fTempArray; } //

  // range of temperature readings/values during the run 
  Int_t ScanTemperatureInfo(); // go through the temperature info
  Double_t GetMinTemp() const { return fMinTemp; } // 
  Double_t GetMaxTemp() const { return fMaxTemp; } // 
  UInt_t GetMinTime() const { return fMinTime; } // 
  UInt_t GetMaxTime() const { return fMaxTime; } // 
  Double_t GetRangeOfTempMeasureInHours() const; // 
  Double_t GetRangeOfTempMeasureInDegrees() const; // 

  // basic calibration info
  Double_t GetTemperature(UInt_t timeStamp) const; // for all sensors, all SuperModules
  Double_t GetTemperatureSM(int imod, UInt_t timeStamp) const; // for all sensors, in a SuperModule
  Double_t GetTemperatureSMSensor(int imod, int isens, UInt_t timeStamp) const; // for a sensor, in a SuperModule
  // End of Temperature Section

  // control parameters
  void SetTemperatureResolution(Double_t d) { fTemperatureResolution = d; } // value for checking at which level we care about temperature differences
  void SetTimeBinsPerHour(Int_t i) { fTimeBinsPerHour = i; } // size of the time-bins we use for corrections
  Double_t GetTemperatureResolution() const { return fTemperatureResolution; } // value for checking at which level we care about temperature differences
  Int_t GetTimeBinsPerHour() const { return fTimeBinsPerHour; } // size of the time-bins we use foc corrections

  void SetHighLowGainFactor(Double_t value) {fHighLowGainFactor = value;}
  Double_t GetHighLowGainFactor() const {return fHighLowGainFactor;}

  // access to other pointers
  AliCaloCalibSignal  * GetCalibSignal() const { return fCalibSignal; } //
  AliEMCALBiasAPD  * GetBiasAPD() const { return fBiasAPD; } //
  AliEMCALCalibMapAPD  * GetCalibMapAPD() const { return fCalibMapAPD; } //
  AliEMCALCalibReference  * GetCalibReference() const { return fCalibReference; } //
  AliEMCALCalibTimeDepCorrection  * GetCalibTimeDepCorrection() const { return fCalibTimeDepCorrection; } //

  // storage and access of the correction info
  Int_t CalcCorrection(); //
  AliEMCALCalibTimeDepCorrection  * GetTimeDepCorrection() 
    const { return fCalibTimeDepCorrection; } //

  Double_t GetTempCoeff(Double_t IDark, Double_t M) const; //

 private:

  void Reset(); //

  void GetTemperatureInfo(); // pick up Preprocessor output
  void GetCalibSignalInfo(); // pick up Preprocessor output
  void GetBiasAPDInfo(); // pick up OCDB info
  void GetCalibMapAPDInfo(); // pick up OCDB info
  void GetCalibReferenceInfo(); // pick up OCDB info

  Int_t CalcLEDCorrection(Int_t nSM, Int_t nBins); // based on LED signals, and reference photodiodes
  Int_t CalcLEDCorrectionStripBasis(Int_t nSM, Int_t nBins); // based on LED signals, and reference photodiodes
  Int_t CalcTemperatureCorrection(Int_t nSM, Int_t nBins); // based on temperetare info

  //
  Int_t fRun;
  UInt_t fStartTime;
  UInt_t fEndTime;
  // temperature stuff
  Double_t fMinTemp;
  Double_t fMaxTemp;
  UInt_t fMinTime;
  UInt_t fMaxTime;
  //
  Double_t fTemperatureResolution; // value for checking at which level we care about temperature differences
  Int_t fTimeBinsPerHour; // size of the time-bins we use for corrections

  Double_t fHighLowGainFactor;     //gain factor to convert between high and low gain
  
  // pointers to the different used classes
  AliEMCALSensorTempArray  *fTempArray;     // CDB class for temperature sensors
  AliCaloCalibSignal *fCalibSignal; //
  AliEMCALBiasAPD *fBiasAPD; //
  AliEMCALCalibMapAPD *fCalibMapAPD; //
  AliEMCALCalibReference *fCalibReference; 
  AliEMCALCalibTimeDepCorrection *fCalibTimeDepCorrection; // 

  //
  ClassDef(AliEMCALCalibTimeDep,3)    // EMCAL time-dep Calibration data
};

#endif
