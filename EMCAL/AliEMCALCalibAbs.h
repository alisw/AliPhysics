#ifndef ALIEMCALCALIBABS_H
#define ALIEMCALCALIBABS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
#include "AliEMCALGeoParams.h"
class TString;

/*
  Objects of this class contain basis for absolute calibrations
  AliEMCALCalibAbs inherits TObject only to use AliLog "functions".
*/

class AliEMCALCalibAbs : public TObject {

public:

  enum kProblemType {kNoLED=-999};// code in possible problems

  AliEMCALCalibAbs();

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // debug purposes
  void ReadCalibAbsInfo(Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  void WriteCalibAbsInfo(const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  virtual ~AliEMCALCalibAbs();

  // total calibration factor is a product of
  // a) overall calibration factor [fAbsoluteGain]
  // b) individual gain factor per tower [fRelativeGain]
  // c) time-dependent correction
  // In this class we store a), b) and the needed static ingredients for c)

  // values per single tower
  struct AliEMCALCalibAbsVal {
    Float_t fRelativeGain; // (ADC>GeV relative gain/conversion), value around 1
    Float_t fHighLowRatio; // value around 16 or so
    Int_t fHighLow; // 0 (low) or 1 (high) gain, used for LEDAmp info
    Float_t fLEDAmp; // LED amplitude
    Float_t fLEDAmpRMS; // RMS
  };

  // 1 SuperModule's worth of info: info on where the different APDs are
  struct AliEMCALSuperModuleCalibAbs {
    // first: overall values for the whole SuperModule
    Int_t fSuperModuleNum; // which SuperModule is this?
    Int_t fCalibMethod; // a la 0=cosmics, 1=pi0, 2=electrons,3=using ecore,
    Int_t fCalibPass; // which analysis iteration is this.. 1,2,..N
    Int_t fCalibTime; // t0, unix timestamp
    Float_t fAbsoluteGain; // (ADC>GeV absolute gain/conversion)

    // second: additional info for LED Reference and SM temperature
    Float_t fLEDRefAmp[AliEMCALGeoParams::fgkEMCALLEDRefs]; // LED amplitude at  t0, low gain equivalent
    Float_t fLEDRefAmpRMS[AliEMCALGeoParams::fgkEMCALLEDRefs]; // RMS
    Float_t fLEDRefHighLowRatio[AliEMCALGeoParams::fgkEMCALLEDRefs]; // value around 16 or so
    Int_t fLEDRefHighLow[AliEMCALGeoParams::fgkEMCALLEDRefs]; // 0 (low) or 1 (high) gain

    Float_t fTemperature[AliEMCALGeoParams::fgkEMCALTempSensors]; // temperature at t0
    Float_t fTemperatureRMS[AliEMCALGeoParams::fgkEMCALTempSensors]; // RMS

    // third: individual info for each tower
    AliEMCALCalibAbsVal fAPDVal[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; // at t0
  };

  // pointer to stored info.

  Int_t GetNSuperModule() const { return fNSuperModule; }; 
  AliEMCALSuperModuleCalibAbs * GetSuperModuleData() const { return fSuperModuleData; };

  // - via the index in the stored array:
  virtual AliEMCALSuperModuleCalibAbs GetSuperModuleCalibAbsId(Int_t smIndex) const;
  // - or via the actual SM number
  virtual AliEMCALSuperModuleCalibAbs GetSuperModuleCalibAbsNum(Int_t smNum) const;

protected:

  Int_t 	  fNSuperModule; // Number of supermodules.
  AliEMCALSuperModuleCalibAbs *fSuperModuleData; // SuperModule data

private:

  AliEMCALCalibAbs(const AliEMCALCalibAbs &);
  AliEMCALCalibAbs &operator = (const AliEMCALCalibAbs &);

  ClassDef(AliEMCALCalibAbs, 1) //CalibAbs data reader
};

#endif
