#ifndef ALIEMCALCALIBABS_H
#define ALIEMCALCALIBABS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
#include "AliEMCALGeoParams.h"
class TString;
class TTree;

/*
  Objects of this class contain basis for absolute calibrations
*/

// total calibration factor is a product of
// a) overall calibration factor [fAbsoluteGain]
// b) individual gain factor per tower [fRelativeGain]
// c) time-dependent correction
// In this class we store a), b) and the needed static ingredients for c)

// ******* internal class definition *************
// values per single tower
class AliEMCALCalibAbsVal : public TObject {

 public:
  AliEMCALCalibAbsVal() : TObject(), // just init values
    fRelativeGain(0),
    fHighLowRatio(0),
    fHighLow(0),
    fLEDAmp(0),
    fLEDAmpRMS(0)
    {
    }
  
  void Init() {
    fRelativeGain = 0;
    fHighLowRatio = 0;
    fHighLow = 0;
    fLEDAmp = 0;
    fLEDAmpRMS = 0;
    return;
  }

 public:
  Float_t fRelativeGain; // (ADC>GeV relative gain/conversion), value around 1
  Float_t fHighLowRatio; // value around 16 or so
  Int_t fHighLow; // 0 (low) or 1 (high) gain, used for LEDAmp info
  Float_t fLEDAmp; // LED amplitude
  Float_t fLEDAmpRMS; // RMS

  ClassDef(AliEMCALCalibAbsVal, 1) // help class
};

// 1 SuperModule's worth of info: info on where the different APDs are
class AliEMCALSuperModuleCalibAbs : public TObject {

 public:
  AliEMCALSuperModuleCalibAbs() : TObject(), // just init values
    fSuperModuleNum(0),
    fCalibMethod(0),
    fCalibPass(0),
    fCalibTime(0),
    fAbsoluteGain(0)
    {
      for (int iref=0; iref<AliEMCALGeoParams::fgkEMCALLEDRefs; iref++) {
	fLEDRefAmp[iref] = 0;
	fLEDRefAmpRMS[iref] = 0;
	fLEDRefHighLowRatio[iref] = 0;
	fLEDRefHighLow[iref] = 0;
      }

      for (int itemp=0; itemp<AliEMCALGeoParams::fgkEMCALTempSensors; itemp++) {
	fTemperature[itemp] = 0;
	fTemperatureRMS[itemp] = 0;
      }

      for (int icol=0; icol<AliEMCALGeoParams::fgkEMCALCols; icol++) {
	for (int irow=0; irow<AliEMCALGeoParams::fgkEMCALRows; irow++) {
	  fAPDVal[icol][irow].Init();
	}
      }
    }

 public:
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

  ClassDef(AliEMCALSuperModuleCalibAbs, 1) // help class
};
// ******* end of internal class definition *************

class AliEMCALCalibAbs : public TObject {

public:

  enum kProblemType {kNoLED=-999};// code in possible problems

  AliEMCALCalibAbs();

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadTextCalibAbsInfo(Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void WriteTextCalibAbsInfo(const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void ReadRootCalibAbsInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void ReadTreeCalibAbsInfo(TTree *tree, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void WriteRootCalibAbsInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  virtual ~AliEMCALCalibAbs();

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

  ClassDef(AliEMCALCalibAbs, 2) //CalibAbs data reader
};

#endif
