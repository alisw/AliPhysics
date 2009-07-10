#ifndef ALIEMCALCALIBMAPAPD_H
#define ALIEMCALCALIBMAPAPD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
#include "AliEMCALGeoParams.h"
#include <cmath>
class TString;

/*
  Objects of this class read txt file with APD data
  AliEMCALCalibMapAPD inherits TObject only to use AliLog "functions".
*/

class AliEMCALCalibMapAPD : public TObject {

public:

  enum kValType {kCalibTemp=25};// 25 deg C used for all APD calibrations

  AliEMCALCalibMapAPD();

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadCalibMapAPDInfo(Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  void WriteCalibMapAPDInfo(const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  virtual ~AliEMCALCalibMapAPD();

  // values per single APD
  struct AliEMCALCalibMapAPDVal {

    Int_t fHardWareId; // HardWareIndex
    // info from APD calibrations
    Int_t fAPDNum;    // assigned APD-PA number; Catania 10000-, Houston: 20000-
    Float_t fV30;      // Catania/Houston Voltage V30 (V) at T = 25 deg C
    Float_t fPar[3];   // fit parameters, p0,p1,p2 - for ADC vs bias measurement
    Float_t fParErr[3]; // error on fit parameters	

    Int_t fBreakDown; // Hamamatsu Breakdown Voltage (V)	
    Float_t fDarkCurrent; // Hamamatsu Dark Current (A)	
  };

  // 1 SuperModule's worth of info: info on where the different APDs are
  struct AliEMCALSuperModuleCalibMapAPD {
    Int_t fSuperModuleNum;
    AliEMCALCalibMapAPDVal fAPDVal[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows];
  };

  // pointer to stored info.

  Int_t GetNSuperModule() const { return fNSuperModule; }; 
  AliEMCALSuperModuleCalibMapAPD * GetSuperModuleData() const { return fSuperModuleData; };

  // - via the index in the stored array:
  virtual AliEMCALSuperModuleCalibMapAPD GetSuperModuleCalibMapAPDId(Int_t smIndex) const;
  // - or via the actual SM number
  virtual AliEMCALSuperModuleCalibMapAPD GetSuperModuleCalibMapAPDNum(Int_t smNum) const;

  // method to calculate gain M from fit parameters, and HV value
  Float_t GetGain(Float_t fitPar[3], Float_t HV) const 
    { return (fitPar[0] + fitPar[1] * exp(fitPar[2]*HV)); };
  Float_t GetGain(Float_t fitPar0, Float_t fitPar1, Float_t fitPar2, Float_t HV) const 
    { return (fitPar0 + fitPar1 * exp(fitPar2*HV)); };

protected:

  Int_t 	  fNSuperModule; // Number of supermodules.
  AliEMCALSuperModuleCalibMapAPD *fSuperModuleData; // SuperModule data

private:

  AliEMCALCalibMapAPD(const AliEMCALCalibMapAPD &);
  AliEMCALCalibMapAPD &operator = (const AliEMCALCalibMapAPD &);

  ClassDef(AliEMCALCalibMapAPD, 1) //MapAPD data reader
};

#endif
