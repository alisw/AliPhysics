#ifndef ALIEMCALCALIBTEMPCOEFF_H
#define ALIEMCALCALIBTEMPCOEFF_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
#include <TObjArray.h>
#include "AliEMCALGeoParams.h"
class TString;
class TTree;

/*
  Objects of this class contain temperature-dependence coefficients
*/

// total calibration factor is a product of
// a) overall calibration factor [fAbsoluteCalib]
// b) individual gain factor per tower [fRelativeCalib]
// c) time-dependent correction
// In this class we store factors needed for c)

// ******* internal class definition *************

// 1 SuperModule's worth of info: info on where the different APDs are
class AliEMCALSuperModuleCalibTempCoeff : public TObject {

 public:
  AliEMCALSuperModuleCalibTempCoeff(const int smNum=0) : TObject(), // just init values
    fSuperModuleNum(smNum)
    {
      for (int icol=0; icol<AliEMCALGeoParams::fgkEMCALCols; icol++) {
	for (int irow=0; irow<AliEMCALGeoParams::fgkEMCALRows; irow++) {
	  fTC[icol][irow] = 1.0;
	  fSrc[icol][irow] = 0;
	}
      }
    }

 public:

  void SetSuperModuleNum(Int_t i) { fSuperModuleNum = i;}; // 
  Int_t GetSuperModuleNum() const { return fSuperModuleNum;}; // 
  void SetTC(int icol, int irow, Float_t f) { fTC[icol][irow] = f; }; //
  Float_t GetTC(int icol, int irow) const { return fTC[icol][irow]; }; //
  void SetSrc(int icol, int irow, Int_t i) { fSrc[icol][irow] = i; }; //
  Int_t GetSrc(int icol, int irow) const { return fSrc[icol][irow]; }; //

 private:
  Int_t fSuperModuleNum; // which SuperModule is this?
  // individual info for each tower
  Float_t fTC[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; // Temperature Coefficient values (nominally around 2% change per deg C)
  Int_t fSrc[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; // from DarkCurrent, or LED or ..

  ClassDef(AliEMCALSuperModuleCalibTempCoeff, 1) // help class
};
// ******* end of internal class definition *************

class AliEMCALCalibTempCoeff : public TObject {

public:

  enum kSrcType {kDarkCurrent=0, kLED=1};// code in possible sources

  AliEMCALCalibTempCoeff(const int nSM = AliEMCALGeoParams::fgkEMCALModules);

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadTextCalibTempCoeffInfo(Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void WriteTextCalibTempCoeffInfo(const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void ReadRootCalibTempCoeffInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void ReadTreeCalibTempCoeffInfo(TTree *tree, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void WriteRootCalibTempCoeffInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  virtual ~AliEMCALCalibTempCoeff();

  // pointer to stored info.
  Int_t GetNSuperModule() const { return fNSuperModule; }; 

  // - via the index in the stored array:
  virtual AliEMCALSuperModuleCalibTempCoeff * GetSuperModuleCalibTempCoeffId(Int_t smIndex) const
   { return (AliEMCALSuperModuleCalibTempCoeff*) fSuperModuleData[smIndex]; };

  // - or via the actual SM number
  virtual AliEMCALSuperModuleCalibTempCoeff * GetSuperModuleCalibTempCoeffNum(Int_t smNum) const;

protected:

  Int_t 	  fNSuperModule; // Number of supermodules.
  TObjArray fSuperModuleData; // SuperModule data

private:

  AliEMCALCalibTempCoeff(const AliEMCALCalibTempCoeff &);
  AliEMCALCalibTempCoeff &operator = (const AliEMCALCalibTempCoeff &);

  ClassDef(AliEMCALCalibTempCoeff, 1) //CalibTempCoeff data info
};

#endif
