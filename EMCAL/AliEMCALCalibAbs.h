#ifndef ALIEMCALCALIBABS_H
#define ALIEMCALCALIBABS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
#include <TObjArray.h>
#include "AliEMCALGeoParams.h"
class TString;
class TTree;

/*
  Objects of this class contain basis for absolute calibrations
*/

// total calibration factor is a product of
// a) overall calibration factor [fAbsoluteCalib]
// b) individual gain factor per tower [fRelativeCalib]
// c) time-dependent correction
// In this class we store a), b) 

// ******* internal class definition *************

// 1 SuperModule's worth of info: info on where the different APDs are
class AliEMCALSuperModuleCalibAbs : public TObject {

 public:
  AliEMCALSuperModuleCalibAbs(const int smNum=0) : TObject(), // just init values
    fSuperModuleNum(smNum),
    fCalibMethod(0),
    fCalibPass(0),
    fAbsoluteCalib(0)
    {
      for (int icol=0; icol<AliEMCALGeoParams::fgkEMCALCols; icol++) {
	for (int irow=0; irow<AliEMCALGeoParams::fgkEMCALRows; irow++) {
	  fRelativeCalib[icol][irow] = 1.0;
	}
      }
    }

 public:
  // first
  void SetSuperModuleNum(Int_t i) { fSuperModuleNum = i;}; // 
  Int_t GetSuperModuleNum() const { return fSuperModuleNum;}; // 
  void SetCalibMethod(Int_t i) { fCalibMethod = i;}; // 
  Int_t GetCalibMethod() const { return fCalibMethod;}; // 
  void SetCalibPass(Int_t i) { fCalibPass = i;}; // 
  Int_t GetCalibPass() const { return fCalibPass;}; // 
  void SetAbsoluteCalib(Float_t f) { fAbsoluteCalib = f;}; // 
  Float_t GetAbsoluteCalib() const { return fAbsoluteCalib;}; // 

  // third
  void SetRelativeCalib(int icol, int irow, Float_t f) { fRelativeCalib[icol][irow] = f; }; //
  Float_t GetRelativeCalib(int icol, int irow) const { return fRelativeCalib[icol][irow]; }; //

 private:
  // first: overall values for the whole SuperModule
  Int_t fSuperModuleNum; // which SuperModule is this?
  Int_t fCalibMethod; // a la 0=cosmics, 1=pi0, 2=electrons,3=using ecore,
  Int_t fCalibPass; // which analysis iteration is this.. 1,2,..N
  Float_t fAbsoluteCalib; // (ADC>GeV absolute gain/conversion)
  
  // third: individual info for each tower
  Float_t fRelativeCalib[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; // values around 1, if gains are well balanced

  ClassDef(AliEMCALSuperModuleCalibAbs, 3) // help class
};
// ******* end of internal class definition *************

class AliEMCALCalibAbs : public TObject {

public:

  enum kProblemType {kNoLED=-999};// code in possible problems

  AliEMCALCalibAbs(const int nSM = AliEMCALGeoParams::fgkEMCALModules);

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

  // - via the index in the stored array:
  virtual AliEMCALSuperModuleCalibAbs * GetSuperModuleCalibAbsId(Int_t smIndex) const
   { return (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[smIndex]; };

  // - or via the actual SM number
  virtual AliEMCALSuperModuleCalibAbs * GetSuperModuleCalibAbsNum(Int_t smNum) const;

protected:

  Int_t 	  fNSuperModule; // Number of supermodules.
  TObjArray fSuperModuleData; // SuperModule data

private:

  AliEMCALCalibAbs(const AliEMCALCalibAbs &);
  AliEMCALCalibAbs &operator = (const AliEMCALCalibAbs &);

  ClassDef(AliEMCALCalibAbs, 4) //CalibAbs data info
};

#endif
