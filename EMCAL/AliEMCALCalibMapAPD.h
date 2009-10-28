#ifndef ALIEMCALCALIBMAPAPD_H
#define ALIEMCALCALIBMAPAPD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
#include <TObjArray.h>
#include "AliEMCALGeoParams.h"
#include <cmath>
class TString;
class TTree;

/*
  Objects of this class contain info on APD calibration and map info,
  such as V30 and other parameters from the tests in Catania/Houston,
  as well as info on which APD is located where.  
*/

// ******* internal class definition *************
// values per single APD
class AliEMCALCalibMapAPDVal : public TObject {

 public:
  AliEMCALCalibMapAPDVal() : TObject(), // just init values
    fHardWareId(0),
    fAPDNum(0),
    fV30(0),
    fBreakDown(0),
    fDarkCurrent(0) 
    {
      Init();
    }
  
  void Init() {
    fHardWareId = 0;
    fAPDNum = 0;
    fV30 = 0;
    fBreakDown = 0;
    fDarkCurrent = 0; 
    for (int i=0; i<3; i++) {
      fPar[i] = 0;
      fParErr[i] = 0;
    }
    return;
  }
  
 public:
  void SetHardWareId(Int_t i) { fHardWareId = i; }; //
  Int_t GetHardWareId() const { return fHardWareId; }; //
  void SetAPDNum(Int_t i) { fAPDNum = i; }; //
  Int_t GetAPDNum() const { return fAPDNum; }; //
  void SetV30(Float_t f) { fV30 = f; }; //
  Float_t GetV30() const { return fV30; }; // 
  void SetPar(int ip, Float_t f) { fPar[ip] = f; }; //
  Float_t GetPar(int ip) const { return fPar[ip]; }; // 
  void SetParErr(int ip, Float_t f) { fParErr[ip] = f; }; //
  Float_t GetParErr(int ip) const { return fParErr[ip]; }; // 
  void SetBreakDown(Int_t i) { fBreakDown = i; }; //
  Int_t GetBreakDown() const { return fBreakDown; }; //
  void SetDarkCurrent(Float_t f) { fDarkCurrent = f; }; //
  Float_t GetDarkCurrent() const { return fDarkCurrent; }; // 

 private:
  Int_t fHardWareId; // HardWareIndex
  // info from APD calibrations
  Int_t fAPDNum;    // assigned APD-PA number; Catania 10000-, Houston: 20000-
  Float_t fV30;      // Catania/Houston Voltage V30 (V) at T = 25 deg C
  Float_t fPar[3];   // fit parameters, p0,p1,p2 - for ADC vs bias measurement
  Float_t fParErr[3]; // error on fit parameters	
  
  Int_t fBreakDown; // Hamamatsu Breakdown Voltage (V)	
  Float_t fDarkCurrent; // Hamamatsu Dark Current (A)	
  
  ClassDef(AliEMCALCalibMapAPDVal, 2) // help class
}; // AliEMCALCalibAPDVal

// 1 SuperModule's worth of info: info on where the different APDs are
class AliEMCALSuperModuleCalibMapAPD : public TObject {
 public:
  AliEMCALSuperModuleCalibMapAPD(const int smNum=0) : TObject(), // just init values
    fSuperModuleNum(smNum)
    {
      for (int icol=0; icol<AliEMCALGeoParams::fgkEMCALCols; icol++) {
	for (int irow=0; irow<AliEMCALGeoParams::fgkEMCALRows; irow++) {
	  fAPDVal[icol][irow].Init();
	}
      }
    }
  
 public:
  void SetSuperModuleNum(Int_t i) { fSuperModuleNum = i;}; // 
  Int_t GetSuperModuleNum() const { return fSuperModuleNum;}; // 
  AliEMCALCalibMapAPDVal * GetAPDVal(int icol, int irow) 
    { return &fAPDVal[icol][irow]; };

 private:
  Int_t fSuperModuleNum; // SuperModule index
  AliEMCALCalibMapAPDVal fAPDVal[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; // APD calibration info
  
  ClassDef(AliEMCALSuperModuleCalibMapAPD, 2) // help class
};
// ******* end of internal class definition *************    
    
class AliEMCALCalibMapAPD : public TObject {

public:

  enum kValType {kCalibTemp=25};// 25 deg C used for all APD calibrations

  AliEMCALCalibMapAPD(const int nSM = AliEMCALGeoParams::fgkEMCALModules);

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadTextCalibMapAPDInfo(Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void WriteTextCalibMapAPDInfo(const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void ReadRootCalibMapAPDInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void ReadTreeCalibMapAPDInfo(TTree *tree, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void WriteRootCalibMapAPDInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  virtual ~AliEMCALCalibMapAPD();

  // pointer to stored info.
  Int_t GetNSuperModule() const { return fNSuperModule; }; 

  // on a SuperModule level
  virtual AliEMCALSuperModuleCalibMapAPD * GetSuperModuleCalibMapAPDId(Int_t smIndex) const
    { return (AliEMCALSuperModuleCalibMapAPD*) fSuperModuleData[smIndex]; }; // - via the index in the stored array:

  virtual AliEMCALSuperModuleCalibMapAPD * GetSuperModuleCalibMapAPDNum(Int_t smNum) const;   // - or via the actual SM number
  
  // method to calculate gain M from fit parameters, and HV value
  Float_t GetGain(Float_t fitPar[3], Float_t HV) const 
    { return (fitPar[0] + fitPar[1] * exp(fitPar[2]*HV)); };
  Float_t GetGain(Float_t fitPar0, Float_t fitPar1, Float_t fitPar2, Float_t HV) const 
    { return (fitPar0 + fitPar1 * exp(fitPar2*HV)); };

protected:

  Int_t 	  fNSuperModule; // Number of supermodules.
  TObjArray fSuperModuleData; // SuperModule data

private:

  AliEMCALCalibMapAPD(const AliEMCALCalibMapAPD &);
  AliEMCALCalibMapAPD &operator = (const AliEMCALCalibMapAPD &);

  ClassDef(AliEMCALCalibMapAPD, 3) //CalibMapAPD data info
};

#endif
