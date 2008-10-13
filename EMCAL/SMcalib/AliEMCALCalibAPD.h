#ifndef ALIEMCALCALIBAPD_H
#define ALIEMCALCALIBAPD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
class TString;

/*
  Objects of this class read txt file with APD data
  AliEMCALCalibAPD inherits TObject only to use AliLog "functions".
*/

class AliEMCALCalibAPD : public TObject {
public:
  AliEMCALCalibAPD();

  void ReadCalibAPDInfo(Int_t nAPD, const TString &txtFileName); // info file is for nSm=1 to 12 SuperModules
  void WriteCalibAPDInfo(const TString &txtFileName); // info file is for nSm=1 to 12 SuperModules

  virtual ~AliEMCALCalibAPD();

  struct AliEMCALCalibAPDData {
    Int_t fAPDNum;    // assigned APD-PA number; Catania 10000-, Houston: 20000-
    UInt_t fSerialNum; // Serial Number; from Hamamatsu	
    Char_t fStatus[10]; // Status info: should be "tested"
    Char_t fLocation[20]; // where was the test done: "Catania" or "Houston"
    Int_t fRunNum; // DATE run at test station
    Int_t fTestPos; // location of APD during test

    Float_t fV30;      // Catania/Houston Voltage V30 (V) at T = 25 deg C
    Float_t fV50;      // Catania/Houston Voltage V50 (V) at T = 25 deg C
    Float_t fVoltCoeff; // 1/M x dM/dV
    Float_t fPar[3];   // fit parameters, p0,p1,p2 - for ADC vs bias measurement
    Float_t fParErr[3]; // error on fit parameters	

    Int_t fBreakDown; // Hamamatsu Breakdown Voltage (V)	
    Float_t fHamV50;       // Hamamatsu Voltage V50 (V)	
    Float_t fDarkCurrent; // Hamamatsu Dark Current (A)	
    Float_t fTestTemp; // Hamamatsu Testing Temperature (deg C)	

  };

  // pointer to stored info.
  Int_t GetNCalibAPD() const { return fNCalibAPD; }; 
  AliEMCALCalibAPDData * GetCalibAPDData() const { return fData; };

  // - via the index in the stored array:
  virtual AliEMCALCalibAPDData GetCalibAPDDataId(Int_t apdIndex) const;
  // - or via the actual APD number
  virtual AliEMCALCalibAPDData GetCalibAPDDataNum(Int_t apdNum) const;

protected:

  Int_t 	  fNCalibAPD; // Number of APDs
  AliEMCALCalibAPDData *fData; // array with the data

private:

  AliEMCALCalibAPD(const AliEMCALCalibAPD &);
  AliEMCALCalibAPD &operator = (const AliEMCALCalibAPD &);

  ClassDef(AliEMCALCalibAPD, 1) //CalibAPD data reader
};

#endif
