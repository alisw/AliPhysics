#ifndef ALIEMCALCALIBTIMEDEPCORRECTION_H
#define ALIEMCALCALIBTIMEDEPCORRECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
#include <TArrayF.h>
#include "AliEMCALGeoParams.h"
class TString;
class TArrayF;

/*
  Objects of this class read txt file with APD data
  AliEMCALCalibTimeDepCorrection inherits TObject only to use AliLog "functions".
*/

class AliEMCALCalibTimeDepCorrection : public TObject {
public:
  AliEMCALCalibTimeDepCorrection();

  // interface methods; getting the whole struct should be more efficient though
  void InitCorrection(Int_t nSM, Int_t nBins, Float_t val=1.0); // assign a certain value to all 
  // use the methods below with caution: take care that your argument ranges are valid
  void SetCorrection(Int_t smIndex, Int_t iCol, Int_t iRow, Int_t iBin, Float_t val=1.0); // assign a certain value to a given bin
  Float_t GetCorrection(Int_t smIndex, Int_t iCol, Int_t iRow, Int_t iBin) const; // assign a certain value to a given bin
 

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadCalibTimeDepCorrectionInfo(Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  void WriteCalibTimeDepCorrectionInfo(const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  virtual ~AliEMCALCalibTimeDepCorrection();

  struct AliEMCALSuperModuleCalibTimeDepCorrection {
    Int_t fSuperModuleNum;
    TArrayF fCorrection[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 
  };

  // pointer to stored info.
  Int_t GetNSuperModule() const { return fNSuperModule; }; 
  AliEMCALSuperModuleCalibTimeDepCorrection * GetSuperModuleData() const { return fSuperModuleData; };

  // - via the index in the stored array:
  virtual AliEMCALSuperModuleCalibTimeDepCorrection GetSuperModuleCalibTimeDepCorrectionId(Int_t smIndex) const;
  // - or via the actual SM number
  virtual AliEMCALSuperModuleCalibTimeDepCorrection GetSuperModuleCalibTimeDepCorrectionNum(Int_t smNum) const;

  void SetStartTime(UInt_t ui) { fStartTime = ui; } //
  void SetNTimeBins(Int_t i) { fNTimeBins = i; } // 
  void SetTimeBinSize(Int_t i) { fTimeBinSize = i; } // 

  Int_t GetStartTime() const { return fStartTime; } //
  Int_t GetNTimeBins() const { return fNTimeBins; } // 
  Int_t GetTimeBinSize() const { return fTimeBinSize; } // 

protected:

  Int_t 	  fNSuperModule; // Number of supermodules.
  AliEMCALSuperModuleCalibTimeDepCorrection *fSuperModuleData; // SuperModule data

private:

  AliEMCALCalibTimeDepCorrection(const AliEMCALCalibTimeDepCorrection &);
  AliEMCALCalibTimeDepCorrection &operator = (const AliEMCALCalibTimeDepCorrection &);

  UInt_t fStartTime; // timestamp for start of run/first bin
  Int_t fNTimeBins; // how many timestamp bins do we have
  Int_t fTimeBinSize; // seconds per time-bin

  ClassDef(AliEMCALCalibTimeDepCorrection, 2) //
};

#endif
