#ifndef ALIEMCALCALIBTIMEDEPCORRECTION_H
#define ALIEMCALCALIBTIMEDEPCORRECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
#include <TObjArray.h>
#include <TArrayF.h>
#include "AliEMCALGeoParams.h"
class TString;
class TTree;

/*
  Objects of this class contain info on time-dependent corrections
*/

// ******* internal class definition *************
// 1 SuperModule's worth of info
class AliEMCALSuperModuleCalibTimeDepCorrection : public TObject {
 public:
  AliEMCALSuperModuleCalibTimeDepCorrection(const int smNum=0) : TObject(), // just init values
    fSuperModuleNum(smNum)
    {
      for (int icol=0; icol<AliEMCALGeoParams::fgkEMCALCols; icol++) {
	for (int irow=0; irow<AliEMCALGeoParams::fgkEMCALRows; irow++) {
	  fCorrection[icol][irow].Reset();
	}
      }
    }

 public:
  void SetSuperModuleNum(Int_t i) { fSuperModuleNum = i;}; // 
  Int_t GetSuperModuleNum() const { return fSuperModuleNum;}; // 
  TArrayF * GetCorrection(int icol, int irow) 
    { return &fCorrection[icol][irow]; };

 private:
  Int_t fSuperModuleNum; // SM id
  TArrayF fCorrection[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; // values

  ClassDef(AliEMCALSuperModuleCalibTimeDepCorrection, 2) // help class
};
// ******* end of internal class definition *************

class AliEMCALCalibTimeDepCorrection : public TObject {
 public:
  AliEMCALCalibTimeDepCorrection(const int nSM = AliEMCALGeoParams::fgkEMCALModules);

  // interface methods; getting the whole struct should be more efficient though
  void InitCorrection(Int_t nSM, Int_t nBins, Float_t val=1.0); // assign a certain value to all 
  // use the methods below with caution: take care that your argument ranges are valid
  void SetCorrection(Int_t smIndex, Int_t iCol, Int_t iRow, Int_t iBin, Float_t val=1.0); // assign a certain value to a given bin
  Float_t GetCorrection(Int_t smIndex, Int_t iCol, Int_t iRow, Int_t iBin) const; // assign a certain value to a given bin
 
  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadTextInfo(Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void WriteTextInfo(const TString &txtFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void ReadRootInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void ReadTreeInfo(TTree *treeGlob, TTree *treeCorr, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules
  void WriteRootInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); // info file is for nSm=1 to 12 SuperModules

  virtual ~AliEMCALCalibTimeDepCorrection();

  // pointer to stored info.
  Int_t GetNSuperModule() const { return fNSuperModule; }; 
 
  // - via the index in the stored array:
  virtual AliEMCALSuperModuleCalibTimeDepCorrection * GetSuperModuleCalibTimeDepCorrectionId(Int_t smIndex) const
    { return (AliEMCALSuperModuleCalibTimeDepCorrection*) fSuperModuleData[smIndex]; };

  // - or via the actual SM number
  virtual AliEMCALSuperModuleCalibTimeDepCorrection * GetSuperModuleCalibTimeDepCorrectionNum(Int_t smNum) const;

  void SetStartTime(UInt_t ui) { fStartTime = ui; } //
  void SetNTimeBins(Int_t i) { fNTimeBins = i; } // 
  void SetTimeBinSize(Int_t i) { fTimeBinSize = i; } // 

  Int_t GetStartTime() const { return fStartTime; } //
  Int_t GetNTimeBins() const { return fNTimeBins; } // 
  Int_t GetTimeBinSize() const { return fTimeBinSize; } // 

  static Int_t GetMaxTimeBins() { return fgkMaxTimeBins; }

protected:

  Int_t 	  fNSuperModule; // Number of supermodules.
  TObjArray fSuperModuleData; // SuperModule data

private:

  AliEMCALCalibTimeDepCorrection(const AliEMCALCalibTimeDepCorrection &);
  AliEMCALCalibTimeDepCorrection &operator = (const AliEMCALCalibTimeDepCorrection &);

  UInt_t fStartTime; // timestamp for start of run/first bin
  Int_t fNTimeBins; // how many timestamp bins do we have
  Int_t fTimeBinSize; // seconds per time-bin

  static const Int_t fgkMaxTimeBins = 50; // we are not going to have more correction time bins than this for a single runnumber.. 

  ClassDef(AliEMCALCalibTimeDepCorrection, 3) //
};

#endif
