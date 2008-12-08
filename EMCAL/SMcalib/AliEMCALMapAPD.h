#ifndef ALIEMCALMAPAPD_H
#define ALIEMCALMAPAPD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

#include <TObject.h>
class TString;

static const int fgkEmCalRows = 24; // number of rows per module for EMCAL
static const int fgkEmCalCols = 48; // number of columns per module for EMCAL

/*
  Objects of this class read txt file with APD data
  AliEMCALMapAPD inherits TObject only to use AliLog "functions".
*/

class AliEMCALMapAPD : public TObject {
public:
  AliEMCALMapAPD();

  void ReadMapAPDInfoStripBasis(Int_t nSM, const TString &txtFileName); // info file is for nSm=1 to 12 SuperModules
  void ReadMapAPDInfoSingleStripBasis(Int_t iSM, Int_t iStrip, const TString &txtFileName); // info file is for one single SuperModule and StripModule
 
  void ReadMapAPDInfo(Int_t nSM, const TString &txtFileName); // info file is for nSm=1 to 12 SuperModules

  void WriteMapAPDInfo(const TString &txtFileName); // info file is for nSm=1 to 12 SuperModules

  void GenerateDummyAPDInfo(Int_t nSM, Int_t * iSM); // for debug purposes 

  int CheckForDuplicates(); // see if the same APD numbers occur more than once

  virtual ~AliEMCALMapAPD();

  struct AliEMCALSuperModuleMapAPD {
    Int_t fSuperModuleNum;
    Int_t fAPDNum[fgkEmCalCols][fgkEmCalRows];
  };

  // pointer to stored info.
  Int_t GetNSuperModule() const { return fNSuperModule; }; 
  AliEMCALSuperModuleMapAPD * GetSuperModuleData() const { return fSuperModuleData; };

  // - via the index in the stored array:
  virtual AliEMCALSuperModuleMapAPD GetSuperModuleMapAPDId(Int_t smIndex) const;
  // - or via the actual SM number
  virtual AliEMCALSuperModuleMapAPD GetSuperModuleMapAPDNum(Int_t smNum) const;

protected:

  Int_t 	  fNSuperModule; // Number of supermodules.
  AliEMCALSuperModuleMapAPD *fSuperModuleData; // SuperModule data

private:

  AliEMCALMapAPD(const AliEMCALMapAPD &);
  AliEMCALMapAPD &operator = (const AliEMCALMapAPD &);

  ClassDef(AliEMCALMapAPD, 1) //MapAPD data reader
};

#endif
