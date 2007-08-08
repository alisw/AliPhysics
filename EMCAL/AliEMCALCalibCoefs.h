#ifndef ALIEMCALCALIBCOEFS_H
#define ALIEMCALCALIBCOEFS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//    Table of Calibration coefficients  
//                  
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

// --- ROOT system ---
#include <TTable.h>

// unit is GeV
struct  calibCoef {
  Int_t    absId; // absolute id of cell 
  Double_t cc;    // Calib. coef
  Double_t eCc;   // Calib. coef. error
};

class TH1F;

class AliEMCALCalibCoefs : public TTable {
 public:
  enum EEmcalCalibType {kMC, kEQUALIZATION, kMIP, kPI0}; // type of EMCAL calibrations 

  void  SetCalibMethod(Int_t var) {fCalibMethod=var;}
  Int_t GetCalibMethod() {return fCalibMethod;}
  calibCoef* GetRow(const int absId);
  // Get initial Calib Data from DB
  static AliEMCALCalibCoefs *GetCalibTableFromDb(const char *tn="CCIN");
  static TH1F *GetHistOfCalibTableFromDb(const char *tn="CCTMP");
  // Menu
  void PrintTable();                 // *MENU*
  void PrintTable(const Int_t i);    // *MENU*
  void PrintRec(calibCoef *r);

 protected:
  Int_t fCalibMethod;  // method of calibration - EEmcalCalibType

  ClassDefTable(AliEMCALCalibCoefs , calibCoef)
  ClassDef(AliEMCALCalibCoefs,1) // Table of Calibration coefficients  
};

#endif // ALIEMCALCalibCoefs_H
