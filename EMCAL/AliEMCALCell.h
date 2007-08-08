#ifndef ALIEMCALCELL_H
#define ALIEMCALCELL_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  EMCAL cell - keep everyrhing for calibration task     
//                  
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

#include <TObjectSet.h>

class TList;
class TH1;
class TF1;
class TNtuple;

class AliEMCALCalibData;
class AliEMCALCalibCoefs;

class AliEMCALCell : public TObjectSet {

 public:
  
  AliEMCALCell(); 
  AliEMCALCell(const Int_t absId, const char* title="EMCAL cell");

  virtual ~AliEMCALCell();

  void SetCCfromDB(AliEMCALCalibData *ccDb);  // obsolete
  void SetCCfromCCTable(AliEMCALCalibCoefs *t);

  TList*   GetHists() {return (TList*)fObj;}
  Int_t    GetAbsId()  const {return fAbsId;}
  Int_t    GetSupMod() const {return fSupMod;}
  Int_t    GetModule() const {return fModule;}

  Double_t GetCcIn()  {return  fCcIn;}
  Double_t GetCcOut() {return  fCcOut;}
  TF1*     GetFunction() {return fFun;}

  void FillEffMass(const Double_t mgg);
  void FillCellNtuple(TNtuple *nt);

  static void FitHist(TH1* h, const char* name="",const char* opt="");
  // Menu
  void FitEffMassHist(const char* opt=""); //*MENU*
  void Print();                            //*MENU*
 protected:
  Int_t fAbsId;   // abs cell id 
  Int_t fSupMod;  // super module number
  Int_t fModule;  // module number inside SM
  Int_t fPhi;     // phi number of cell inside module  
  Int_t fEta;     // eta number of cell inside module  
  Int_t fPhiCell; // phi number of cell SM  
  Int_t fEtaCell; // eta number of cell SM  
  // CC staf
  Double_t fCcIn;  // input  cc in GeV (from Db or table
  Double_t fCcOut; // output cc  in GeV (from fir now)

  TF1*   fFun;     //! fitting function - gaus + pol2
  //
  TList* BookHists();

  ClassDef(AliEMCALCell,1) // EMCAL cell
    
};

#endif // ALIEMCALCELL_H
