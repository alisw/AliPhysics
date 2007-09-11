#ifndef ALIEMCALSUPERMODULE_H
#define ALIEMCALSUPERMODULE_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Emcal Super Module     
//                  
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA)
//  Super Module folder
//  Initial version was created with TDataSet staf
//  TObjectSet -> TFolder; Sep 6, 2007


#include <TFolder.h>

class TList;
class AliEMCALCell;

class AliEMCALSuperModule : public TFolder {

 public:
  
  AliEMCALSuperModule(); 
  AliEMCALSuperModule(const Int_t m, const char* title="Emcal Super Module");

  virtual ~AliEMCALSuperModule();

  void Init();
  void   AddCellToEtaRow(AliEMCALCell *cell, const Int_t etaRow);
  TList*   GetHists()  {return fLh;}
  TObject* GetParent() {return fParent;}
  void     SetParent(TObject *parent) {fParent=parent;}
  // MENU
  void FitForAllCells(); //*MENU*  
  void FitEffMassHist(); //*MENU*  
  void PrintInfo();      //*MENU* 
  void DrawCC(int iopt=1); //*MENU* 
  //
  Int_t GetNumberOfCells();
  Int_t GetSMNumber() const {return fSMNumber;}
 protected:
  TList* BookHists();
  //
  TObject* fParent; // parent
  TList*   fLh;     // List of hists
  Int_t    fSMNumber;

  ClassDef(AliEMCALSuperModule,2) // EMCAL SuperModule
    
};

#endif // ALIEMCALSUPERMODULE_H
