#ifndef ALIEMCALSUPERMODULE_H
#define ALIEMCALSUPERMODULE_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Emcal Super Module     
//                  
//*-- Author: Aleksei Pavlinov (WSU, Detroit, USA) 

// --- ROOT system ---

#include <TObjectSet.h>

class TList;
class AliEMCALCell;

class AliEMCALSuperModule : public TObjectSet {

 public:
  
  AliEMCALSuperModule(); 
  AliEMCALSuperModule(const Int_t m, const char* title="Emcal Super Module");

  virtual ~AliEMCALSuperModule();

  void Init();
  void   AddCellToEtaRow(AliEMCALCell *cell, const Int_t etaRow);
  TList* GetHists() {return (TList*)fObj;}
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
  Int_t  fSMNumber;

  ClassDef(AliEMCALSuperModule,1) // EMCAL SuperModule
    
};

#endif // ALIEMCALSUPERMODULE_H
