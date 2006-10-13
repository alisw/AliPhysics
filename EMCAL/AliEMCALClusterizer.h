#ifndef ALIEMCALCLUSTERIZER_H
#define ALIEMCALCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
// --- ROOT system ---

#include "TTask.h" 
#include "AliConfig.h"

// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALClusterizer : public TTask {

public:

  AliEMCALClusterizer() ;        // default ctor
  AliEMCALClusterizer(const TString alirunFileName, const TString eventFolderName = AliConfig::GetDefaultEventFolderName()) ;
  AliEMCALClusterizer(const AliEMCALClusterizer &); //copy ctor
  virtual ~AliEMCALClusterizer() ; // dtorEM

  virtual Float_t GetTowerClusteringThreshold()const {Warning("GetTowerClusteringThreshold", "Not Defined") ; return 0. ; }
  virtual Float_t GetTowerLocalMaxCut()const {Warning("GetTowerLocalMaxCut", "Not Defined") ; return 0. ; }
  virtual Float_t GetTowerLogWeight()const {Warning("GetTowerLogWeight", "Not Defined") ; return 0. ; }
  virtual Float_t GetTimeGate() const {Warning("GetTimeGate", "Not Defined") ; return 0. ; }
  virtual const char *  GetRecPointsBranch() const {Warning("GetRecPointsBranch", "Not Defined") ; return 0 ; }
  virtual Int_t GetRecPointsInRun()  const {Warning("GetRecPointsInRun", "Not Defined") ; return 0 ; }
  virtual const char *  GetDigitsBranch() const  {Warning("GetDigitsBranch", "Not Defined") ; return 0 ; }

  virtual void MakeClusters() = 0;

  virtual void SetECAClusteringThreshold(Float_t) = 0;
  virtual void SetECALocalMaxCut(Float_t)         = 0;
  virtual void SetECALogWeight(Float_t)           = 0;
  virtual void SetTimeGate(Float_t)               = 0;
  virtual void SetUnfolding(Bool_t)               = 0;
  void SetEventRange(Int_t first=0, Int_t last=-1) {fFirstEvent=first; fLastEvent=last; }
  void SetEventFolderName(TString name) { fEventFolderName = name ; }

  AliEMCALClusterizer & operator = (const AliEMCALClusterizer & /*rvalue*/)  {return *this ;} 

  virtual const char * Version() const {Warning("Version", "Not Defined") ; return 0 ; } 

protected:
  TString fEventFolderName ;  // event folder name
  Int_t   fFirstEvent;        // first event to process
  Int_t   fLastEvent;         // last  event to process

  ClassDef(AliEMCALClusterizer,4)  // Clusterization algorithm class 

} ;

#endif // AliEMCALCLUSTERIZER_H
