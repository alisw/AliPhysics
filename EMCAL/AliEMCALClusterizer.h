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

// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALClusterizer : public TTask {

public:

  AliEMCALClusterizer() ;        // default ctor
  AliEMCALClusterizer(const TString alirunFileName, const TString eventFolderName = AliConfig::fgkDefaultEventFolderName) ;
  virtual ~AliEMCALClusterizer() ; // dtor

  virtual Float_t GetTowerClusteringThreshold()const {Warning("GetTowerClusteringThreshold", "Not Defined") ; return 0. ; }
  virtual Float_t GetTowerLocalMaxCut()const {Warning("GetTowerLocalMaxCut", "Not Defined") ; return 0. ; }
  virtual Float_t GetTowerLogWeight()const {Warning("GetTowerLogWeight", "Not Defined") ; return 0. ; }
  virtual Float_t GetTimeGate() const {Warning("GetTimeGate", "Not Defined") ; return 0. ; }
  virtual const char *  GetRecPointsBranch() const {Warning("GetRecPointsBranch", "Not Defined") ; return 0 ; }
  virtual const Int_t GetRecPointsInRun()  const {Warning("GetRecPointsInRun", "Not Defined") ; return 0 ; }
  virtual const char *  GetDigitsBranch() const  {Warning("GetDigitsBranch", "Not Defined") ; return 0 ; }

  virtual void MakeClusters() const {Warning("MakeClusters", "Not Defined") ; }
  virtual void Print(Option_t * /*option*/)const {Warning("Print", "Not Defined") ; }

  virtual void SetECAClusteringThreshold(Float_t) = 0;
  virtual void SetECALocalMaxCut(Float_t)         = 0;
  virtual void SetECALogWeight(Float_t)           = 0;
  virtual void SetTimeGate(Float_t)               = 0;
  virtual void SetUnfolding(Bool_t)               = 0;
  virtual const char * Version() const {Warning("Version", "Not Defined") ; return 0 ; } 

protected:
  TString fEventFolderName ;  // event folder name

  ClassDef(AliEMCALClusterizer,4)  // Clusterization algorithm class 

} ;

#endif // AliEMCALCLUSTERIZER_H
