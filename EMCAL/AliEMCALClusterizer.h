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

//#include "AliEMCALDigit.h"

class AliEMCALClusterizer : public TTask {

public:

  AliEMCALClusterizer() ;        // default ctor
  AliEMCALClusterizer(const char * headerFile, const char * name, const Bool_t toSplit) ;
  virtual ~AliEMCALClusterizer() ; // dtor

  virtual Float_t GetTowerClusteringThreshold()const {Warning("GetTowerClusteringThreshold", "Not Defined") ; return 0. ; }
  virtual Float_t GetTowerLocalMaxCut()const {Warning("GetTowerLocalMaxCut", "Not Defined") ; return 0. ; }
  virtual Float_t GetTowerLogWeight()const {Warning("GetTowerLogWeight", "Not Defined") ; return 0. ; }
  virtual Float_t GetTimeGate() const {Warning("GetTimeGate", "Not Defined") ; return 0. ; }
  virtual Float_t GetPreShoClusteringThreshold()const {Warning("GetPreShoClusteringThreshold", "Not Defined") ; return 0. ; }
  virtual Float_t GetPreShoLocalMaxCut()const {Warning("GetPreShoLocalMaxCut", "Not Defined") ; return 0. ; }
  virtual Float_t GetPreShoLogWeight()const {Warning("GetPreShoLogWeight", "Not Defined") ; return 0. ; }
  virtual const char *  GetRecPointsBranch() const {Warning("GetRecPointsBranch", "Not Defined") ; return 0 ; }
  virtual const Int_t GetRecPointsInRun()  const {Warning("GetRecPointsInRun", "Not Defined") ; return 0 ; }
  virtual const char *  GetDigitsBranch() const  {Warning("GetDigitsBranch", "Not Defined") ; return 0 ; }

  virtual void MakeClusters() {Warning("MakeClusters", "Not Defined") ; }
  virtual void Print(Option_t * option)const {Warning("Print", "Not Defined") ; }

  virtual void SetECClusteringThreshold(Float_t cluth) {Warning("SetECClusteringThreshold", "Not Defined") ; }
  virtual void SetECLocalMaxCut(Float_t cut) {Warning("SetECLocalMaxCut", "Not Defined") ; }
  virtual void SetECLogWeight(Float_t w) {Warning("SetECLogWeight", "Not Defined") ; }
  virtual void SetHCClusteringThreshold(Float_t cluth) {Warning("SetHCLogWeight", "Not Defined") ; }
  virtual void SetHCLocalMaxCut(Float_t cut)            {Warning("SetHCLogWeight", "Not Defined") ; }
  virtual void SetHCLogWeight(Float_t w)                {Warning("SetHCLogWeight", "Not Defined") ; }
  virtual void SetTimeGate(Float_t gate) {Warning("SetTimeGate", "Not Defined") ; }
  virtual void SetPREClusteringThreshold(Float_t cluth) {Warning("SetPREClusteringThreshold", "Not Defined") ; }
  virtual void SetPRELocalMaxCut(Float_t cut) {Warning("SetPRELocalMaxCut", "Not Defined") ; }
  virtual void SetPRELogWeight(Float_t w) {Warning("SetPRELogWeight", "Not Defined") ; }
  virtual void SetDigitsBranch(const char * title) {Warning("SetDigitsBranch", "Not Defined") ; }
  virtual void SetRecPointsBranch(const char *title) {Warning("SetRecPointsBranch", "Not Defined") ; } 
  virtual void SetUnfolding(Bool_t toUnfold ) {Warning("SetUnfolding", "Not Defined") ; }
  virtual const char * Version() const {Warning("Version", "Not Defined") ; return 0 ; } 

protected:
  
  TFile * fSplitFile ;             //! file in which RecPoints will eventually be stored
  Bool_t  fToSplit ;               //! Should we write to splitted file

  ClassDef(AliEMCALClusterizer,2)  // Clusterization algorithm class 

} ;

#endif // AliEMCALCLUSTERIZER_H
