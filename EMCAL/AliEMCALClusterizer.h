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
#include <Riostream.h> 

// --- AliRoot header files ---

//#include "AliEMCALDigit.h"

class AliEMCALClusterizer : public TTask {

public:

  AliEMCALClusterizer() ;        // default ctor
  AliEMCALClusterizer(const char * headerFile, const char * name, const Bool_t toSplit) ;
  virtual ~AliEMCALClusterizer() ; // dtor

  virtual Float_t GetTowerClusteringThreshold()const {cout << "Not Defined" << endl ; return 0. ; }
  virtual Float_t GetTowerLocalMaxCut()const {cout << "Not Defined" << endl ; return 0. ; }
  virtual Float_t GetTowerLogWeight()const {cout << "Not Defined" << endl ; return 0. ; }
  virtual Float_t GetTimeGate() const {cout << "Not Defined" << endl ; return 0. ; }
  virtual Float_t GetPreShoClusteringThreshold()const {cout << "Not Defined" << endl ; return 0. ; }
  virtual Float_t GetPreShoLocalMaxCut()const {cout << "Not Defined" << endl ; return 0. ; }
  virtual Float_t GetPreShoLogWeight()const {cout << "Not Defined" << endl ; return 0. ; }
  virtual const char *  GetRecPointsBranch() const {cout << "Not Defined" << endl ; return 0 ; }
  virtual const Int_t GetRecPointsInRun()  const {cout << "Not Defined" << endl ; return 0 ; }
  virtual const char *  GetDigitsBranch() const  {cout << "Not Defined" << endl ; return 0 ; }

  virtual void MakeClusters() {cout << "Not Defined" << endl ; }
  virtual void Print(Option_t * option)const {cout << "Not Defined" << endl ; }

  virtual void SetTowerClusteringThreshold(Float_t cluth) {cout << "Not Defined" << endl ; }
  virtual void SetTowerLocalMaxCut(Float_t cut) {cout << "Not Defined" << endl ; }
  virtual void SetTowerLogWeight(Float_t w) {cout << "Not Defined" << endl ; }
  virtual void SetTimeGate(Float_t gate) {cout << "Not Defined" << endl ; }
  virtual void SetPreShoClusteringThreshold(Float_t cluth) {cout << "Not Defined" << endl ; }
  virtual void SetPreShoLocalMaxCut(Float_t cut) {cout << "Not Defined" << endl ; }
  virtual void SetPreShoLogWeight(Float_t w) {cout << "Not Defined" << endl ; }
  virtual void SetDigitsBranch(const char * title) {cout << "Not Defined" << endl ; }
  virtual void SetRecPointsBranch(const char *title) {cout << "Not Defined" << endl ; } 
  virtual void SetUnfolding(Bool_t toUnfold ) {cout << "Not Defined" << endl ; }
  virtual const char * Version() const {cout << "Not Defined" << endl ; return 0 ; } 

protected:
  
  TFile * fSplitFile ;             //! file in which RecPoints will eventually be stored
  Bool_t  fToSplit ;               //! Should we write to splitted file

  ClassDef(AliEMCALClusterizer,2)  // Clusterization algorithm class 

} ;

#endif // AliEMCALCLUSTERIZER_H
