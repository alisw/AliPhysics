#ifndef ALIEMCALCLUSTERIZER_H
#define ALIEMCALCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
// --- ROOT system ---

#include "TTask.h" 
class TFile ; 
// --- Standard library ---

// --- AliRoot header files ---

//#include "AliEMCALDigit.h"

class AliEMCALClusterizer : public TTask {

public:

  AliEMCALClusterizer() ;        // default ctor
  AliEMCALClusterizer(const char * headerFile, const char * name) ;
  virtual ~AliEMCALClusterizer() ; // dtor

  const TString GetHitsFileName() const { return fHitsFileName ; }
  const TString GetSDigitsFileName() const { return fSDigitsFileName ; }
  const TString GetDigitsFileName() const { return fDigitsFileName ; }

  virtual Float_t GetEmcClusteringThreshold()const = 0 ; 
  virtual Float_t GetEmcLocalMaxCut()const = 0 ; 
  virtual Float_t GetEmcLogWeight()const = 0 ; 
  virtual Float_t GetTimeGate() const = 0 ;
  virtual Float_t GetCpvClusteringThreshold()const = 0 ; 
  virtual Float_t GetCpvLocalMaxCut()const = 0 ; 
  virtual Float_t GetCpvLogWeight()const = 0 ; 
  virtual char *  GetRecPointsBranch() const = 0 ;
  virtual const Int_t GetRecPointsInRun()  const = 0 ; 
  virtual char *  GetDigitsBranch() const = 0 ;

  virtual void MakeClusters() = 0 ; 
  virtual void Print(Option_t * option)const = 0;

  virtual void SetTowerClusteringThreshold(Float_t cluth) = 0 ; 
  virtual void SetTowerLocalMaxCut(Float_t cut) = 0 ; 
  virtual void SetTowerLogWeight(Float_t w) = 0 ; 
  virtual void SetTimeGate(Float_t gate) = 0 ;
  virtual void SetPreShoClusteringThreshold(Float_t cluth) = 0 ; 
  virtual void SetPreShoLocalMaxCut(Float_t cut) = 0 ; 
  virtual void SetPreShoLogWeight(Float_t w) = 0 ; 
  virtual void SetDigitsBranch(const char * title) = 0 ;
  virtual void SetRecPointsBranch(const char *title) = 0 ;
  void SetSplitFile(const TString splitFileName = "EMCAL.RecPoints.root") ;
  virtual void SetUnfolding(Bool_t toUnfold ) = 0 ;
  virtual const char * Version() const = 0 ;  

protected:
  
  TString fHitsFileName ;          // file name that contains the original hits
  TString fSDigitsFileName ;       // file name that contains the original SDigits
  TString fDigitsFileName ;        // file name that contains the original Digits
  TFile * fSplitFile ;             //! file in which RecPoints will eventually be stored

  ClassDef(AliEMCALClusterizer,1)  // Clusterization algorithm class 

} ;

#endif // AliEMCALCLUSTERIZER_H
