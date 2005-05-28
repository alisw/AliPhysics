#ifndef ALIPHOSCLUSTERIZER_H
#define ALIPHOSCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

//_________________________________________________________________________
//  Base class for the clusterization algorithm 
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
// --- ROOT system ---

#include "TTask.h" 
#include "AliConfig.h"

class TFile ; 

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSClusterizer : public TTask {

public:

  AliPHOSClusterizer() ;        // default ctor
  AliPHOSClusterizer(const TString alirunFileName, const TString eventFolderName = AliConfig::GetDefaultEventFolderName()) ;
  AliPHOSClusterizer(const AliPHOSClusterizer & clusterizer) : TTask(clusterizer) { ; }
  virtual ~AliPHOSClusterizer() ; // dtor
  virtual Float_t GetEmcClusteringThreshold()const {Warning("GetEmcClusteringThreshold", "Not Defined" ) ; return 0. ; }  
  virtual Float_t GetEmcLocalMaxCut()const {Warning("GetEmcLocalMaxCut", "Not Defined" ) ; return 0. ; } 
  virtual Float_t GetEmcLogWeight()const {Warning("GetEmcLogWeight", "Not Defined" ) ; return 0. ; } 
  virtual Float_t GetEmcTimeGate() const {Warning("GetEmcTimeGate", "Not Defined" ) ; return 0. ; }  ;
  virtual Float_t GetCpvClusteringThreshold()const {Warning("GetCpvClusteringThreshold", "Not Defined" ) ; return 0. ; } 
  virtual Float_t GetCpvLocalMaxCut()const {Warning("GetCpvLocalMaxCut", "Not Defined" ) ; return 0. ; } 
  virtual Float_t GetCpvLogWeight()const {Warning("GetCpvLogWeight", "Not Defined" ) ; return 0. ; } 
  virtual Int_t GetRecPointsInRun()  const {Warning("GetRecPointsInRun", "Not Defined" ) ; return 0 ; } 

  virtual void MakeClusters() {Warning("MakeClusters", "Not Defined" ) ; } 
  virtual void Print(const Option_t * = "")const {Warning("Print", "Not Defined" ) ; } 

  virtual void SetEmcClusteringThreshold(Float_t) = 0;
  virtual void SetEmcLocalMaxCut(Float_t )        = 0;
    
  virtual void SetEmcLogWeight(Float_t)           = 0;
  virtual void SetEmcTimeGate(Float_t)            = 0;
  virtual void SetCpvClusteringThreshold(Float_t) = 0;
  virtual void SetCpvLocalMaxCut(Float_t)         = 0;
  virtual void SetCpvLogWeight(Float_t)           = 0;
  virtual void SetUnfolding(Bool_t)               = 0;
  void SetEventRange(Int_t first=0, Int_t last=-1) {fFirstEvent=first; fLastEvent=last; }
  void SetEventFolderName(TString name) { fEventFolderName = name ; }

  AliPHOSClusterizer & operator = (const AliPHOSClusterizer & /*rvalue*/)  {return *this ;} 
 
  virtual const char * Version() const {Warning("Version", "Not Defined" ) ; return 0 ; }  

protected:
  TString fEventFolderName ;  // event folder name
  Int_t   fFirstEvent;        // first event to process
  Int_t   fLastEvent;         // last  event to process

  ClassDef(AliPHOSClusterizer,4)  // Clusterization algorithm class 

} ;

#endif // AliPHOSCLUSTERIZER_H
