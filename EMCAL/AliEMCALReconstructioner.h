#ifndef ALIEMCALRECONSTRUCTIONER_H
#define ALIEMCALRECONSTRUCTIONER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Wrapping class for reconstruction
//*--
//*-- Author: Yves Schutz (SUBATECH) 
//*--         Dmitri Peressounko (SUBATECH & Kurchatov Institute)


// --- ROOT system ---

#include "TTask.h"
class AliEMCALDigitizer ;
class AliEMCALClusterizer ;
class AliEMCALPID ;
class AliEMCALSDigitizer ;
class AliESD ;

// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALReconstructioner : public TTask {

public:

  AliEMCALReconstructioner() ; //ctor            
  AliEMCALReconstructioner(const char * headerFile, const char * branchName = "Default");
  AliEMCALReconstructioner(const AliEMCALReconstructioner & rec) : TTask(rec) {
    // cpy ctor: 
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
  }
   
  virtual ~AliEMCALReconstructioner() ;

  virtual void Exec(Option_t *) ;
  void Clusters2Tracks(Int_t ievent, AliESD *event);

  AliEMCALClusterizer       * GetClusterizer()const { return fClusterizer ; }
  AliEMCALPID               * GetPID()        const { return fPID;          }
  void SetEventRange(Int_t first=0, Int_t last=-1) ; 

  void Print()const ;

  AliEMCALReconstructioner & operator = (const AliEMCALReconstructioner & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  
private:
  void Init() ;  

private:
  
  TString  fRecPointBranch ;    // Title of RecPoints branch   
  TString  fRecPartBranch ;     // Title of RecParticles branch 


  AliEMCALClusterizer       * fClusterizer ; //! Pointer to AliEMCALClusterizer
  AliEMCALPID               * fPID ;         //! Pointer to AliEMCALPID
  Bool_t  fIsInitialized ; // kTRUE if reconstructioner is initialized
  Int_t   fFirstEvent;        // first event to process
  Int_t   fLastEvent;         // last  event to process

  ClassDef(AliEMCALReconstructioner,1)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIEMCALRECONSTRUCTIONER_H
