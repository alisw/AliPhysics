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

  //#include <stdlib.h>

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
  AliEMCALReconstructioner(const char * headreFile, const char * branchName = "Default");
  AliEMCALReconstructioner(const AliEMCALReconstructioner & rec) : TTask(rec) {
    // cpy ctor: 
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
  }
   
  virtual ~AliEMCALReconstructioner() ;

  virtual void Exec(Option_t *) ;
  void Clusters2Tracks(Int_t ievent, AliESD *event);

  AliEMCALDigitizer         * GetDigitizer()  const { return fDigitizer   ; }
  AliEMCALClusterizer       * GetClusterizer()const { return fClusterizer ; }
  AliEMCALPID               * GetPID()        const { return fPID;          }
  AliEMCALSDigitizer        * GetSDigitizer() const { return fSDigitizer  ; }

  void Print()const ;
  
  //  void SetBranchTitle(const char* branch,const char * title) ;
  //            // Sets the branch titles to separate different reconstruction flows 
  //
  //  void StartFrom(char * module = "SDigitizer",char * title = "Default") ;
  //            // From wich step reconstruction begins, 
  //            // title to be set to all reconstructed branches

  AliEMCALReconstructioner & operator = (const AliEMCALReconstructioner & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  

private:
  void Init() ;  

private:
  
  TString  fDigitsBranch ;      // Title of digits branch
  TString  fRecPointBranch ;    // Title of RecPoints branch   
  TString  fRecPartBranch ;     // Title of RecParticles branch 
  TString  fSDigitsBranch ;     // Title of SDigits branch      


  AliEMCALDigitizer         * fDigitizer ;   //! Pointer to AliEMCALDigitizer
  AliEMCALClusterizer       * fClusterizer ; //! Pointer to AliEMCALClusterizer
  AliEMCALPID               * fPID ;         //! Pointer to AliEMCALPID
  AliEMCALSDigitizer        * fSDigitizer ;  //! Pointer to AliEMCALSDigitizer

  Bool_t   fIsInitialized ; // kTRUE if reconstructioner is initialized
 
  ClassDef(AliEMCALReconstructioner,1)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIEMCALRECONSTRUCTIONER_H
