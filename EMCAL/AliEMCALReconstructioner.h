#ifndef ALIEMCALRECONSTRUCTIONER_H
#define ALIEMCALRECONSTRUCTIONER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Wrapping class for reconstruction
//*--
//*-- Author: Gines Martinez & Yves Schutz (SUBATECH) 
//*--         Dmitri Peressounko (SUBATECH & Kurchatov Institute)

#include <stdlib.h>

// --- ROOT system ---

#include "TTask.h"
class AliEMCALDigitizer ;
class AliEMCALClusterizerv1 ;
//class AliEMCALTrackSegmentMaker ;
class AliEMCALPID ;
class AliEMCALSDigitizer ;

// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALReconstructioner : public TTask {

public:

  AliEMCALReconstructioner() ; //ctor            
  AliEMCALReconstructioner(const char * headreFile, const char * branchName = "Default",Bool_t toSplit = kFALSE) ;
  AliEMCALReconstructioner(const AliEMCALReconstructioner & rec):TTask(rec) {
    // cpy ctor: 
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;  
  }
   
  virtual ~AliEMCALReconstructioner() ;

  virtual void Exec(Option_t * option) ;

  AliEMCALDigitizer         * GetDigitizer()  const { return fDigitizer   ; }
  AliEMCALClusterizerv1       * GetClusterizer()const { return fClusterizer ; }
  //AliEMCALPID               * GetPID()        const { return fPID;          }
  //AliEMCALTrackSegmentMaker * GetTSMaker()    const { return fTSMaker ;     }
  AliEMCALSDigitizer        * GetSDigitizer() const { return fSDigitizer  ; }

  void Print(Option_t * option)const ;
  
  AliEMCALReconstructioner & operator = (const AliEMCALReconstructioner & rvalue)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;  
    return *this ; 
  }
  

private:
  void Init() ;  

private:
  
  Bool_t   fToSplit ; 
  TString  fHeaderFileName ;    // File with headers and gAlice
  TString  fDigitsBranch ;      // Title of digits branch
  TString  fRecPointBranch ;    // Title of RecPoints branch   
  TString  fTSBranch  ;         // Title of TrackSegments branch
  TString  fRecPartBranch ;     // Title of RecParticles branch 
  TString  fSDigitsBranch ;     // Title of SDigits branch      


  AliEMCALDigitizer         * fDigitizer ;   //! Pointer to AliEMCALDigitizer
  AliEMCALClusterizerv1       * fClusterizer ; //! Pointer to AliEMCALClusterizer
  //AliEMCALPID               * fPID ;         //! Pointer to AliEMCALPID
  //AliEMCALTrackSegmentMaker * fTSMaker ;     //! Pointer to AliEMCALTrackSegmentMaker
  AliEMCALSDigitizer        * fSDigitizer ;  //! Pointer to AliEMCALSDigitizer

  Bool_t   fIsInitialized ; // kTRUE if reconstructioner is initialized
 
  ClassDef(AliEMCALReconstructioner,3)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIEMCALRECONSTRUCTIONER_H
