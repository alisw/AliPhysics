#ifndef ALIPHOSRECONSTRUCTIONER_H
#define ALIPHOSRECONSTRUCTIONER_H
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
class AliPHOSDigitizer ;
class AliPHOSClusterizer ;
class AliPHOSTrackSegmentMaker ;
class AliPHOSPID ;
class AliPHOSSDigitizer ;

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSReconstructioner : public TTask {

public:

  AliPHOSReconstructioner() ; //ctor            
  AliPHOSReconstructioner(const char * headreFile, const char * branchName = "Default");
  AliPHOSReconstructioner(const AliPHOSReconstructioner & rec) : TTask(rec) {
    // cpy ctor: 
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
  }
   
  virtual ~AliPHOSReconstructioner() ;

  virtual void Exec(Option_t) ;

  AliPHOSDigitizer         * GetDigitizer()  const { return fDigitizer   ; }
  AliPHOSClusterizer       * GetClusterizer()const { return fClusterizer ; }
  AliPHOSPID               * GetPID()        const { return fPID;          }
  AliPHOSTrackSegmentMaker * GetTSMaker()    const { return fTSMaker ;     }
  AliPHOSSDigitizer        * GetSDigitizer() const { return fSDigitizer  ; }

  void Print()const ;
  
  //  void SetBranchTitle(const char* branch,const char * title) ;
  //            // Sets the branch titles to separate different reconstruction flows 
  //
  //  void StartFrom(char * module = "SDigitizer",char * title = "Default") ;
  //            // From wich step reconstruction begins, 
  //            // title to be set to all reconstructed branches

  AliPHOSReconstructioner & operator = (const AliPHOSReconstructioner & rvalue)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implementeyd") ;
    return *this ; 
  }
  

private:
  void Init() ;  

private:
  
  TString  fDigitsBranch ;      // Title of digits branch
  TString  fRecPointBranch ;    // Title of RecPoints branch   
  TString  fTSBranch  ;         // Title of TrackSegments branch
  TString  fRecPartBranch ;     // Title of RecParticles branch 
  TString  fSDigitsBranch ;     // Title of SDigits branch      


  AliPHOSDigitizer         * fDigitizer ;   //! Pointer to AliPHOSDigitizer
  AliPHOSClusterizer       * fClusterizer ; //! Pointer to AliPHOSClusterizer
  AliPHOSPID               * fPID ;         //! Pointer to AliPHOSPID
  AliPHOSTrackSegmentMaker * fTSMaker ;     //! Pointer to AliPHOSTrackSegmentMaker
  AliPHOSSDigitizer        * fSDigitizer ;  //! Pointer to AliPHOSSDigitizer

  Bool_t   fIsInitialized ; // kTRUE if reconstructioner is initialized
 
  ClassDef(AliPHOSReconstructioner,1)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIPHOSRECONSTRUCTIONER_H
