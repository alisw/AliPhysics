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
  AliPHOSReconstructioner(const char * headreFile) ;
  AliPHOSReconstructioner(const AliPHOSReconstructioner & rec) {
    // cpy ctor: 
    // requested by the Coding Convention
    abort() ; 
  }
   
  virtual ~AliPHOSReconstructioner() ;

  virtual void Exec(Option_t * option) ;

  AliPHOSDigitizer         * GetDigitizer()  { return fDigitizer   ; }
  AliPHOSClusterizer       * GetClusterizer(){ return fClusterizer ; }
  AliPHOSPID               * GetPID()        { return fPID;          }
  AliPHOSTrackSegmentMaker * GetTSMaker()    { return fTSMaker ;     }
  AliPHOSSDigitizer        * GetSDigitizer() { return fSDigitizer  ; }

  void Print(Option_t * option)const ;
  
  void SetBranchTitle(const char* branch,const char * title) ;
            // Sets the branch titles to separate different reconstruction flows 

  void StartFrom(char * module = "SDigitizer",char * title = 0) ;
            // From wich step reconstruction begins, 
            // title to be set to all reconstructed branches

  AliPHOSReconstructioner & operator = (const AliPHOSReconstructioner & rvalue)  {
    // assignement operator requested by coding convention but not needed
    abort() ;
    return *this ; 
  }
  

private:
  void Init() ;  

private:
  
  TString  fHeaderFileName ;    // File with headers and gAlice
  TString  fDigitsBranch ;      // Title of digits branch
  TString  fRecPointBranch ;    // -"-      RecPoints     -"-
  TString  fTSBranch  ;         // -"-      TrackSegments -"-
  TString  fRecPartBranch ;     // -"-      RecParticles  -"-
  TString  fSDigitsBranch ;     // -"-      SDigits       -"-


  AliPHOSDigitizer         * fDigitizer ;
  AliPHOSClusterizer       * fClusterizer ;
  AliPHOSPID               * fPID ;
  AliPHOSTrackSegmentMaker * fTSMaker ;
  AliPHOSSDigitizer        * fSDigitizer ;

  Bool_t   fIsInitialized ;
 
  ClassDef(AliPHOSReconstructioner,1)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIPHOSRECONSTRUCTIONER_H
