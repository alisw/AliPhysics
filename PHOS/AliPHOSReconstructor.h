#ifndef ALIPHOSRECONSTRUCTOR_H
#define ALIPHOSRECONSTRUCTOR_H
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
class AliESD ;

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSReconstructor : public TTask {

public:

  AliPHOSReconstructor() ; //ctor            
  AliPHOSReconstructor(const char * headerFile, const char * branchName = "Default",const TString taskName="CTP");
  AliPHOSReconstructor(const AliPHOSReconstructor & rec) : TTask(rec) {
    // cpy ctor: 
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
  }
   
  virtual ~AliPHOSReconstructor() ;

  virtual void Exec(Option_t *) ;
  void Clusters2Tracks(Int_t ievent, AliESD *event);

  AliPHOSClusterizer       * GetClusterizer()const { return fClusterizer ; }
  AliPHOSPID               * GetPID()        const { return fPID;          }
  AliPHOSTrackSegmentMaker * GetTSMaker()    const { return fTSMaker ;     }
  void SetEventRange(Int_t first=0, Int_t last=-1) ; 

  void Print()const ;

  AliPHOSReconstructor & operator = (const AliPHOSReconstructor & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  

private:
  void Init() ;  

private:
  
  TString  fRecPointBranch ;    // Title of RecPoints branch   
  TString  fTSBranch  ;         // Title of TrackSegments branch
  TString  fRecPartBranch ;     // Title of RecParticles branch 


  AliPHOSClusterizer       * fClusterizer ; //! Pointer to AliPHOSClusterizer
  AliPHOSPID               * fPID ;         //! Pointer to AliPHOSPID
  AliPHOSTrackSegmentMaker * fTSMaker ;     //! Pointer to AliPHOSTrackSegmentMaker
  Bool_t  fIsInitialized ; // kTRUE if reconstructioner is initialized
  Int_t   fFirstEvent;        // first event to process
  Int_t   fLastEvent;         // last  event to process
 
  ClassDef(AliPHOSReconstructor,1)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIPHOSRECONSTRUCTOR_H
