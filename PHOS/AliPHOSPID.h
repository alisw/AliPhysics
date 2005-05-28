#ifndef ALIPHOSPID_H
#define ALIPHOSPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

//_________________________________________________________________________
//  Algorithm class for the identification of particles detected in PHOS        
//  base  class                             
//  of identified particles                
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TTask.h"
#include "AliConfig.h"

class TFormula ;
class TClonesArray ;

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSGeometry ;
class AliPHOSClusterizer ;
class AliPHOSTrackSegmentMaker ;

class AliPHOSPID : public TTask {

 public:

  AliPHOSPID() ;          // ctor            
  AliPHOSPID (const TString alirunFileName, const TString eventFolderName = AliConfig::GetDefaultEventFolderName()) ;
  AliPHOSPID(const AliPHOSPID & pid) : TTask(pid) {;} 
  virtual ~AliPHOSPID() ; // dtor

  virtual Int_t GetRecParticlesInRun()  const { Warning("GetRecParticlesInRun", "not defined" ) ; return 0 ;} 
  virtual void Print(const Option_t * = "") const { Warning("Print", "not defined" ) ;}
  void SetEventRange(Int_t first=0, Int_t last=-1) {fFirstEvent=first; fLastEvent=last; }
  void SetEventFolderName(TString name) { fEventFolderName = name ; }
  virtual const char * Version() const { Warning("Version", "not defined" ) ; return 0 ; }  
  virtual void WriteRecParticles() = 0;

protected:
  TString fEventFolderName ;  // event folder name
  Int_t   fFirstEvent;        // first event to process
  Int_t   fLastEvent;         // last  event to process

private: 
  virtual void Init() { Warning("Init", "not defined" ) ; } 

  ClassDef(AliPHOSPID,4)  // Particle Identifier algorithm (base class)

} ;

#endif // ALIPHOSPID_H
