#ifndef ALIPHOSPID_H
#define ALIPHOSPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

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
  AliPHOSPID (const TString alirunFileName, const TString eventFolderName = AliConfig::fgkDefaultEventFolderName) ;
  AliPHOSPID(const AliPHOSPID & pid) : TTask(pid) {;} 
  virtual ~AliPHOSPID() ; // dtor

  virtual void Exec(Option_t *) = 0;
  virtual const Int_t GetRecParticlesInRun()  const { Warning("GetRecParticlesInRun", "not defined" ) ; return 0 ;} 
  virtual void Print() const { Warning("Print", "not defined" ) ;}
  void   SetEventFolderName(TString name) { fEventFolderName = name ; }
  virtual const char * Version() const { Warning("Version", "not defined" ) ; return 0 ; }  
  virtual void WriteRecParticles() = 0;

private: 
  virtual void Init() { Warning("Init", "not defined" ) ; } 

protected:
  TString fEventFolderName ;  // event folder name

  ClassDef(AliPHOSPID,3)  // Particle Identifier algorithm (base class)

} ;

#endif // ALIPHOSPID_H
