#ifndef ALIEMCALPID_H
#define ALIEMCALPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

//_________________________________________________________________________
//  Algorithm class for the identification of particles detected in EMCAL        
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

class AliEMCALGeometry ;
class AliEMCALClusterizer ;
class AliEMCALTrackSegmentMaker ;

class AliEMCALPID : public TTask {

 public:

  AliEMCALPID() ;          // ctor            
  AliEMCALPID(const TString alirunFileName, const TString eventFolderName = AliConfig::GetDefaultEventFolderName()) ;
  AliEMCALPID(const AliEMCALPID & pid):TTask(pid) {;} 
  virtual ~AliEMCALPID() ; // dtor

  virtual void Exec(Option_t *) = 0;
  virtual Int_t GetRecParticlesInRun()  const { Warning("GetRecParticlesInRun", "not defined" ) ; return 0 ;} 
  virtual void Print() const { Warning("Print", "not defined" ) ;}
  void SetEventRange(Int_t first=0, Int_t last=-1) {fFirstEvent=first; fLastEvent=last; }
  void SetEventFolderName(TString name) { fEventFolderName = name ; }
  virtual void SetPREtoECADistanceCut(Float_t, TString, Float_t) { Warning("SetCpvtoEmcDistanceCut", "not defined" ) ;}
  virtual void SetTimeGate(Float_t /*Cluster_En*/, TString /*Eff_Pur*/, Float_t /*gate*/) { Warning("SetTimeGate", "not defined" ) ; }
  virtual const char * Version() const { Warning("Version", "not defined" ) ; return 0 ; }  
  virtual void WriteRecParticles() = 0;

private: 
  virtual void Init() { Warning("Init", "not defined" ) ; } 

protected:
  TString fEventFolderName ;  // event folder name
  Int_t   fFirstEvent;        // first event to process
  Int_t   fLastEvent;         // last  event to process
 
  ClassDef(AliEMCALPID,3)  // Particle Identifier algorithm (base class)

} ;

#endif // ALIEMCALPID_H
