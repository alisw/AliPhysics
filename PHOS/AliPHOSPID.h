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
  AliPHOSPID(const char* headerFile,const char * tsBranch = 0) ;
  virtual ~AliPHOSPID() ; // dtor

  virtual void Exec(Option_t * option) = 0 ;
  virtual char * GetRecParticlesBranch()const = 0 ;
  virtual char * GetTrackSegmentsBranch()const = 0 ;
  virtual void Init()= 0 ;

  virtual void Print(Option_t * option) const = 0 ; 
  virtual void PlotDispersionCuts()const = 0;
  virtual Bool_t ReadTrackSegments()= 0 ;

  virtual void SetIdentificationMethod(char * option = "CPV DISP" ) = 0 ;

  virtual void SetShowerProfileCut(char *  formula) = 0  ; 
  virtual void SetDispersionCut(Float_t cut) = 0  ;   
  virtual void SetCpvtoEmcDistanceCut(Float_t cut ) = 0;

  virtual void SetTrackSegmentsBranch(const char* title) = 0 ;
  virtual void SetRecParticlesBranch (const char* title) = 0 ;

  virtual void WriteRecParticles()= 0 ; 

protected:

  ClassDef(AliPHOSPID,1)  // Particle Identifier algorithm (base class)

} ;

#endif // ALIPHOSPID_H
