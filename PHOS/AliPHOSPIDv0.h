#ifndef ALIPHOSPIDV0_H
#define ALIPHOSPIDV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//_________________________________________________________________________
// Implementation version v0 of the PHOS particle identifier 
// Identification is based on information from PPSD and EMC
// Oh yeah                 
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
class TFormula ;
class TVector3 ;

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSEmcRecPoint ;
class AliPHOSRecPoint ;

#include "AliPHOSPID.h"

class  AliPHOSPIDv0 : public AliPHOSPID {

public:

  AliPHOSPIDv0() ;          // ctor            
  AliPHOSPIDv0(const char* evFolderName, const char * tsBranch = "Default");
  virtual ~AliPHOSPIDv0() ; // dtor

  virtual void Exec(Option_t * option);
  //  virtual char * GetRecParticlesBranch()const {return (char*) fRecParticlesTitle.Data() ;}      
  //  virtual char * GetTrackSegmentsBranch()const{return (char*) fTrackSegmentsTitle.Data(); }
  virtual Int_t GetRecParticlesInRun() const  {return fRecParticlesInRun ;}  

  virtual void PlotDispersionCuts()const ;
  virtual void Print()const ; 
  virtual void SetIdentificationMethod(const char * option = "CPV DISP" ){fIDOptions = option ;} 
  virtual void SetShowerProfileCut(const char * formula = "0.35*0.35 - (x-1.386)*(x-1.386) - 1.707*1.707*(y-1.008)*(y-1.008)") ;
  virtual void SetDispersionCut(Float_t cut){fDispersion = cut ; } 
  virtual void SetCpvtoEmcDistanceCut(Float_t cut )      {fCpvEmcDistance = cut ;}
  virtual void SetTimeGate(Float_t gate)                 {fTimeGate = gate ;}
  //  virtual void SetTrackSegmentsBranch(const char* title) { fTrackSegmentsTitle = title;}
  //  virtual void SetRecParticlesBranch (const char* title) { fRecParticlesTitle = title;} 
  virtual const char * Version() const { return "pid-v0" ; }  
                     
 private:
  
  virtual void Init() ;
  void     MakeRecParticles(void ) ;
  Float_t  GetDistance(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, Option_t * Axis)const ; // Relative Distance CPV-EMC
  TVector3 GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv)const ;
  void     PrintRecParticles(Option_t * option) ;
  virtual void WriteRecParticles(); 

 private:

  TString                fTrackSegmentsTitle; // branch name with track segments
  TString                fRecPointsTitle ;    // branch name with rec points
  TString                fRecParticlesTitle ; // branch name with rec particles
  TString                fIDOptions ;         // PID option
  Int_t                  fNEvent ;            // current event number

  AliPHOSClusterizer   * fClusterizer ;       // !
  AliPHOSTrackSegmentMaker * fTSMaker ;       // !

  TFormula             * fFormula ;           // formula to define cut on the shower elips axis
  Float_t                fDispersion ;        // dispersion cut
  Float_t                fCpvEmcDistance ;    // Max EMC-CPV distance
  Float_t                fTimeGate ;          // Time of the latest EmcRecPoint accepted as EM
  Int_t                  fRecParticlesInRun ; //! Total number of recparticles in one run

  ClassDef( AliPHOSPIDv0,1)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV0_H
