#ifndef ALIPHOSPIDV1_H
#define ALIPHOSPIDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//_________________________________________________________________________
// Implementation version v1 of the PHOS particle identifier 
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

class  AliPHOSPIDv1 : public AliPHOSPID {

public:

  AliPHOSPIDv1() ;          // ctor            
  AliPHOSPIDv1(const char* headerFile, const char * tsBranch = "Default") ;
  virtual ~AliPHOSPIDv1() ; // dtor

  virtual void Exec(Option_t * option);
  virtual char * GetRecParticlesBranch()const {return (char*) fRecParticlesTitle.Data() ;}      
  virtual char * GetTrackSegmentsBranch()const{return (char*) fTrackSegmentsTitle.Data(); }
  virtual const Int_t GetRecParticlesInRun() const  {return fRecParticlesInRun ;}  
  
  virtual void Init() ;
  virtual void PlotDispersionCuts()const ;
  virtual void Print(Option_t * option)const ; 
  virtual void SetIdentificationMethod(char * option = "CPV DISP" ){fIDOptions = option ;} 
  virtual void SetShowerProfileCut(char * formula = "0.35*0.35 - (x-1.386)*(x-1.386) - 1.707*1.707*(y-1.008)*(y-1.008)") ;
  virtual void SetDispersionCut(Float_t cut){fDispersion = cut ; } 
  virtual void SetCpvtoEmcDistanceCut(Float_t cut ) {fCpvEmcDistance = cut ;}
  virtual void SetTrackSegmentsBranch(const char* title) { fTrackSegmentsTitle = title;}
  virtual void SetRecParticlesBranch (const char* title) { fRecParticlesTitle = title;} 
  virtual const char * Version() const { return "pid-v1" ; }  
                     
 private:

  void     MakeRecParticles(void ) ;
  Float_t  GetDistance(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, Option_t * Axis)const ; // Relative Distance PPSD-EMC
  TVector3 GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, AliPHOSRecPoint * ppsd)const ;
  void     PrintRecParticles(Option_t * option) ;
  virtual Bool_t ReadTrackSegments(Int_t event) ;
  virtual void WriteRecParticles(Int_t event) ; 

 private:

  TString                fHeaderFileName ;    // file name with event header
  TString                fTrackSegmentsTitle; // branch name with track segments
  TString                fRecPointsTitle ;    // branch name with rec points
  TString                fRecParticlesTitle ; // branch name with rec particles
  TString                fIDOptions ;         // PID option
  Int_t                  fNEvent ;            // current event number
  TObjArray            * fEmcRecPoints ;      // ! initial EMC RecPoints
  TObjArray            * fCpvRecPoints ;      // ! initial CPV RecPoints
  TClonesArray         * fTrackSegments;      // ! initial list of TrackSegments
  TClonesArray         * fRecParticles ;      // ! output

  AliPHOSClusterizer   * fClusterizer ;       // !
  AliPHOSTrackSegmentMaker * fTSMaker ;       // !

  TFormula             * fFormula ;           // formula to define cut on the shouer elips axis
  Float_t                fDispersion ;        // dispersion cut
  Float_t                fCpvEmcDistance ;    // Max EMC-CPV distance
  Int_t                  fRecParticlesInRun ; //! Total number of recparticles in one run
  ClassDef( AliPHOSPIDv1,1)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV1_H
