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
  AliPHOSPIDv1(const char* headerFile, const char * tsBranch = 0) ;
  virtual ~AliPHOSPIDv1() ; // dtor

  virtual void Exec(Option_t * option);
  virtual char * GetRecParticlesBranch()const {return (char*) fRecparticlesTitle.Data() ;}      
  virtual char * GetTrackSegmentsBranch()const{return (char*) fTSTitle.Data(); }

  virtual void Init() ;

  virtual void PlotDispersionCuts()const ;

  virtual void Print(Option_t * option)const ; 
  
  virtual Bool_t ReadTrackSegments() ;

  virtual void SetIdentificationMethod(char * option = "CPV DISP" ){fIDOptions = option ;} 

  virtual void SetShowerProfileCut(char * formula = 
				   "0.35*0.35 - (x-1.386)*(x-1.386) - 1.707*1.707*(y-1.008)*(y-1.008)") ;

  virtual void SetDispersionCut(Float_t cut){fDispersion = cut ; } 
  virtual void SetCpvtoEmcDistanceCut(Float_t cut ) {fCpvEmcDistance = cut ;}
  virtual void SetTrackSegmentsBranch(const char* title) { fTSTitle = title;}
  virtual void SetRecParticlesBranch (const char* title) { fRecparticlesTitle = title;} 

  virtual void WriteRecParticles() ; 
                     

 private:
  void     MakeRecParticles(void ) ;
  Float_t  GetDistance(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, Option_t * Axis)const ; 
                                     // Relative Distance PPSD-EMC
  TVector3 GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, AliPHOSRecPoint * ppsd)const ;

  void     PrintRecParticles(Option_t * option) ;

 private:

  TString                fHeaderFileName ;
  TString                fTSTitle;
  TString                fRecPointsTitle ;
  TString                fRecparticlesTitle ;

  TString                fIDOptions ;

  Int_t                  fNEvent ;
  TObjArray            * fEmcRecPoints ;  // ! initial EMC RecPoints
  TObjArray            * fCpvRecPoints ;  // ! initial CPV RecPoints
  TClonesArray         * fTrackSegments;  // ! initial list of TrackSegments
  TClonesArray         * fRecParticles ;  // ! output

  AliPHOSClusterizer   * fClusterizer ;    // !
  AliPHOSTrackSegmentMaker * fTSMaker ;    // !

  AliPHOSGeometry      * fGeom ;           // !pointer to PHOS geometry  
  TFormula             * fFormula ;        // formula to define cut on the shouer elips axis
  Float_t                fDispersion ;     // dispersion cut
  Float_t                fCpvEmcDistance ; 

  Bool_t                 fIsInitialized ;




  ClassDef( AliPHOSPIDv1,1)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV1_H
