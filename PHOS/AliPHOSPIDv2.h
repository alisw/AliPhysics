#ifndef ALIPHOSPIDV2_H
#define ALIPHOSPIDV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//_________________________________________________________________________
// Implementation version v2 of the PHOS particle identifier 
// Identification is based on information from PPSD and EMC
// Oh yeah                 
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
//class TFormula ;
class TVector3 ;
class TEllipse ;
class TPrincipal ;

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSEmcRecPoint ;
class AliPHOSRecPoint ;

#include "AliPHOSPID.h"

class  AliPHOSPIDv2 : public AliPHOSPID {

public:

  AliPHOSPIDv2() ;          // ctor            
  AliPHOSPIDv2(const char* headerFile, const char * tsBranch = "Default", const char * from = 0) ;
  virtual ~AliPHOSPIDv2() ; // dtor

  virtual void Exec(Option_t * option);
  virtual char * GetRecParticlesBranch()const {return (char*) fRecParticlesTitle.Data() ;}      
  virtual char * GetTrackSegmentsBranch()const{return (char*) fTrackSegmentsTitle.Data(); }
  virtual const Int_t GetRecParticlesInRun() const  {return fRecParticlesInRun ;}  
 
  virtual void Print(Option_t * option)const ; 
  Float_t GetCpvtoEmcDistanceCut() const {return fCpvEmcDistance ;}
  Float_t GetTimeGate()            const {return fTimeGate ;}
  virtual void SetCpvtoEmcDistanceCut(Float_t cut )      {fCpvEmcDistance = cut ;} 
  virtual void SetTimeGate(Float_t gate)                 {fTimeGate = gate ;}

  void SetEllipseXCenter(Float_t x)              {fX_center = x ;}
  void SetEllipseYCenter(Float_t y)              {fY_center = y ;}
  void SetEllipseAParameter(Float_t a)           {fA = a  ;} 
  void SetEllipseBParameter(Float_t b)           {fB = b ;}
  void SetEllipseAngle(Float_t angle)            {fAngle = angle ;}
  void SetEllipseParameters(Float_t x, Float_t y, Float_t a, Float_t b, Float_t angle);
 
  virtual void SetTrackSegmentsBranch(const char* title) { fTrackSegmentsTitle = title;}
  virtual void SetRecParticlesBranch (const char* title) { fRecParticlesTitle = title;} 
  virtual const char * Version() const { return "pid-v2" ; }  
                     
 private:

  virtual void Init() ;
  void     MakeRecParticles(void ) ;
  Float_t  GetDistance(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv, Option_t * Axis)const ; // Relative Distance CPV-EMC
  Int_t    GetPrincipalSign(Double_t* P )const ; //Principal cut
  TVector3 GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv)const ;
  void     PrintRecParticles(Option_t * option) ;
  virtual void WriteRecParticles(Int_t event) ; 

 private:

  TString                fFrom ;              // name of Recpoints and TrackSegments 
  TString                fHeaderFileName ;    // file name with event header
  TString                fTrackSegmentsTitle; // branch name with track segments
  TString                fRecPointsTitle ;    // branch name with rec points
  TString                fRecParticlesTitle ; // branch name with rec particles
 
  Int_t                  fNEvent ;            // current event number

  AliPHOSClusterizer   * fClusterizer ;       //!
  AliPHOSTrackSegmentMaker * fTSMaker ;       //!
  TPrincipal           * fPrincipal ;         //!
  
  Float_t                fCpvEmcDistance ;    // Max EMC-CPV distance
  Float_t                fTimeGate ;          // Time of the latest EmcRecPoint accepted as EM
  Int_t                  fRecParticlesInRun ; //! Total number of recparticles in one run

  Double_t*              fX ; //! Principal data 
  Double_t*              fP ; //! Principal eigenvalues

  Double_t               fX_center  ; 
  Double_t               fY_center ; 
  Double_t               fA  ; 
  Double_t               fB  ; 
  Double_t               fAngle ; 

  TString                fFileName;
  
  ClassDef( AliPHOSPIDv2,1)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV2_H
