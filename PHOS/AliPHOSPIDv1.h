#ifndef ALIPHOSPIDV1_H
#define ALIPHOSPIDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//_________________________________________________________________________
// Implementation version v1 of the PHOS particle identifier 
// Identification is based on information from CPV and EMC
// Oh yeah                 
//*-- Author: Yves Schutz (SUBATECH), Gustavo Conesa.

// --- ROOT system ---
class TVector3 ;
class TMatrix ;
class TPrincipal ;

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSEmcRecPoint ;
class AliPHOSCpvRecPoint ;

#include "AliPHOSPID.h"

class  AliPHOSPIDv1 : public AliPHOSPID {
  
public:
  
  AliPHOSPIDv1() ;          // ctor   
  AliPHOSPIDv1(const TString alirunFileNameFile, const TString eventFolderName = AliConfig::GetDefaultEventFolderName()) ;
  AliPHOSPIDv1(const AliPHOSPIDv1 & pid) ;          // cpy ctor            
  
  virtual ~AliPHOSPIDv1() ; // dtor
  
  virtual void Exec(Option_t *option);  // Does the job

  //Get file name that contain the PCA
  const TString GetFileNamePrincipal(TString particle) const;

  //Get file name that contain PID parameters
  const TString GetFileNameParameters()      const {return fFileNameParameters ;}

  // Get number of rec.particles in this run
  virtual Int_t GetRecParticlesInRun() const {return fRecParticlesInRun ;}  

  // Get PID parameters as they are defined in fParameters
  Float_t GetParameterCalibration    (Int_t i)               const;
  Float_t GetParameterCpv2Emc        (Int_t i, TString axis) const;
  Float_t GetParameterTimeGate       (Int_t i)               const;
  Float_t GetParameterToCalculateEllipse(TString particle, TString param, Int_t i) const  ;     
  Float_t GetParameterPhotonBoundary (Int_t i)               const;
  Float_t GetParameterPi0Boundary    (Int_t i)               const;

  // Get energy-dependent PID parameters
  Float_t GetCalibratedEnergy    (Float_t e)                 const;
  Float_t GetCpv2EmcDistanceCut  (TString axis, Float_t e)   const ;
  Float_t GetEllipseParameter    (TString particle, TString param, Float_t e) const;

  // Set PID parameters to change appropriate element of fParameters
  void SetParameterCalibration   (Int_t i, Float_t param);
  void SetParameterCpv2Emc       (Int_t i, TString axis, Float_t cut)  ; 
  void SetParameterTimeGate      (Int_t i, Float_t gate)  ; 
  void SetParameterToCalculateEllipse(TString particle, TString param, Int_t i, Float_t value) ;
  void SetParameterPhotonBoundary(Int_t i, Float_t param);
  void SetParameterPi0Boundary   (Int_t i, Float_t param);

  void Print() const ; 

  virtual const char * Version() const { return "pid-v1" ; }  

  AliPHOSPIDv1 & operator = (const AliPHOSPIDv1 & /*pid*/) { return *this ;} 
  
private:
  
  const TString BranchName() const ; 
  virtual void  Init() ;
  virtual void  InitParameters() ;
  void          MakeRecParticles(void ) ;
  void          MakePID(void ) ;
 // Relative Distance CPV-EMC
  Float_t GetDistance     (AliPHOSEmcRecPoint * emc, AliPHOSCpvRecPoint * cpv, Option_t * axis)const ; 
  Int_t   GetCPVBit       (AliPHOSEmcRecPoint * emc, AliPHOSCpvRecPoint * cpv, Int_t EffPur, Float_t e) const;
  Int_t   GetPrincipalBit (TString particle, const Double_t* P, Int_t EffPur, Float_t e)const ; //Principal cut
  Int_t   GetHardPhotonBit(AliPHOSEmcRecPoint * emc) const;
  Int_t   GetHardPi0Bit   (AliPHOSEmcRecPoint * emc) const;
  TVector3      GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSCpvRecPoint * cpv)const ;
  void          PrintRecParticles(Option_t * option) ;
  virtual void  WriteRecParticles() ; 
  void          SetParameters() ; //Fills the matrix of parameters
  void Unload(); 

private:

  Bool_t      fDefaultInit;              //! kTRUE if the task was created by defaut ctor (only parameters are initialized)
  Int_t       fNEvent ;                  //! current event number
  TString     fFileNamePrincipalPhoton ; //  File name of the photon principals
  TString     fFileNamePrincipalPi0 ;    //  File name of the pi0 principals
  TString     fFileNameParameters ;      //  File name with PID parameters
  TPrincipal *fPrincipalPhoton ;         //! TPrincipal from photon pca file 
  TPrincipal *fPrincipalPi0 ;            //! TPrincipal from pi0 pca file 
  Double_t   *fX ;                       //! Shower shape for the principal data 
  Double_t   *fPPhoton ;                 //! Principal photon eigenvalues
  Double_t   *fPPi0 ;                    //! Principal pi0 eigenvalues
  Int_t       fRecParticlesInRun ;       //! Total number of recparticles in one run
  TMatrix    *fParameters;               //! Matrix of identification Parameters
  // response function parameters
  // ToF
  Double_t fTphoton[3] ;                 // gaussian response for photon
  TFormula * fTFphoton ;                 // the formula   
  Double_t fTelectron[3] ;               // gaussian response for electrons
  TFormula * fTFelectron ;               // the formula   
  Double_t fTchargedhadron[3] ;          // landau   response for charged hadrons
  TFormula * fTFchargedhadron ;          // the formula   
  Double_t fTneutralhadron[3] ;          // landau   response for neutral hadrons
  TFormula * fTFneutralhadron ;          // the formula   


  ClassDef( AliPHOSPIDv1,10)  // Particle identifier implementation version 1

};

#endif // AliPHOSPIDV1_H
