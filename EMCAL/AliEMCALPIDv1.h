#ifndef ALIEMCALPIDV1_H
#define ALIEMCALPIDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//_________________________________________________________________________
// Implementation version v1 of the EMCAL particle identifier 
// Identification is based on information from CPV and EMC
// Oh yeah                 
//*-- Author: Yves Schutz (SUBATECH), Gustavo Conesa.

// --- ROOT system ---
//class TFormula ;
class TVector3 ;
class TMatrix ;
class TPrincipal ;

// --- Standard library ---

// --- AliRoot header files ---
class AliEMCALRecPoint ;

#include "AliEMCALPID.h"

class  AliEMCALPIDv1 : public AliEMCALPID {
  
 public:
  
  AliEMCALPIDv1() ;          // ctor            
  AliEMCALPIDv1(const TString alirunFileNameFile, const TString eventFolderName = AliConfig::GetDefaultEventFolderName()) ;
  AliEMCALPIDv1(const AliEMCALPIDv1 & pid) ;          // cpy ctor            

  virtual ~AliEMCALPIDv1() ; // dtor
  
  virtual void Exec(Option_t * option) ;

  virtual Int_t GetRecParticlesInRun() const  {return fRecParticlesInRun ;}  
  Float_t  GetCalibratedEnergy    (Float_t e)                 const;
  TVector3 GetMomentumDirection(AliEMCALRecPoint * emc)const ;

  virtual void Print(Option_t * option) const ;
 
  //To turn on or off the Pi0 analysis
  Bool_t GetPi0Analysis() const {return fPi0Analysis;}
  void  SetPi0Analysis(Bool_t turnonoff){ fPi0Analysis = turnonoff; }
   
 private:
  
  const TString BranchName() const ; 
  virtual void Init() ;
  virtual void InitParameters() ;
  void     MakeRecParticles(void ) ;
  void     PrintRecParticles(Option_t * option) ;
  virtual  void WriteRecParticles() ; 
  void Unload(); 

  void Print() const ; 

  virtual const char * Version() const { return "pid-v1" ; }  

  AliEMCALPIDv1 & operator = (const AliEMCALPIDv1 & /*pid*/) { return *this ;} 

 private:

  Bool_t  fDefaultInit;        //! Says if the task was created by defaut ctor (only parameters are initialized)
  Int_t   fNEvent ;            //! current event number
  TString fFileName ;          // File that contains the Principal file for analysis
  TString fFileNamePar ;       // File that contains the parameters for analysis
  Bool_t  fPi0Analysis;        //! Pi0 analysis on or off  
  Int_t   fRecParticlesInRun ; //! Total number of recparticles in one run

  ClassDef( AliEMCALPIDv1,8)  // Particle identifier implementation version 1

};

#endif // AliEMCALPIDV1_H
