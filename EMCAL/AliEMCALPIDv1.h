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
class TMatrixD ;
class TPrincipal ;

// --- Standard library ---

// --- AliRoot header files ---
class AliEMCALEmcRecPoint ;
class AliEMCALRecPoint ;

#include "AliEMCALPID.h"

class  AliEMCALPIDv1 : public AliEMCALPID {
  
 public:
  
  AliEMCALPIDv1() ;          // ctor            
  AliEMCALPIDv1(const char* headerFile, const char * tsBranch = "Default", const Bool_t toSplit=kFALSE) ;
  
  virtual ~AliEMCALPIDv1() ; // dtor
  
  virtual void Exec(Option_t * option) ;
  //  virtual char * GetRecParticlesBranch()const {return (char*) fRecParticlesTitle.Data() ;}      
  //  virtual char * GetTrackSegmentsBranch()const{return (char*) fTrackSegmentsTitle.Data(); }
  virtual const Int_t GetRecParticlesInRun() const  {return fRecParticlesInRun ;}  
  
  virtual void Print(Option_t * option) const {}
  void Print() ; 
  
 
  //To turn on or off the Pi0 analysis
  const Bool_t GetPi0Analysis(){return fPi0Analysis;}
  void  SetPi0Analysis(Bool_t turnonoff){ fPi0Analysis = turnonoff; }
 
  virtual const char * Version() const { return "pid-v1" ; }  
  
 private:
  
  const TString BranchName() const ; 
  virtual void Init() ;
  virtual void InitParameters() ;
  void     MakeRecParticles(void ) ;
  void     PrintRecParticles(Option_t * option) ;
  virtual  void WriteRecParticles(Int_t event) ; 
  
 

 private:

  Bool_t                 fDefaultInit;        //! Says if the task was created by defaut ctor (only parameters are initialized)
  TString    fFileName ;  // File that contains the Principal file for analysis
  TString    fFileNamePar ;// File that contains the parameters for analysis
  Int_t                      fNEvent ;            //! current event number
  Bool_t                     fPi0Analysis;        //! Pi0 analysis on or off  
  Int_t                      fRecParticlesInRun ; //! Total number of recparticles in one run

  ClassDef( AliEMCALPIDv1,6)  // Particle identifier implementation version 1

};

#endif // AliEMCALPIDV1_H
