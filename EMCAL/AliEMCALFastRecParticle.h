#ifndef ALIEMCALFASTRECPARTICLE_H
#define ALIEMCALFASTRECPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A  Particle modified by EMCAL response and produced by AliEMCALvFast
//  To become a general class of AliRoot ?    
//               
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

class TClonesArray;
#include "TParticle.h"

// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALFastRecParticle : public TParticle {
  
 public:
  
  AliEMCALFastRecParticle() ;
  
  AliEMCALFastRecParticle(const AliEMCALFastRecParticle & rp) ;  // ctor
  AliEMCALFastRecParticle(const TParticle & p) ;  // ctor
  virtual ~AliEMCALFastRecParticle(){
    // dtor
  }
  virtual Int_t DistancetoPrimitive(Int_t px, Int_t py) ; 
  virtual void Draw(Option_t *option) ;  
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py) ;
  Int_t GetIndexInList() const { 
    // returns the index of this in the list
    return fIndexInList ; 
  } 
  virtual const Int_t GetNPrimaries() const {return 0 ;}
  virtual const TParticle * GetPrimary(Int_t /*index=0*/) const  {return 0 ;} 
  const Int_t GetType() const { 
    // returns the type of the particle
    return fType ; 
  } 
 
  void SetPIDBit(UInt_t fSet) {
    // Set PID bit number fSet
    fType |= (1<<fSet) ; 
  } 
  
  Bool_t TestPIDBit(UInt_t fTest) const {
    // Check PID bit number fTest
    if (fType & (1<<fTest) ) return  kTRUE ;	
    else return kFALSE ;
  }
  
  Bool_t IsPhoton           (TString purity = "low") const;
  Bool_t IsPi0              (TString purity = "low") const;
  Bool_t IsElectron         (TString purity = "low") const;
  Bool_t IsHardPhoton       () const;
  Bool_t IsHardPi0          () const;
  Bool_t IsHadron           () const;
  Bool_t IsChargedHadron    () const;
  Bool_t IsNeutralHadron    () const;
  Bool_t IsFastChargedHadron() const;
  Bool_t IsSlowChargedHadron() const;
  Bool_t IsFastNeutralHadron() const;
  Bool_t IsSlowNeutralHadron() const;
  Bool_t IsEleCon(TString purity = "low") const;  

  TString Name() const ; 
  virtual void Paint(Option_t * option="");
  virtual void Print(Option_t * option = "") const ; 
  
  void SetType(Int_t type) ;
  
  void SetIndexInList(Int_t val) { 
    // sets the value of the index in the list 
    fIndexInList = val ; 
  }
  //This has to disappear
  enum EParticleType { kUNDEFINED=-1, 
		       kNEUTRALEMFAST, kNEUTRALHAFAST,  kNEUTRALEMSLOW, kNEUTRALHASLOW, 
		       kCHARGEDEMFAST, kCHARGEDHAFAST,  kCHARGEDEMSLOW, kCHARGEDHASLOW } ; 
  
  
  typedef TClonesArray  FastRecParticlesList ; 
  
 protected:

  Int_t fIndexInList ; // the index of this RecParticle in the list stored in TreeR (to be set by analysis)
  Int_t fType ;        // particle type obtained by "virtual" reconstruction
 private:


  ClassDef(AliEMCALFastRecParticle,2)  // Reconstructed Particle produced by the fast simulation 

};

#endif // AliEMCALFASTRECPARTICLE_H
