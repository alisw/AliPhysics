#ifndef ALIPHOSFASTRECPARTICLE_H
#define ALIPHOSFASTRECPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A  Particle modified by PHOS response and produced by AliPHOSvFast
//  To become a general class of AliRoot ?    
//               
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TParticle.h"
#include "TVector3.h"

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSFastRecParticle : public TParticle {
  
 public:
  
  AliPHOSFastRecParticle() ;
  
  AliPHOSFastRecParticle(const AliPHOSFastRecParticle & rp) ;  // ctor
  AliPHOSFastRecParticle(const TParticle & p) ;  // ctor
  virtual ~AliPHOSFastRecParticle(){
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
  virtual const TParticle * GetPrimary(Int_t index=0) const  {return 0 ;} 
  const Int_t GetType() const { 
    // returns the type of the particle
    return fType ; 
  } 
  
  void SetPIDBit(UInt_t fSet)
    {
      fType |= (1<<fSet) ; 
    } 
  
  Bool_t TestPIDBit(UInt_t fTest){
    if (fType & (1<<fTest) )
      return  kTRUE ;	
    else
      return kFALSE ;
  }
  
  Bool_t IsPhotonHiPu_LoEf()  {
    Bool_t pid=kFALSE ;
    if(TestPIDBit(8)&& //PCA
       TestPIDBit(5)&& //TOF
       TestPIDBit(2))  //RCPV
      pid = kTRUE;
    return pid ;
  }
  
  Bool_t IsPhotonMed_Pu_Ef(){
    Bool_t pid=kFALSE ;
    if(TestPIDBit(7)&& //PCA
       TestPIDBit(4)&& //TOF
       TestPIDBit(1))  //RCPV
      pid = kTRUE ;
    return pid ;
  }
  
  Bool_t IsPhotonHiEf_LoPu()  {
    Bool_t pid=kFALSE ;
    if(TestPIDBit(6)&& //PCA
       TestPIDBit(3)&& //TOF
       TestPIDBit(0))  //RCPV
      pid = kTRUE ;
    return pid ;
  }
  Bool_t IsPhoton()  {
    Bool_t pid=kFALSE ;
    if(IsPhotonHiEf_LoPu()) pid = kTRUE ;
    return pid ;
  }
  
  Bool_t IsFastChargedHadron()  {
    Bool_t pid=kFALSE ;
    if(TestPIDBit(5)&&TestPIDBit(4)&&TestPIDBit(3)) //TOF
      pid = kTRUE ;
    return pid ;
  }
  Bool_t IsSlowChargedHadron()  {
    Bool_t pid=kFALSE ;
    if(TestPIDBit(1)||TestPIDBit(0)) //CPV
      pid = kTRUE ;
    return pid ;
  }
  Bool_t IsFastNeutralHadron()  {
    Bool_t pid=kFALSE ;
    if(TestPIDBit(5)&&TestPIDBit(4)&&TestPIDBit(3)&& //TOF
       TestPIDBit(2)&&TestPIDBit(1)&&TestPIDBit(0))//RCPV
      pid = kTRUE ;
    return pid ;
  }
  Bool_t IsSlowNeutralHadron()  {
    Bool_t pid=kFALSE ;
    if(TestPIDBit(2)&&TestPIDBit(1)&&TestPIDBit(0))//RCPV
      pid = kTRUE ;
    return pid ;
  }

  Bool_t IsFastChargedEM()  {
    Bool_t pid=kFALSE ;
    if((TestPIDBit(8)||TestPIDBit(7)||TestPIDBit(6))&&
       TestPIDBit(5)&&TestPIDBit(4)&&TestPIDBit(3))//TOF
      pid = kTRUE ;
    return pid ;
  }
  
  Bool_t IsSlowChargedEM()  {
    Bool_t pid=kFALSE ;
    if(TestPIDBit(8)||TestPIDBit(7)||TestPIDBit(6))
      pid = kTRUE ;
    return pid ;
  }
  
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


  ClassDef(AliPHOSFastRecParticle,2)  // Reconstructed Particle produced by the fast simulation 

};

#endif // AliPHOSFASTRECPARTICLE_H
