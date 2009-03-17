#ifndef ALIANALYSISTASKPHOSEXAMPLE_H
#define ALIANALYSISTASKPHOSEXAMPLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// A basic analysis task to analyse photon detected by PHOS
// A basic analysis task to analyse photon detected by PHOS
// Adapted for AliAnalysisTaskSE and AOD production 
// by Gustavo Conesa
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

class TTree ; 
#include "AliAnalysisTaskSE.h"  

class AliESDEvent ; 
class AliAODEvent ; 
class TNtuple ;
class TH1D ; 
class TH1I ; 

class AliAnalysisTaskPHOSExample : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskPHOSExample() ;
  AliAnalysisTaskPHOSExample(const char *name) ;
  AliAnalysisTaskPHOSExample(const AliAnalysisTaskPHOSExample& ap) ;   
  AliAnalysisTaskPHOSExample& operator = (const AliAnalysisTaskPHOSExample& ap) ;
  virtual ~AliAnalysisTaskPHOSExample() ;
   
  virtual void UserCreateOutputObjects(); 
  virtual void Init() ; 	
  virtual void LocalInit() { Init() ; }
  virtual void UserExec(Option_t * opt = "") ;
  Float_t  GetPhotonId() const { return fPhotonId ; }
  void SetDebugLevel(Int_t level) { fDebug = level ; }
  void SetPhotonId(Float_t threshold) { fPhotonId = threshold ; }
  virtual void Terminate(Option_t * opt = "") ;

private:
  // input and output
  Int_t          fDebug ;         // Debug flag
  TClonesArray * fAODPhotons ;    //! reconstructed photons
  Int_t          fPhotonsInPhos ; //! number of photons found
  // task parameters
  Float_t   fPhotonId ;  // threshold for photon identification 

  // Histograms
  TList   * fOutputList ;	//! output data list
  TNtuple * fhPHOSPos ;		//! PHOS (x,y)
  TNtuple * fhPHOS ;		//! all PHOS parameters
  TH1D    * fhPHOSEnergy ;	//! PHOS energy 
  TH1I    * fhPHOSDigits ;	//! PHOS numer of SDigits 
  TH1D    * fhPHOSRecParticles ;//! PHOS number of RecParticles
  TH1I    * fhPHOSPhotons ;	//! PHOS number of photons
  TH1D    * fhPHOSInvariantMass ; //! PHOS invariant mass
  TH1I    * fhPHOSDigitsEvent ;	  //! PHOS numbet of Sdigits per event	
   
  ClassDef(AliAnalysisTaskPHOSExample, 1); // a PHOS photon analysis task 
};
#endif // ALIANALYSISTASKPHOSEXAMPLE_H
