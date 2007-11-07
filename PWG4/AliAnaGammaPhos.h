#ifndef ALIANAGAMMAPHOS_H
#define ALIANAGAMMAPHOS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// A basic analysis task to analyse photon detected by PHOS
// A basic analysis task to analyse photon detected by PHOS
// A basic analysis task to analyse photon detected by PHOS
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include "AliAnalysisTask.h"  

class AliESDEvent ; 
class AliAOD ; 
class AliAODEvent ; 
class TNtuple ;
class TH1D ; 
class TH1I ; 
class TTree ; 

class AliAnaGammaPhos : public AliAnalysisTask {

public:
  AliAnaGammaPhos() ;
  AliAnaGammaPhos(const char *name) ;
  AliAnaGammaPhos(const AliAnaGammaPhos& ap) ;   
  AliAnaGammaPhos& operator = (const AliAnaGammaPhos& ap) ;
  virtual ~AliAnaGammaPhos() ;
   
  virtual void ConnectInputData(Option_t * = "");
  virtual void CreateOutputObjects(); 
  virtual void Init() ; 	
  virtual void LocalInit() { Init() ; }
  virtual void Exec(Option_t * opt = "") ;
  const Float_t  GetPhotonId() const { return fPhotonId ; }
  void SetDebugLevel(Int_t level) { fDebug = level ; }
  void SetPhotonId(Float_t threshold) { fPhotonId = threshold ; }
  virtual void Terminate(Option_t * opt = "") ;

private:
  // input and output
  TTree        * fChain ;         //!pointer to the analyzed TTree or TChain
  Int_t          fDebug ;         // Debug flag
  AliESDEvent       * fESD ;           //! ESD
  AliAODEvent  * fAOD ;           //! AOD
  TClonesArray * fAODPhotons ;    //! reconstructed photons
  Int_t          fPhotonsInPhos ; //! number of photons found
  TTree        * fTreeA ;         // tree of identified photons 
  // task parameters
  Float_t   fPhotonId ;  // threshold for photon identification 

  // Histograms
  TList   * fOutputList ;	//! output data list
  TNtuple * fhPHOSPos ;		// PHOS (x,y)
  TNtuple * fhPHOS ;		// all PHOS parameters
  TH1D    * fhPHOSEnergy ;	// PHOS energy 
  TH1I    * fhPHOSDigits ;	// PHOS numer of SDigits 
  TH1D    * fhPHOSRecParticles ;// PHOS number of RecParticles
  TH1I    * fhPHOSPhotons ;		// PHOS number of photons
  TH1D    * fhPHOSInvariantMass ;// PHOS invariant mass
  TH1I    * fhPHOSDigitsEvent ;		//PHOS numbet of Sdigits per event	
   
  ClassDef(AliAnaGammaPhos, 1); // a PHOS photon analysis task 
};
#endif // ALIANAGAMMAPHOS_H
