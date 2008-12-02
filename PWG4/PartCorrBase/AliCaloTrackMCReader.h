#ifndef ALICALOTRACKMCREADER_H
#define ALICALOTRACKMCREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Class for reading data (Kinematics) in order to do prompt gamma or other particle  correlations
//

//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---

// --- AliRoot system ---
#include "AliCaloTrackReader.h" 
class AliAODCaloCluster ;
class AliAODTrack ;

class AliCaloTrackMCReader : public AliCaloTrackReader {
  
 public: 
  
  AliCaloTrackMCReader() ; // ctor
  AliCaloTrackMCReader(const AliCaloTrackMCReader & g) ; // cpy ctor
  AliCaloTrackMCReader & operator = (const AliCaloTrackMCReader & g) ;//cpy assignment
  virtual ~AliCaloTrackMCReader() ;//virtual dtor
  
  enum clonesType {kTParticle, kAliAOD};
  
  void InitParameters();
  
  void Print(const Option_t * opt) const; 
  
  void SwitchOnPi0Decay()  { fDecayPi0 = kTRUE ; } 
  void SwitchOffPi0Decay() { fDecayPi0 = kFALSE ; } 
  Int_t IsPi0DecaySwitchedOn() const { return fDecayPi0 ; } 
  
  void SetClonesArrayType(Int_t type){fClonesArrayType = type ;} 
  Bool_t GetClonesArrayType() const {return fClonesArrayType ;} 
    
  void AddNeutralParticlesArray(TArrayI & array)  
  { fNeutralParticlesArray   = new TArrayI(array) ; }
  TArrayI * GetNeutralParticlesArray() const   {return  fNeutralParticlesArray;}
  Bool_t SkipNeutralParticles(Int_t pdg) const ;
  
  void AddChargedParticlesArray(TArrayI & array)  
  { fChargedParticlesArray   = new TArrayI(array) ; }
  TArrayI * GetChargedParticlesArray() const   {return  fChargedParticlesArray;}
  Bool_t KeepChargedParticles(Int_t pdg) const ;

  void AddStatusArray(TArrayI & array)  
  { fStatusArray   = new TArrayI(array) ; }
  TArrayI * GetStatusArray() const   {return  fStatusArray;}
  
  void SwitchOnStatusSelection()  {fKeepAllStatus = kFALSE;}
  void SwitchOffStatusSelection() {fKeepAllStatus = kTRUE;}
  Bool_t KeepParticleWithStatus(Int_t status) const ;
  
  void GetVertex(Double_t v[3]) const ;

  void FillInputEvent() ;
  AliVEvent*  GetInputEvent() const {return GetMC();}
  void SetInputEvent(TObject* esd, TObject* aod, TObject* mc) ;
  
  void SetCaloClusterPID(const Int_t pdgCode, AliAODCaloCluster *calo) const ;
  void SetTrackChargeAndPID(const Int_t pdgCode, AliAODTrack *track) const ;
  
 private:
  
  void MakePi0Decay(TLorentzVector &p0, TLorentzVector &p1, TLorentzVector &p2) const ;//, Double_t &angle); 
  void FillCalorimeters(const Int_t iParticle, TParticle* particle, TLorentzVector momentum,   
			Int_t &indexPHOS, Int_t &indexEMCAL) ;
  
  private:
  Bool_t      fDecayPi0; //If not decayed, decay pi0 by hand
  TArrayI * fNeutralParticlesArray ; //Do not keep neutral particles of this list in calorimeter.
  TArrayI * fChargedParticlesArray ; //Keep charged particles of this list in calorimeter.
  TArrayI * fStatusArray ; //Keep particles with status of the list.
  Bool_t fKeepAllStatus ; //Do or do not select particles depending on their status code.
  Int_t fClonesArrayType; //Analysis with TParticles or AliAODCaloCluster/Track?
  
  ClassDef(AliCaloTrackMCReader,1)
    } ;


#endif //ALICALOTRACKMCREADER_H



