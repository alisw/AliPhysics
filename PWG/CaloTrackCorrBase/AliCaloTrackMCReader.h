#ifndef ALICALOTRACKMCREADER_H
#define ALICALOTRACKMCREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloTrackMCReader
/// \ingroup CaloTrackCorrelationsBase
/// \brief Class for filtering generated MC particles and prepare them as input for the analysis.
///
/// Class for filtering generated MC particles and prepare them as input for the analysis.
///
/// Separates generated particles into charged (CTS) 
/// and neutral (PHOS or EMCAL acceptance)
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
class TArrayI   ;

// --- AliRoot system ---
#include "AliCaloTrackReader.h" 
class AliVCluster ;
class AliAODTrack ;
class AliAODEvent ;
class AliMCEvent  ;
class AliVEvent   ;

class AliCaloTrackMCReader : public AliCaloTrackReader {
  
 public: 
  
  AliCaloTrackMCReader() ;          // ctor
  virtual ~AliCaloTrackMCReader() ; // virtual dtor

  // Main methods in source file
  
  void   CheckOverlap(Float_t anglethres, Int_t imom, Int_t & iPrimary, Int_t & index, Int_t & pdg);

  void   FillCalorimeters(Int_t & iParticle, Int_t motherIndex, Int_t pdg) ;

  Bool_t FillInputEvent(Int_t iEntry, const char * currentFileName) ;
  
  void   InitParameters();
  
  void   MakePi0Decay() ;//, Double_t &angle);

  void   Print(const Option_t * opt) const; 
  
  void   SetCaloClusterPID   (Int_t pdgCode, AliVCluster *calo ) const ;
  
  void   SetTrackChargeAndPID(Int_t pdgCode, AliAODTrack *track) const ;
  
  void   SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ;

  // Data members setters and getters
  
  AliVEvent*  GetInputEvent()               const { return (AliVEvent *) GetMC()   ; }
  
  void      GetVertex(Double_t v[3]) const ;
  Double_t* GetVertex(Int_t evtIndex) const {return fVertex[evtIndex];}
  void      GetVertex(Double_t vertex[3], Int_t evtIndex) const 
  { vertex[0]=fVertex[evtIndex][0];  vertex[1]=fVertex[evtIndex][1];  vertex[2]=fVertex[evtIndex][2]; }  
  
  // Particle type, status, arrays 
  void      AddNeutralParticlesArray(TArrayI & array)  
      { fNeutralParticlesArray   = new TArrayI(array) ; }
  TArrayI * GetNeutralParticlesArray()      const { return  fNeutralParticlesArray ; }
  Bool_t    SkipNeutralParticles(Int_t pdg) const ;
  
  void      AddChargedParticlesArray(TArrayI & array)  
      { fChargedParticlesArray   = new TArrayI(array) ; }
  TArrayI * GetChargedParticlesArray()      const { return  fChargedParticlesArray ; }
  Bool_t    KeepChargedParticles(Int_t pdg) const ;

  void      AddStatusArray(TArrayI & array)  
      { fStatusArray   = new TArrayI(array) ; }
  TArrayI * GetStatusArray()                const { return  fStatusArray           ; }
  
  void      SwitchOnStatusSelection()             { fKeepAllStatus = kFALSE        ; }
  void      SwitchOffStatusSelection()            { fKeepAllStatus = kTRUE         ; }
  Bool_t    KeepParticleWithStatus(Int_t status) const ;
  
  void      SwitchOnOnlyGeneratorParticles()      { fOnlyGeneratorParticles = kTRUE  ; }
  void      SwitchOffOnlyGeneratorParticles()     { fOnlyGeneratorParticles = kFALSE ; }
  
  // Pi0 Overlapps, decays 
  
  void      SwitchOnPi0Decay()                    { fDecayPi0 = kTRUE              ; } 
  void      SwitchOffPi0Decay()                   { fDecayPi0 = kFALSE             ; } 
  Int_t     IsPi0DecaySwitchedOn()          const { return fDecayPi0               ; }   
  
  void      SwitchOnOverlapCheck()                { fCheckOverlap = kTRUE          ; }
  void      SwitchOffOverlapCheck()               { fCheckOverlap = kFALSE         ; }

  Float_t   GetEMCALOverlapAngle()          const { return fEMCALOverlapAngle      ; }
  void      SetEMCALOverlapAngle(Float_t angle)   { fEMCALOverlapAngle = angle     ; }
  
  Float_t   GetPHOSOverlapAngle()           const { return fPHOSOverlapAngle       ; }
  void      SetPHOSOverlapAngle (Float_t angle)   { fPHOSOverlapAngle  = angle     ; }
    
private:

  Bool_t    fDecayPi0 ;              ///< If not decayed, decay pi0 by hand.
  TArrayI * fNeutralParticlesArray ; ///< Do not keep neutral particles of this list in calorimeter.
  TArrayI * fChargedParticlesArray ; ///< Keep charged particles of this list in calorimeter.
  TArrayI * fStatusArray ;           ///< Keep particles with status of the list.
  Bool_t    fKeepAllStatus ;         ///< Do or do not select particles depending on their status code.
  Bool_t    fCheckOverlap;           ///< Check of overlapped photons from pi0 enter the calorimeter.
  Float_t   fEMCALOverlapAngle;      ///< Aperture angle of photons from decay that is not resolved by EMCAL, in radians.
  Float_t   fPHOSOverlapAngle;       ///< Aperture angle of photons from decay that is not resolved by PHOS, in radians.
  Int_t     fIndex2ndPhoton;         ///< Check overlap of first decay photon already done, internal use.
  Bool_t    fOnlyGeneratorParticles; ///< Use particles only generated by PYTHIA/HERWIG/... and not by the MC tranport G3/G4/FLUKA ...
  
  TLorentzVector fMomentum;          //!<! Momentum
  TLorentzVector fPi0Momentum;       //!<! Pi0 momentum
  TLorentzVector fGamDecayMom1;      //!<! Gamma decay 1 momentum
  TLorentzVector fGamDecayMom2;      //!<! Gamma decay 2 momentum
  
  /// Copy constructor not implemented.
  AliCaloTrackMCReader(              const AliCaloTrackMCReader & r) ; 
  
  /// Assignment operator not implemented.
  AliCaloTrackMCReader & operator = (const AliCaloTrackMCReader & r) ;
  
  /// \cond CLASSIMP
  ClassDef(AliCaloTrackMCReader,5) ; 
  /// \endcond 
  
} ;

#endif //ALICALOTRACKMCREADER_H



