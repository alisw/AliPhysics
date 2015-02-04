#ifndef ALIANAGENERATORKINE_H
#define ALIANAGENERATORKINE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//___________________________________________________________________________
// Do direct photon/decay photon (eta, pi0, other)/pi0/eta isolation
// and correlation with partons/jets/hadrons analysis at the generator level.
// For MC kinematics at ESD and AOD level
// Jets only considered in the case of Pythia, check what to do with others.//
//
//-- Author: Gustavo Conesa (LPSC-CNRS-Grenoble)

// --- ROOT ---
class TH2F ;
class TParticle ;
class AliStack ;
class TLorentzVector ;

// --- ANALYSIS ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaGeneratorKine : public AliAnaCaloTrackCorrBaseClass {
       
public:
  
  AliAnaGeneratorKine() ; // default ctor

  virtual ~AliAnaGeneratorKine() { delete fFidCutTrigger ; } //virtual dtor
  
  enum mcPrimTypes { kmcPrimPhoton = 0, kmcPrimPi0Decay = 1, kmcPrimEtaDecay  = 2, kmcPrimOtherDecay  = 3,
                     kmcPrimPi0    = 4, kmcPrimEta      = 5                                               } ;
  
  static const Int_t fgkNmcPrimTypes = 6;
  static const Int_t fgkNLead = 2;
  static const Int_t fgkNIso  = 4;

  Bool_t CorrelateWithPartonOrJet(Int_t   indexTrig,
                                  Int_t   pdgTrig,
                                  Bool_t  leading [fgkNIso],
                                  Bool_t  isolated[fgkNIso],
                                  Int_t & iparton) ; 
  
  TList * GetCreateOutputObjects() ;
  
  void    GetPartonsAndJets() ;
    
  void    GetXE(Int_t   indexTrig,
                Int_t   pdgTrig,
                Bool_t  leading [fgkNIso],
                Bool_t  isolated[fgkNIso],
                Int_t   iparton) ;
  
  void    InitParameters() ;
  
  void    IsLeadingAndIsolated(Int_t  indexTrig,
                               Int_t  pdgTrig,
                               Bool_t leading [fgkNIso],     
                               Bool_t isolated[fgkNIso]) ;
    
  void    MakeAnalysisFillHistograms() ;
  
  void    SetTriggerDetector( TString & det ) ;
  void    SetTriggerDetector( Int_t  det )    ;
  
  void    SetMinChargedPt   ( Float_t pt )   { fMinChargedPt    = pt   ; }
  void    SetMinNeutralPt   ( Float_t pt )   { fMinNeutralPt    = pt   ; }
  
  // Detector for trigger particles acceptance
  AliFiducialCut * GetFiducialCutForTrigger()
  { if(!fFidCutTrigger)  fFidCutTrigger  = new AliFiducialCut(); return  fFidCutTrigger  ; }
  virtual void     SetFiducialCut(AliFiducialCut * fc)
  { delete fFidCutTrigger;  fFidCutTrigger  = fc      ; }
  
private:

  Int_t            fTriggerDetector ;       // Detector : EMCAL, PHOS, CTS
  TString          fTriggerDetectorString ; // Detector : EMCAL, PHOS, CTS

  AliFiducialCut * fFidCutTrigger;          //! fiducial cut for the trigger detector
  
  Float_t          fMinChargedPt;           //! Minimum energy for charged particles in correlation
  Float_t          fMinNeutralPt;           //! Minimum energy for neutral particles in correlation
  
  AliStack       * fStack;                  //! access ESD stack
  TClonesArray   * fAODMCparticles ;        //! access AOD stack

//  TParticle      * fParton2;              //! Initial state Parton
//  TParticle      * fParton3;              //! Initial state Parton
  
  TLorentzVector   fParton6;                //! Final state Parton
  TLorentzVector   fParton7;                //! Final state Parton
  
  Int_t            fParton6PDG;             //! Final state Parton PDG
  Int_t            fParton7PDG;             //! Final state Parton PDG
  
  TLorentzVector   fJet6;                   //! Pythia jet close to parton in position 6
  TLorentzVector   fJet7;                   //! Pythia jet close to parton in position 7

  TLorentzVector   fTrigger;                //! Trigger momentum, avoid generating TLorentzVectors per event
  TLorentzVector   fLVTmp;                  //! momentum, avoid generating TLorentzVectors per event
  
  Int_t            fNPrimaries;             //! N primaries
  Float_t          fPtHard;                 //! Generated pT hard
  
  // Histograms
  
  TH1F      * fhPtHard;                     //! pt of parton
  TH1F      * fhPtParton;                   //! pt of parton  
  TH1F      * fhPtJet;                      //! pt of jet 
  
  TH2F      * fhPtPartonPtHard;             //! pt of parton divided to pt hard, trigger is photon 
  TH2F      * fhPtJetPtHard;                //! pt of jet divided to pt hard, trigger is photon 
  TH2F      * fhPtJetPtParton;              //! pt of parton divided to pt parton, trigger is photon 

  TH1F      * fhPt      [fgkNmcPrimTypes];  //! Input particle
  TH2F      * fhPtOrigin[fgkNmcPrimTypes];  //! Input particle vs particle type originating it (if meson decay, its mother)
  TH2F      * fhPtOtherDecayMesonId;        //! Decay photons, not originating from pi0 or eta, ID of the particle
  
  TH2F      * fhPtPi0Status;                //! Input pi0 status
  TH2F      * fhPtEtaStatus;                //! Input eta status
  TH2F      * fhPtPi0DecayStatus;           //! Input pi0 decay meson status
  TH2F      * fhPtEtaDecayStatus;           //! Input eta decay meson status
  TH2F      * fhPtOtherDecayStatus;         //! Input other decay particle status

  TH1F      * fhPtPi0Not2Gamma;             //! Input pi0 not 2 gamma decay
  TH1F      * fhPtEtaNot2Gamma;             //! Input eta not 2 gamma decay
  TH1F      * fhPtGammaFromPi0Not2Gamma;    //! Input gamma from pi0 not 2 gamma decay
  TH1F      * fhPtGammaFromEtaNot2Gamma;    //! Input gamma from eta not 2 gamma decay
  
  // Histograms arrays for 4 isolation options and 2 options on leading or not leading particle

  TH2F      * fhPtAcceptedGammaJet                      [fgkNLead][fgkNIso]; //! gamma-jet pair in acceptance (jet in good eta window)

  
  TH1F      * fhPtLeading               [fgkNmcPrimTypes]          [fgkNIso]; //! pT

  TH2F      * fhPtLeadingSumPt          [fgkNmcPrimTypes]          [fgkNIso]; //! pT vs sum in cone
  
  TH1F      * fhPtLeadingIsolated       [fgkNmcPrimTypes]          [fgkNIso]; //! isolated

  TH2F      * fhPtPartonTypeNear        [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! particle pt versus originating parton type
  
  TH2F      * fhPtPartonTypeNearIsolated[fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! pt versus originating parton type

  TH2F      * fhPtPartonTypeAway        [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! pt versus away side parton type

  TH2F      * fhPtPartonTypeAwayIsolated[fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! isolated, particle pt versus away side parton type
  
  TH2F      * fhZHard                   [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! zHard
 
  TH2F      * fhZHardIsolated           [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! isolated, zHard
  
  TH2F      * fhZParton                 [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! zHard

  TH2F      * fhZPartonIsolated         [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! isolated, zHard

  TH2F      * fhZJet                    [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! zHard

  TH2F      * fhZJetIsolated            [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! isolated, zHard
  
  TH2F      * fhXE                      [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! xE away side

  TH2F      * fhXEIsolated              [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! xE away side
  
  TH2F      * fhXEUE                    [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! xE away side

  TH2F      * fhXEUEIsolated            [fgkNmcPrimTypes][fgkNLead][fgkNIso]; //! xE away side
  
  
  AliAnaGeneratorKine              (const AliAnaGeneratorKine & gk) ; // cpy ctor
  AliAnaGeneratorKine & operator = (const AliAnaGeneratorKine & gk) ; // cpy assignment
  
  ClassDef(AliAnaGeneratorKine,7)
  
} ;


#endif //ALIANAGENERATORKINE_H



