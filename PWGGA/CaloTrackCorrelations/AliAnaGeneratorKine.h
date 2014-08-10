#ifndef ALIANAGENERATORKINE_H
#define ALIANAGENERATORKINE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//___________________________________________________________________________
// Do photon/pi0 analysis for isolation and correlation
// at the generator level. Only for kine stack (ESDs)
//
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
  virtual ~AliAnaGeneratorKine() { ; } //virtual dtor              
  
  Bool_t CorrelateWithPartonOrJet(TLorentzVector trigger,
                                  Int_t   indexTrig,
                                  Int_t   pdgTrig,
                                  Bool_t  leading[4],
                                  Bool_t  isolated[4],
                                  Int_t & iparton) ; 
  
  TList * GetCreateOutputObjects() ;
  
  void    GetPartonsAndJets() ;
    
  void    GetXE(TLorentzVector trigger,
                Int_t   indexTrig,
                Int_t   pdgTrig,
                Bool_t  leading[4],
                Bool_t  isolated[4],
                Int_t   iparton) ;
  
  void    InitParameters() ;
  
  void    IsLeadingAndIsolated(TLorentzVector trigger,
                               Int_t  indexTrig,
                               Int_t  pdgTrig,
                               Bool_t leading[4],     
                               Bool_t isolated[4]) ;
    
  void    MakeAnalysisFillHistograms() ;
  
  void    SetTriggerDetector( TString name ) { fTriggerDetector = name ; }
  void    SetCalorimeter    ( TString name ) { fCalorimeter     = name ; }
  
  void    SetMinChargedPt   ( Float_t pt )   { fMinChargedPt    = pt   ; }
  void    SetMinNeutralPt   ( Float_t pt )   { fMinNeutralPt    = pt   ; }
  
    
private:
  
  TString     fTriggerDetector;             //! trigger detector, for fiducial region
  TString     fCalorimeter;                 //! detector neutral particles, for fiducial region
  
  Float_t     fMinChargedPt;                //! Minimum energy for charged particles in correlation
  Float_t     fMinNeutralPt;                //! Minimum energy for neutral particles in correlation
  
  AliStack  * fStack;                       //! access stack
  
  TParticle * fParton2;                     //! Initial state Parton
  TParticle * fParton3;                     //! Initial state Parton
  
  TParticle * fParton6;                     //! Final state Parton
  TParticle * fParton7;                     //! Final state Parton
  
  TLorentzVector fJet6;                     //! Pythia jet close to parton in position 6
  TLorentzVector fJet7;                     //! Pythia jet close to parton in position 7

  Float_t     fPtHard;                      //! Generated pT hard
  
  TH1F      * fhPtHard;                     //! pt of parton 
  TH1F      * fhPtParton;                   //! pt of parton  
  TH1F      * fhPtJet;                      //! pt of jet 
  
  TH2F      * fhPtPartonPtHard;             //! pt of parton divided to pt hard, trigger is photon 
  TH2F      * fhPtJetPtHard;                //! pt of jet divided to pt hard, trigger is photon 
  TH2F      * fhPtJetPtParton;              //! pt of parton divided to pt parton, trigger is photon 

  TH1F      * fhPtPhoton;                   //! Input photon
  TH1F      * fhPtPi0;                      //! Input pi0
  
  // Histograms arrays for 4 isolation options and 2 options on leading or not leading particle
  
  TH1F      * fhPtPhotonLeading[4];         //! Leading photon pT
  TH1F      * fhPtPi0Leading[4];            //! Leading pi0 pT

  TH2F      * fhPtPhotonLeadingSumPt[4];    //! Leading photon pT vs sum in cone
  TH2F      * fhPtPi0LeadingSumPt[4];       //! Leading pi0 pT vs sum in cone
  
  TH1F      * fhPtPhotonLeadingIsolated[4]; //! Leading photon, isolated
  TH1F      * fhPtPi0LeadingIsolated[4];    //! Leading pi0, isolated

  TH2F      * fhPtPartonTypeNearPhoton[2][4];           //! Leading photon, particle pt versus originating parton type
  TH2F      * fhPtPartonTypeNearPi0[2][4];              //! Leading pi0, particle pt versus originating parton type
  TH2F      * fhPtPartonTypeNearPhotonIsolated[2][4];   //! Leading photon, particle pt versus originating parton type
  TH2F      * fhPtPartonTypeNearPi0Isolated[2][4];      //! Leading pi0, particle pt versus originating parton type
  
  TH2F      * fhPtPartonTypeAwayPhoton[2][4];           //! Leading photon, particle pt versus away side parton type
  TH2F      * fhPtPartonTypeAwayPi0[2][4];              //! Leading pi0, particle pt versus away side parton type
  TH2F      * fhPtPartonTypeAwayPhotonIsolated[2][4];   //! Leading photon, isolated, particle pt versus away side parton type 
  TH2F      * fhPtPartonTypeAwayPi0Isolated[2][4];      //! Leading pi0, isolated, particle pt versus away side parton type
  
  TH2F      * fhZHardPhoton[2][4];           //! Leading photon, zHard
  TH2F      * fhZHardPi0[2][4];              //! Leading pi0, zHard
  TH2F      * fhZHardPhotonIsolated[2][4];   //! Leading photon, isolated, zHard
  TH2F      * fhZHardPi0Isolated[2][4];      //! Leading pi0, isolated, zHard
  
  TH2F      * fhZPartonPhoton[2][4];         //! Leading photon, zHard
  TH2F      * fhZPartonPi0[2][4];            //! Leading pi0, zHard
  TH2F      * fhZPartonPhotonIsolated[2][4]; //! Leading photon, isolated, zHard
  TH2F      * fhZPartonPi0Isolated[2][4];    //! Leading pi0, isolated, zHard

  TH2F      * fhZJetPhoton[2][4];            //! Leading photon, zHard
  TH2F      * fhZJetPi0[2][4];               //! Leading pi0, zHard
  TH2F      * fhZJetPhotonIsolated[2][4];    //! Leading photon, isolated, zHard
  TH2F      * fhZJetPi0Isolated[2][4];       //! Leading pi0, isolated, zHard
  
  TH2F      * fhXEPhoton[2][4];              //! Leading photon, xE away side
  TH2F      * fhXEPi0[2][4];                 //! Leading pi0, xE away side
  TH2F      * fhXEPhotonIsolated[2][4];      //! Leading photon, xE away side
  TH2F      * fhXEPi0Isolated[2][4];         //! Leading pi0, isolated, xE away side
  
  TH2F      * fhXEUEPhoton[2][4];              //! Leading photon, xE away side
  TH2F      * fhXEUEPi0[2][4];                 //! Leading pi0, xE away side
  TH2F      * fhXEUEPhotonIsolated[2][4];      //! Leading photon, xE away side
  TH2F      * fhXEUEPi0Isolated[2][4];         //! Leading pi0, isolated, xE away side
  
  AliAnaGeneratorKine              (const AliAnaGeneratorKine & gk) ; // cpy ctor
  AliAnaGeneratorKine & operator = (const AliAnaGeneratorKine & gk) ; // cpy assignment
  
  ClassDef(AliAnaGeneratorKine,3)
  
} ;


#endif //ALIANAGENERATORKINE_H



