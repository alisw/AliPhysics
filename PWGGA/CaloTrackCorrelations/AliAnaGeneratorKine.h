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
  
  Bool_t CorrelateWithPartonOrJet(const TLorentzVector trigger,  
                                  const Int_t   indexTrig,                     
                                  const Int_t   pdgTrig, 
                                  const Bool_t  leading[4], 
                                  const Bool_t  isolated[4], 
                                  Int_t & iparton) ; 
  
  TList * GetCreateOutputObjects() ;
  
  void    GetPartonsAndJets() ;
    
  void    GetXE(const TLorentzVector trigger,  
                const Int_t   indexTrig,                     
                const Int_t   pdgTrig, 
                const Bool_t  leading[4], 
                const Bool_t  isolated[4], 
                const Int_t   iparton) ;    
  
  void    InitParameters() ;
  
  void    IsLeadingAndIsolated(const TLorentzVector trigger, 
                               const Int_t   indexTrig,                     
                               const Int_t   pdgTrig, 
                               Bool_t leading[4],     
                               Bool_t isolated[4]) ;
  
  void    MakeAnalysisFillAOD()  { ; }
  
  void    MakeAnalysisFillHistograms() ; 
    
private:
  
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
  
  TH1F      * fhPtPhotonLeading[4];         //! Leading photon
  TH1F      * fhPtPi0Leading[4];            //! Leading pi0
  
  TH1F      * fhPtPhotonLeadingIsolated[4]; //! Leading photon, isolated
  TH1F      * fhPtPi0LeadingIsolated[4];    //! Leading pi0, isolated

  TH2F      * fhPtPartonTypeNearPhotonLeading[4];           //! Leading photon, particle pt versus originating parton type
  TH2F      * fhPtPartonTypeNearPi0Leading[4];              //! Leading pi0, particle pt versus originating parton type
  TH2F      * fhPtPartonTypeNearPhotonLeadingIsolated[4];   //! Leading photon, particle pt versus originating parton type
  TH2F      * fhPtPartonTypeNearPi0LeadingIsolated[4];      //! Leading pi0, particle pt versus originating parton type
  
  TH2F      * fhPtPartonTypeAwayPhotonLeading[4];           //! Leading photon, particle pt versus away side parton type
  TH2F      * fhPtPartonTypeAwayPi0Leading[4];              //! Leading pi0, particle pt versus away side parton type
  TH2F      * fhPtPartonTypeAwayPhotonLeadingIsolated[4];   //! Leading photon, isolated, particle pt versus away side parton type 
  TH2F      * fhPtPartonTypeAwayPi0LeadingIsolated[4];      //! Leading pi0, isolated, particle pt versus away side parton type
  
  TH2F      * fhZHardPhotonLeading[4];           //! Leading photon, zHard
  TH2F      * fhZHardPi0Leading[4];              //! Leading pi0, zHard
  TH2F      * fhZHardPhotonLeadingIsolated[4];   //! Leading photon, isolated, zHard
  TH2F      * fhZHardPi0LeadingIsolated[4];      //! Leading pi0, isolated, zHard
  
  TH2F      * fhZPartonPhotonLeading[4];         //! Leading photon, zHard
  TH2F      * fhZPartonPi0Leading[4];            //! Leading pi0, zHard
  TH2F      * fhZPartonPhotonLeadingIsolated[4]; //! Leading photon, isolated, zHard
  TH2F      * fhZPartonPi0LeadingIsolated[4];    //! Leading pi0, isolated, zHard

  TH2F      * fhZJetPhotonLeading[4];            //! Leading photon, zHard
  TH2F      * fhZJetPi0Leading[4];               //! Leading pi0, zHard
  TH2F      * fhZJetPhotonLeadingIsolated[4];    //! Leading photon, isolated, zHard
  TH2F      * fhZJetPi0LeadingIsolated[4];       //! Leading pi0, isolated, zHard
  
  TH2F      * fhXEPhotonLeading[4];              //! Leading photon, xE away side
  TH2F      * fhXEPi0Leading[4];                 //! Leading pi0, xE away side
  TH2F      * fhXEPhotonLeadingIsolated[4];      //! Leading photon, xE away side
  TH2F      * fhXEPi0LeadingIsolated[4];         //! Leading pi0, isolated, xE away side
  
  TH2F      * fhXEUEPhotonLeading[4];              //! Leading photon, xE away side
  TH2F      * fhXEUEPi0Leading[4];                 //! Leading pi0, xE away side
  TH2F      * fhXEUEPhotonLeadingIsolated[4];      //! Leading photon, xE away side
  TH2F      * fhXEUEPi0LeadingIsolated[4];         //! Leading pi0, isolated, xE away side
  
  AliAnaGeneratorKine              (const AliAnaGeneratorKine & gk) ; // cpy ctor
  AliAnaGeneratorKine & operator = (const AliAnaGeneratorKine & gk) ; // cpy assignment
  
  ClassDef(AliAnaGeneratorKine,1)
  
} ;


#endif //ALIANAGENERATORKINE_H



