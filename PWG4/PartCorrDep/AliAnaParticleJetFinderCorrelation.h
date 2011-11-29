#ifndef ALIANAPARTICLEJETFINDERCORRELATION_H
#define ALIANAPARTICLEJETFINDERCORRELATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: AliAnaParticleJetFinderCorrelation.h 21839 2007-10-29 13:49:42Z gustavo $ */


//_________________________________________________________________________
// Class that contains the algorithm for the analysis of particle (direct gamma) - jet 
// (standard jet found with JETAN) correlation 
// Particle and jet for correlation found by independent algorithms.
// For Example direct isolated photon found in AliAnaGammaDirect and the jet with JETAN
//
//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TH2F;

//---- Analysis system ----
#include "AliAnaPartCorrBaseClass.h"

class AliAnaParticleJetFinderCorrelation : public AliAnaPartCorrBaseClass {
       
 public:   
  AliAnaParticleJetFinderCorrelation() ;            // default ctor
  virtual ~AliAnaParticleJetFinderCorrelation() {;} // virtual dtor

  // General methods
  
  void     InitParameters();
  
  TList *  GetCreateOutputObjects();

  void     MakeAnalysisFillAOD() ;
  
  void     MakeAnalysisFillHistograms() ;
  
  Int_t    SelectJet(AliAODPWG4Particle * particle, const AliAODEvent * event) const ;

  void     Print(const Option_t * opt)           const;
  
  // Settings
  
  Bool_t   OnlyIsolated()                        const { return fSelectIsolated              ; }
  void     SelectIsolated(Bool_t select)               { fSelectIsolated = select            ; }
  
  Float_t  GetConeSize()                         const { return fConeSize                    ; }
  Float_t  GetPtThresholdInCone()                const { return fPtThresholdInCone           ; }	   
  Double_t GetDeltaPhiMaxCut()                   const { return fDeltaPhiMaxCut              ; }
  Double_t GetDeltaPhiMinCut()                   const { return fDeltaPhiMinCut              ; }
  Double_t GetRatioMaxCut()                      const { return fRatioMaxCut                 ; }
  Double_t GetRatioMinCut()                      const { return fRatioMinCut                 ; }  	   
  Bool_t   AreJetRefTracks()                     const { return fUseJetRefTracks             ; }
  Bool_t   IsCorrelationMadeInHistoMaker()       const { return fMakeCorrelationInHistoMaker ; } 
  
  void     SetConeSize(Float_t cone)                   { fConeSize = cone                    ; }
  void     SetPtThresholdInCone(Float_t pt)            { fPtThresholdInCone = pt             ; }	   
  void     SetDeltaPhiCutRange(Double_t phimin, Double_t phimax)
            { fDeltaPhiMaxCut =phimax;  fDeltaPhiMinCut =phimin                              ; }
  void     SetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
            { fRatioMaxCut =ratiomax;  fRatioMinCut = ratiomin                               ; }
  void     UseJetRefTracks(Bool_t use)                 { fUseJetRefTracks = use              ; }	
  void     SetMakeCorrelationInHistoMaker(Bool_t make) { fMakeCorrelationInHistoMaker = make ; }	
    
private:

  //selection parameters  
  Double_t   fDeltaPhiMaxCut ;    //! Minimum Delta Phi Gamma-Leading
  Double_t   fDeltaPhiMinCut ;    //!  Maximum Delta Phi Gamma-Leading
  Double_t   fRatioMaxCut ;       //! Jet/ particle Ratio cut maximum
  Double_t   fRatioMinCut ;       //! Jet/particle Ratio cut minimum
  
  Double_t   fConeSize  ;         //! Jet cone size 
  Double_t   fPtThresholdInCone ; //! Jet pT threshold in jet cone
  Bool_t     fUseJetRefTracks ;   //! Use track references from JETAN not the AOD tracks
  Bool_t	   fMakeCorrelationInHistoMaker ; //!Make particle-jet correlation in histogram maker
  Bool_t     fSelectIsolated ;    // Select only trigger particles isolated
  
  // Histograms
  TH2F *     fhDeltaEta;          //! Difference of jet eta and trigger particle eta as function of trigger particle pT
  TH2F *     fhDeltaPhi;          //! Difference of jet phi and trigger particle phi as function of trigger particle pT
  TH2F *     fhDeltaPt;           //! Difference of jet pT and trigger particle pT as function of trigger particle pT
  TH2F *     fhPtRatio;           //! Ratio of jet pT and trigger particle pT as function of trigger particle pT
  TH2F *     fhPt;                //! jet pT vs trigger particle pT 
  
  TH2F *     fhFFz ;              //! Accepted reconstructed jet fragmentation function, z=ptjet/pttrig
  TH2F *     fhFFxi;              //! Accepted reconstructed jet fragmentation function, xsi = ln(pttrig/ptjet)
  TH2F *     fhFFpt;              //! Jet particle pt distribution in cone
  TH2F *     fhNTracksInCone;     //! jet multiplicity in cone
  
  AliAnaParticleJetFinderCorrelation(const AliAnaParticleJetFinderCorrelation & g) ;               // cpy ctor
  AliAnaParticleJetFinderCorrelation & operator = (const AliAnaParticleJetFinderCorrelation & g) ; // cpy assignment
  
  ClassDef(AliAnaParticleJetFinderCorrelation,2)
  
 } ;

#endif //ALIANAPARTICLEJETFINDERCORRELATION_H



