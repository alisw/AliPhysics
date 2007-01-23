#ifndef ALIANAGAMMAJET_H
#define ALIANAGAMMAJET_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 *
 */

//_________________________________________________________________________
//  Class for the analysis of gamma-jet correlations.     
//  Basically it seaches for a prompt photon in the Calorimeters acceptance, 
//  if so we construct a jet around the highest pt particle in the opposite 
//  side in azimuth. This jet has to fullfill several conditions to be 
//  accepted. Then the fragmentation function of this jet is constructed 
//  Class created from old AliPHOSGammaJet

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TROOT.h>
#include <TChain.h>
#include "TTask.h"
#include "TArrayD.h"
#include "TChain.h"
#include <TH2F.h>
#include <TTree.h> 
#include "AliAnaGammaDirect.h" 

class AliESD ; 
 
class AliAnaGammaJet : public AliAnaGammaDirect {

public: 

  AliAnaGammaJet(const char *name) ; // default ctor
  AliAnaGammaJet(const AliAnaGammaJet & gj) ; // cpy ctor
  virtual ~AliAnaGammaJet() ; //virtual dtor
  virtual void Exec(Option_t * opt = "") ;
  virtual void Init(Option_t * opt = "");
  virtual void Terminate(Option_t * opt = "");

  Bool_t     AreJetOnlyInCTS() const {return fJetsOnlyInCTS ; } 
  Double_t GetAngleMaxParam(Int_t i) const {return fAngleMaxParam.At(i) ; }
  Double_t GetEtaEMCALCut() const {return fEtaEMCALCut;}
  Double_t GetPhiEMCALCut(Int_t i) const {return fPhiEMCALCut[i];}
  Double_t GetInvMassMaxCut() const {return fInvMassMaxCut ; }
  Double_t GetInvMassMinCut() const {return fInvMassMinCut ; }
  Double_t GetPhiMaxCut() const {return fPhiMaxCut ; }
  Double_t GetPhiMinCut() const {return fPhiMinCut ; }
  Double_t GetPtJetSelectionCut() const {return fPtJetSelectionCut ; }
  Double_t GetJetRatioMaxCut() const {return fJetRatioMaxCut ; }
  Double_t GetJetRatioMinCut() const {return fJetRatioMinCut ; }
  Double_t GetRatioMaxCut() const {return fRatioMaxCut ; }
  Double_t GetRatioMinCut() const {return fRatioMinCut ; }
  Int_t       GetNCones() const {return fNCone ; }
  Int_t       GetNPtThres() const {return fNPt ; }
  Float_t    GetCone() const {return fCone ; }
  Float_t    GetPtThreshold() const {return fPtThreshold ; }
  Float_t    GetCones(Int_t i) const {return fCones[i] ; }
  Float_t    GetPtThreshold(Int_t i) const {return fPtThres[i] ; }
  TString   GetConeName(Int_t i) const {return fNameCones[i] ; }
  TString   GetPtThresName(Int_t i) const {return fNamePtThres[i] ; }
  
  Bool_t   AreSeveralConeAndPtCuts() const {return fSeveralConeAndPtCuts ; }
  Bool_t   IsPbPb() const {return fPbPb ; }


  void Print(const Option_t * opt)const;
  
  void SetAngleMaxParam(Int_t i, Double_t par)
  {fAngleMaxParam.AddAt(par,i) ; }
  void SetSeveralConeAndPtCuts(Bool_t several){fSeveralConeAndPtCuts = several ;}
  void SetEtaEMCALCut(Double_t etacut) {fEtaEMCALCut = etacut ; }
  void SetPhiEMCALCut(Double_t phi, Int_t i){fPhiEMCALCut[i] = phi; }
  void SetPbPb(Bool_t opt){fPbPb = opt; }
  void SetNCones(Int_t n){fNCone = n ; }
  void SetNPtThresholds(Int_t n){fNPt = n ; }
  void SetCones(Int_t i, Float_t cone, TString sc)
    {fCones[i] = cone ; fNameCones[i] = sc; };
  void SetCone(Float_t cone)
    {fCone = cone; }
  void SetPtThreshold(Float_t pt){fPtThreshold = pt; };
  void SetPtThresholds(Int_t i,Float_t pt, TString spt){fPtThres[i] = pt ; 
  fNamePtThres[i] = spt; };
  void SetInvMassCutRange(Double_t invmassmin, Double_t invmassmax)
    {fInvMassMaxCut =invmassmax;  fInvMassMinCut =invmassmin;}	
  void SetJetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
    {fJetRatioMaxCut =ratiomax;  fJetRatioMinCut = ratiomin ; }
  void SetJetsOnlyInCTS(Bool_t opt){fJetsOnlyInCTS = opt; }
  void SetJetCTSRatioCutRange(Double_t ratiomin, Double_t ratiomax)
    {fJetCTSRatioMaxCut =ratiomax;  fJetCTSRatioMinCut = ratiomin ; }
  void SetPhiCutRange(Double_t phimin, Double_t phimax)
  {fPhiMaxCut =phimax;  fPhiMinCut =phimin;}
  void SetPtJetSelectionCut(Double_t cut){fPtJetSelectionCut = cut; }
  void SetJetSelection(UInt_t select){ fSelect= select ; }
  void SetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
  {fRatioMaxCut = ratiomax;  fRatioMinCut = ratiomin;}


  private:

  Double_t CalculateJetRatioLimit(const Double_t ptg, const Double_t *param, 
				  const Double_t *x);
  
  void FillJetHistos(TClonesArray * pl, Double_t ptg, Double_t ptl,TString type, TString lastname);

  Bool_t IsAngleInWindow(const Float_t angle, const Float_t e);
  Bool_t IsJetSelected(const Double_t ptg, const Double_t ptjet);

  void MakeJet(TClonesArray * particleList, TParticle * pGamma, TParticle* pLeading, TString lastname); 

  void GetLeadingCharge(TClonesArray * pl, TParticle *pGamma, TParticle * pLeading) const ;
  void GetLeadingPi0   (TClonesArray * pl, TParticle *pGamma, TParticle * pLeading)  ;
  Bool_t GetLeadingParticle(TClonesArray * plCTS, TClonesArray * plNe, TParticle *pGamma, TParticle * pLeading) ;
  void MakeHistos() ;
 
  void SetJet(TParticle * part, Bool_t & b, Float_t cone, Double_t eta, 
	      Double_t phi);

 private:

  //  TTree       *fChain ;   //!pointer to the analyzed TTree or TChain
  //AliESD       *fESD ;     //! Declaration of leave types
 
  Bool_t       fSeveralConeAndPtCuts;     //  To play with the jet cone size and pt th.
  //Bool_t       fPrintInfo ;       //Print most interesting information on screen
		               // and give interesting information

  Bool_t       fPbPb;          // PbPb event
  Bool_t       fJetsOnlyInCTS ;    // Jets measured only in TPC+ITS.
  Double_t   fEtaEMCALCut ;         // Eta EMCAL acceptance
  Double_t   fPhiEMCALCut[2] ; // Phi cut maximum
  Double_t   fPhiMaxCut ;      // Phi EMCAL maximum acceptance
  Double_t   fPhiMinCut ;      // Phi EMCAL minimum acceptance
  // Double_t   fGammaPtCut ;  // Min pt in Calorimeter
  Double_t   fInvMassMaxCut ;  // Invariant Mass cut maximum
  Double_t   fInvMassMinCut ;  // Invariant Masscut minimun
  Double_t   fRatioMaxCut ;    // Leading particle/gamma Ratio cut maximum
  Double_t   fRatioMinCut ;    // Leading particle/gamma Ratio cut minimum

  //Jet selection parameters
  //Fixed cuts (old)
  Double_t   fJetCTSRatioMaxCut ; // Leading particle/gamma Ratio cut maximum
  Double_t   fJetCTSRatioMinCut ; // Leading particle/gamma Ratio cut minimum
  Double_t   fJetRatioMaxCut ; // Jet/gamma Ratio cut maximum
  Double_t   fJetRatioMinCut ; // Jet/gamma Ratio cut minimum

  //Cuts depending on jet pt
  Double_t fJetE1[2];    //Rec. jet energy parameters
  Double_t fJetE2[2];    //Rec. jet energy parameters
  Double_t fJetSigma1[2];//Rec. sigma of jet energy  parameters
  Double_t fJetSigma2[2];//Rec. sigma of jet energy  parameters
  Double_t fBkgMean[6];  //Background mean energy 
  Double_t fBkgRMS[6];   //Background RMS
  Double_t fJetXMin1[6]; //X Factor to set jet min limit for pp
  Double_t fJetXMin2[6]; //X Factor to set jet min limit for PbPb
  Double_t fJetXMax1[6]; //X Factor to set jet max limit for pp
  Double_t fJetXMax2[6]; //X Factor to set jet max limit for PbPb

  Int_t         fNCone ;            // Number of jet cones sizes
  Int_t         fNPt   ;            // Number of jet particle pT threshold
  Double_t   fCone  ;            // Jet cone sizes under study (!fSeveralConeAndPtCuts)
  Double_t   fCones[10];         // Jet cone sizes under study (fSeveralConeAndPtCuts)
  TString     fNameCones[10];     // String name of cone to append to histos
  Double_t   fPtThreshold;       // Jet pT threshold under study(!fSeveralConeAndPtCuts)
  Double_t   fPtThres[10];       // Jet pT threshold under study(fSeveralConeAndPtCuts)
  Double_t   fPtJetSelectionCut; // Jet pt to change to low pt jets analysis
  TObjArray  *fOutputContainer ; //! output data container
  TString     fNamePtThres[10];   // String name of pt th to append to histos
  TArrayD    fAngleMaxParam ; //Max opening angle selection parameters
  UInt_t       fSelect  ;   //kTRUE: Selects all jets, no limits.
  //TString     fCalorimeter ; //PHOS or EMCAL detects Gamma
  //Bool_t      fEMCALPID ; //Fill EMCAL particle lists with particles with corresponding pid
  //Bool_t      fPHOSPID;   //Fill PHOS particle lists with particles with corresponding pid

  //Histograms

  TH2F * fhChargeRatio  ; 
  TH2F * fhPi0Ratio   ; 
  TH2F * fhDeltaPhiCharge  ;  
  TH2F * fhDeltaPhiPi0   ; 
  TH2F * fhDeltaEtaCharge  ; 
  TH2F * fhDeltaEtaPi0  ; 
  TH2F * fhAnglePair  ; 
  TH2F * fhAnglePairAccepted  ; 
  TH2F * fhAnglePairNoCut  ; 
  TH2F * fhAnglePairLeadingCut  ; 
  TH2F * fhAnglePairAngleCut   ; 
  TH2F * fhAnglePairAllCut   ; 
  TH2F * fhAnglePairLeading  ; 
  TH2F * fhInvMassPairNoCut    ; 
  TH2F * fhInvMassPairLeadingCut  ; 
  TH2F * fhInvMassPairAngleCut  ; 
  TH2F * fhInvMassPairAllCut   ; 
  TH2F * fhInvMassPairLeading  ; 
  TH1F * fhNBkg   ; 
  TH2F * fhNLeading  ; 
  TH1F * fhNJet  ; 
  TH2F * fhJetRatio  ; 
  TH2F * fhJetPt   ; 
  TH2F * fhBkgRatio   ; 
  TH2F * fhBkgPt  ; 
  TH2F * fhJetFragment  ; 
  TH2F * fhBkgFragment  ; 
  TH2F * fhJetPtDist  ; 
  TH2F * fhBkgPtDist  ; 

  TH2F * fhJetRatios[5][5];
  TH2F * fhJetPts[5][5];
  TH2F * fhBkgRatios[5][5];  
  TH2F * fhBkgPts[5][5];
  
  TH2F * fhNLeadings[5][5];
  TH1F * fhNJets[5][5];
  TH1F * fhNBkgs[5][5];
  
  TH2F * fhJetFragments[5][5];
  TH2F * fhBkgFragments[5][5];
  TH2F * fhJetPtDists[5][5];
  TH2F * fhBkgPtDists[5][5];
  
  ClassDef(AliAnaGammaJet,0)
} ;
 

#endif //ALIANAGAMMAJET_H



