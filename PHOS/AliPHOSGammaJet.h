#ifndef ALIPHOSGammaJet_H
#define ALIPHOSGammaJet_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

//_________________________________________________________________________
//  Class for the analysis of gamma-jet correlations     
//                  
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN)

// --- ROOT system ---
#include "TTask.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrix.h"
#include "TList.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSFastGlobalReconstruction.h"
#include "../PYTHIA6/AliGenPythia.h"
// --- AliRoot header files ---

class AliPHOSGammaJet : public TTask {

public: 

  AliPHOSGammaJet() ; // default ctor
  AliPHOSGammaJet(const TString inputfilename) ; //ctor 
  AliPHOSGammaJet(const AliPHOSGammaJet & gj) ; // cpy ctor
  ~AliPHOSGammaJet() ; // dtor
  virtual void   Exec(Option_t *option); 
  void List() const; 
  Double_t GetAngleMaxParam(Int_t i){return fAngleMaxParam.At(i) ; }
  Double_t GetEtaCut(){return fEtaCut;}
  Double_t GetPhiEMCALCut(Int_t i){return fPhiEMCALCut[i];}
  TString  GetHIJINGFileName(){return fHIJINGFileName ; }
  TString  GetHistosFileName(){return fOutputFileName ; }
  Double_t GetInvMassMaxCut(){return fInvMassMaxCut ; }
  Double_t GetInvMassMinCut(){return fInvMassMinCut ; }
  Double_t GetPhiMaxCut(){return fPhiMaxCut ; }
  Double_t GetPhiMinCut(){return fPhiMinCut ; }
  Double_t GetPtCut(){return fPtCut ; }
  Double_t GetNeutralPtCut(){return fNeutralPtCut ; }
  Double_t GetChargedPtCut(){return fChargedPtCut ; }
  Double_t GetPtJetSelectionCut(){return fPtJetSelectionCut ; }
  Double_t GetMinDistance(){return fMinDistance ; }
  Double_t GetJetRatioMaxCut(){return fJetRatioMaxCut ; }
  Double_t GetJetRatioMinCut(){return fJetRatioMinCut ; }
  Double_t GetRatioMaxCut(){return fRatioMaxCut ; }
  Double_t GetRatioMinCut(){return fRatioMinCut ; }
  Int_t    GetNEvent(){return fNEvent ; }
  Int_t    GetNCones(){return fNCone ; }
  Int_t    GetNPtThres(){return fNPt ; }
  Float_t  GetCone(){return fCone ; }
  Float_t  GetPtThreshold(){return fPtThreshold ; }
  Float_t  GetCones(Int_t i){return fCones[i] ; }
  Float_t  GetPtThreshold(Int_t i){return fPtThres[i] ; }
  TString  GetConeName(Int_t i){return fNameCones[i] ; }
  TString  GetPtThresName(Int_t i){return fNamePtThres[i] ; }
  Bool_t   GetTPCCutsLikeEMCAL(){return fTPCCutsLikeEMCAL ; }
 

  Bool_t   IsAnyConeOrPt(){return fAnyConeOrPt ; }
  Bool_t   IsFastReconstruction(){return fOptFast ; }
  Bool_t   IsHIJING(){return fHIJING ; }
  Bool_t   IsOnlyCharged(){return fOnlyCharged ; }

  void Plot(TString what="all", Option_t *option="") const;
  void Print(char * opt);

  void SetAngleMaxParam(Int_t i, Double_t par)
  {fAngleMaxParam.AddAt(par,i) ; }
  void SetAnyConeOrPt(Bool_t any){fAnyConeOrPt = any ;}
  void SetEtaCut(Double_t etacut) {fEtaCut = etacut ; }
  void SetPhiEMCALCut(Double_t phi, Int_t i){fPhiEMCALCut[i] = phi; }
  void SetFastReconstruction(Bool_t fast){fOptFast = fast ; }
  void SetHIJING(Bool_t opt){fHIJING = opt; }
  void SetHIJINGFileName(TString file){fHIJINGFileName = file ; }
  void SetNCones(Int_t n){fNCone = n ; }
  void SetNPtThresholds(Int_t n){fNPt = n ; }
  void SetCones(Int_t i, Float_t cone, TString sc)
    {fCones[i] = cone ; fNameCones[i] = sc; };
  void SetCone(Float_t cone)
    {fCone = cone; }
  void SetPtThreshold(Float_t pt){fPtThreshold = pt; };
  void SetPtThresholds(Int_t i,Float_t pt, TString spt){fPtThres[i] = pt ; 
  fNamePtThres[i] = spt; };
  void SetOnlyCharged(Bool_t opt){fOnlyCharged = opt; }
  void SetPythiaFileName(TString file){fInputFileName = file ; }
  void SetHistosFileName(TString file){fOutputFileName = file ; }
  void SetInvMassCutRange(Double_t invmassmin, Double_t invmassmax)
    {fInvMassMaxCut =invmassmax;  fInvMassMinCut =invmassmin;}	
  void SetJetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
    {fJetRatioMaxCut =ratiomax;  fJetRatioMinCut = ratiomin ; }
 void SetJetTPCRatioCutRange(Double_t ratiomin, Double_t ratiomax)
    {fJetTPCRatioMaxCut =ratiomax;  fJetTPCRatioMinCut = ratiomin ; }
  void SetNEvent(Int_t n){fNEvent  = n ; }
  void SetMinDistance(Double_t min){fMinDistance  = min ; }
  void SetPhiCutRange(Double_t phimin, Double_t phimax)
  {fPhiMaxCut =phimax;  fPhiMinCut =phimin;}
  void SetPtCut(Double_t ptcut)
  {fPtCut =ptcut;}
  void SetNeutralPtCut(Double_t ptcut)
  {fNeutralPtCut =ptcut;}
  void SetChargedPtCut(Double_t ptcut)
  {fChargedPtCut =ptcut;}
  void SetPtJetSelectionCut(Double_t cut){fPtJetSelectionCut = cut; }
  void SetJetSelection(Bool_t select){ fSelect= select ; }
  void SetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
  {fRatioMaxCut = ratiomax;  fRatioMinCut = ratiomin;}
  void SetTPCCutsLikeEMCAL(Bool_t b){ fTPCCutsLikeEMCAL= b ; }
 

 private:
//   void AddHIJINGToList(TList & particleList, TList & particleListCh, 
// 		       TList & particleListNe, const Int_t iEvent, 
// 		       const TLorentzVector gamma, Double_t & rot ); 
 
  void AddHIJINGToList(Int_t iEvent, TClonesArray * particleList, 
		       TClonesArray * plCh, TClonesArray * plNe, 
		       TClonesArray * plNePHOS, const AliPHOSGeometry * geom); 


  Double_t CalculateJetRatioLimit(const Double_t ptg, const Double_t *param, 
				  const Double_t *x);

  void CreateParticleList(Int_t iEvent, TClonesArray * particleList, 
			  TClonesArray * plCh, TClonesArray * plNe, 
			  TClonesArray * plNePHOS, 
			  const AliPHOSGeometry * geom );
  
  void FillJetHistos(TClonesArray * pl, Double_t ptg, TString conf, TString type);

  void FillJetHistosAnyConeOrPt( TClonesArray * pl, Double_t ptg, TString conf, 
				 TString type, TString cone, TString ptcut);
  Bool_t IsAngleInWindow(const Float_t angle, const Float_t e);
  Bool_t IsJetSelected(const Double_t ptg, const Double_t ptjet, 
		       const TString type);

  void MakeJet(TClonesArray * particleList, 
	       Double_t ptg, Double_t phig,
	       Double_t ptl, Double_t phil, Double_t etal, 
	       TString  type, TLorentzVector & jet); 
  void MakeJetAnyConeOrPt(TClonesArray * particleList, Double_t ptg, 
			  Double_t phil, Double_t ptl, Double_t phil, 
			  Double_t etal, TString  type); 
  void GetGammaJet(TClonesArray * pl,  Double_t &pt, 
		   Double_t &phi, Double_t &eta, Bool_t &Is) ;

  void GetLeadingCharge(TClonesArray * pl, 
			Double_t ptg,  Double_t phig, 
			Double_t &pt, Double_t &eta, Double_t &phi) ;
  void GetLeadingPi0   (TClonesArray * pl, 
			Double_t ptg, Double_t phig, 
			Double_t &pt, Double_t &eta, Double_t &phi) ;

  void InitParameters();
  Double_t MakeEnergy(const Double_t energy) ;
  void MakeHistos() ;
  void MakePhoton(TLorentzVector & particle) ; 
  TVector3 MakePosition(const Double_t energy, const TVector3 pos) ;
 
  void Pi0Decay(Double_t mPi0, TLorentzVector &p0, 
		TLorentzVector &p1, TLorentzVector &p2, Double_t &angle) ;
  Double_t SigmaE(Double_t energy) ;
  Double_t SigmaP(Double_t energy) ;

  void SetJet(TParticle * part, Bool_t & b, Float_t cone, Double_t eta, 
	      Double_t phi);

 private: 
  Bool_t     fAnyConeOrPt; 
  Option_t * fOption ;         //! Fill most interesting histograms 
		               // and give interesting information
  TFile *    fOutputFile ;     //! Output file
  TString    fOutputFileName;  //! Output file Name
  TString    fInputFileName;   //!
  TString    fHIJINGFileName;  //!
  Bool_t     fHIJING;
  Double_t   fEtaCut ;         // Eta cut
  Bool_t     fOnlyCharged ;    // Only jets of charged particles
  Double_t   fPhiEMCALCut[2] ; // Phi cut maximum
  Double_t   fPhiMaxCut ;      // Phi cut maximum
  Double_t   fPhiMinCut ;      // Phi cut minimun
  Double_t   fPtCut ;          // Min pt in PHOS
  Double_t   fNeutralPtCut ;   // Min pt detected in PHOS
  Double_t   fChargedPtCut ;   // Min pt detected in TPC
  Double_t   fInvMassMaxCut ;  // Invariant Mass cut maximum
  Double_t   fInvMassMinCut ;  // Invariant Masscut minimun
  Double_t   fMinDistance ;    // Minimal distance to resolve gamma decay.
  Double_t   fRatioMaxCut ;    // Leading particle/gamma Ratio cut maximum
  Double_t   fRatioMinCut ;    // Leading particle/gamma Ratio cut minimum
  Bool_t     fTPCCutsLikeEMCAL ; //Same jet energy ratio limits for both conf.

  //Jet selection parameters
  //Fixed cuts (old)
  Double_t   fJetTPCRatioMaxCut ; // Leading particle/gamma Ratio cut maximum
  Double_t   fJetTPCRatioMinCut ; // Leading particle/gamma Ratio cut minimum
  Double_t   fJetRatioMaxCut ; // Jet/gamma Ratio cut maximum
  Double_t   fJetRatioMinCut ; // Jet/gamma Ratio cut minimum

  //Cuts depending on jet pt
  Double_t fJetE1[2];
  Double_t fJetE2[2];
  Double_t fJetSigma1[2];
  Double_t fJetSigma2[2];
  Double_t fBkgMean[6];
  Double_t fBkgRMS[6];
  Double_t fJetXMin1[6];
  Double_t fJetXMin2[6];
  Double_t fJetXMax1[6];
  Double_t fJetXMax2[6];

  Int_t      fNEvent ;         // Number of events to analyze
  Int_t      fNCone ;          // Number of jet cones sizes
  Int_t      fNPt   ;          // Number of jet particle pT threshold
  Double_t   fCone  ;        // Jet cone sizes under study
  Double_t   fCones[10];        // Jet cone sizes under study
  TString    fNameCones[10];
  Double_t   fPtThreshold;
  Double_t   fPtThres[10];     // Jet pT threshold under study
  Double_t   fPtJetSelectionCut;     // Jet pT threshold under study
  TString    fNamePtThres[10]; 
  TObjArray * fListHistos ;    //! list of Histograms
  AliPHOSFastGlobalReconstruction * fFastRec;
  Bool_t     fOptFast;
  TRandom    fRan ;       //! random number generator
  //Energy and position parameters
  Double_t   fResPara1 ;  // parameter for the energy resolution dependence  
  Double_t   fResPara2 ;  // parameter for the energy resolution dependence  
  Double_t   fResPara3 ;  // parameter for the energy resolution dependence 
  Double_t   fPosParaA ;  // parameter for the position resolution
  Double_t   fPosParaB ;  // parameter for the position resolution 
  TArrayD    fAngleMaxParam ;
  Bool_t fSelect  ;  //Select jet within limits

  ClassDef(AliPHOSGammaJet,2)
} ;
 

#endif //ALIPHOSGammaJet_H



