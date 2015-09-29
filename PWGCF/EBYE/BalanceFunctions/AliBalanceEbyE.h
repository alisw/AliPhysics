#ifndef ALIBALANCEEBYE_H
#define ALIBALANCEEBYE_H
/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                          Class AliBalanceEbyE
//   This is the class for the Balance Function analysis on an EbyE basis
//
//    Origin: m.weber@cern.ch
//    Based on AliBalancePsi
//-------------------------------------------------------------------------
#include "TObject.h"
#include "TString.h"

class TH2D;
class TH3D;

class AliBalanceEbyE : public TObject {
 public:

  AliBalanceEbyE();
  AliBalanceEbyE(const AliBalanceEbyE& balance);
  ~AliBalanceEbyE();

   
  void SetAnalysisLevel(const char* analysisLevel) {
    fAnalysisLevel = analysisLevel;}
  void SetShuffle(Bool_t shuffle) {fShuffle = shuffle;}
  void SetDeltaEtaMax(Double_t receivedDeltaEtaMax){ fDeltaEtaMax = receivedDeltaEtaMax; }

  void InitHistograms(void);

  const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  Int_t GetNumberOfAnalyzedEvent() {return fAnalyzedEvents;}

  void CalculateBalance(Double_t gReactionPlane, 
			TObjArray* particles,
			TObjArray* particlesMixed,
			Float_t bSign,
			Double_t kMultorCent = -100,
			Double_t vertexZ = 0);

  //++++++++++++++++++//
  TH3F *GetHistBF() {return fHistBFSum;}  
  TH2D *GetQAHistHBTbefore() {return fHistHBTbefore;}
  TH2D *GetQAHistHBTafter() {return fHistHBTafter;}
  TH3D *GetQAHistConversionbefore() {return fHistConversionbefore;}
  TH3D *GetQAHistConversionafter() {return fHistConversionafter;}
  TH2D *GetQAHistPsiMinusPhi() {return fHistPsiMinusPhi;}
  TH3D *GetQAHistResonancesBefore() {return fHistResonancesBefore;}
  TH3D *GetQAHistResonancesRho() {return fHistResonancesRho;}
  TH3D *GetQAHistResonancesK0() {return fHistResonancesK0;}
  TH3D *GetQAHistResonancesLambda() {return fHistResonancesLambda;}
  TH3D *GetQAHistQbefore() {return fHistQbefore;}
  TH3D *GetQAHistQafter() {return fHistQafter;}

  void UseResonancesCut() {fResonancesCut = kTRUE;}
  void UseHBTCut(Double_t setHBTCutValue = 0.02) {
    fHBTCut = kTRUE; fHBTCutValue = setHBTCutValue;}
  void UseConversionCut(Double_t setInvMassCutConversion = 0.04) {
    fConversionCut = kTRUE; fInvMassCutConversion = setInvMassCutConversion; }
  void UseMomentumDifferenceCut(Double_t gDeltaPtCutMin) {
    fQCut = kTRUE; fDeltaPtMin = gDeltaPtCutMin;}


 private:
  Float_t   GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign); 

  Bool_t fShuffle; //shuffled balance function object
  TString fAnalysisLevel; //ESD, AOD or MC
  Int_t fAnalyzedEvents; //number of events that have been analyzed

  //EbyE histograms
  TH2F *fHistPN;//! ebye delta eta - delta phi histogram +-
  TH2F *fHistNP;//! ebye delta eta - delta phi histogram -+
  TH2F *fHistPP;//! ebye delta eta - delta phi histogram ++
  TH2F *fHistNN;//! ebye delta eta - delta phi histogram --
  TH2F *fHistBF;//! ebye delta eta - delta phi histogram BF
  TH3F *fHistBFSum;//! centrality - delta eta - delta phi histogram BF for all events

  //QA histograms
  TH2D *fHistHBTbefore; // Delta Eta vs. Delta Phi before HBT inspired cuts
  TH2D *fHistHBTafter; // Delta Eta vs. Delta Phi after HBT inspired cuts
  TH3D *fHistConversionbefore; // 3D histogram (Deta,Dphi,Invmass) before Conversion cuts
  TH3D *fHistConversionafter; // 3D histogram (Deta,Dphi,Invmass) before Conversion cuts
  TH2D *fHistPsiMinusPhi;// psi - phi QA histogram
  TH3D *fHistResonancesBefore; // 3D histogram (Deta,Dphi,Invmass) before resonance cuts
  TH3D *fHistResonancesRho;    // 3D histogram (Deta,Dphi,Invmass) after removing rho 
  TH3D *fHistResonancesK0;     // 3D histogram (Deta,Dphi,Invmass) after removing rho, K0 
  TH3D *fHistResonancesLambda; // 3D histogram (Deta,Dphi,Invmass) after removing rho, K0, and Lambda
  TH3D *fHistQbefore; // Delta Eta vs. Delta Phi before cut on momentum difference
  TH3D *fHistQafter; // Delta Eta vs. Delta Phi after cut on momentum difference

  Double_t fPsiInterval;// interval in Psi-phi1
  Double_t fDeltaEtaMax;// maximum delta eta for output THnSparse

  Bool_t fResonancesCut;//resonances cut
  Bool_t fHBTCut;//cut for two-track efficiency (like HBT group)
  Double_t fHBTCutValue;// value for two-track efficiency cut (default = 0.02 from dphicorrelations)
  Bool_t fConversionCut;//conversion cut
  Double_t fInvMassCutConversion;//invariant mass for conversion cut
  Bool_t fQCut;//cut on momentum difference to suppress femtoscopic effect correlations
  Double_t fDeltaPtMin;//delta pt cut: minimum value

  AliBalanceEbyE & operator=(const AliBalanceEbyE & ) {return *this;}

  ClassDef(AliBalanceEbyE, 1)
};

#endif
