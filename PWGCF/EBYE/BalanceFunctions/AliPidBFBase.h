#ifndef ALIPIDBFBASE_H
#define ALIPIDBFBASE_H

#include <vector>
#include <TObject.h>
#include "TString.h"
#include "TH2D.h"

#include "AliTHn.h"

using std::vector;

#define ANALYSIS_TYPES	7
#define MAXIMUM_NUMBER_OF_STEPS	1024
#define MAXIMUM_STEPS_IN_PSI 360

//-------------------------------------------------------------------------
// Base class for filling the mixed and same events histogram for BF (++, +- , -- ,-+)
// Noor Alam : vecc, Kolkata, noor1989phyalam@gmail.com
// Thanks to Panos : 
//-------------------------------------------------------------------------

class TH1D;
class TH2D;
class TH3D;


class AliPidBFBase : public TObject {
 public:

  AliPidBFBase();
  AliPidBFBase(const AliPidBFBase& balance);
  ~AliPidBFBase();

  void SetCentralityIdentifier(const char* centralityId) {
    fCentralityId = centralityId;}
  
  void SetAnalysisLevel(const char* analysisLevel) {
    fAnalysisLevel = analysisLevel;}
  void SetCentralityInterval(Double_t cStart, Double_t cStop)  { fCentStart = cStart; fCentStop = cStop;};
  void SetEventClass(TString receivedEventClass){ fEventClass = receivedEventClass; } 
  void SetDeltaEtaMax(Double_t receivedDeltaEtaMax){ fDeltaEtaMax = receivedDeltaEtaMax; }
  void SetVertexZBinning(Bool_t receivedVertexBinning=kTRUE){ fVertexBinning = receivedVertexBinning; }
  void SetCustomBinning(TString receivedCustomBinning) { fCustomBinning = receivedCustomBinning; }
  void SetTrackVariableSingle(Int_t SingleVariable) {kTrackVariablesSingle=SingleVariable;}
  void SetTrackVariablePair(Int_t PairVariable) {kTrackVariablesPair=PairVariable;}


  void InitHistograms(void);

  const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  Int_t GetNumberOfAnalyzedEvent() {return fAnalyzedEvents;}

  void CalculateBalance(Double_t gReactionPlane, 
			TObjArray* particles,
			TObjArray* particlesMixed,
			Float_t bSign,
			Double_t kMultorCent = -100,
			Double_t vertexZ = 0);


  AliTHn *GetHistNp() {return fHistP;}
  AliTHn *GetHistNn() {return fHistN;}
  AliTHn *GetHistNpn() {return fHistPN;}
  AliTHn *GetHistNnp() {return fHistNP;}
  AliTHn *GetHistNpp() {return fHistPP;}
  AliTHn *GetHistNnn() {return fHistNN;}

  void SetHistNp(AliTHn *gHist) {
    fHistP = gHist; }//fHistP->FillParent(); fHistP->DeleteContainers();}
  void SetHistNn(AliTHn *gHist) {
    fHistN = gHist; }//fHistN->FillParent(); fHistN->DeleteContainers();}
  void SetHistNpn(AliTHn *gHist) {
    fHistPN = gHist; }//fHistPN->FillParent(); fHistPN->DeleteContainers();}
  void SetHistNnp(AliTHn *gHist) {
    fHistNP = gHist; }//fHistNP->FillParent(); fHistNP->DeleteContainers();}
  void SetHistNpp(AliTHn *gHist) {
    fHistPP = gHist; }//fHistPP->FillParent(); fHistPP->DeleteContainers();}
  void SetHistNnn(AliTHn *gHist) {
    fHistNN = gHist; }//fHistNN->FillParent(); fHistNN->DeleteContainers();}

  TH2D *GetQAHistHBTbefore() {return fHistHBTbefore;}
  TH2D *GetQAHistHBTafter() {return fHistHBTafter;}
  TH2D *GetQAHistPhiStarHBTbefore() {return fHistPhiStarHBTbefore;}
  TH2D *GetQAHistPhiStarHBTafter() {return fHistPhiStarHBTafter;}
  TH3D *GetQAHistConversionbefore() {return fHistConversionbefore;}
  TH3D *GetQAHistConversionafter() {return fHistConversionafter;}
  TH2D *GetQAHistPsiMinusPhi() {return fHistPsiMinusPhi;}
  TH3D *GetQAHistResonancesBefore() {return fHistResonancesBefore;}
  TH3D *GetQAHistResonancesRho() {return fHistResonancesRho;}
  TH3D *GetQAHistResonancesK0() {return fHistResonancesK0;}
  TH3D *GetQAHistResonancesLambda() {return fHistResonancesLambda;}
  TH3D *GetQAHistQbefore() {return fHistQbefore;}
  TH3D *GetQAHistQafter() {return fHistQafter;}

  void UseMomentumOrdering(Bool_t momentumOrdering = kTRUE) {fMomentumOrdering = momentumOrdering;}
  void UseResonancesCut() {fResonancesCut = kTRUE;}
  void UseHBTCut(Double_t setHBTCutValue = 0.02) {
    fHBTCut = kTRUE; fHBTCutValue = setHBTCutValue;}
  void UseConversionCut(Double_t setInvMassCutConversion = 0.04) {
    fConversionCut = kTRUE; fInvMassCutConversion = setInvMassCutConversion; }
  void UseMomentumDifferenceCut(Double_t gDeltaPtCutMin) {
    fQCut = kTRUE; fDeltaPtMin = gDeltaPtCutMin;}

  // related to customized binning of output AliTHn
  Bool_t    IsUseVertexBinning() { return fVertexBinning; }
  TString   GetBinningString()   { return fBinningString; }
  Double_t* GetBinning(const char* configuration, const char* tag, Int_t& nBins);

 private:
  Float_t   GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign); 

  TString fAnalysisLevel; //ESD, AOD or MC
  Int_t fAnalyzedEvents; //number of events that have been analyzed

  TString fCentralityId;//Centrality identifier to be used for the histo naming

  Double_t fCentStart;
  Double_t fCentStop;

  AliTHn *fHistP; //N+
  AliTHn *fHistN; //N-
  AliTHn *fHistPN; //N+-
  AliTHn *fHistNP; //N-+
  AliTHn *fHistPP; //N++
  AliTHn *fHistNN; //N--

  //QA histograms
  TH2D *fHistHBTbefore; // Delta Eta vs. Delta Phi before HBT inspired cuts
  TH2D *fHistHBTafter; // Delta Eta vs. Delta Phi after HBT inspired cuts  
  TH2D *fHistPhiStarHBTbefore; // Delta Eta vs. Delta Phi* before HBT inspired cuts
  TH2D *fHistPhiStarHBTafter; // Delta Eta vs. Delta Phi* after HBT inspired cuts
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

  Bool_t fMomentumOrdering;//use momentum ordering pT,trig > pT,assoc (default = kTRUE)
  Bool_t fResonancesCut;//resonances cut
  Bool_t fHBTCut;//cut for two-track efficiency (like HBT group)
  Double_t fHBTCutValue;// value for two-track efficiency cut (default = 0.02 from dphicorrelations)
  Bool_t fConversionCut;//conversion cut
  Double_t fInvMassCutConversion;//invariant mass for conversion cut
  Bool_t fQCut;//cut on momentum difference to suppress femtoscopic effect correlations
  Double_t fDeltaPtMin;//delta pt cut: minimum value
  Bool_t fVertexBinning;//use vertex z binning in AliTHn
  TString fCustomBinning;//for setting customized binning
  TString fBinningString;//final binning string

  Int_t kTrackVariablesSingle; // For binning AliTHn (single)
  Int_t kTrackVariablesPair;   // For binning AliTHn (Pair)

  TString fEventClass;

  AliPidBFBase & operator=(const AliPidBFBase & ) {return *this;}

  ClassDef(AliPidBFBase, 2)
};

#endif
