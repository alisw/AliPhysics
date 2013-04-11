#ifndef ALIBALANCEPSI_H
#define ALIBALANCEPSI_H
/*  See cxx source for full Copyright notice */


/* $Id: AliBalancePsi.h 54125 2012-01-24 21:07:41Z miweber $ */

//-------------------------------------------------------------------------
//                          Class AliBalancePsi
//   This is the class for the Balance Function writ Psi analysis
//
//    Origin: Panos Christakoglou, Nikhef, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <vector>
#include <TObject.h>
#include "TString.h"
#include "TH2D.h"

#include "AliTHn.h"

using std::vector;

#define ANALYSIS_TYPES	7
#define MAXIMUM_NUMBER_OF_STEPS	1024
#define MAXIMUM_STEPS_IN_PSI 360

class TH1D;
class TH2D;
class TH3D;

const Int_t kTrackVariablesSingle = 3;       // track variables in histogram (event class, pTtrig, vertexZ)
const Int_t kTrackVariablesPair   = 6;       // track variables in histogram (event class, dEta, dPhi, pTtrig, ptAssociated, vertexZ)
const TString gBFPsiAnalysisType[ANALYSIS_TYPES] = {"y","eta","qlong","qout","qside","qinv","phi"};

class AliBalancePsi : public TObject {
 public:
  enum EAnalysisType {
    kRapidity = 0, 
    kEta      = 1, 
    kQlong    = 2, 
    kQout     = 3, 
    kQside    = 4, 
    kQinv     = 5, 
    kPhi      = 6 
  };

  AliBalancePsi();
  AliBalancePsi(const AliBalancePsi& balance);
  ~AliBalancePsi();

  void SetCentralityIdentifier(const char* centralityId) {
    fCentralityId = centralityId;}
  
  void SetAnalysisLevel(const char* analysisLevel) {
    fAnalysisLevel = analysisLevel;}
  void SetShuffle(Bool_t shuffle) {fShuffle = shuffle;}
  void SetCentralityInterval(Double_t cStart, Double_t cStop)  { fCentStart = cStart; fCentStop = cStop;};
  void SetEventClass(TString receivedEventClass){ fEventClass = receivedEventClass; } 
  void SetDeltaEtaMax(Double_t receivedDeltaEtaMax){ fDeltaEtaMax = receivedDeltaEtaMax; }
  void SetVertexZBinning(Bool_t receivedVertexBinning=kTRUE){ fVertexBinning = receivedVertexBinning; }

  void InitHistograms(void);

  const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  Int_t GetNumberOfAnalyzedEvent() {return fAnalyzedEvents;}

  void CalculateBalance(Double_t gReactionPlane, 
			TObjArray* particles,
			TObjArray* particlesMixed,
			Float_t bSign,
			Double_t kMultorCent = -100,
			Double_t vertexZ = 0);
  
  TH2D   *GetCorrelationFunction(TString type,
				 Double_t psiMin, Double_t psiMax,
				 Double_t vertexZMin=-1,
				 Double_t vertexZMax=-1,
				 Double_t ptTriggerMin=-1.,
				 Double_t ptTriggerMax=-1.,
				 Double_t ptAssociatedMin=-1.,
				 Double_t ptAssociatedMax=-1,
				 AliBalancePsi *bMixed=NULL);
  TH2D   *GetCorrelationFunctionPN(Double_t psiMin, Double_t psiMax,
				   Double_t vertexZMin=-1,
				   Double_t vertexZMax=-1,
				   Double_t ptTriggerMin=-1.,
				   Double_t ptTriggerMax=-1.,
				   Double_t ptAssociatedMin=-1.,
				   Double_t ptAssociatedMax=-1);
  TH2D   *GetCorrelationFunctionNP(Double_t psiMin, Double_t psiMax,
				   Double_t vertexZMin=-1,
				   Double_t vertexZMax=-1,
				   Double_t ptTriggerMin=-1.,
				   Double_t ptTriggerMax=-1.,
				   Double_t ptAssociatedMin=-1.,
				   Double_t ptAssociatedMax=-1);
  TH2D   *GetCorrelationFunctionPP(Double_t psiMin, Double_t psiMax,
				   Double_t vertexZMin=-1,
				   Double_t vertexZMax=-1,
				   Double_t ptTriggerMin=-1.,
				   Double_t ptTriggerMax=-1.,
				   Double_t ptAssociatedMin=-1.,
				   Double_t ptAssociatedMax=-1);
  TH2D   *GetCorrelationFunctionNN(Double_t psiMin, Double_t psiMax,
				   Double_t vertexZMin=-1,
				   Double_t vertexZMax=-1,
				   Double_t ptTriggerMin=-1.,
				   Double_t ptTriggerMax=-1.,
				   Double_t ptAssociatedMin=-1.,
				   Double_t ptAssociatedMax=-1);  
  TH2D   *GetCorrelationFunctionChargeIndependent(Double_t psiMin, Double_t psiMax,
						  Double_t vertexZMin=-1,
						  Double_t vertexZMax=-1,
						  Double_t ptTriggerMin=-1.,
						  Double_t ptTriggerMax=-1.,
						  Double_t ptAssociatedMin=-1.,
						  Double_t ptAssociatedMax=-1);

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

  TH1D *GetBalanceFunctionHistogram(Int_t iVariableSingle,
				    Int_t iVariablePair,
				    Double_t psiMin, Double_t psiMax,
				    Double_t vertexZMin=-1,
				    Double_t vertexZMax=-1,
				    Double_t ptTriggerMin=-1.,
				    Double_t ptTriggerMax=-1.,
				    Double_t ptAssociatedMin=-1.,
				    Double_t ptAssociatedMax=-1);   //

  TH1D *GetBalanceFunctionHistogram2pMethod(Int_t iVariableSingle,
					    Int_t iVariablePair,
					    Double_t psiMin, Double_t psiMax,
					    Double_t vertexZMin=-1,
					    Double_t vertexZMax=-1,
					    Double_t ptTriggerMin=-1.,
					    Double_t ptTriggerMax=-1.,
					    Double_t ptAssociatedMin=-1.,
					    Double_t ptAssociatedMax=-1,
					    AliBalancePsi *bfMix=NULL);

  TH2D *GetBalanceFunctionDeltaEtaDeltaPhi(Double_t psiMin, Double_t psiMax,
					   Double_t vertexZMin=-1,
					   Double_t vertexZMax=-1,
					   Double_t ptTriggerMin=-1.,
					   Double_t ptTriggerMax=-1.,
					   Double_t ptAssociatedMin=-1.,
					   Double_t ptAssociatedMax=-1);   
  
  TH2D *GetBalanceFunctionDeltaEtaDeltaPhi2pMethod(Double_t psiMin, Double_t psiMax,
						   Double_t vertexZMin=-1,
						   Double_t vertexZMax=-1,
						   Double_t ptTriggerMin=-1.,
						   Double_t ptTriggerMax=-1.,
						   Double_t ptAssociatedMin=-1.,
						   Double_t ptAssociatedMax=-1.,
						   AliBalancePsi *bfMix=NULL);

  TH1D *GetBalanceFunction1DFrom2D2pMethod(Bool_t bPhi,
					   Double_t psiMin, Double_t psiMax,
					   Double_t vertexZMin=-1,
					   Double_t vertexZMax=-1,
					   Double_t ptTriggerMin=-1.,
					   Double_t ptTriggerMax=-1.,
					   Double_t ptAssociatedMin=-1.,
					   Double_t ptAssociatedMax=-1.,
					   AliBalancePsi *bfMix=NULL);

  Bool_t GetMomentsAnalytical(Int_t fVariable, TH1D* gHist,
			      Double_t &mean, Double_t &meanError,
			      Double_t &sigma, Double_t &sigmaError,
			      Double_t &skewness, Double_t &skewnessError,
			      Double_t &kurtosis, Double_t &kurtosisError);
  
  TH2D *GetQAHistHBTbefore() {return fHistHBTbefore;}
  TH2D *GetQAHistHBTafter() {return fHistHBTafter;}
  TH2D *GetQAHistConversionbefore() {return fHistConversionbefore;}
  TH2D *GetQAHistConversionafter() {return fHistConversionafter;}
  TH2D *GetQAHistPsiMinusPhi() {return fHistPsiMinusPhi;}
  TH3D *GetQAHistResonancesBefore() {return fHistResonancesBefore;}
  TH3D *GetQAHistResonancesRho() {return fHistResonancesRho;}
  TH3D *GetQAHistResonancesK0() {return fHistResonancesK0;}
  TH3D *GetQAHistResonancesLambda() {return fHistResonancesLambda;}
  TH3D *GetQAHistQbefore() {return fHistQbefore;}
  TH3D *GetQAHistQafter() {return fHistQafter;}

  void UseResonancesCut() {fResonancesCut = kTRUE;}
  void UseHBTCut() {fHBTCut = kTRUE;}
  void UseConversionCut() {fConversionCut = kTRUE;}
  void UseMomentumDifferenceCut() {fQCut = kTRUE;}

 private:
  Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign); 

  Bool_t fShuffle; //shuffled balance function object
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
  TH2D *fHistConversionbefore; // Delta Eta vs. Delta Phi before Conversion cuts
  TH2D *fHistConversionafter; // Delta Eta vs. Delta Phi before Conversion cuts
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
  Bool_t fConversionCut;//conversion cut
  Bool_t fQCut;//cut on momentum difference to suppress femtoscopic effect correlations
  Bool_t fVertexBinning;//use vertex z binning in AliTHn

  TString fEventClass;

  AliBalancePsi & operator=(const AliBalancePsi & ) {return *this;}

  ClassDef(AliBalancePsi, 1)
};

#endif
