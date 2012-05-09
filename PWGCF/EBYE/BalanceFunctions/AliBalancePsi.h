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

#include <TObject.h>
#include "TString.h"

#define ANALYSIS_TYPES	7
#define MAXIMUM_NUMBER_OF_STEPS	1024
#define MAXIMUM_STEPS_IN_PSI 360

class TGraphErrors;
class TH1D;
class TH2D;
class TH3D;
class AliTHn;

const Int_t nTrackVariablesSingle = 4;       // track variables in histogram (eta, phi, pTtrig, centrality)
const Int_t nTrackVariablesPair   = 4;       // track variables in histogram (dEta, dPhi, pT, pTtrig, centrality)
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
  void SetShuffle(Bool_t shuffle) {bShuffle = shuffle;}
  void SetCentralityInterval(Double_t cStart, Double_t cStop)  { fCentStart = cStart; fCentStop = cStop;};

  void InitHistograms(void);

  const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  Int_t GetNumberOfAnalyzedEvent() {return fAnalyzedEvents;}

  void CalculateBalance(Float_t fCentrality, Double_t gReactionPlane, vector<Double_t> **chargeVector);
  
  AliTHn *GetHistNp() {return fHistP;}
  AliTHn *GetHistNn() {return fHistN;}
  AliTHn *GetHistNpn() {return fHistPN;}
  AliTHn *GetHistNnp() {return fHistNP;}
  AliTHn *GetHistNpp() {return fHistPP;}
  AliTHn *GetHistNnn() {return fHistNN;}

  void SetHistNp(AliTHn *gHist) {fHistP = gHist;}
  void SetHistNn(AliTHn *gHist) {fHistN = gHist;}
  void SetHistNpn(AliTHn *gHist) {fHistPN = gHist;}
  void SetHistNnp(AliTHn *gHist) {fHistNP = gHist;}
  void SetHistNpp(AliTHn *gHist) {fHistPP = gHist;}
  void SetHistNnn(AliTHn *gHist) {fHistNN = gHist;}

  TH1D *GetBalanceFunctionHistogram(Int_t iVariableSingle,
				    Int_t iVariablePair,
				    Double_t centrMin, 
				    Double_t centrMax, 
				    Double_t psiMin, Double_t psiMax);

 private:
  Bool_t bShuffle; //shuffled balance function object
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

  Double_t fPsiInterval;// interval in Psi-phi1

  AliBalancePsi & operator=(const AliBalancePsi & ) {return *this;}

  ClassDef(AliBalancePsi, 1)
};

#endif
