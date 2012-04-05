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

  void SetPsiInterval(Double_t psiInterval) {
    fPsiInterval = psiInterval; 
    fPsiNumberOfBins = (Int_t)(360./fPsiInterval);}
  void SetPsiNumberOfBins(Int_t psiNumberOfBins) {
    fPsiNumberOfBins = psiNumberOfBins;
    fPsiInterval = 360./fPsiNumberOfBins;}

  void SetCentralityIdentifier(const char* centralityId) {
    fCentralityId = centralityId;}
  
  void SetAnalysisLevel(const char* analysisLevel) {
    fAnalysisLevel = analysisLevel;}
  void SetShuffle(Bool_t shuffle) {bShuffle = shuffle;}
  void SetInterval(Int_t iAnalysisType, Double_t p1Start, Double_t p1Stop,
		   Int_t ibins, Double_t p2Start, Double_t p2Stop,
		   Double_t psiInterval);
  void SetCentralityInterval(Double_t cStart, Double_t cStop)  { fCentStart = cStart; fCentStop = cStop;};
  void SetNp(Int_t analysisType, Int_t psiBin, Double_t NpSet)   { fNp[analysisType][psiBin] = NpSet; }
  void SetNn(Int_t analysisType, Int_t psiBin, Double_t NnSet)   { fNn[analysisType][psiBin] = NnSet; }
  void SetNpp(Int_t analysisType, Int_t psiBin, Int_t ibin, Double_t NppSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNpp[analysisType][psiBin][ibin] = NppSet; }
  void SetNpn(Int_t analysisType, Int_t psiBin, Int_t ibin, Double_t NpnSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNpn[analysisType][psiBin][ibin] = NpnSet; }
  void SetNnp(Int_t analysisType, Int_t psiBin, Int_t ibin, Double_t NnpSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNnp[analysisType][psiBin][ibin] = NnpSet; }
  void SetNnn(Int_t analysisType, Int_t psiBin, Int_t ibin, Double_t NnnSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNnn[analysisType][psiBin][ibin] = NnnSet; }

  void InitHistograms(void);

  const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  Int_t GetNumberOfAnalyzedEvent() {return fAnalyzedEvents;}

  Int_t GetNumberOfBins(Int_t ibin) {return fNumberOfBins[ibin];}
  Double_t GetP1Start(Int_t ibin){return fP1Start[ibin];}
  Double_t GetP1Stop(Int_t ibin){return fP1Stop[ibin];}   
  Double_t GetP2Start(Int_t ibin){return fP2Start[ibin];}
  Double_t GetP2Stop(Int_t ibin){return fP2Stop[ibin];}    
 
  Double_t GetNp(Int_t analysisType, Int_t psiBin) const { return 1.0*fNp[analysisType][psiBin]; }
  Double_t GetNn(Int_t analysisType, Int_t psiBin) const { return 1.0*fNn[analysisType][psiBin]; }
  Double_t GetNnn(Int_t analysisType, Int_t psiBin, Int_t p2) const { 
    return 1.0*fNnn[analysisType][psiBin][p2]; }
  Double_t GetNpp(Int_t analysisType, Int_t psiBin, Int_t p2) const { 
    return 1.0*fNpp[analysisType][psiBin][p2]; }
  Double_t GetNpn(Int_t analysisType, Int_t psiBin, Int_t p2) const { 
    return 1.0*fNpn[analysisType][psiBin][p2]; }  
  Double_t GetNnp(Int_t analysisType, Int_t psiBin, Int_t p2) const { 
    return 1.0*fNnp[analysisType][psiBin][p2]; }

  void CalculateBalance(Float_t fCentrality, Double_t gReactionPlane, vector<Double_t> **chargeVector);
  
  Double_t GetBalance(Int_t iAnalysisType, Int_t psiBin, Int_t p2);
  Double_t GetError(Int_t iAnalysisType, Int_t psiBin, Int_t p2);

  TH3D *GetHistNp(Int_t iAnalysisType) { return fHistP[iAnalysisType];}
  TH3D *GetHistNn(Int_t iAnalysisType) { return fHistN[iAnalysisType];}
  TH3D *GetHistNpn(Int_t iAnalysisType) { return fHistPN[iAnalysisType];}
  TH3D *GetHistNnp(Int_t iAnalysisType) { return fHistNP[iAnalysisType];}
  TH3D *GetHistNpp(Int_t iAnalysisType) { return fHistPP[iAnalysisType];}
  TH3D *GetHistNnn(Int_t iAnalysisType) { return fHistNN[iAnalysisType];}

  void PrintAnalysisSettings();
  TGraphErrors *DrawBalance(Int_t fAnalysisType, Int_t psiBin);

  void SetHistNp(Int_t iAnalysisType, TH3D *gHist) { 
    fHistP[iAnalysisType] = gHist;}
  void SetHistNn(Int_t iAnalysisType, TH3D *gHist) { 
    fHistN[iAnalysisType] = gHist;}
  void SetHistNpn(Int_t iAnalysisType, TH3D *gHist) { 
    fHistPN[iAnalysisType] = gHist;}
  void SetHistNnp(Int_t iAnalysisType, TH3D *gHist) { 
    fHistNP[iAnalysisType] = gHist;}
  void SetHistNpp(Int_t iAnalysisType, TH3D *gHist) { 
    fHistPP[iAnalysisType] = gHist;}
  void SetHistNnn(Int_t iAnalysisType, TH3D *gHist) { 
    fHistNN[iAnalysisType] = gHist;}

  TH1D *GetBalanceFunctionHistogram(Int_t iAnalysisType,
				    Double_t centrMin, 
				    Double_t centrMax, 
				    Double_t psiMin, Double_t psiMax);
  void PrintResults(Int_t iAnalysisType, TH1D *gHist);

 private:
  Bool_t bShuffle; //shuffled balance function object
  TString fAnalysisLevel; //ESD, AOD or MC
  Int_t fAnalyzedEvents; //number of events that have been analyzed

  TString fCentralityId;//Centrality identifier to be used for the histo naming

  Int_t fNumberOfBins[ANALYSIS_TYPES];//number of bins of the analyzed interval
  Double_t fP1Start[ANALYSIS_TYPES];
  Double_t fP1Stop[ANALYSIS_TYPES];
  Double_t fP2Start[ANALYSIS_TYPES];
  Double_t fP2Stop[ANALYSIS_TYPES];
  Double_t fP2Step[ANALYSIS_TYPES]; 
  Double_t fCentStart;
  Double_t fCentStop;

 	
  Double_t fNnn[ANALYSIS_TYPES][MAXIMUM_STEPS_IN_PSI][MAXIMUM_NUMBER_OF_STEPS]; //N(--)
  Double_t fNpp[ANALYSIS_TYPES][MAXIMUM_STEPS_IN_PSI][MAXIMUM_NUMBER_OF_STEPS]; //N(++)
  Double_t fNpn[ANALYSIS_TYPES][MAXIMUM_STEPS_IN_PSI][MAXIMUM_NUMBER_OF_STEPS]; //N(+-)
  Double_t fNnp[ANALYSIS_TYPES][MAXIMUM_STEPS_IN_PSI][MAXIMUM_NUMBER_OF_STEPS]; //N(-+)
  Double_t fNn[ANALYSIS_TYPES][MAXIMUM_STEPS_IN_PSI], fNp[ANALYSIS_TYPES][MAXIMUM_STEPS_IN_PSI]; //number of pos./neg. inside the analyzed interval
  
  Double_t fB[ANALYSIS_TYPES][MAXIMUM_STEPS_IN_PSI][MAXIMUM_NUMBER_OF_STEPS]; //BF matrix
  Double_t ferror[ANALYSIS_TYPES][MAXIMUM_STEPS_IN_PSI][MAXIMUM_NUMBER_OF_STEPS]; //error of the BF
  
  TH3D *fHistP[ANALYSIS_TYPES]; //N+
  TH3D *fHistN[ANALYSIS_TYPES]; //N-
  TH3D *fHistPN[ANALYSIS_TYPES]; //N+-
  TH3D *fHistNP[ANALYSIS_TYPES]; //N-+
  TH3D *fHistPP[ANALYSIS_TYPES]; //N++
  TH3D *fHistNN[ANALYSIS_TYPES]; //N--

  Double_t fPsiInterval;// interval in Psi-phi1
  Int_t fPsiNumberOfBins;// number of bins

  AliBalancePsi & operator=(const AliBalancePsi & ) {return *this;}

  ClassDef(AliBalancePsi, 1)
};

#endif
