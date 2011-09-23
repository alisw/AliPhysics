#ifndef ALIBALANCE_H
#define ALIBALANCE_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliBalance
//   This is the class for the Balance Function analysis
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>
#include "TString.h"

#define ANALYSIS_TYPES	7
#define MAXIMUM_NUMBER_OF_STEPS	1024

class TGraphErrors;
class TObjArray;
class TH1D;

const TString gBFAnalysisType[ANALYSIS_TYPES] = {"y","eta","qlong","qout","qside","qinv","phi"};

class AliBalance : public TObject {
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

  AliBalance();
  AliBalance(const AliBalance& balance);
  ~AliBalance();
  
  //void SetNumberOfBins(Int_t ibin, Int_t ibins);
  void SetAnalysisLevel(const char* analysisLevel) {
    fAnalysisLevel = analysisLevel;}
  void SetShuffle(Bool_t shuffle) {bShuffle = shuffle;}
  void SetInterval(Int_t iAnalysisType, Double_t p1Start, Double_t p1Stop,
		   Int_t ibins, Double_t p2Start, Double_t p2Stop);
  void SetNp(Int_t analysisType, Double_t NpSet)   { fNp[analysisType] = NpSet; }
  void SetNn(Int_t analysisType, Double_t NnSet)   { fNn[analysisType] = NnSet; }
  void SetNpp(Int_t analysisType, Int_t ibin, Double_t NppSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNpp[analysisType][ibin] = NppSet; }
  void SetNpn(Int_t analysisType, Int_t ibin, Double_t NpnSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNpn[analysisType][ibin] = NpnSet; }
  void SetNnp(Int_t analysisType, Int_t ibin, Double_t NnpSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNnp[analysisType][ibin] = NnpSet; }
  void SetNnn(Int_t analysisType, Int_t ibin, Double_t NnnSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNnn[analysisType][ibin] = NnnSet; }

  void InitHistograms(void);

  const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  Int_t GetNumberOfAnalyzedEvent() {return fAnalyzedEvents;}

  Int_t GetNumberOfBins(Int_t ibin) {return fNumberOfBins[ibin];}
  Double_t GetP1Start(Int_t ibin){return fP1Start[ibin];}
  Double_t GetP1Stop(Int_t ibin){return fP1Stop[ibin];}   
  Double_t GetP2Start(Int_t ibin){return fP2Start[ibin];}
  Double_t GetP2Stop(Int_t ibin){return fP2Stop[ibin];}    
 
  Double_t GetNp(Int_t analysisType) const { return 1.0*fNp[analysisType]; }
  Double_t GetNn(Int_t analysisType) const { return 1.0*fNn[analysisType]; }
  Double_t GetNnn(Int_t analysisType, Int_t p2) const { 
    return 1.0*fNnn[analysisType][p2]; }
  Double_t GetNpp(Int_t analysisType, Int_t p2) const { 
    return 1.0*fNpp[analysisType][p2]; }
  Double_t GetNpn(Int_t analysisType, Int_t p2) const { 
    return 1.0*fNpn[analysisType][p2]; }  
  Double_t GetNnp(Int_t analysisType, Int_t p2) const { 
    return 1.0*fNnp[analysisType][p2]; }

  void CalculateBalance(TObjArray *gTrackArray);
  
  Double_t GetBalance(Int_t a, Int_t p2);
  Double_t GetError(Int_t a, Int_t p2);

  TH1D *GetHistNp(Int_t iAnalysisType) { return fHistP[iAnalysisType];}
  TH1D *GetHistNn(Int_t iAnalysisType) { return fHistN[iAnalysisType];}
  TH1D *GetHistNpn(Int_t iAnalysisType) { return fHistPN[iAnalysisType];}
  TH1D *GetHistNnp(Int_t iAnalysisType) { return fHistNP[iAnalysisType];}
  TH1D *GetHistNpp(Int_t iAnalysisType) { return fHistPP[iAnalysisType];}
  TH1D *GetHistNnn(Int_t iAnalysisType) { return fHistNN[iAnalysisType];}

  void PrintAnalysisSettings();
  TGraphErrors *DrawBalance(Int_t fAnalysisType);

  void SetHistNp(Int_t iAnalysisType, TH1D *gHist) { 
    fHistP[iAnalysisType] = gHist;}
  void SetHistNn(Int_t iAnalysisType, TH1D *gHist) { 
    fHistN[iAnalysisType] = gHist;}
  void SetHistNpn(Int_t iAnalysisType, TH1D *gHist) { 
    fHistPN[iAnalysisType] = gHist;}
  void SetHistNnp(Int_t iAnalysisType, TH1D *gHist) { 
    fHistNP[iAnalysisType] = gHist;}
  void SetHistNpp(Int_t iAnalysisType, TH1D *gHist) { 
    fHistPP[iAnalysisType] = gHist;}
  void SetHistNnn(Int_t iAnalysisType, TH1D *gHist) { 
    fHistNN[iAnalysisType] = gHist;}

  TH1D *GetBalanceFunctionHistogram(Int_t iAnalysisType);
  void PrintResults(Int_t iAnalysisType, TH1D *gHist);

 private:

  Bool_t bShuffle; //shuffled balance function object
  TString fAnalysisLevel; //ESD, AOD or MC
  Int_t fAnalyzedEvents; //number of events that have been analyzed

  Int_t fNumberOfBins[ANALYSIS_TYPES]; //number of bins of the analyzed interval
  Double_t fP1Start[ANALYSIS_TYPES];
  Double_t fP1Stop[ANALYSIS_TYPES];
  Double_t fP2Start[ANALYSIS_TYPES];
  Double_t fP2Stop[ANALYSIS_TYPES];
  Double_t fP2Step[ANALYSIS_TYPES]; 
 	
  Double_t fNnn[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //N(--)
  Double_t fNpp[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //N(++)
  Double_t fNpn[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //N(+-)
  Double_t fNnp[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //N(-+)
  Double_t fNn[ANALYSIS_TYPES], fNp[ANALYSIS_TYPES]; //number of pos./neg. inside the analyzed interval
  
  Double_t fB[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //BF matrix
  Double_t ferror[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //error of the BF
  
  TH1D *fHistP[ANALYSIS_TYPES]; //N+
  TH1D *fHistN[ANALYSIS_TYPES]; //N-
  TH1D *fHistPN[ANALYSIS_TYPES]; //N+-
  TH1D *fHistNP[ANALYSIS_TYPES]; //N-+
  TH1D *fHistPP[ANALYSIS_TYPES]; //N++
  TH1D *fHistNN[ANALYSIS_TYPES]; //N--

  AliBalance & operator=(const AliBalance & ) {return *this;}

  ClassDef(AliBalance, 3)
};

#endif
