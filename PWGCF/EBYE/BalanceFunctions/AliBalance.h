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

#include <vector>
#include <TObject.h>
#include "TString.h"

using std::vector;

#define ANALYSIS_TYPES	7
#define MAXIMUM_NUMBER_OF_STEPS	1024

class TGraphErrors;
class TH1D;
class TH2D;

const TString kBFAnalysisType[ANALYSIS_TYPES] = {"y","eta","qlong","qout","qside","qinv","phi"};

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

  void SetCentralityIdentifier(const char* centralityId) {
    fCentralityId = centralityId;}
  
  void SetAnalysisLevel(const char* analysisLevel) {
    fAnalysisLevel = analysisLevel;}
  void SetShuffle(Bool_t shuffle) {fShuffle = shuffle;}
  void SetHBTcut(Bool_t HBTcut) {fHBTcut = HBTcut;}
  void SetConversionCut(Bool_t ConversionCut) {fConversionCut = ConversionCut;}
  void SetInterval(Int_t iAnalysisType, Double_t p1Start, Double_t p1Stop,
		   Int_t ibins, Double_t p2Start, Double_t p2Stop);
  void SetCentralityInterval(Double_t cStart, Double_t cStop)  { fCentStart = cStart; fCentStop = cStop;};
  void SetNp(Int_t analysisType, Double_t NpSet)   { fNp[analysisType] = NpSet; }
  void SetNn(Int_t analysisType, Double_t NnSet)   { fNn[analysisType] = NnSet; }
  void SetNpp(Int_t analysisType, Int_t ibin, Double_t NppSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNpp[analysisType][ibin] = NppSet; }
  void SetNpn(Int_t analysisType, Int_t ibin, Double_t NpnSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNpn[analysisType][ibin] = NpnSet; }
  void SetNnp(Int_t analysisType, Int_t ibin, Double_t NnpSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNnp[analysisType][ibin] = NnpSet; }
  void SetNnn(Int_t analysisType, Int_t ibin, Double_t NnnSet) { if(ibin > -1 && ibin < MAXIMUM_NUMBER_OF_STEPS) fNnn[analysisType][ibin] = NnnSet; }

  void InitHistograms(void);

  const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  Int_t GetNumberOfAnalyzedEvent()  const {return fAnalyzedEvents;}

  Int_t GetNumberOfBins(Int_t ibin) const {return fNumberOfBins[ibin];}
  Double_t GetP1Start(Int_t ibin)   const {return fP1Start[ibin];}
  Double_t GetP1Stop(Int_t ibin)    const {return fP1Stop[ibin];}   
  Double_t GetP2Start(Int_t ibin)   const {return fP2Start[ibin];}
  Double_t GetP2Stop(Int_t ibin)    const {return fP2Stop[ibin];}    
 
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

  void CalculateBalance(Float_t fCentrality, vector<Double_t> **chargeVector,Float_t bSign = 0.);
  
  Double_t GetBalance(Int_t a, Int_t p2);
  Double_t GetError(Int_t a, Int_t p2);

  TH2D *GetHistNp(Int_t iAnalysisType)  const { return fHistP[iAnalysisType];}
  TH2D *GetHistNn(Int_t iAnalysisType)  const { return fHistN[iAnalysisType];}
  TH2D *GetHistNpn(Int_t iAnalysisType) const { return fHistPN[iAnalysisType];}
  TH2D *GetHistNnp(Int_t iAnalysisType) const { return fHistNP[iAnalysisType];}
  TH2D *GetHistNpp(Int_t iAnalysisType) const { return fHistPP[iAnalysisType];}
  TH2D *GetHistNnn(Int_t iAnalysisType) const { return fHistNN[iAnalysisType];}

  TH2D *GetQAHistHBTbefore()         {return fHistHBTbefore;};
  TH2D *GetQAHistHBTafter()          {return fHistHBTafter;};
  TH2D *GetQAHistConversionbefore()  {return fHistConversionbefore;};
  TH2D *GetQAHistConversionafter()   {return fHistConversionafter;};

  void PrintAnalysisSettings();
  TGraphErrors *DrawBalance(Int_t fAnalysisType);

  void SetHistNp(Int_t iAnalysisType, TH2D *gHist) { 
    fHistP[iAnalysisType] = gHist;}
  void SetHistNn(Int_t iAnalysisType, TH2D *gHist) { 
    fHistN[iAnalysisType] = gHist;}
  void SetHistNpn(Int_t iAnalysisType, TH2D *gHist) { 
    fHistPN[iAnalysisType] = gHist;}
  void SetHistNnp(Int_t iAnalysisType, TH2D *gHist) { 
    fHistNP[iAnalysisType] = gHist;}
  void SetHistNpp(Int_t iAnalysisType, TH2D *gHist) { 
    fHistPP[iAnalysisType] = gHist;}
  void SetHistNnn(Int_t iAnalysisType, TH2D *gHist) { 
    fHistNN[iAnalysisType] = gHist;}

  TH1D *GetBalanceFunctionHistogram(Int_t iAnalysisType,Double_t centrMin, Double_t centrMax, Double_t etaWindow = -1, Bool_t correctWithEfficiency = kFALSE);
  void PrintResults(Int_t iAnalysisType, TH1D *gHist);

 private:
  inline Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign); 

  Bool_t fShuffle; // shuffled balance function object
  Bool_t fHBTcut;  // apply HBT like cuts
  Bool_t fConversionCut;  // apply conversion cuts

  TString fAnalysisLevel; //ESD, AOD or MC
  Int_t fAnalyzedEvents; //number of events that have been analyzed

  TString fCentralityId;//Centrality identifier to be used for the histo naming

  Int_t fNumberOfBins[ANALYSIS_TYPES];//number of bins of the analyzed interval
  Double_t fP1Start[ANALYSIS_TYPES];//lower boundaries for single particle histograms 
  Double_t fP1Stop[ANALYSIS_TYPES];//upper boundaries for single particle histograms 
  Double_t fP2Start[ANALYSIS_TYPES];//lower boundaries for pair histograms 
  Double_t fP2Stop[ANALYSIS_TYPES];//upper boundaries for pair histograms 
  Double_t fP2Step[ANALYSIS_TYPES];//bin size for pair histograms 
  Double_t fCentStart;//lower boundary for centrality
  Double_t fCentStop;//upper boundary for centrality

 	
  Double_t fNnn[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //N(--)
  Double_t fNpp[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //N(++)
  Double_t fNpn[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //N(+-)
  Double_t fNnp[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //N(-+)
  Double_t fNn[ANALYSIS_TYPES], fNp[ANALYSIS_TYPES]; //number of pos./neg. inside the analyzed interval
  
  Double_t fB[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //BF matrix
  Double_t ferror[ANALYSIS_TYPES][MAXIMUM_NUMBER_OF_STEPS]; //error of the BF
  
  TH2D *fHistP[ANALYSIS_TYPES]; //N+
  TH2D *fHistN[ANALYSIS_TYPES]; //N-
  TH2D *fHistPN[ANALYSIS_TYPES]; //N+-
  TH2D *fHistNP[ANALYSIS_TYPES]; //N-+
  TH2D *fHistPP[ANALYSIS_TYPES]; //N++
  TH2D *fHistNN[ANALYSIS_TYPES]; //N--

  //QA histograms
  TH2D *fHistHBTbefore; // Delta Eta vs. Delta Phi before HBT inspired cuts
  TH2D *fHistHBTafter; // Delta Eta vs. Delta Phi after HBT inspired cuts
  TH2D *fHistConversionbefore; // Delta Eta vs. Delta Phi before Conversion cuts
  TH2D *fHistConversionafter; // Delta Eta vs. Delta Phi before Conversion cuts



  AliBalance & operator=(const AliBalance & ) {return *this;}

  ClassDef(AliBalance, 3)
};

Float_t AliBalance::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign)
{ 
  //
  // calculates dphistar
  //
  
  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
  static const Double_t kPi = TMath::Pi();
  
  // circularity
//   if (dphistar > 2 * kPi)
//     dphistar -= 2 * kPi;
//   if (dphistar < -2 * kPi)
//     dphistar += 2 * kPi;
  
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;
  
  return dphistar;
}

#endif
