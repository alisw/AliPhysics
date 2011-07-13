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

#define NUMBER_OF_ANALYSES	7
#define MAXIMUM_NUMBER_OF_STEPS	1024

class TGraphErrors;
class TObjArray;
class TH1F;

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
  
  void SetNumberOfBins(Int_t ibin, Int_t ibins);
  void SetAnalysisLevel(const char* analysisLevel) {fAnalysisLevel = analysisLevel;}
  void SetInterval(Int_t ibin, Double_t p2Start, Double_t p2Stop);

  Int_t GetNumberOfBins(Int_t ibin) {return fNumberOfBins[ibin];}
  Double_t GetP2Start(Int_t ibin){return fP2Start[ibin];}
  Double_t GetP2Stop(Int_t ibin){return fP2Stop[ibin];}    
  const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  Int_t GetNumberOfAnalyzedEvent() {return fAnalyzedEvents;}
 
  Double_t GetNp() const { return 1.0*fNp; }
  Double_t GetNn() const { return 1.0*fNn; }
  Double_t GetNnn(Int_t a, Int_t p2) const { return 1.0*fNnn[a][p2]; }
  Double_t GetNpp(Int_t a, Int_t p2) const { return 1.0*fNpp[a][p2]; }
  Double_t GetNpn(Int_t a, Int_t p2) const { return 1.0*fNpn[a][p2]; }

  void CalculateBalance(TObjArray *gTrackArray);
  
  Double_t GetBalance(Int_t a, Int_t p2);
  Double_t GetError(Int_t a, Int_t p2);

  void PrintAnalysisSettings();

 private:
  TString fAnalysisLevel; //ESD, AOD or MC
  Int_t fAnalyzedEvents; //number of events that have been analyzed
 
  Int_t fNumberOfBins[NUMBER_OF_ANALYSES]; //number of bins of the analyzed interval
  Double_t fP2Start[NUMBER_OF_ANALYSES];
  Double_t fP2Stop[NUMBER_OF_ANALYSES];
  Double_t fP2Step[NUMBER_OF_ANALYSES]; 
 	
  Double_t fNnn[NUMBER_OF_ANALYSES][MAXIMUM_NUMBER_OF_STEPS]; //N(--)
  Double_t fNpp[NUMBER_OF_ANALYSES][MAXIMUM_NUMBER_OF_STEPS]; //N(++)
  Double_t fNpn[NUMBER_OF_ANALYSES][MAXIMUM_NUMBER_OF_STEPS]; //N(+-)
  Double_t fNn, fNp; //number of pos./neg. inside the analyzed interval
  
  Double_t fB[NUMBER_OF_ANALYSES][MAXIMUM_NUMBER_OF_STEPS]; //BF matrix
  Double_t ferror[NUMBER_OF_ANALYSES][MAXIMUM_NUMBER_OF_STEPS]; //error of the BF
  
  AliBalance & operator=(const AliBalance & ) {return *this;}

  ClassDef(AliBalance, 3)
};

#endif
