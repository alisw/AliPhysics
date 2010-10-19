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

#define MAXIMUM_NUMBER_OF_STEPS	1024

class TGraphErrors;
class TObjArray;
class TH1F;

class AliBalance : public TObject {
 public:
  AliBalance();
  AliBalance(Double_t p2Start, Double_t p2Stop, Int_t p2Steps);
  AliBalance(const AliBalance& balance);
  ~AliBalance();
  
  void SetNumberOfBins(Int_t ibins);
  void SetAnalysisLevel(const char* analysisLevel) {fAnalysisLevel = analysisLevel;}
  void SetAnalysisType(Int_t iType);
  void SetInterval(Double_t p2Start, Double_t p2Stop);

  Int_t GetNumberOfBins() {return fNumberOfBins;}
  const char* GetAnalysisLevel() {return fAnalysisLevel.Data();}
  Int_t GetNumberOfAnalyzedEvent() {return fAnalyzedEvents;}

  void SetNp(Int_t np) {fNp = np;}
  void SetNn(Int_t nn) {fNn = nn;}
  void SetNnn(Double_t *nn);
  void SetNpp(Double_t *pp);
  void SetNpn(Double_t *pn);
 
  Double_t GetNp() const { return 1.0*fNp; }
  Double_t GetNn() const { return 1.0*fNn; }
  Double_t GetNnn(Int_t p2) const { return 1.0*fNnn[p2]; }
  Double_t GetNpp(Int_t p2) const { return 1.0*fNpp[p2]; }
  Double_t GetNpn(Int_t p2) const { return 1.0*fNpn[p2]; }
 
  TH1F *GetHistNnn();
  TH1F *GetHistNpp();
  TH1F *GetHistNpn();

  void CalculateBalance(TObjArray *gTrackArray);
  
  Double_t GetBalance(Int_t p2);
  Double_t GetError(Int_t p2);

  TGraphErrors *DrawBalance();
  void PrintResults();
  void PrintAnalysisSettings();

  void Merge(AliBalance *b);

 private:
  TString fAnalysisLevel; //ESD, AOD or MC
  Int_t fNumberOfBins; //number of bins of the analyzed interval
  Int_t fAnalysisType; //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
  Int_t fAnalyzedEvents; //number of events that have been analyzed
 
  Double_t fP2Start, fP2Stop, fP2Step; //inerval info
 	
  Double_t fNnn[MAXIMUM_NUMBER_OF_STEPS]; //N(--)
  Double_t fNpp[MAXIMUM_NUMBER_OF_STEPS]; //N(++)
  Double_t fNpn[MAXIMUM_NUMBER_OF_STEPS]; //N(+-)
  Double_t fNn, fNp; //number of pos./neg. inside the analyzed interval
  
  Double_t fB[MAXIMUM_NUMBER_OF_STEPS]; //BF matrix
  Double_t ferror[MAXIMUM_NUMBER_OF_STEPS]; //error of the BF
  
  TH1F *fHistfNnn; //N(--) in a histo
  TH1F *fHistfNpp; //N(++) in a histo
  TH1F *fHistfNpn; //N(+-) in a histo

  AliBalance & operator=(const AliBalance & ) {return *this;}

  ClassDef(AliBalance, 1)
};

#endif
