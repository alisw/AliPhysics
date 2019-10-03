#ifndef ALIEVENTPLANERESOLUTIONHANDLER_H
#define ALIEVENTPLANERESOLUTIONHANDLER_H


/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to compute event plane resolution starting from root    //
// file with event plane correlation histos                      //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TObject.h"

class TFile;
class TH1F;

class AliEventPlaneResolutionHandler : public TObject{
 public:
  AliEventPlaneResolutionHandler();
  AliEventPlaneResolutionHandler(TString filename);
  virtual ~AliEventPlaneResolutionHandler(){};

  void SetRootInputFile(TString filn){
    fRootFileName=filn.Data();
  }
  void SetCentralityRange(Double_t minC, Double_t maxC){
    fMinCent=minC;
    fMaxCent=maxC;
  }
  void SetUseNoWeights(){fWeight=0;}
  void SetUseNcollWeights(){fWeight=1;}
  void SetUseLowPtDYieldWeights(){fWeight=2;}
  void SetUseMidPtDYieldWeights(){fWeight=3;}
  void SetUseHighPtDYieldWeights(){fWeight=4;}
  void SetEventPlane(Int_t ep){fEventPlane=ep;}
  void SetTPCFullEtaEventPlane(){fEventPlane=kTPCFullEta;}
  void SetTPCPositiveEtaEventPlane(){fEventPlane=kTPCPosEta;}
  void SetVZEROFullEventPlane(){fEventPlane=kVZERO;}
  void SetVZEROAEventPlane(){fEventPlane=kVZEROA;}
  void SetVZEROCEventPlane(){fEventPlane=kVZEROC;}

  void SetResolutionOption(Int_t res){fResolOption=res;}
  void SetResolutionFromTwoRandomSubevents(){fResolOption=kTwoRandSub;}
  void SetResolutionFromTwoChargeSubevents(){fResolOption=kTwoChargeSub;}
  void SetResolutionFromTwoEtaSubevents(){fResolOption=kTwoEtaSub;}
  void SetResolutionFromThreeSubevents(){fResolOption=kThreeSub;}
  void SetResolutionFromThreeSubeventsTPCGap(){fResolOption=kThreeSubTPCGap;}

  Double_t GetEventPlaneResolution(){
    return GetEventPlaneResolution(fMinCent,fMaxCent);
  }
  Double_t GetEventPlaneResolution(Double_t minCent, Double_t maxCent);

  TH1F* GetSubEvABCorrelHisto() const{
    return fHistoAB;
  }
  TH1F* GetSubEvBCCorrelHisto() const{
    return fHistoBC;
  }
  TH1F* GetSubEvACCorrelHisto() const{
    return fHistoAC;
  }


  enum EEventPlane {kTPCFullEta, kTPCPosEta,kVZERO,kVZEROA,kVZEROC};
  enum EResolOption {kTwoRandSub,kTwoChargeSub,kTwoEtaSub,kThreeSub,kThreeSubTPCGap};

 private:
  Bool_t SetHistoNames();
  void InitializeWeights();
  AliEventPlaneResolutionHandler(const AliEventPlaneResolutionHandler &source);
  AliEventPlaneResolutionHandler& operator=(const AliEventPlaneResolutionHandler& source); 
 


  Int_t    fEventPlane;      // event plane determination (see enum)
  Int_t    fResolOption;     // option to get the resolution (see enum)
  Double_t fMinCent;         // minimum centrality
  Double_t fMaxCent;         // maximum centrality
  TString  fCorrHistoName1;  // name of 1st histogram
  TString  fCorrHistoName2;  // name of 2nd histogram
  TString  fCorrHistoName3;  // name of 3rd histogram
  Int_t    fNsubevents;      // number of subevents to be used
  Bool_t   fExtrapToFull;    // flag for extrapolating to full event
  Int_t    fWeight; // flag use/not use Ncoll to weight centrality sub-ranges
  Double_t fNcoll[20];       // values of Ncoll in the subranges
  Double_t fYield24[20];     // values of D0 yield 2<pt<4 in the subranges
  Double_t fYield46[20];     // values of D0 yield 4<pt<6 in the subranges
  Double_t fYield612[20];    // values of D0 yield 6<pt<12 in the subranges
  TH1F*    fHistoAB;         // histo of subevents A-B   
  TH1F*    fHistoBC;         // histo of subevents B-C   
  TH1F*    fHistoAC;         // histo of subevents A-C   
  TString  fRootFileName;    // file with histos of correlations

  ClassDef(AliEventPlaneResolutionHandler,0) 
};
#endif
