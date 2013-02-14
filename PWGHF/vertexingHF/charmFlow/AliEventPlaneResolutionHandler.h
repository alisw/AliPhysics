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

  void SetCentralityRange(Double_t minC, Double_t maxC){
    fMinCent=minC;
    fMaxCent=maxC;
  }
  void SetUseNcollWeights(Bool_t opt){fUseNcollWeights=opt;}
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



  enum EEventPlane {kTPCFullEta, kTPCPosEta,kVZERO,kVZEROA,kVZEROC};
  enum EResolOption {kTwoRandSub,kTwoChargeSub,kTwoEtaSub,kThreeSub,kThreeSubTPCGap};

 private:
  Bool_t SetHistoNames();
  void InitializeNcoll();
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
  Bool_t   fUseNcollWeights; // flag use/not use Ncoll to weight centrality sub-ranges
  Double_t fNcoll[20];       // values of Ncoll in the subranges
  TString  fRootFileName;    // file with histos of correlations

  ClassDef(AliEventPlaneResolutionHandler,0) 
};
#endif
