#ifndef ALITRDQABLACKEVENTS_H
#define ALITRDQABLACKEVENTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaBlackEvents.h 23387 2008-01-17 17:25:16Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  QA of black events                                                    //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TH1D;
class TH2D;
class TH2S;
class TH3F;
class TGraph;
class TObjArray;
class AliTRDrawStreamTB;
class AliRawReader;

class AliTRDqaBlackEvents : public TObject {

 public:
  
  AliTRDqaBlackEvents();
  AliTRDqaBlackEvents(const AliTRDqaBlackEvents &qa);
  ~AliTRDqaBlackEvents() {}
  AliTRDqaBlackEvents& operator = (const AliTRDqaBlackEvents& /*qa*/) { return *this; };

  void Init();
  void Reset();
  //Int_t AddEvent(AliTRDrawStreamTB *data, AliRawReader *reader);

  void StartEvent();
  void AddBuffer(AliTRDrawStreamTB *data, AliRawReader *reader);
  void FinishEvent();

  void Process(const char* filename);
  
  //TH2D *GetChamberPedestal(Int_t sm, Int_t layer, Int_t stack) {return 0;}
  TH2D *GetChamberPedestal(Int_t det) {return fChPed[det];}
  
  //TH2D *GetChamberNoise(Int_t sm, Int_t layer, Int_t stack) {return 0;}
  TH2D *GetChamberNoise(Int_t det) {return fChNoise[det];}
  
  void SetNoiseLevel(Double_t min, Double_t max) {fMinNoise = min; fMaxNoise = max;}
  void SetFitMethod(Int_t fit) {fFitType = fit;} 

  void SetRefFile(const char *filename);

  void DrawChamber(const char *filename, Int_t det, Int_t w=700, Int_t h=400);
  //void ScanChamber(const char *filename, Int_t first, Int_t last);
  void DrawSm(const char *filename, Int_t sm, Int_t w=900, Int_t h=700);

  void SetFullView(Int_t v, Int_t thresh, Int_t count) {
    fCreateFull = v;
    fThresh = thresh;
    fCount = count;
  }

 private:
  
  Int_t fnEvents;         // number of events processed  
  Int_t fCreateFull;      // flag if to create a full view
  Int_t fThresh;          // threshold to analyze MCM data
  Int_t fCount;           // minimum number of entries above threshold
  
  Int_t fRefEv;           // reference event number

  Char_t fRefFileName[256];  // name of the file with reference distributions


  // geometry constants 
  enum {
    kDET = 540,
    kROB = 8,
    kMCM = 16,
    kADC = 21,
    kTB  = 30,
    kROW = 16,
    kPAD = 144,
    kSM  = 18,
    kCH  = 50
  };

  // histograms per detector

  TH1D *fOccupancy;       // how many times is a pad present in data
  TH2D *fDetRob;          // detector -- read out board

  // histograms per chamber

  TH1D *fPed[kDET];        // reconstructed pedestals distribution (on hist per chamber)
  TH1D *fNoise[kDET];      // reconstructed noise distribution (on hist per chamber)
  TH1D *fChPP[kDET];       // peak to peak for each chamber
  TH1D *fNPointDist[kDET]; // distributin of the number of points
  TH2D *fChPed[kDET];      // Some histograms
  TH2D *fChNoise[kDET];    // Some histograms
  TH2D *fNPoint[kDET];     // number of data points
  TH3F *fData[kDET];       // Some histograms
  TH1D *fSignal[kDET];     // Some histograms
  TH2D *fnEntriesRM[kDET];     // number of entries for ROB - MCM
  TH1D *fnEntriesRMDist[kDET]; // distribtion of number of entries per ROB-MCM

  // direct access to data
  Float_t  fDataDirect[kDET][kROW][kPAD][kCH];
  Double_t fSignalDirect[kDET][kCH]; 

  // after reference subtraction
  TH2D *fChPedRes[kDET];    // histograms after reference subtraction
  TH2D *fChNoiseRes[kDET];  // histograms after reference subtraction

  TH2D *fTBEvent;    // coherent noise

  TH2D *fRefHistPed;        // reference distributions
  TH2D *fRefHistNoise;      // reference distributions

  TH2S *fFullSignal[kDET*kROB*kMCM];     // one histogram per MCM  
  Short_t fFullCounter[kDET*kROB*kMCM];  // counts a number of entries with high signal
  
  // error codes
  TH1D *fErrorHC;          // number of errors HC
  TH1D *fErrorMCM;         // number of errors MCM
  TH1D *fErrorADC;         // number of errors ADC
  
  TH1D *fErrorSMHC;        // number of errors in HC per SM
  TH1D *fErrorSMMCM;       // number of errors in MCM per SM
  TH1D *fErrorSMADC;       // number of errors in ADC per SM

  TH2D *fErrorLocHC[kDET];       // location of errors
  TH2D *fErrorLocMCM[kDET];      // location
  TH2D *fErrorLocADC[kDET];      // errors in ADC

  // error fraction
  TGraph *fErrorGraphHC;
  TGraph *fErrorGraphMCM;
  TGraph *fErrorGraphADC;

  TGraph *fGraphMCM;         // number of strange MCMs detected 
  TGraph *fGraphPP[3];
  

  // mcm trackles
  TObjArray *fMcmTracks;

  // problematic MCMs
  TH2D *fMapMCM;
  TH1D *fFracMCM;
  
  // full detector view
  TH2D *fSMHCped;
  TH2D *fSMHCerr;
  TH2D *fSMLink[3];
  TGraph *fGrLink[3];
  
  //TH1D *fZSsize;
  
  
  // number of fired ADC channels in total and per SM
  TGraph *fNumberADC[kSM+1];
  
  //Int_t fChkDe
  
  TH1D *fNoiseTotal;
  TH1D *fPP;
  
  TH1D *fSmNoiseRms[kSM];
  TH1D *fSmNoiseFit[kSM];
  TH1D *fSmPP[kSM];    


  TH1D *fEvNoDist[1000];

  //
  Double_t fMinNoise;   // Minimum noise
  Double_t fMaxNoise;   // Maximum noise
  Int_t fFitType;
  
  // variables keeping info in one event
  Int_t fnErrorHC[2];  // 0 good, 1 error
  Int_t fnErrorMCM[2]; //
  Int_t fnErrorADC[2]; // 
    
  Int_t fppThresh[3];      // thersholds for storing pp
  Int_t fnPP[3];           // number of entries above the thershold
  Int_t fnLink[3];         // links present, beaf-beaf, good
  Int_t fnADCinSM[kSM+1];  // number of ADC channels in a SuperModule
  //



  // private function
  void  ReadRefHists(Int_t det);
  Int_t CheckMCM(Int_t index);
  
  Int_t FillBits(TH1D *hist, Int_t code, Int_t offset);


  ClassDef(AliTRDqaBlackEvents,0) // QA for black events  

};
#endif
