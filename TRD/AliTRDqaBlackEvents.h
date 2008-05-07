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
class AliTRDrawStreamTB;

class AliTRDqaBlackEvents : public TObject {

 public:
  
  AliTRDqaBlackEvents();
  AliTRDqaBlackEvents(const AliTRDqaBlackEvents &qa);
  ~AliTRDqaBlackEvents() {}
  AliTRDqaBlackEvents& operator = (const AliTRDqaBlackEvents& /*qa*/) { return *this; };

  void Init();
  void Reset();
  Int_t AddEvent(AliTRDrawStreamTB *data);
  void Process(const char* filename);
  
  //TH2D *GetChamberPedestal(Int_t sm, Int_t layer, Int_t stack) {return 0;}
  TH2D *GetChamberPedestal(Int_t det) {return fChPed[det];}
  
  //TH2D *GetChamberNoise(Int_t sm, Int_t layer, Int_t stack) {return 0;}
  TH2D *GetChamberNoise(Int_t det) {return fChNoise[det];}
  
  void SetNoiseLevel(Double_t min, Double_t max) {fMinNoise = min; fMaxNoise = max;}
  void SetFitMethod(Int_t fit) {fFitType = fit;} 

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
  
  // geometry constants 
  enum {
    kDET = 540,
    kROB = 8,
    kMCM = 16,
    kADC = 21,
    kTB  = 30,
    kCOL = 16,
    kPAD = 144
  };

  // histograms per detector

  TH1D *fOccupancy;       // how many times is a pad present in data
  TH2D *fDetRob;          // detector -- read out board

  // histograms per chamber

  TH1D *fPed[kDET];        // reconstructed pedestals distribution (on hist per chamber)
  TH1D *fNoise[kDET];      // reconstructed noise distribution (on hist per chamber)
  TH1D *fNPointDist[kDET]; // distributin of the number of points
  TH2D *fChPed[kDET];      // Some histograms
  TH2D *fChNoise[kDET];    // Some histograms
  TH2D *fNPoint[kDET];     // number of data points
  TH3F *fData[kDET];       // Some histograms
  TH1D *fSignal[kDET];     // Some histograms
  TH2D *fnEntriesRM[kDET];     // number of entries for ROB - MCM
  TH1D *fnEntriesRMDist[kDET]; // distribtion of number of entries per ROB-MCM

  TH2S *fFullSignal[kDET*kROB*kMCM];     // one histogram per MCM  
  Short_t fFullCounter[kDET*kROB*kMCM];  // counts a number of entries with high signal
  
  TH2D *fTBEvent;    // coherent noise


  Int_t fFitType;

  Double_t fMinNoise;   // Minimum noise
  Double_t fMaxNoise;   // Maximum noise

  ClassDef(AliTRDqaBlackEvents,0) // QA for black events  

};
#endif
