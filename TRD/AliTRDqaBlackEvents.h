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
class TH3D;
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

 private:
  
  Int_t fnEvents;         // number of events processed
  
  TH1D *fOccupancy;       // how many times is a pad present in data

  TH1D *fPed[540];        // reconstructed pedestals distribution (on hist per chamber)
  TH1D *fNoise[540];      // reconstructed noise distribution (on hist per chamber)
  TH1D *fNPointDist[540]; // distributin of the number of points
  TH2D *fChPed[540];      // Some histograms
  TH2D *fChNoise[540];    // Some histograms
  TH2D *fNPoint[540];     // number of data points
  TH3D *fData[540];       // Some histograms
  TH1D *fSignal[540];     // Some histograms

  Int_t fFitType;
  Double_t fMinNoise;   // Minimum noise
  Double_t fMaxNoise;   // Maximum noise

  ClassDef(AliTRDqaBlackEvents,0) // QA for black events  

};
#endif
