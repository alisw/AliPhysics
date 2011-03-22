#ifndef ALIANAVZEROQA_H
#define ALIANAVZEROQA_H

//------------------------------
// Analysis task for quality-assurance
// of VZERO ESD
//
// 05/12/2009 cvetan.cheshkov@cern.ch
//------------------------------

class TH1F;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AliAnaVZEROQA : public AliAnalysisTaskSE
{
 public:	 
  AliAnaVZEROQA();
  AliAnaVZEROQA(const char *name);
  virtual ~AliAnaVZEROQA() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  TH1F* CreateHisto1D(const char* name, const char* title, Int_t nBins, Double_t xMin, Double_t xMax,
		      const char* xLabel = NULL, const char* yLabel = NULL);
  TH2F* CreateHisto2D(const char* name, const char* title, Int_t nBinsX, Double_t xMin, Double_t xMax,
		      Int_t nBinsY, Double_t yMin, Double_t yMax,
		      const char* xLabel = NULL, const char* yLabel = NULL);

 private:

  TList* fListOfHistos;   // List of output histos

  TH1F *fhAdcNoTimeA;     // ADC spectra (no time measurement) for V0A
  TH1F *fhAdcWithTimeA;   // ADC spectra (with time measurement) for V0A
  TH1F *fhAdcNoTimeC;     // ADC spectra (no time measurement) for V0C
  TH1F *fhAdcWithTimeC;   // ADC spectra (with time measurement) for V0C

  TH2F *fhAdcPMTNoTime;   // ADC spectra per PMT (no time measurement)
  TH2F *fhAdcPMTWithTime; // ADC spectra per PMT (with time measurement)
 
  TH1F *fhTimeA;          // Time spectra for V0A
  TH1F *fhTimeC;          // Time spectra for V0C

  TH1F *fhWidthA;         // Signal width for V0A
  TH1F *fhWidthC;         // Signal width for V0C

  TH2F *fhTimePMT;        // Time spectra per PMT
  TH2F *fhWidthPMT;       // Signal width per PMT

  TH2F *fhAdcWidthA;      // ADC vs Signal width for V0A
  TH2F *fhAdcWidthC;      // ADC vs Signal width for V0C

  TH2F *fhTimeCorr;       // Corrected mean time V0C vs V0A

  TH2F *fhAdcTimeA;       // ADC vs Time for V0A
  TH2F *fhAdcTimeC;       // ADC vs Time for V0C

  TH1F *fV0a;             // Number of fired PMTs in V0A
  TH1F *fV0c;             // Number of fired PMTs in V0C
  TH1F *fV0multA;         // Mutiplicity in V0A
  TH1F *fV0multC;         // Mutiplicity in V0C
  TH1F *fV0ampl;          // ADC spectra for both rings

  TH2F *fhEvents;         // Event statistics histogram

  TH2F *fhVtxXYBB;        // XY vertex for beam-beam events
  TH1F *fhVtxZBB;         // Z vertex for beam-beam events
  TH2F *fhVtxXYBGA;       // XY vertex for beam-gas (A side) events
  TH1F *fhVtxZBGA;        // Z vertex for beam-gas (A side) events
  TH2F *fhVtxXYBGC;       // XY vertex for beam-gas (C side) events
  TH1F *fhVtxZBGC;        // Z vertex for beam-gas (C side) events

  AliAnaVZEROQA(const AliAnaVZEROQA&); // not implemented
  AliAnaVZEROQA& operator=(const AliAnaVZEROQA&); // not implemented

  ClassDef(AliAnaVZEROQA, 1) // VZERO QA task
};

#endif
