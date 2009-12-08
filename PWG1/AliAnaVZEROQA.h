#ifndef ALIANAVZEROQA_H
#define ALIANAVZEROQA_H

//------------------------------
// Analysis task for quality-assurance
// of VZERO ESD
//
// 05/12/2009 cvetan.cheshkov@cern.ch
//------------------------------

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

  Float_t CorrectLeadingTime(Int_t i, Float_t time, Float_t adc);
  
 private:

  TList* fListOfHistos;

  TH1F *fhAdcNoTimeA;
  TH1F *fhAdcWithTimeA;
  TH1F *fhAdcNoTimeC;
  TH1F *fhAdcWithTimeC;

  TH2F *fhAdcPMTNoTime;
  TH2F *fhAdcPMTWithTime;
 
  TH1F *fhTimeA;
  TH1F *fhTimeC;

  TH1F *fhWidthA;
  TH1F *fhWidthC;

  TH2F *fhTimePMT;
  TH2F *fhWidthPMT;

  TH2F *fhAdcWidthA;
  TH2F *fhAdcWidthC;

  TH2F *fhTimeCorr;

  TH2F *fhAdcTimeA;
  TH2F *fhAdcTimeC;

  TH1F *fV0a;
  TH1F *fV0c;
  TH1F *fV0multA;
  TH1F *fV0multC;
  TH1F *fV0ampl;

  TH2F *fhTimePMTCorr;
  TH2F *fhEvents;

  TH2F *fhVtxXYBB;
  TH1F *fhVtxZBB;
  TH2F *fhVtxXYBGA;
  TH1F *fhVtxZBGA;
  TH2F *fhVtxXYBGC;
  TH1F *fhVtxZBGC;

  AliAnaVZEROQA(const AliAnaVZEROQA&); // not implemented
  AliAnaVZEROQA& operator=(const AliAnaVZEROQA&); // not implemented

  ClassDef(AliAnaVZEROQA, 1) // VZERO QA task
};

#endif
