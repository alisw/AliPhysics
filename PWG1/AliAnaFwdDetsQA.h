#ifndef ALIANAFWDDETSQA_H
#define ALIANAFWDDETSQA_H

//------------------------------
// Analysis task for quality-assurance
// of forward detectors ESD
//
// 12/06/2009 cvetan.cheshkov@cern.ch
//------------------------------


class TH1;
class TH1F;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AliAnaFwdDetsQA : public AliAnalysisTaskSE
{
 public:	 
  AliAnaFwdDetsQA();
  AliAnaFwdDetsQA(const char *name);
  virtual ~AliAnaFwdDetsQA() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  TH1F* CreateHisto(const char* name, const char* title, Int_t nBins, Double_t xMin, 
			Double_t xMax, const char* xLabel = NULL, const char* yLabel = NULL);
  TH1F* CreateEffHisto(const TH1F* hGen, const TH1F* hRec);
  Bool_t FitHisto(TH1* histo, Double_t& res, Double_t& resError);
  
 private:

  TList* fListOfHistos;

  TH1F* fT0vtxRec;
  TH2F* fT0vtxRecGen;
  TH1F* fT0time;
  TH1F* fT0time2;
  TH1F* fT0mult;
  TH1F* fT0vtxRes;
  TH1F* fT0ampl;

  TH1F* fV0a;
  TH1F* fV0c;
  TH1F* fV0multA;
  TH1F* fV0multC;
  TH2F* fV0multAcorr;
  TH2F* fV0multCcorr;
  TH2F* fV0Acorr;
  TH2F* fV0Ccorr;
  TH1F* fV0ampl;

  AliAnaFwdDetsQA(const AliAnaFwdDetsQA&); // not implemented
  AliAnaFwdDetsQA& operator=(const AliAnaFwdDetsQA&); // not implemented

  ClassDef(AliAnaFwdDetsQA, 1) // example of analysis
};

#endif
