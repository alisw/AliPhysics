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

  TList* fListOfHistos; // Container for output histos

  TH1F* fT0vtxRec;      // T0 reconstructed z vertex
  TH2F* fT0vtxRecGen;   // T0 reconstructed vs generate z vertex
  TH1F* fT0time;        // T0 time0
  TH1F* fT0time2;       // T0 time0
  TH1F* fT0mult;        // T0 multiplicity
  TH1F* fT0vtxRes;      // T0 z vertex resolution
  TH1F* fT0ampl;        // T0 signals amplitude

  TH1F* fV0a;           // V0 number of fired PMs A side
  TH1F* fV0c;           // V0 number of fired PMs C side
  TH1F* fV0multA;       // V0 multiplicity on A side
  TH1F* fV0multC;       // V0 multiplicity on C side
  TH2F* fV0multAcorr;   // V0 reconstructed vs generated multiplicity on A side
  TH2F* fV0multCcorr;   // V0 reconstructed vs generated multiplicity on C side
  TH2F* fV0Acorr;       // V0 number of fired PMs (reco vs gen) A side
  TH2F* fV0Ccorr;       // V0 number of fired PMs (reco vs gen) C side
  TH1F* fV0ampl;        // V0 multiplicity in single channel

  AliAnaFwdDetsQA(const AliAnaFwdDetsQA&); // not implemented
  AliAnaFwdDetsQA& operator=(const AliAnaFwdDetsQA&); // not implemented

  ClassDef(AliAnaFwdDetsQA, 1) // example of analysis
};

#endif
