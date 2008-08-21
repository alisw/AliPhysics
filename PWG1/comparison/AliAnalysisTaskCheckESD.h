#ifndef AliAnalysisTaskCheckESD_H
#define AliAnalysisTaskCheckESD_H

class TH1F;
class TArrayI;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckESD : public AliAnalysisTaskSE
{
 public:	 
  AliAnalysisTaskCheckESD();
  AliAnalysisTaskCheckESD(const char *name);
  virtual ~AliAnalysisTaskCheckESD() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  TH1F* CreateHisto(const char* name, const char* title, Int_t nBins, Double_t xMin, 
			Double_t xMax, const char* xLabel = NULL, const char* yLabel = NULL);
  TH1F* CreateEffHisto(TH1F* hGen, TH1F* hRec);
  Bool_t FitHisto(TH1* histo, Double_t& res, Double_t& resError);
  
 private:

  TList* fListOfHistos;

  TH1F * hGen;
  TH1F* hRec;
  TH1F* hResPtInv;
  TH1F* hResPhi;
  TH1F* hResTheta;
  TH2F* hDEdxRight;
  TH2F* hDEdxWrong;
  TH1F* hResTOFRight;
  TH1F* hResTOFWrong;
  TH1F* hEPHOS;
  TH1F* hEEMCAL;
  TH1F* hPtMUON;
  TH1F* hMassK0;
  TH1F* hMassLambda;
  TH1F* hMassLambdaBar;
  TH1F* hMassXi;
  TH1F* hMassOmega;
  TH1F* hScalars;
  TH1F* hArrayHist;  
   
  AliAnalysisTaskCheckESD(const AliAnalysisTaskCheckESD&); // not implemented
  AliAnalysisTaskCheckESD& operator=(const AliAnalysisTaskCheckESD&); // not implemented

  ClassDef(AliAnalysisTaskCheckESD, 1); // example of analysis
};

#endif
