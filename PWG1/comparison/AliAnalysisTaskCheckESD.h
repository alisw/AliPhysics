#ifndef ALIANALYSISTASKCHECKESD_H
#define ALIANALYSISTASKCHECKESD_H

//------------------------------
// Proof-enabled 
// version of CheckESD.C
//------------------------------

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
  TH1F* CreateEffHisto(const TH1F* hGen, const TH1F* hRec);
  Bool_t FitHisto(TH1* histo, Double_t& res, Double_t& resError);
  
 private:

  TList* fListOfHistos;

  TH1F * fGen;
  TH1F* fRec;
  TH1F* fResPtInv;
  TH1F* fResPhi;
  TH1F* fResTheta;
  TH2F* fDEdxRight;
  TH2F* fDEdxWrong;
  TH1F* fResTOFRight;
  TH1F* fResTOFWrong;
  TH1F* fEPHOS;
  TH1F* fEEMCAL;
  TH1F* fPtMUON;
  TH1F* fMassK0;
  TH1F* fMassLambda;
  TH1F* fMassLambdaBar;
  TH1F* fMassXi;
  TH1F* fMassOmega;
  TH1F* fScalars;
  TH1F* fArrayHist;  
   
  AliAnalysisTaskCheckESD(const AliAnalysisTaskCheckESD&); // not implemented
  AliAnalysisTaskCheckESD& operator=(const AliAnalysisTaskCheckESD&); // not implemented

  ClassDef(AliAnalysisTaskCheckESD, 1); // example of analysis
};

#endif
