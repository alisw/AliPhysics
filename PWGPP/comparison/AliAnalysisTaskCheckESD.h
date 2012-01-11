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

  TList* fListOfHistos;   // list of the histos

  TH1F * fGen;            // generated tracks
  TH1F* fRec;             // reconstructed tracks
  TH1F* fResPtInv;        // 1/pt resolution
  TH1F* fResPhi;          // phi resolution
  TH1F* fResTheta;        // theta resolution
  TH2F* fDEdxRight;       // dE/dx in TPC
  TH2F* fDEdxWrong;       // 
  TH1F* fResTOFRight;     // TOF resolution
  TH1F* fResTOFWrong;     //
  TH1F* fEPHOS;           // PHOS E
  TH1F* fEEMCAL;          // EMCAL E
  TH1F* fPtMUON;          // MUON pt 
  TH1F* fMassK0;          // K0s mass
  TH1F* fMassLambda;      // Lambda mass
  TH1F* fMassLambdaBar;   // anti-Lambda mass
  TH1F* fMassXi;          // Xi mass
  TH1F* fMassOmega;       // Omega mass 
  TH1F* fScalars;         // container of scalars
  TH1F* fArrayHist;       // container of arrays
   
  AliAnalysisTaskCheckESD(const AliAnalysisTaskCheckESD&); // not implemented
  AliAnalysisTaskCheckESD& operator=(const AliAnalysisTaskCheckESD&); // not implemented

  ClassDef(AliAnalysisTaskCheckESD, 1); // example of analysis
};

#endif
