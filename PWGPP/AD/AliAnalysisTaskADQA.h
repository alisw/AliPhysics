#ifndef ALIANALYSISTASKADQA_H
#define ALIANALYSISTASKADQA_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskADQA class
//            This task is for QAing the AD data from ESD/AOD
//              Origin: Michal Broz
//-----------------------------------------------------------------

class TString;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TF1;

#include "AliAnalysisTaskSE.h"

class AliESDEvent;


class AliAnalysisTaskADQA : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskADQA();
  AliAnalysisTaskADQA(const char *name);
 ~AliAnalysisTaskADQA();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  TH1F* CreateHist1D(const char* name, const char* title, Int_t nBins, Double_t xMin, Double_t xMax,
		      const char* xLabel = NULL, const char* yLabel = NULL);
  TH2F* CreateHist2D(const char* name, const char* title, Int_t nBinsX, Double_t xMin, Double_t xMax,
		      Int_t nBinsY, Double_t yMin, Double_t yMax,
		      const char* xLabel = NULL, const char* yLabel = NULL);

private:

  TList       *fListHist;                       //! List of histograms
 
  TH1F        *fHistTotalChargePerEventADA;
  TH1F        *fHistTotalChargePerEventADC;
  TH2F        *fHistChargePerPM;
  TH2F        *fHistTimePerPM;
  TH2F        *fHistWidthPerPM;
  TH2F	      *fHistTimeVsCharge;
  TH2F	      *fHistWidthVsCharge;
  TH1F        *fHistNBBflagsADA;
  TH1F        *fHistNBBflagsADC;
  TH2F	      *fHistNBBflagsADAVsADC;
  TH1F	      *fHistNCoincidencesADA;
  TH1F	      *fHistNCoincidencesADC;
  TH2F	      *fHistNCoincidencesADAVsADC;
  TH1F	      *fHistChargeNoFlag;
  TH1F	      *fHistTimeNoFlag;
  TH1F	      *fHistChargeNoTime;
  
  
   
  AliAnalysisTaskADQA(const AliAnalysisTaskADQA&);            // not implemented
  AliAnalysisTaskADQA& operator=(const AliAnalysisTaskADQA&); // not implemented
  
  ClassDef(AliAnalysisTaskADQA, 1);
};

#endif
