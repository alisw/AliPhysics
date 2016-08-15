#ifndef ALIANALYSISTASKADPILOT_H
#define ALIANALYSISTASKADPILOT_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskADPilot class
//            This task is for QAing the AD data from ESD/AOD
//              Origin: Michal Broz
//-----------------------------------------------------------------

class TString;
class TList;
class TTree;
class TH1F;
class TH2F;
class TH3F;
class TF1;
class AliADCalibData;
class AliAnalysisUtils;
class AliTriggerAnalysis; 

#include "AliAnalysisTaskSE.h"

class AliESDEvent;


class AliAnalysisTaskADPilot : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskADPilot();
  AliAnalysisTaskADPilot(const char *name);
 ~AliAnalysisTaskADPilot();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  TH1F* CreateHist1D(const char* name, const char* title, Int_t nBins, Double_t xMin, Double_t xMax,
		      const char* xLabel = NULL, const char* yLabel = NULL);
  TH2F* CreateHist2D(const char* name, const char* title, Int_t nBinsX, Double_t xMin, Double_t xMax,
		      Int_t nBinsY, Double_t yMin, Double_t yMax,
		      const char* xLabel = NULL, const char* yLabel = NULL, const char* zLabel = NULL);
  TH3F* CreateHist3D(const char* name, const char* title, Int_t nBinsX, Double_t xMin, Double_t xMax,
		      Int_t nBinsY, Double_t yMin, Double_t yMax,
		      Int_t nBinsZ, Double_t zMin, Double_t zMax,
		      const char* xLabel = NULL, const char* yLabel = NULL, const char* zLabel = NULL);

private:

  TList       *fListHist;                       //! List of histograms
  /*From ESD*/
TH1F        *fHistTotalChargePerEventADA;
  TH1F        *fHistTotalChargePerEventADC;
  TH2F        *fHistChargePerPM_All;
  TH2F        *fHistChargePerPM_BB;
  TH2F        *fHistChargePerPM_BG;
  TH2F        *fHistChargePerPM_Time;
  TH2F        *fHistTimePerPM_Corr;
  TH2F	      *fHistTimeVsChargeADA_Corr;
  TH2F	      *fHistTimeVsChargeADC_Corr;
  TH2F	      *fHistTimeVsChargeADA_Cut;
  TH2F	      *fHistTimeVsChargeADC_Cut;
  TH2F        *fHistWidthPerPM;
  TH2F	      *fHistWidthVsCharge;
  TH1F        *fHistNBBflagsADA;
  TH1F        *fHistNBBflagsADC;
  TH2F	      *fHistNBBflagsADAVsADC;
  TH1F        *fHistNBGflagsADA;
  TH1F        *fHistNBGflagsADC;
  TH2F	      *fHistNBGflagsADAVsADC;
  TH1F	      *fHistNBBCoincidencesADA;
  TH1F	      *fHistNBBCoincidencesADC;
  TH2F	      *fHistNBBCoincidencesADAVsADC;
  TH1F	      *fHistNBGCoincidencesADA;
  TH1F	      *fHistNBGCoincidencesADC;
  TH2F	      *fHistNBGCoincidencesADAVsADC;
  TH1F	      *fHistChargeNoFlag;
  TH2F	      *fHistTimeNoFlag;
  TH2F	      *fHistChargeNoTime;
  TH1F	      *fHistFlagNoTime;	
  TH2F	      *fHistChargePerCoincidence;
  
  TH1F	      *fHistMeanTimeADA;
  TH1F	      *fHistMeanTimeADC;
  TH1F	      *fHistMeanTimeDifference;
  TH2F	      *fHistMeanTimeCorrelation;
  TH2F	      *fHistMeanTimeSumDiff;
  TH2F	      *fHistDecision;
  
  TH1F	      *fHistTriggerMasked;
  TH1F	      *fHistTriggerUnMasked;
  TH1F	      *fHistTriggerOthers;
  
  /*From ESD friend*/
  TH2F	      *fHistChargeVsClockInt0;
  TH2F	      *fHistChargeVsClockInt1;
  TH2F	      *fHistBBFlagVsClock;
  TH2F	      *fHistBGFlagVsClock;
  TH2F	      *fHistBBFlagPerChannel;
  TH2F	      *fHistBGFlagPerChannel;
  TH2F	      *fHistMaxChargeClock;
  TH2F	      *fHistMaxChargeValueInt0;
  TH2F	      *fHistMaxChargeValueInt1;
  TH2F        *fHistTimePerPM_UnCorr;
  TH2F	      *fHistTimeVsChargeADA_UnCorr;
  TH2F	      *fHistTimeVsChargeADC_UnCorr;
  TH3F	      *fHistTimeVsChargePerPM_UnCorr;
  
  /*Aging monitoring*/
  TH3F	      *fHistTriggerChargePerPMPerV0Flag;
  TH3F	      *fHistTailChargePerPMPerV0Flag;
  TH3F	      *fHistIntegratedChargePerPMPerV0Flag;
  
  TH2F	      *fHistTriggerChargePerChannel;
  TH2F	      *fHistTriggerChargePerChannel_PF;
  TH2F	      *fHistTriggerChargePerChannel_TVX;
  TH2F	      *fHistTriggerChargePerChannel_PF_TVX;
  
  TH2F	      *fHistTailChargePerChannel;
  TH2F	      *fHistTailChargePerChannel_PF;
  TH2F	      *fHistTailChargePerChannel_TVX;
  TH2F	      *fHistTailChargePerChannel_PF_TVX;  
  
  TH2F        *fHistIntegratedChargePerChannel; 
  TH2F        *fHistIntegratedChargePerChannel_PF;
  TH2F        *fHistIntegratedChargePerChannel_TVX;  
  TH2F        *fHistIntegratedChargePerChannel_PF_TVX;  
    
  TH1F	      *fHistADMOn; 
  TH1F	      *fHistADMOff;
   
  AliAnalysisUtils *fAnalysisUtils;
  AliTriggerAnalysis *fTriggerAnalysis;
  

  AliAnalysisTaskADPilot(const AliAnalysisTaskADPilot&);            // not implemented
  AliAnalysisTaskADPilot& operator=(const AliAnalysisTaskADPilot&); // not implemented
  
  ClassDef(AliAnalysisTaskADPilot, 5);
};

#endif
