/**************************************************************************
 * Author: Michal Broz                                               *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//                 AliAnalysisTaskADPilot class
//            This task is for QAing the AD data from ESD/AOD
//              Origin:michal.broz@cern.ch
//-----------------------------------------------------------------
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "TSpline.h"

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"
#include "AliESDVZERO.h"
#include "AliADCalibData.h"

#include "AliAnalysisTaskADPilot.h"

ClassImp(AliAnalysisTaskADPilot)

//________________________________________________________________________
AliAnalysisTaskADPilot::AliAnalysisTaskADPilot() 
  : AliAnalysisTaskSE(),fListHist(0),fHistTotalChargePerEventADA(0),fHistTotalChargePerEventADC(0),
    fHistChargePerPM_All(0), fHistChargePerPM_BB(0), fHistChargePerPM_BG(0), fHistChargePerPM_Time(0), fHistTimePerPM_Corr(0), fHistTimePerPM_UnCorr(0),
    fHistWidthPerPM(0),fHistTimeVsChargeADA_Corr(0), fHistTimeVsChargeADC_Corr(0),fHistTimeVsChargeADA_UnCorr(0),fHistTimeVsChargeADC_UnCorr(0),fHistWidthVsCharge(0),
    fHistTimeVsChargeADA_Cut(0), fHistTimeVsChargeADC_Cut(0),
    fHistNBBflagsADA(0),fHistNBBflagsADC(0),fHistNBBflagsADAVsADC(0),
    fHistNBBCoincidencesADA(0),fHistNBBCoincidencesADC(0),fHistNBBCoincidencesADAVsADC(0),
    fHistNBGflagsADA(0),fHistNBGflagsADC(0),fHistNBGflagsADAVsADC(0),
    fHistNBGCoincidencesADA(0),fHistNBGCoincidencesADC(0),fHistNBGCoincidencesADAVsADC(0),
    fHistChargeNoFlag(0),fHistTimeNoFlag(0), fHistChargeNoTime(0),fHistFlagNoTime(0),fHistChargePerCoincidence(0),
    fHistMeanTimeADA(0),fHistMeanTimeADC(0),fHistMeanTimeDifference(0),fHistMeanTimeCorrelation(0),fHistMeanTimeSumDiff(0),fHistDecision(0),fHistDecisionBasic(0),fHistDecisionRobust(0),
    fHistTriggerMasked(0),fHistTriggerUnMasked(0),fHistTriggerOthers(0),
    fHistChargeVsClockInt0(0),fHistChargeVsClockInt1(0),fHistBBFlagVsClock(0),fHistBGFlagVsClock(0),fHistBBFlagPerChannel(0),fHistBGFlagPerChannel(0),fHistMaxChargeClock(0),
    fHistMaxChargeValueInt0(0),fHistMaxChargeValueInt1(0),
    fHistTimeVsChargePerPM_UnCorr(0),
    fHistChargeTriggerPerChannel(0),fHistChargeTriggerPerChannel_ADAND(0),fHistChargeTriggerPerChannel_PF(0),fHistChargeTriggerPerChannel_ADANDPF(0),
    fHistChargeTriggerADA(0),fHistChargeTriggerADA_ADAND(0),fHistChargeTriggerADA_PF(0),fHistChargeTriggerADA_ADANDPF(0),
    fHistChargeTriggerADC(0),fHistChargeTriggerADC_ADAND(0),fHistChargeTriggerADC_PF(0),fHistChargeTriggerADC_ADANDPF(0),
    fHistMedianTimeADA(0),fHistMedianTimeADC(0),fHistNTimesMedianADA(0),fHistNTimesMedianADC(0),fHistRobustTimeADA(0),fHistRobustTimeADC(0),fHistNTimesRobustADA(0),fHistNTimesRobustADC(0),
    fHistMedianIndDiffVsChargeADA(0),fHistMedianIndDiffVsChargeADC(0),
    fHistTimePairSumDiffADA_NoCut(0),fHistTimePairSumDiffADC_NoCut(0),fHistTimePairSumDiffADA_Cut(0),fHistTimePairSumDiffADC_Cut(0),
    fHistTimeVsChargeADA_Ex(0), fHistTimeVsChargeADC_Ex(0),
    fRun(0),fOldRun(0),fCalibData(0)
{
  // Dummy constructor
}
//________________________________________________________________________
AliAnalysisTaskADPilot::AliAnalysisTaskADPilot(const char *name) 
  : AliAnalysisTaskSE(name),fListHist(0),fHistTotalChargePerEventADA(0),fHistTotalChargePerEventADC(0),
    fHistChargePerPM_All(0), fHistChargePerPM_BB(0), fHistChargePerPM_BG(0), fHistChargePerPM_Time(0), fHistTimePerPM_Corr(0), fHistTimePerPM_UnCorr(0),
    fHistWidthPerPM(0),fHistTimeVsChargeADA_Corr(0), fHistTimeVsChargeADC_Corr(0),fHistTimeVsChargeADA_UnCorr(0),fHistTimeVsChargeADC_UnCorr(0),fHistWidthVsCharge(0),
    fHistTimeVsChargeADA_Cut(0), fHistTimeVsChargeADC_Cut(0),
    fHistNBBflagsADA(0),fHistNBBflagsADC(0),fHistNBBflagsADAVsADC(0),
    fHistNBBCoincidencesADA(0),fHistNBBCoincidencesADC(0),fHistNBBCoincidencesADAVsADC(0),
    fHistNBGflagsADA(0),fHistNBGflagsADC(0),fHistNBGflagsADAVsADC(0),
    fHistNBGCoincidencesADA(0),fHistNBGCoincidencesADC(0),fHistNBGCoincidencesADAVsADC(0),
    fHistChargeNoFlag(0),fHistTimeNoFlag(0), fHistChargeNoTime(0),fHistFlagNoTime(0),fHistChargePerCoincidence(0),
    fHistMeanTimeADA(0),fHistMeanTimeADC(0),fHistMeanTimeDifference(0),fHistMeanTimeCorrelation(0),fHistMeanTimeSumDiff(0),fHistDecision(0),fHistDecisionBasic(0),fHistDecisionRobust(0),
    fHistTriggerMasked(0),fHistTriggerUnMasked(0),fHistTriggerOthers(0),
    fHistChargeVsClockInt0(0),fHistChargeVsClockInt1(0),fHistBBFlagVsClock(0),fHistBGFlagVsClock(0),fHistBBFlagPerChannel(0),fHistBGFlagPerChannel(0),fHistMaxChargeClock(0),
    fHistMaxChargeValueInt0(0),fHistMaxChargeValueInt1(0),
    fHistTimeVsChargePerPM_UnCorr(0),
    fHistChargeTriggerPerChannel(0),fHistChargeTriggerPerChannel_ADAND(0),fHistChargeTriggerPerChannel_PF(0),fHistChargeTriggerPerChannel_ADANDPF(0),
    fHistChargeTriggerADA(0),fHistChargeTriggerADA_ADAND(0),fHistChargeTriggerADA_PF(0),fHistChargeTriggerADA_ADANDPF(0),
    fHistChargeTriggerADC(0),fHistChargeTriggerADC_ADAND(0),fHistChargeTriggerADC_PF(0),fHistChargeTriggerADC_ADANDPF(0),
    fHistMedianTimeADA(0),fHistMedianTimeADC(0),fHistNTimesMedianADA(0),fHistNTimesMedianADC(0),fHistRobustTimeADA(0),fHistRobustTimeADC(0),fHistNTimesRobustADA(0),fHistNTimesRobustADC(0),
    fHistMedianIndDiffVsChargeADA(0),fHistMedianIndDiffVsChargeADC(0),
    fHistTimePairSumDiffADA_NoCut(0),fHistTimePairSumDiffADC_NoCut(0),fHistTimePairSumDiffADA_Cut(0),fHistTimePairSumDiffADC_Cut(0),
    fHistTimeVsChargeADA_Ex(0), fHistTimeVsChargeADC_Ex(0),
    fRun(0),fOldRun(0),fCalibData(0)
{
  // Constructor
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskADPilot::~AliAnalysisTaskADPilot(){
  // Destructor
  if (fListHist) { delete fListHist; fListHist = 0x0; }
}
//________________________________________________________________________
void AliAnalysisTaskADPilot::UserCreateOutputObjects()
{

  // Create histograms
  // Called once

  const Int_t kNRawTimeBins = 3062;
  const Float_t kRawTimeMin = 0.976562;
  const Float_t kRawTimeMax = 300;
  
  const Int_t kNCorrTimeBins = 1638;
  const Float_t kCorrTimeMin = -79.980469; 
  const Float_t kCorrTimeMax = 79.980469;
  
  const Int_t kNTdcWidthBins =  153;
  const Float_t kTdcWidthMin = 2.343750;
  const Float_t kTdcWidthMax = 121.875000;
  
  const Int_t kNChargeChannelBins = 5500;
  const Int_t kChargeChannelMin = 0;
  const Int_t kChargeChannelMax = 5500;
  
  const Int_t kNChargeSideBins = 40000;
  const Int_t kChargeSideMin = 0;
  const Int_t kChargeSideMax = 40000;
    
  const Int_t kNChannelBins  =   16;
  const Float_t kChannelMin    =    -0.5;
  const Float_t kChannelMax    =   15.5;
  

  fListHist = new TList();
  fListHist->SetOwner(); 

if (!fHistTotalChargePerEventADA) {
    fHistTotalChargePerEventADA = CreateHist1D("fHistTotalChargePerEventADA","Total Charge in ADA per event",kNChargeSideBins,kChargeSideMin,kChargeSideMax,"ADC counts","Entries");
    fListHist->Add(fHistTotalChargePerEventADA);
  }
if (!fHistTotalChargePerEventADC) {
    fHistTotalChargePerEventADC = CreateHist1D("fHistTotalChargePerEventADC","Total Charge in ADC per event",kNChargeSideBins,kChargeSideMin,kChargeSideMax,"ADC counts","Entries");
    fListHist->Add(fHistTotalChargePerEventADC);
  }        
if (!fHistChargePerPM_All) {
    fHistChargePerPM_All = CreateHist2D("fHistChargePerPM_All","Charge per PM all events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
    fListHist->Add(fHistChargePerPM_All);
  }
if (!fHistChargePerPM_BB) {
    fHistChargePerPM_BB = CreateHist2D("fHistChargePerPM_BB","Charge per PM BB events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
    fListHist->Add(fHistChargePerPM_BB);
  } 
if (!fHistChargePerPM_BG) {
    fHistChargePerPM_BG = CreateHist2D("fHistChargePerPM_BG","Charge per PM BG events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
    fListHist->Add(fHistChargePerPM_BG);
  }
if (!fHistChargePerPM_Time) {
    fHistChargePerPM_Time = CreateHist2D("fHistChargePerPM_Time","Charge per PM Time events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
    fListHist->Add(fHistChargePerPM_Time);
  }    
if (!fHistTimePerPM_Corr) {
    fHistTimePerPM_Corr = CreateHist2D("fHistTimePerPM_Corr","Corrected Time per PM",kNChannelBins, kChannelMin, kChannelMax, kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"PM number","Leading time [ns]");
    fListHist->Add(fHistTimePerPM_Corr);
  } 
if (!fHistTimePerPM_UnCorr) {
    fHistTimePerPM_UnCorr = CreateHist2D("fHistTimePerPM_UnCorr","Raw Time per PM",kNChannelBins, kChannelMin, kChannelMax, kNRawTimeBins, kRawTimeMin, kRawTimeMax,"PM number","Leading time [ns]");
    fListHist->Add(fHistTimePerPM_UnCorr);
  }   
if (!fHistWidthPerPM) {
    fHistWidthPerPM = CreateHist2D("fHistWidthPerPM","Width per PM",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax,"PM number","Time width [ns]");
    fListHist->Add(fHistWidthPerPM);
  } 
if (!fHistTimeVsChargeADA_Cut) {
    fHistTimeVsChargeADA_Cut = CreateHist2D("fHistTimeVsChargeADA_Cut","Cutted Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsChargeADA_Cut);
  } 
if (!fHistTimeVsChargeADA_Corr) {
    fHistTimeVsChargeADA_Corr = CreateHist2D("fHistTimeVsChargeADA_Corr","Corrected Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsChargeADA_Corr);
  }
if (!fHistTimeVsChargeADA_UnCorr) {
    fHistTimeVsChargeADA_UnCorr = CreateHist2D("fHistTimeVsChargeADA_UnCorr","Raw Time vs Charge",kNRawTimeBins, kRawTimeMin, kRawTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsChargeADA_UnCorr);
  }
if (!fHistTimeVsChargeADC_Cut) {
    fHistTimeVsChargeADC_Cut = CreateHist2D("fHistTimeVsChargeADC_Cut","Cutted Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsChargeADC_Cut);
  }
if (!fHistTimeVsChargeADC_Corr) {
    fHistTimeVsChargeADC_Corr = CreateHist2D("fHistTimeVsChargeADC_Corr","Corrected Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsChargeADC_Corr);
  }
if (!fHistTimeVsChargeADC_UnCorr) {
    fHistTimeVsChargeADC_UnCorr = CreateHist2D("fHistTimeVsChargeADC_UnCorr","Raw Time vs Charge",kNRawTimeBins, kRawTimeMin, kRawTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsChargeADC_UnCorr);
  }  
if (!fHistWidthVsCharge) {
    fHistWidthVsCharge = CreateHist2D("fHistWidthVsCharge","Width vs Charge",kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax,kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Time width [ns]","ADC counts");
    fListHist->Add(fHistWidthVsCharge);
  } 
if (!fHistNBBflagsADA) {
    fHistNBBflagsADA = CreateHist1D("fHistNBBflagsADA","Number of BB flags in ADA",9,-0.5,8.5,"Number of BB flags in ADA per event","Entries");
    fListHist->Add(fHistNBBflagsADA);
  } 
if (!fHistNBBflagsADC) {
    fHistNBBflagsADC = CreateHist1D("fHistNBBflagsADC","Number of BB flags in ADC",9,-0.5,8.5,"Number of BB flags in ADC per event","Entries");
    fListHist->Add(fHistNBBflagsADC);
  }
if (!fHistNBBflagsADAVsADC) {
    fHistNBBflagsADAVsADC = CreateHist2D("fHistNBBflagsADAVsADC","Number of BB flags",9,-0.5,8.5,9,-0.5,8.5,"Number of BB flags in ADA","Number of BB flags in ADC");
    fListHist->Add(fHistNBBflagsADAVsADC);
  }

if (!fHistNBGflagsADA) {
    fHistNBGflagsADA = CreateHist1D("fHistNBGflagsADA","Number of BG flags in ADA",9,-0.5,8.5,"Number of BG flags in ADA per event","Entries");
    fListHist->Add(fHistNBGflagsADA);
  } 
if (!fHistNBGflagsADC) {
    fHistNBGflagsADC = CreateHist1D("fHistNBGflagsADC","Number of BG flags in ADC",9,-0.5,8.5,"Number of BG flags in ADC per event","Entries");
    fListHist->Add(fHistNBGflagsADC);
  }
if (!fHistNBGflagsADAVsADC) {
    fHistNBGflagsADAVsADC = CreateHist2D("fHistNBGflagsADAVsADC","Number of BG flags",9,-0.5,8.5,9,-0.5,8.5,"Number of BG flags in ADA","Number of BG flags in ADC");
    fListHist->Add(fHistNBGflagsADAVsADC);
  }
  
if (!fHistNBBCoincidencesADA) {
    fHistNBBCoincidencesADA = CreateHist1D("fHistNBBCoincidencesADA","Number of BB coincidences in ADA",5,-0.5,4.5,"Number of BB coincidences in ADA per event","Entries");
    fListHist->Add(fHistNBBCoincidencesADA);
  } 
if (!fHistNBBCoincidencesADC) {
    fHistNBBCoincidencesADC = CreateHist1D("fHistNBBCoincidencesADC","Number of BB coincidences in ADC",5,-0.5,4.5,"Number of BB coincidences in ADC per event","Entries");
    fListHist->Add(fHistNBBCoincidencesADC);
  } 
if (!fHistNBBCoincidencesADAVsADC) {
    fHistNBBCoincidencesADAVsADC = CreateHist2D("fHistNBBCoincidencesADAVsADC","Number of BB coincidences",5,-0.5,4.5, 5,-0.5,4.5,"Number of BB coincidences in ADA","Number of BB coincidences in ADC");
    fListHist->Add(fHistNBBCoincidencesADAVsADC);
  }
if (!fHistNBGCoincidencesADA) {
    fHistNBGCoincidencesADA = CreateHist1D("fHistNBGCoincidencesADA","Number of BG coincidences in ADA",5,-0.5,4.5,"Number of BG coincidences in ADA per event","Entries");
    fListHist->Add(fHistNBGCoincidencesADA);
  } 
if (!fHistNBGCoincidencesADC) {
    fHistNBGCoincidencesADC = CreateHist1D("fHistNBGCoincidencesADC","Number of BG coincidences in ADC",5,-0.5,4.5,"Number of BG coincidences in ADC per event","Entries");
    fListHist->Add(fHistNBGCoincidencesADC);
  } 
if (!fHistNBGCoincidencesADAVsADC) {
    fHistNBGCoincidencesADAVsADC = CreateHist2D("fHistNBGCoincidencesADAVsADC","Number of BG coincidences",5,-0.5,4.5, 5,-0.5,4.5,"Number of BG coincidences in ADA","Number of BG coincidences in ADC");
    fListHist->Add(fHistNBGCoincidencesADAVsADC);
  }      
if (!fHistChargeNoFlag) {
    fHistChargeNoFlag = CreateHist1D("fHistChargeNoFlag","Charge in PM without BB flag",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"Charge","Entries");
    fListHist->Add(fHistChargeNoFlag);
  } 
if (!fHistTimeNoFlag) {
    fHistTimeNoFlag = CreateHist2D("fHistTimeNoFlag","Time in PM without BB flag",kNChannelBins, kChannelMin, kChannelMax,kNRawTimeBins, kRawTimeMin, kRawTimeMax,"Channel","Time","Entries");
    fListHist->Add(fHistTimeNoFlag);
  }   
if (!fHistChargeNoTime) {
    fHistChargeNoTime = CreateHist2D("fHistChargeNoTime","Charge in PM without time measurement",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"Channel","Charge","Entries");
    fListHist->Add(fHistChargeNoTime);
  }
if (!fHistFlagNoTime) {
    fHistFlagNoTime = CreateHist1D("fHistFlagNoTime","events with BB/BG flag but no time",kNChannelBins, kChannelMin, kChannelMax,"Channel","Entries");
    fListHist->Add(fHistFlagNoTime);
  } 
if (!fHistChargePerCoincidence) {
     fHistChargePerCoincidence = CreateHist2D("fHistChargePerCoincidence","Charge in PM per coincidence",5,-0.5,4.5,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"Number of BB coincidences","Charge","Entries");
     fListHist->Add(fHistChargePerCoincidence);
     }

if (!fHistMeanTimeADA) {
    fHistMeanTimeADA = CreateHist1D("fHistMeanTimeADA","Mean Time in ADA",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time","Entries");
    fListHist->Add(fHistMeanTimeADA);
  }   

if (!fHistMeanTimeADC) {
    fHistMeanTimeADC = CreateHist1D("fHistMeanTimeADC","Mean Time in ADC",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time","Entries");
    fListHist->Add(fHistMeanTimeADC);
  }   

if (!fHistMeanTimeDifference) {
    fHistMeanTimeDifference = CreateHist1D("fHistMeanTimeDifference","Mean Time difference",1024,-150,150,"AD Mean time t_{A} - t_{C} [ns]","Entries");
    fListHist->Add(fHistMeanTimeDifference);
  }  
   
if (!fHistMeanTimeCorrelation) {
    fHistMeanTimeCorrelation = CreateHist2D("fHistMeanTimeCorrelation","Mean Time in ADA-ADC",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time ADA","Time ADC");
    fListHist->Add(fHistMeanTimeCorrelation);
  }   

if (!fHistMeanTimeSumDiff) {
    fHistMeanTimeSumDiff = CreateHist2D("fHistMeanTimeSumDiff","Mean Time in ADA-ADC",307, -150.000000, 149.804688, 410, 0.000000, 400.390625,"AD Mean time t_{A} - t_{C} [ns]","AD Mean time t_{A} + t_{C} [ns]");
    fListHist->Add(fHistMeanTimeSumDiff);
  }
if (!fHistDecision) {
    fHistDecision = CreateHist2D("fHistDecision","Offline decision in ADA-ADC based on  time",4,0 ,4,4,0,4,"ADA","ADC");
    fHistDecision->SetOption("coltext");
    fHistDecision->GetXaxis()->SetLabelSize(0.06);
    fHistDecision->GetYaxis()->SetLabelSize(0.06);
    fHistDecision->GetXaxis()->SetNdivisions(808,kFALSE);
    fHistDecision->GetYaxis()->SetNdivisions(808,kFALSE);
    fHistDecision->GetXaxis()->SetBinLabel(1, "Empty");
    fHistDecision->GetXaxis()->SetBinLabel(2, "BB");
    fHistDecision->GetXaxis()->SetBinLabel(3, "BG");
    fHistDecision->GetXaxis()->SetBinLabel(4, "Fake");
    fHistDecision->GetYaxis()->SetBinLabel(1, "Empty");
    fHistDecision->GetYaxis()->SetBinLabel(2, "BB");
    fHistDecision->GetYaxis()->SetBinLabel(3, "BG");
    fHistDecision->GetYaxis()->SetBinLabel(4, "Fake");
    fListHist->Add(fHistDecision);
 }
     
if (!fHistDecisionBasic) {
    fHistDecisionBasic = CreateHist2D("fHistDecisionBasic","Offline decision in ADA-ADC based on basic time",4,0 ,4,4,0,4,"ADA","ADC");
    fHistDecisionBasic->SetOption("coltext");
    fHistDecisionBasic->GetXaxis()->SetLabelSize(0.06);
    fHistDecisionBasic->GetYaxis()->SetLabelSize(0.06);
    fHistDecisionBasic->GetXaxis()->SetNdivisions(808,kFALSE);
    fHistDecisionBasic->GetYaxis()->SetNdivisions(808,kFALSE);
    fHistDecisionBasic->GetXaxis()->SetBinLabel(1, "Empty");
    fHistDecisionBasic->GetXaxis()->SetBinLabel(2, "BB");
    fHistDecisionBasic->GetXaxis()->SetBinLabel(3, "BG");
    fHistDecisionBasic->GetXaxis()->SetBinLabel(4, "Fake");
    fHistDecisionBasic->GetYaxis()->SetBinLabel(1, "Empty");
    fHistDecisionBasic->GetYaxis()->SetBinLabel(2, "BB");
    fHistDecisionBasic->GetYaxis()->SetBinLabel(3, "BG");
    fHistDecisionBasic->GetYaxis()->SetBinLabel(4, "Fake");
    fListHist->Add(fHistDecisionBasic);
 }

if (!fHistDecisionRobust) {
    fHistDecisionRobust = CreateHist2D("fHistDecisionRobust","Offline decision in ADA-ADC based on Robust time",4,0 ,4,4,0,4,"ADA","ADC");
    fHistDecisionRobust->SetOption("coltext");
    fHistDecisionRobust->GetXaxis()->SetLabelSize(0.06);
    fHistDecisionRobust->GetYaxis()->SetLabelSize(0.06);
    fHistDecisionRobust->GetXaxis()->SetNdivisions(808,kFALSE);
    fHistDecisionRobust->GetYaxis()->SetNdivisions(808,kFALSE);
    fHistDecisionRobust->GetXaxis()->SetBinLabel(1, "Empty");
    fHistDecisionRobust->GetXaxis()->SetBinLabel(2, "BB");
    fHistDecisionRobust->GetXaxis()->SetBinLabel(3, "BG");
    fHistDecisionRobust->GetXaxis()->SetBinLabel(4, "Fake");
    fHistDecisionRobust->GetYaxis()->SetBinLabel(1, "Empty");
    fHistDecisionRobust->GetYaxis()->SetBinLabel(2, "BB");
    fHistDecisionRobust->GetYaxis()->SetBinLabel(3, "BG");
    fHistDecisionRobust->GetYaxis()->SetBinLabel(4, "Fake");
    fListHist->Add(fHistDecisionRobust);
 }
    
if(!fHistTriggerMasked) {
    fHistTriggerMasked = CreateHist1D("fHistTriggerMasked","Trigger inputs, from FEE (BC masked)",10,0 ,10,"AD0 Trigger Type","Counts");
    fHistTriggerMasked->SetFillColor(kAzure-8);
    fHistTriggerMasked->SetLineWidth(2);
    fHistTriggerMasked->GetXaxis()->SetLabelSize(0.04);
    fHistTriggerMasked->GetXaxis()->SetNdivisions(808,kFALSE);
    fHistTriggerMasked->GetXaxis()->SetBinLabel(1, "UBA & UBC");
    fHistTriggerMasked->GetXaxis()->SetBinLabel(2, "UBA || UBC");
    fHistTriggerMasked->GetXaxis()->SetBinLabel(3, "UGA & UBC");
    fHistTriggerMasked->GetXaxis()->SetBinLabel(4, "UGA");
    fHistTriggerMasked->GetXaxis()->SetBinLabel(5, "UGC & UBA");
    fHistTriggerMasked->GetXaxis()->SetBinLabel(6, "UGC");
    fHistTriggerMasked->GetXaxis()->SetBinLabel(7, "UBA");
    fHistTriggerMasked->GetXaxis()->SetBinLabel(8, "UBC");
    fHistTriggerMasked->GetXaxis()->SetBinLabel(9, "UGA || UGC");
    fHistTriggerMasked->GetXaxis()->SetBinLabel(10, "(UGA & UBC) || (UGC & UBA)");
    fListHist->Add(fHistTriggerMasked);
    }
 
 if(!fHistTriggerUnMasked) {
    fHistTriggerUnMasked = CreateHist1D("fHistTriggerUnMasked","Trigger inputs, from FEE (BC masked)",10,0 ,10,"AD0 Trigger Type","Counts");
    fHistTriggerUnMasked->SetFillColor(kAzure-8);
    fHistTriggerUnMasked->SetLineWidth(2);
    fHistTriggerUnMasked->GetXaxis()->SetLabelSize(0.04);
    fHistTriggerUnMasked->GetXaxis()->SetNdivisions(808,kFALSE);
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(1, "UBA & UBC");
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(2, "UBA || UBC");
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(3, "UGA & UBC");
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(4, "UGA");
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(5, "UGC & UBA");
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(6, "UGC");
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(7, "UBA");
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(8, "UBC");
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(9, "UGA || UGC");
    fHistTriggerUnMasked->GetXaxis()->SetBinLabel(10, "(UGA & UBC) || (UGC & UBA)");
    fListHist->Add(fHistTriggerUnMasked);
    }
    
if(!fHistTriggerOthers) {
    fHistTriggerOthers = CreateHist1D("fHistTriggerOthers","Trigger inputs, from other detectors",2,0 ,2,"Trigger inputs","Counts");
    fHistTriggerOthers->SetFillColor(kAzure-8);
    fHistTriggerOthers->SetLineWidth(2);
    fHistTriggerOthers->GetXaxis()->SetLabelSize(0.04);
    fHistTriggerOthers->GetXaxis()->SetNdivisions(808,kFALSE);
    fHistTriggerOthers->GetXaxis()->SetBinLabel(1, "VBAND");
    fHistTriggerOthers->GetXaxis()->SetBinLabel(2, "VBOR");

    fListHist->Add(fHistTriggerOthers);
    }

     
if (!fHistChargeVsClockInt0) {
    fHistChargeVsClockInt0 = CreateHist2D("fHistChargeVsClockInt0","Charge Vs Clock (Int0)",16,-0.5, 15.5, 21, -10.5, 10.5,"Channel","LHC clock");
    fListHist->Add(fHistChargeVsClockInt0);
  }  
if (!fHistChargeVsClockInt1) {
    fHistChargeVsClockInt1 = CreateHist2D("fHistChargeVsClockInt1","Charge Vs Clock (Int1)",16,-0.5, 15.5, 21, -10.5, 10.5,"Channel","LHC clock");
    fListHist->Add(fHistChargeVsClockInt1);
  }  
if (!fHistBBFlagVsClock) {
    fHistBBFlagVsClock = CreateHist2D("fHistBBFlagVsClock","BB Flag Vs Clock",16,-0.5, 15.5, 21, -10.5, 10.5,"Channel","LHC clock");
    fListHist->Add(fHistBBFlagVsClock);
  } 
if (!fHistBGFlagVsClock) {
    fHistBGFlagVsClock = CreateHist2D("fHistBGFlagVsClock","BG Flag Vs Clock",16,-0.5, 15.5, 21, -10.5, 10.5,"Channel","LHC clock");
    fListHist->Add(fHistBGFlagVsClock);
  } 
if (!fHistBBFlagPerChannel) {
    fHistBBFlagPerChannel = CreateHist2D("fHistBBFlagPerChannel","Number of BB flags per channel",16,-0.5, 15.5, 25, -0.5, 21.5,"Channel","BB Flags Count");
    fListHist->Add(fHistBBFlagPerChannel);
  } 
if (!fHistBGFlagPerChannel) {
    fHistBGFlagPerChannel = CreateHist2D("fHistBGFlagPerChannel","Number of BG flags per channel",16,-0.5, 15.5, 25, -0.5, 21.5,"Channel","BG Flags Count");
    fListHist->Add(fHistBGFlagPerChannel);
  } 
if (!fHistMaxChargeClock) {
    fHistMaxChargeClock = CreateHist2D("fHistMaxChargeClock","Clock with maximum charge per channel",16,-0.5, 15.5, 21, -10.5, 10.5,"Channel","LHC clock");
    fListHist->Add(fHistMaxChargeClock);
  }
if (!fHistMaxChargeValueInt0) {
    fHistMaxChargeValueInt0 = CreateHist2D("fHistMaxChargeValueInt0","Maximum charge value per PM Int0",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"PM number","ADC counts");
    fListHist->Add(fHistMaxChargeValueInt0);
  }
if (!fHistMaxChargeValueInt1) {
    fHistMaxChargeValueInt1 = CreateHist2D("fHistMaxChargeValueInt1","Maximum charge value per PM Int1",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"PM number","ADC counts");
    fListHist->Add(fHistMaxChargeValueInt1);
  }
if(!fHistTimeVsChargePerPM_UnCorr) {
    fHistTimeVsChargePerPM_UnCorr = CreateHist3D("fHistTimeVsChargePerPM_UnCorr","Raw Time vs Charge per PM",2600, 400, 3000, 
    					200,-4,0,
					kNChannelBins, kChannelMin, kChannelMax,"Leading time [ns]","ADC counts","Channel");
    fListHist->Add(fHistTimeVsChargePerPM_UnCorr);
  } 
//Robust time testing
if (!fHistMedianTimeADA) {
    fHistMedianTimeADA = CreateHist1D("fHistMedianTimeADA","Median Time in ADA",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time","Entries");
    fListHist->Add(fHistMedianTimeADA);
  }   

if (!fHistMedianTimeADC) {
    fHistMedianTimeADC = CreateHist1D("fHistMedianTimeADC","Median Time in ADC",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time","Entries");
    fListHist->Add(fHistMedianTimeADC);
  }   
if (!fHistRobustTimeADA) {
    fHistRobustTimeADA = CreateHist1D("fHistRobustTimeADA","Robust Time in ADA",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time","Entries");
    fListHist->Add(fHistRobustTimeADA);
  }   

if (!fHistRobustTimeADC) {
    fHistRobustTimeADC = CreateHist1D("fHistRobustTimeADC","Robust Time in ADC",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time","Entries");
    fListHist->Add(fHistRobustTimeADC);
  }   
if (!fHistNTimesMedianADA) {
    fHistNTimesMedianADA = CreateHist1D("fHistNTimesMedianADA","Number of channels with time for median in ADA",9,-0.5,8.5,"Number of times in median in ADA per event","Entries");
    fListHist->Add(fHistNTimesMedianADA);
  } 
if (!fHistNTimesMedianADC) {
    fHistNTimesMedianADC = CreateHist1D("fHistNTimesMedianADC","Number of channels with time for median in ADC",9,-0.5,8.5,"Number of times in median in ADC per event","Entries");
    fListHist->Add(fHistNTimesMedianADC);
  } 

if (!fHistNTimesRobustADA) {
    fHistNTimesRobustADA = CreateHist1D("fHistNTimesRobustADA","Number of channels with time for Robust in ADA",9,-0.5,8.5,"Number of times in Robust in ADA per event","Entries");
    fListHist->Add(fHistNTimesRobustADA);
  } 
if (!fHistNTimesRobustADC) {
    fHistNTimesRobustADC = CreateHist1D("fHistNTimesRobustADC","Number of channels with time for Robust in ADC",9,-0.5,8.5,"Number of times in Robust in ADC per event","Entries");
    fListHist->Add(fHistNTimesRobustADC);
  } 
if (!fHistMedianIndDiffVsChargeADA) {
    fHistMedianIndDiffVsChargeADA = CreateHist2D("fHistMedianIndDiffVsChargeADA","Corrected Time vs Charge",1024,-50,50, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Time difference [ns]","ADC counts");
    fListHist->Add(fHistMedianIndDiffVsChargeADA);
  }
if (!fHistMedianIndDiffVsChargeADC) {
    fHistMedianIndDiffVsChargeADC = CreateHist2D("fHistMedianIndDiffVsChargeADC","Corrected Time vs Charge",1024,-50,50, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Time difference [ns]","ADC counts");
    fListHist->Add(fHistMedianIndDiffVsChargeADC);
  }
  
if(!fHistTimePairSumDiffADA_NoCut) {
    fHistTimePairSumDiffADA_NoCut = CreateHist2D("fHistTimePairSumDiffADA_NoCut","Time Sum vs Diff  ADA", 410, 0.000000, 400.390625,410, 0.000000, 40.039062,"t1 + t2","t1 - t2");
    fListHist->Add(fHistTimePairSumDiffADA_NoCut);
  } 

if(!fHistTimePairSumDiffADC_NoCut) {
    fHistTimePairSumDiffADC_NoCut = CreateHist2D("fHistTimePairSumDiffADC_NoCut","Time Sum vs Diff  ADC", 410, 0.000000, 400.390625,410, 0.000000, 40.039062,"t1 + t2","t1 - t2");
    fListHist->Add(fHistTimePairSumDiffADC_NoCut);
  } 
if(!fHistTimePairSumDiffADA_Cut) {
    fHistTimePairSumDiffADA_Cut = CreateHist2D("fHistTimePairSumDiffADA_Cut","Time Sum vs Diff  ADA", 410, 0.000000, 400.390625,410, 0.000000, 40.039062,"t1 + t2","t1 - t2");
    fListHist->Add(fHistTimePairSumDiffADA_Cut);
  } 

if(!fHistTimePairSumDiffADC_Cut) {
    fHistTimePairSumDiffADC_Cut = CreateHist2D("fHistTimePairSumDiffADC_Cut","Time Sum vs Diff  ADC", 410, 0.000000, 400.390625,410, 0.000000, 40.039062,"t1 + t2","t1 - t2");
    fListHist->Add(fHistTimePairSumDiffADC_Cut);
  } 
if (!fHistTimeVsChargeADA_Ex) {
    fHistTimeVsChargeADA_Ex = CreateHist2D("fHistTimeVsChargeADA_Ex","Excluded Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsChargeADA_Ex);
  }
if (!fHistTimeVsChargeADC_Ex) {
    fHistTimeVsChargeADC_Ex = CreateHist2D("fHistTimeVsChargeADC_Ex","Excluded Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsChargeADC_Ex);
  }
  
//---------------------------------------------  
if (!fHistChargeTriggerPerChannel) {
    fHistChargeTriggerPerChannel = CreateHist2D("fHistChargeTriggerPerChannel","Trigger charge per chanel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
    fListHist->Add(fHistChargeTriggerPerChannel);
  }
if (!fHistChargeTriggerPerChannel_ADAND) {
    fHistChargeTriggerPerChannel_ADAND = CreateHist2D("fHistChargeTriggerPerChannel_ADAND","Trigger charge per chanel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
    fListHist->Add(fHistChargeTriggerPerChannel_ADAND);
  }
if (!fHistChargeTriggerPerChannel_PF) {
    fHistChargeTriggerPerChannel_PF = CreateHist2D("fHistChargeTriggerPerChannel_PF","Trigger charge per chanel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
    fListHist->Add(fHistChargeTriggerPerChannel_PF);
  }
if (!fHistChargeTriggerPerChannel_ADANDPF) {
    fHistChargeTriggerPerChannel_ADANDPF = CreateHist2D("fHistChargeTriggerPerChannel_ADANDPF","Trigger charge per chanel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
    fListHist->Add(fHistChargeTriggerPerChannel_ADANDPF);
  }
if (!fHistChargeTriggerADA) {
    fHistChargeTriggerADA = CreateHist1D("fHistChargeTriggerADA","Trigger charge ADA",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
    fListHist->Add(fHistChargeTriggerADA);
  }
if (!fHistChargeTriggerADA_ADAND) {
    fHistChargeTriggerADA_ADAND = CreateHist1D("fHistChargeTriggerADA_ADAND","Trigger charge ADA, ADAND events",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
    fListHist->Add(fHistChargeTriggerADA_ADAND);
  }
if (!fHistChargeTriggerADA_PF) {
    fHistChargeTriggerADA_PF = CreateHist1D("fHistChargeTriggerADA_PF","Trigger charge ADA, PF protection",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
    fListHist->Add(fHistChargeTriggerADA_PF);
  }
if (!fHistChargeTriggerADA_ADANDPF) {
    fHistChargeTriggerADA_ADANDPF = CreateHist1D("fHistChargeTriggerADA_ADANDPF","Trigger charge ADA, ADAND events",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
    fListHist->Add(fHistChargeTriggerADA_ADANDPF);
  }
if (!fHistChargeTriggerADC) {
    fHistChargeTriggerADC = CreateHist1D("fHistChargeTriggerADC","Trigger charge ADC",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
    fListHist->Add(fHistChargeTriggerADC);
  }
if (!fHistChargeTriggerADC_ADAND) {
    fHistChargeTriggerADC_ADAND = CreateHist1D("fHistChargeTriggerADC_ADAND","Trigger charge ADC, ADAND events",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
    fListHist->Add(fHistChargeTriggerADC_ADAND);
  }
if (!fHistChargeTriggerADC_PF) {
    fHistChargeTriggerADC_PF = CreateHist1D("fHistChargeTriggerADC_PF","Trigger charge ADC, PF protection",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
    fListHist->Add(fHistChargeTriggerADC_PF);
  }
if (!fHistChargeTriggerADC_ADANDPF) {
    fHistChargeTriggerADC_ADANDPF = CreateHist1D("fHistChargeTriggerADC_ADANDPF","Trigger charge ADC, ADAND events",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
    fListHist->Add(fHistChargeTriggerADC_ADANDPF);
  }


  // Post output data.
  PostData(1, fListHist);
}
//________________________________________________________________________
void AliAnalysisTaskADPilot::SetTimeSlewing()
{
 
    AliCDBManager *man = AliCDBManager::Instance();
    man->SetDefaultStorage("alien://Folder=/alice/data/2015/OCDB");
    man->SetRun(fRun);

    AliCDBEntry *ent = man->Get("AD/Calib/TimeSlewing");
    TList *fListSplines = (TList*)ent->GetObject();
    
    for(Int_t i=0; i<16; i++) fTimeSlewingSpline[i] = dynamic_cast<TSpline3*> (fListSplines->At(i));
}

//________________________________________________________________________
void AliAnalysisTaskADPilot::SetCalibData()
{
 
    AliCDBManager *man = AliCDBManager::Instance();
    man->SetDefaultStorage("alien://Folder=/alice/data/2015/OCDB");
    man->SetRun(fRun);

    AliCDBEntry *ent = man->Get("AD/Calib/Data");
    fCalibData = (AliADCalibData*)ent->GetObject();
    
}

//____________________________________________________________________________
Float_t AliAnalysisTaskADPilot::CorrectLeadingTime(Int_t i, Float_t time, Float_t adc)
{
  // Correct the leading time
  // for slewing effect and
  // misalignment of the channels
  const Double_t fTOF[4] = {65.2418, 65.1417, 56.6459, 56.7459};
  const Double_t fRes = 25.0/256.0;
   
  if (time < 1e-6) return -1024;

  // In case of pathological signals
  if (adc < 1) return time;

  // Slewing and offset correction
  time -= fTimeSlewingSpline[i]->Eval(TMath::Log10(1.0/adc))*fRes;
  time += fTOF[i/4];
  
  return time;
}

//________________________________________________________________________
void AliAnalysisTaskADPilot::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  AliVEvent* fEvent = InputEvent();
  if (!fEvent) {
    Printf("ERROR: Event not available");
    return;
  }
  AliESDEvent* fESD = (AliESDEvent*)fEvent;
  if (!fESD) {
    Printf("ERROR: ESD not available");
    return;
  }
  AliESDAD* esdAD = fESD->GetADData();
  if (!esdAD) {
    Printf("ERROR: No ESD AD");
    return;
  }
  AliESDfriend *fESDfriend = fESD->FindFriend();
  if (!fESDfriend) {
    Printf("ERROR: No ESD friend");
    return;
  }
  AliESDADfriend* esdADfriend = fESDfriend->GetADfriend();
  
  TString trigger = fESD->GetFiredTriggerClasses();
  if(!trigger.Contains("C0TVX-B") && !trigger.Contains("CINT10-B") && !trigger.Contains("CINT5-B")) return;
  
  fRun=fEvent->GetRunNumber();
  
  if (fRun!=fOldRun){
    SetTimeSlewing();
    SetCalibData();
    fOldRun=fRun;
  }

  Float_t totChargeADA = 0;
  Float_t totChargeADC = 0;
  Int_t nBBflagsADA = 0;
  Int_t nBBflagsADC = 0;
  Int_t nBGflagsADA = 0;
  Int_t nBGflagsADC = 0;
  Float_t fCharges[16];
  Bool_t globalPF = kTRUE;
  Float_t chargeADA   = 0.;
  Float_t chargeADC   = 0.;
  
  for(Int_t i=0; i<16; i++){ 
  	if(i<8)totChargeADC += esdAD->GetAdc(i);
	if(i>7)totChargeADA += esdAD->GetAdc(i);
	if(i<8 && esdAD->GetBBFlag(i)) nBBflagsADC++;
	if(i>7 && esdAD->GetBBFlag(i)) nBBflagsADA++;
	if(i<8 && esdAD->GetBGFlag(i)) nBGflagsADC++;
	if(i>7 && esdAD->GetBGFlag(i)) nBGflagsADA++;
	
	fHistChargePerPM_All->Fill(i,esdAD->GetAdc(i));
	if(esdAD->GetBBFlag(i))fHistChargePerPM_BB->Fill(i,esdAD->GetAdc(i));
	if(esdAD->GetBGFlag(i))fHistChargePerPM_BG->Fill(i,esdAD->GetAdc(i));
	if(esdAD->GetTime(i)>1e-6)fHistChargePerPM_Time->Fill(i,esdAD->GetAdc(i));
	
	fHistWidthPerPM->Fill(i,esdAD->GetWidth(i));
	fHistWidthVsCharge->Fill(esdAD->GetWidth(i),esdAD->GetAdc(i));
	
	if(esdADfriend){
		Float_t corrTime = CorrectLeadingTime(i,esdADfriend->GetTime(i),esdAD->GetAdc(i));
		fHistTimePerPM_Corr->Fill(i,corrTime);
	
		if(i<8)fHistTimeVsChargeADC_Corr->Fill(corrTime,esdAD->GetAdc(i));
        	if(i>7)fHistTimeVsChargeADA_Corr->Fill(corrTime,esdAD->GetAdc(i));
		
		if(!esdAD->GetBBFlag(i)) {
			fHistChargeNoFlag->Fill(esdAD->GetAdc(i));
			fHistTimeNoFlag->Fill(i,esdADfriend->GetTime(i));
			}
		if(esdADfriend->GetTime(i) < 1e-6) fHistChargeNoTime->Fill(i,esdAD->GetAdc(i));
		if(esdADfriend->GetTime(i) < 1e-6) fHistChargeNoTime->Fill(i,esdAD->GetAdc(i));
		if(esdADfriend->GetTime(i) < 1e-6 && (esdAD->GetBBFlag(i)||esdAD->GetBGFlag(i)))fHistFlagNoTime->Fill(i);
	
		fHistTimePerPM_UnCorr->Fill(i,esdADfriend->GetTime(i));
 		if(i<8) fHistTimeVsChargeADC_UnCorr->Fill(esdADfriend->GetTime(i),esdAD->GetAdc(i));
 		if(i>7) fHistTimeVsChargeADA_UnCorr->Fill(esdADfriend->GetTime(i),esdAD->GetAdc(i));
		if(esdAD->GetAdc(i)>0)fHistTimeVsChargePerPM_UnCorr->Fill(TMath::Nint(esdADfriend->GetTime(i)/(25./256.)),TMath::Log10(1.0/esdAD->GetAdc(i)),i);
		Int_t nbbFlag = 0;
      		Int_t nbgFlag = 0;
		Int_t charge[21];
		Bool_t localPF = kTRUE;
		
		for(Int_t iClock=0; iClock<21; iClock++){
			charge[iClock] = esdADfriend->GetPedestal(i,iClock);
			Bool_t intgr = esdADfriend->GetIntegratorFlag(i,iClock);
			Bool_t bbFlag = esdADfriend->GetBBFlag(i,iClock);
			Bool_t bgFlag = esdADfriend->GetBGFlag(i,iClock);
			if(bbFlag) nbbFlag++;
			if(bgFlag) nbgFlag++;
			if((bbFlag || bgFlag) && iClock < 10){globalPF = kFALSE; localPF = kFALSE;}
	
			if(!intgr)fHistChargeVsClockInt0->Fill(i,iClock-10,charge[iClock]);
			if(intgr)fHistChargeVsClockInt1->Fill(i,iClock-10,charge[iClock]);
			fHistBBFlagVsClock->Fill(i,iClock-10, bbFlag);
			fHistBGFlagVsClock->Fill(i,iClock-10, bgFlag);
			}
		//Gain monitoring
	  	Int_t k = i + 16*esdADfriend->GetIntegratorFlag(i,11);
		fCharges[i] = esdADfriend->GetPedestal(i,11) - fCalibData->GetPedestal(k);
		fHistChargeTriggerPerChannel->Fill(i,fCharges[i]);
		if(localPF)fHistChargeTriggerPerChannel_PF->Fill(i,fCharges[i]);
		UShort_t fTriggerBC = esdAD->GetTriggerBits();
  		if(fTriggerBC & (1 << 0)){
			fHistChargeTriggerPerChannel_ADAND->Fill(i,fCharges[i]);
			if(localPF)fHistChargeTriggerPerChannel_ADANDPF->Fill(i,fCharges[i]);
			}	

		fHistBBFlagPerChannel->Fill(i,nbbFlag);
		fHistBGFlagPerChannel->Fill(i,nbgFlag);
		Int_t maxClock = TMath::LocMax(21,charge);
		fHistMaxChargeClock->Fill(i,maxClock-10);
		if(!esdADfriend->GetIntegratorFlag(i,maxClock))fHistMaxChargeValueInt0->Fill(i,charge[maxClock]);
		if( esdADfriend->GetIntegratorFlag(i,maxClock))fHistMaxChargeValueInt1->Fill(i,charge[maxClock]);
		}
  }
  if(esdADfriend){
  //Gain monitoring
  for(Int_t i = 0; i<8; i++)chargeADC += fCharges[i];
  for(Int_t i = 8; i<16; i++)chargeADA += fCharges[i];
  fHistChargeTriggerADA->Fill(chargeADA);
  fHistChargeTriggerADC->Fill(chargeADC);
  if(globalPF){
  	  fHistChargeTriggerADA_PF->Fill(chargeADA);
  	  fHistChargeTriggerADC_PF->Fill(chargeADC);
  	  }
  UShort_t fTriggerBC = esdAD->GetTriggerBits();
  if(fTriggerBC & (1 << 0) ? kTRUE : kFALSE){
	fHistChargeTriggerADA_ADAND->Fill(chargeADA);
  	fHistChargeTriggerADC_ADAND->Fill(chargeADC);
	if(globalPF){
		fHistChargeTriggerADA_ADANDPF->Fill(chargeADA);
  		fHistChargeTriggerADC_ADANDPF->Fill(chargeADC);
		}
	}
  }	
  Int_t nBBCoincidencesADA = 0;
  Int_t nBBCoincidencesADC = 0;
  Int_t nBGCoincidencesADA = 0;
  Int_t nBGCoincidencesADC = 0;
  for(Int_t i=0; i<4; i++){
  	if(esdAD->GetBBFlag(i) && esdAD->GetBBFlag(i+4)) nBBCoincidencesADC++;
	if(esdAD->GetBBFlag(i+8) && esdAD->GetBBFlag(i+12)) nBBCoincidencesADA++;
	if(esdAD->GetBGFlag(i) && esdAD->GetBGFlag(i+4)) nBGCoincidencesADC++;
	if(esdAD->GetBGFlag(i+8) && esdAD->GetBGFlag(i+12)) nBGCoincidencesADA++;
  	}
  for(Int_t i=0; i<8; i++) fHistChargePerCoincidence->Fill(nBBCoincidencesADC,esdAD->GetAdc(i));
  for(Int_t i=8; i<16; i++) fHistChargePerCoincidence->Fill(nBBCoincidencesADA,esdAD->GetAdc(i));
	
  fHistTotalChargePerEventADA->Fill(totChargeADA);
  fHistTotalChargePerEventADC->Fill(totChargeADC);
  fHistNBBflagsADA->Fill(nBBflagsADA);
  fHistNBBflagsADC->Fill(nBBflagsADC);
  fHistNBBflagsADAVsADC->Fill(nBBflagsADA,nBBflagsADC);
  fHistNBBCoincidencesADA->Fill(nBBCoincidencesADA);
  fHistNBBCoincidencesADC->Fill(nBBCoincidencesADC);
  fHistNBBCoincidencesADAVsADC->Fill(nBBCoincidencesADA,nBBCoincidencesADC);
  fHistNBGCoincidencesADA->Fill(nBGCoincidencesADA);
  fHistNBGCoincidencesADC->Fill(nBGCoincidencesADC);
  fHistNBGCoincidencesADAVsADC->Fill(nBGCoincidencesADA,nBGCoincidencesADC);
    
  //Triggers
  Bool_t UBA = kFALSE;
  Bool_t UBC = kFALSE;
  Bool_t UGA = kFALSE;
  Bool_t UGC = kFALSE;

  if(nBBCoincidencesADA>=1) UBA = kTRUE;
  if(nBBCoincidencesADC>=1) UBC = kTRUE;
  if(nBGCoincidencesADA>=1) UGA = kTRUE;
  if(nBGCoincidencesADC>=1) UGC = kTRUE;
  
  if(UBA && UBC) fHistTriggerMasked->Fill(0);
  if(UBA || UBC) fHistTriggerMasked->Fill(1);
  if(UGA && UBC) fHistTriggerMasked->Fill(2);
  if(UGA)fHistTriggerMasked->Fill(3);
  if(UGC && UBA) fHistTriggerMasked->Fill(4);
  if(UGC) fHistTriggerMasked->Fill(5);
  if(UBA)fHistTriggerMasked->Fill(6);
  if(UBC) fHistTriggerMasked->Fill(7); 
  if(UGA || UGC) fHistTriggerMasked->Fill(8);
  if((UGA && UBC) || (UGC && UBA))fHistTriggerMasked->Fill(9);
	
  
  if(esdADfriend){
  	UShort_t fTriggerUnBC = esdADfriend->GetTriggerInputs();
  	for(Int_t i = 0; i<6; i++) if(fTriggerUnBC & (1 << i) ? kTRUE : kFALSE) fHistTriggerUnMasked->Fill(i);
	for(Int_t i = 12; i<16; i++) if(fTriggerUnBC & (1 << i) ? kTRUE : kFALSE) fHistTriggerUnMasked->Fill(i-6);	
  }
  
  //Robusts time testing
  const Int_t minAdc = 1;
  
  Double_t timeBasicADA=0, timeBasicADC=0;
  Double_t weightBasicADA =0., weightBasicADC = 0.;
  
  Double_t timesADA[8], timesADC[8];
  Double_t weightADA[8], weightADC[8];
  UInt_t   ntimeADA=0, ntimeADC=0;
  UInt_t   itimeADA[8], itimeADC[8];
  
  for (Int_t i = 0; i < 16; ++i) {
    Float_t adc = esdAD->GetAdc(i);
    if (adc > minAdc) {
      Float_t time = esdAD->GetTime(i);
	if(time > (-1024 + 1.e-6)){
		Float_t timeErr = 1/adc;
		if (i<8) {
			itimeADC[ntimeADC] = i;
	    		timesADC[ntimeADC] = time;
	    		weightADC[ntimeADC] = 1.0/(timeErr*timeErr);
			timeBasicADC += time/(timeErr*timeErr);
			weightBasicADC += 1.0/(timeErr*timeErr);
			ntimeADC++;
	  		}
		else{
			itimeADA[ntimeADA] = i-8;
	    		timesADA[ntimeADA] = time;
	    		weightADA[ntimeADA] = 1.0/(timeErr*timeErr);
			timeBasicADA += time/(timeErr*timeErr);
			weightBasicADA += 1.0/(timeErr*timeErr);
			ntimeADA++;
	  		}
      		}
    	}
  } // end of loop over channels
  if(weightBasicADA > 1) timeBasicADA /= weightBasicADA; 
  else timeBasicADA = -1024.;
  if(weightBasicADC > 1) timeBasicADC /= weightBasicADC;
  else timeBasicADC = -1024.;
  
  Double_t medianTimeADA = 0;
  if (ntimeADA > 0) medianTimeADA = TMath::Median(ntimeADA,timesADA,weightADA);
  else medianTimeADA = -1024;
  Double_t medianTimeADC = 0;
  if (ntimeADC > 0) medianTimeADC = TMath::Median(ntimeADC,timesADC,weightADC);
  else medianTimeADC = -1024;
  
  const Float_t fADADist = 56.69;
  const Float_t fADCDist = 65.19;
  TF1 *fEarlyHitCutShape = new TF1("fEarlyHitCutShape", " [0]+(x>[2])*[1]*(x-[2])**2");
  fEarlyHitCutShape->SetParameter(0,1.5);
  fEarlyHitCutShape->SetParameter(1,0.02);
  
  Double_t timeRobustADA=0, timeRobustADC=0;
  Double_t weightRobustADA =0., weightRobustADC = 0.;
  UInt_t   ntimeRobustADA=0, ntimeRobustADC=0;
  UInt_t   itimeRobustADA[8], itimeRobustADC[8];
  
  fEarlyHitCutShape->SetParameter(2,2*fADCDist -10);
  for (Int_t i = 0; i < 4; ++i) {
    Float_t adc1 = esdAD->GetAdc(i);
    Float_t adc2 = esdAD->GetAdc(i+4);
    if (adc1 > minAdc && adc2 > minAdc) {
      Float_t time1 = esdAD->GetTime(i);
      Float_t time2 = esdAD->GetTime(i+4);
	if(time1 > (-1024+1.e-6) && time2 > (-1024+1.e-6)){
		Float_t timeErr1 = 1/adc1;
		Float_t timeErr2 = 1/adc2;
		Float_t timeDiff = TMath::Abs(time1-time2);
		Float_t timeSum = time1+time2;
		fHistTimePairSumDiffADC_NoCut->Fill(timeSum,timeDiff);
		Float_t earlyHitCut = 1000;
		if(TMath::Abs(timeSum - 2*fADCDist) < 20) earlyHitCut = fEarlyHitCutShape->Eval(timeSum);
		if(timeDiff < earlyHitCut){
			itimeRobustADC[ntimeRobustADC] = i;
			ntimeRobustADC++;
			timeRobustADC += time1/(timeErr1*timeErr1);
			weightRobustADC += 1./(timeErr1*timeErr1);
			
			itimeRobustADC[ntimeRobustADC] = i+4;
			ntimeRobustADC++;
			timeRobustADC += time2/(timeErr2*timeErr2);
			weightRobustADC += 1./(timeErr2*timeErr2);
			fHistTimePairSumDiffADC_Cut->Fill(timeSum,timeDiff);
			fHistTimeVsChargeADC_Cut->Fill(time1,adc1);
			fHistTimeVsChargeADC_Cut->Fill(time2,adc2);
        		
			}
		}
	 else{
	 	if(time1 > (-1024+1.e-6)) fHistTimeVsChargeADC_Ex->Fill(time1,adc1);
		if(time2 > (-1024+1.e-6)) fHistTimeVsChargeADC_Ex->Fill(time2,adc2);
	 	}
	}
		
  }
  fEarlyHitCutShape->SetParameter(2,2*fADADist -10);
  for (Int_t i = 0; i < 4; ++i) {
    Float_t adc1 = esdAD->GetAdc(i+8);
    Float_t adc2 = esdAD->GetAdc(i+12);
    if (adc1 > minAdc && adc2 > minAdc) {
      Float_t time1 = esdAD->GetTime(i+8);
      Float_t time2 = esdAD->GetTime(i+12);
	if(time1 > (-1024+1.e-6) && time2 > (-1024+1.e-6)){
		Float_t timeErr1 = 1/adc1;
		Float_t timeErr2 = 1/adc2;
		Float_t timeDiff = TMath::Abs(time1-time2);
		Float_t timeSum = time1+time2;
		fHistTimePairSumDiffADA_NoCut->Fill(timeSum,timeDiff);
		Float_t earlyHitCut = 1000;
		if(TMath::Abs(timeSum - 2*fADADist) < 20) earlyHitCut = fEarlyHitCutShape->Eval(timeSum);
		if(timeDiff < earlyHitCut){
			itimeRobustADA[ntimeRobustADA] = i;
			ntimeRobustADA++;
			timeRobustADA += time1/(timeErr1*timeErr1);
			weightRobustADA += 1./(timeErr1*timeErr1);
			
			itimeRobustADA[ntimeRobustADA] = i+4;
			ntimeRobustADA++;
			timeRobustADA += time2/(timeErr2*timeErr2);
			weightRobustADA += 1./(timeErr2*timeErr2);
			fHistTimePairSumDiffADA_Cut->Fill(timeSum,timeDiff);
			fHistTimeVsChargeADA_Cut->Fill(time1,adc1);
			fHistTimeVsChargeADA_Cut->Fill(time2,adc2);
			}
		}
	else{
	 	if(time1 > (-1024+1.e-6)) fHistTimeVsChargeADA_Ex->Fill(time1,adc1);
		if(time2 > (-1024+1.e-6)) fHistTimeVsChargeADA_Ex->Fill(time2,adc2);
	 	}
	}
		
  }
  delete fEarlyHitCutShape;
  
  if(weightRobustADA > 1) timeRobustADA /= weightRobustADA; 
  else timeRobustADA = -1024.;
  if(weightRobustADC > 1) timeRobustADC /= weightRobustADC;
  else timeRobustADC = -1024.;
  
  fHistMeanTimeADA->Fill(timeBasicADA);
  fHistMeanTimeADC->Fill(timeBasicADC);
  if(timeBasicADA!= -1024 && timeBasicADC!= -1024){
  	fHistMeanTimeDifference->Fill(timeBasicADA-timeBasicADC);
	fHistMeanTimeCorrelation->Fill(timeBasicADA,timeBasicADC);
	fHistMeanTimeSumDiff->Fill(timeBasicADA-timeBasicADC,timeBasicADA+timeBasicADC);
	}
	
  //Decisions
  Int_t ADADecision=0;
  Int_t ADCDecision=0; 
  const Float_t TimeWindowBBALow = -2.0;
  const Float_t TimeWindowBBAUp = 2.0;
  const Float_t TimeWindowBBCLow = -1.0;
  const Float_t TimeWindowBBCUp = 1.0;
  const Float_t TimeWindowBGALow = -4.0;
  const Float_t TimeWindowBGAUp = 4.0;
  const Float_t TimeWindowBGCLow = -2.0;
  const Float_t TimeWindowBGCUp = 2.0;

  if(timeBasicADA > (fADADist + TimeWindowBBALow) && timeBasicADA < (fADADist + TimeWindowBBAUp)) ADADecision=1;
  else if(timeBasicADA > (-fADADist + TimeWindowBGALow) && timeBasicADA < (-fADADist + TimeWindowBGAUp)) ADADecision=2;
  else if(timeBasicADA>(-1024.+1.e-6)) ADADecision=3;
  else ADADecision=0;
  
  if(timeBasicADC > (fADCDist + TimeWindowBBCLow) && timeBasicADC < (fADCDist + TimeWindowBBCUp)) ADCDecision=1;
  else if(timeBasicADC > (-fADCDist + TimeWindowBGCLow) && timeBasicADC < (-fADCDist + TimeWindowBGCUp)) ADCDecision=2;
  else if(timeBasicADC>(-1024.+1.e-6)) ADCDecision=3;
  else ADCDecision=0;
	
  fHistDecisionBasic->Fill(ADADecision,ADCDecision);
  
  if(timeRobustADA > (fADADist + TimeWindowBBALow) && timeRobustADA < (fADADist + TimeWindowBBAUp)) ADADecision=1;
  else if(timeRobustADA > (-fADADist + TimeWindowBGALow) && timeRobustADA < (-fADADist + TimeWindowBGAUp)) ADADecision=2;
  else if(timeRobustADA>(-1024.+1.e-6)) ADADecision=3;
  else ADADecision=0;
  
  if(timeRobustADC > (fADCDist + TimeWindowBBCLow) && timeRobustADC < (fADCDist + TimeWindowBBCUp)) ADCDecision=1;
  else if(timeRobustADC > (-fADCDist + TimeWindowBGCLow) && timeRobustADC < (-fADCDist + TimeWindowBGCUp)) ADCDecision=2;
  else if(timeRobustADC>(-1024.+1.e-6)) ADCDecision=3;
  else ADCDecision=0;
	
  fHistDecisionRobust->Fill(ADADecision,ADCDecision);
  
  fHistDecision->Fill(esdAD->GetADADecision(),esdAD->GetADCDecision());
  
  fHistMedianTimeADA->Fill(medianTimeADA);
  fHistMedianTimeADC->Fill(medianTimeADC);
  fHistNTimesMedianADA->Fill(ntimeADA);
  fHistNTimesMedianADC->Fill(ntimeADC);
  
  fHistRobustTimeADA->Fill(timeRobustADA);
  fHistRobustTimeADC->Fill(timeRobustADC);
  fHistNTimesRobustADA->Fill(ntimeRobustADA);
  fHistNTimesRobustADC->Fill(ntimeRobustADC);


  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskADPilot::Terminate(Option_t *) 
{
 // Draw result to the screen
  // Called once at the end of the query

}
//________________________________________________________________________
TH1F * AliAnalysisTaskADPilot::CreateHist1D(const char* name, const char* title,Int_t nBins, 
				    Double_t xMin, Double_t xMax,
				    const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH1F* result = new TH1F(name, title, nBins, xMin, xMax);
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  return result;
}
//________________________________________________________________________
TH2F * AliAnalysisTaskADPilot::CreateHist2D(const char* name, const char* title,Int_t nBinsX, 
				    Double_t xMin, Double_t xMax,
				    Int_t nBinsY,
				    Double_t yMin, Double_t yMax,
				    const char* xLabel, const char* yLabel, const char* zLabel)
{
  // create a histogram
  TH2F* result = new TH2F(name, title, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
  result->SetOption("COLZ");
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  if (zLabel) result->GetYaxis()->SetTitle(zLabel);
  return result;
}

//________________________________________________________________________
TH3F * AliAnalysisTaskADPilot::CreateHist3D(const char* name, const char* title,Int_t nBinsX, 
				    Double_t xMin, Double_t xMax,
				    Int_t nBinsY,
				    Double_t yMin, Double_t yMax,
				    Int_t nBinsZ,
				    Double_t zMin, Double_t zMax,
				    const char* xLabel, const char* yLabel, const char* zLabel)
{
  // create a histogram
  TH3F* result = new TH3F(name, title, nBinsX, xMin, xMax, nBinsY, yMin, yMax, nBinsZ, zMin, zMax);
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  if (zLabel) result->GetYaxis()->SetTitle(zLabel);
  return result;
}

