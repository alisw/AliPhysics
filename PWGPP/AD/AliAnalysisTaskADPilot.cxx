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
#include <Riostream.h>
#include <bitset>

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
#include "AliAnalysisUtils.h"

#include "AliAnalysisTaskADPilot.h"

ClassImp(AliAnalysisTaskADPilot)

//________________________________________________________________________
AliAnalysisTaskADPilot::AliAnalysisTaskADPilot() 
  : AliAnalysisTaskSE(),
    	fListHist(0),			
	fHistTotalChargePerEventADA(0),
	fHistTotalChargePerEventADC(0),
	fHistChargePerPM_All(0),
	fHistChargePerPM_BG(0),
	fHistTimePerPM_Corr(0),
	fHistTimeVsChargeADA_Corr(0),
	fHistTimeVsChargeADC_Corr(0),
	fHistTimeVsChargeADA_Cut(0),
	fHistTimeVsChargeADC_Cut(0),
	fHistWidthPerPM(0),
	fHistWidthVsCharge(0),
	fHistNBBflagsADA(0),
	fHistNBBflagsADC(0),
	fHistNBBflagsADAVsADC(0),
	fHistNBGflagsADA(0),
	fHistNBGflagsADC(0),
	fHistNBGflagsADAVsADC(0),
	fHistNBBCoincidencesADA(0),
	fHistNBBCoincidencesADC(0),
	fHistNBBCoincidencesADAVsADC(0),
	fHistNBGCoincidencesADA(0),
	fHistNBGCoincidencesADC(0),
	fHistNBGCoincidencesADAVsADC(0),
	fHistChargeNoFlag(0),
	fHistTimeNoFlag(0),
	fHistChargeNoTime(0),
	fHistFlagNoTime(0), 
	fHistChargePerCoincidence(0),
	fHistMeanTimeADA(0),
	fHistMeanTimeADC(0),
	fHistMeanTimeDifference(0),
	fHistMeanTimeCorrelation(0),
	fHistMeanTimeSumDiff(0),
	fHistDecision(0),
	fHistTriggerMasked(0),
	fHistTriggerUnMasked(0),
	fHistTriggerOthers(0),
	fHistChargeVsClockInt0(0),
	fHistChargeVsClockInt1(0),
	fHistBBFlagVsClock(0),
	fHistBGFlagVsClock(0),
	fHistBBFlagPerChannel(0),
	fHistBGFlagPerChannel(0),
	fHistMaxChargeClock(0),
	fHistMaxChargeValueInt0(0),
	fHistMaxChargeValueInt1(0),
	fHistTimePerPM_UnCorr(0),
	fHistTimeVsChargeADA_UnCorr(0),
	fHistTimeVsChargeADC_UnCorr(0),
	fHistTimeVsChargePerPM_UnCorr(0),
	//fHistTimeVsChargePerPM_UnCorr_Cut(0),
	//fHistTimePairCorrelation(0),
	//fHistTimePairCorrelation_Cut(0),
	fTimeSlewingTree(0),
	fHistChargeTriggerPerPMPerV0Flag(0),
	fHistChargeTailPerPMPerV0Flag(0),
	fHistChargeBBPerPMPerV0Flag(0), 
	fHistChargeTriggerPerChannel(0),
	fHistChargeTriggerPerChannel_PF(0),
	fHistChargeTriggerPerChannel_TVX(0),
	fHistChargeTriggerPerChannel_PF_TVX(0),
	fHistChargeTailPerChannel(0),
	fHistChargeTailPerChannel_PF(0),
	fHistChargeTailPerChannel_TVX(0),
	fHistChargeTailPerChannel_PF_TVX(0), 
	fHistChargePerPM_BB(0), 
	fHistChargePerPM_BB_PF(0),
	fHistChargePerPM_BB_TVX(0),  
	fHistChargePerPM_BB_PF_TVX(0),
	fHistChargePerPM_Time(0), 
	fHistChargePerPM_Time_PF(0),
	fHistChargePerPM_Time_TVX(0),  
	fHistChargePerPM_Time_PF_TVX(0),
    	fRun(0),fOldRun(0),fCalibData(0),fAnalysisUtils(0)
{
  // Dummy constructor
}
//________________________________________________________________________
AliAnalysisTaskADPilot::AliAnalysisTaskADPilot(const char *name) 
  : AliAnalysisTaskSE(name),
  	fListHist(0),			
	fHistTotalChargePerEventADA(0),
	fHistTotalChargePerEventADC(0),
	fHistChargePerPM_All(0),
	fHistChargePerPM_BG(0),
	fHistTimePerPM_Corr(0),
	fHistTimeVsChargeADA_Corr(0),
	fHistTimeVsChargeADC_Corr(0),
	fHistTimeVsChargeADA_Cut(0),
	fHistTimeVsChargeADC_Cut(0),
	fHistWidthPerPM(0),
	fHistWidthVsCharge(0),
	fHistNBBflagsADA(0),
	fHistNBBflagsADC(0),
	fHistNBBflagsADAVsADC(0),
	fHistNBGflagsADA(0),
	fHistNBGflagsADC(0),
	fHistNBGflagsADAVsADC(0),
	fHistNBBCoincidencesADA(0),
	fHistNBBCoincidencesADC(0),
	fHistNBBCoincidencesADAVsADC(0),
	fHistNBGCoincidencesADA(0),
	fHistNBGCoincidencesADC(0),
	fHistNBGCoincidencesADAVsADC(0),
	fHistChargeNoFlag(0),
	fHistTimeNoFlag(0),
	fHistChargeNoTime(0),
	fHistFlagNoTime(0), 
	fHistChargePerCoincidence(0),
	fHistMeanTimeADA(0),
	fHistMeanTimeADC(0),
	fHistMeanTimeDifference(0),
	fHistMeanTimeCorrelation(0),
	fHistMeanTimeSumDiff(0),
	fHistDecision(0),
	fHistTriggerMasked(0),
	fHistTriggerUnMasked(0),
	fHistTriggerOthers(0),
	fHistChargeVsClockInt0(0),
	fHistChargeVsClockInt1(0),
	fHistBBFlagVsClock(0),
	fHistBGFlagVsClock(0),
	fHistBBFlagPerChannel(0),
	fHistBGFlagPerChannel(0),
	fHistMaxChargeClock(0),
	fHistMaxChargeValueInt0(0),
	fHistMaxChargeValueInt1(0),
	fHistTimePerPM_UnCorr(0),
	fHistTimeVsChargeADA_UnCorr(0),
	fHistTimeVsChargeADC_UnCorr(0),
	fHistTimeVsChargePerPM_UnCorr(0),
	//fHistTimeVsChargePerPM_UnCorr_Cut(0),
	//fHistTimePairCorrelation(0),
	//fHistTimePairCorrelation_Cut(0),
	fTimeSlewingTree(0),
	fHistChargeTriggerPerPMPerV0Flag(0),
	fHistChargeTailPerPMPerV0Flag(0),
	fHistChargeBBPerPMPerV0Flag(0),
	fHistChargeTriggerPerChannel(0),
	fHistChargeTriggerPerChannel_PF(0),
	fHistChargeTriggerPerChannel_TVX(0),
	fHistChargeTriggerPerChannel_PF_TVX(0), 
	fHistChargeTailPerChannel(0),
	fHistChargeTailPerChannel_PF(0),
	fHistChargeTailPerChannel_TVX(0),
	fHistChargeTailPerChannel_PF_TVX(0),
	fHistChargePerPM_BB(0), 
	fHistChargePerPM_BB_PF(0),
	fHistChargePerPM_BB_TVX(0),  
	fHistChargePerPM_BB_PF_TVX(0),
	fHistChargePerPM_Time(0), 
	fHistChargePerPM_Time_PF(0),
	fHistChargePerPM_Time_TVX(0),  
	fHistChargePerPM_Time_PF_TVX(0),
    	fRun(0),fOldRun(0),fCalibData(0),fAnalysisUtils(new AliAnalysisUtils)
{
  // Constructor
  // Output slot #1 writes into a TList container
  for(Int_t i=0; i<16; i++){ fAdc[i] = 0; fTdc[i] = 0; fAdcTail[i] = 0;}
  
  fAnalysisUtils->SetUseOutOfBunchPileUp(kTRUE);
 
  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}
//________________________________________________________________________
AliAnalysisTaskADPilot::~AliAnalysisTaskADPilot(){
  // Destructor
  if (fListHist) { delete fListHist; fListHist = 0x0; }
  if(fTimeSlewingTree){ delete fTimeSlewingTree; fTimeSlewingTree = 0x0; }
}
//________________________________________________________________________
void AliAnalysisTaskADPilot::UserCreateOutputObjects()
{

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
  
  fTimeSlewingTree = new TTree("fTimeSlewingTree", "fTimeSlewingTree");
  fTimeSlewingTree ->Branch("fAdc", &fAdc[0], "fAdc[16]/s");
  fTimeSlewingTree ->Branch("fAdcTail", &fAdcTail[0], "fAdcTail[16]/s");
  fTimeSlewingTree ->Branch("fTdc", &fTdc[0], "fTdc[16]/s");
  

  fListHist = new TList();
  fListHist->SetOwner();
  
  fHistChargeTriggerPerPMPerV0Flag = CreateHist3D("fHistChargeTriggerPerPMPerV0Flag","Charge per PM per VZERO flag",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,
					65,-0.5,64.5,"ADC counts","Channel","Number of VZERO BB flags");
  fListHist->Add(fHistChargeTriggerPerPMPerV0Flag);
  
  fHistChargeTailPerPMPerV0Flag = CreateHist3D("fHistChargeTailPerPMPerV0Flag","Charge per PM per VZERO flag",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,
					65,-0.5,64.5,"ADC counts","Channel","Number of VZERO BB flags");
  fListHist->Add(fHistChargeTailPerPMPerV0Flag);
  
  fHistChargeBBPerPMPerV0Flag = CreateHist3D("fHistChargeBBPerPMPerV0Flag","Charge per PM per VZERO flag",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,
					65,-0.5,64.5,"ADC counts","Channel","Number of VZERO BB flags");
  fListHist->Add(fHistChargeBBPerPMPerV0Flag);
  
  fHistChargeTriggerPerChannel = CreateHist2D("fHistChargeTriggerPerChannel","Trigger charge per channel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  fListHist->Add(fHistChargeTriggerPerChannel);
  
  fHistChargeTriggerPerChannel_PF = CreateHist2D("fHistChargeTriggerPerChannel_PF","Trigger charge per channel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  fListHist->Add(fHistChargeTriggerPerChannel_PF); 
  
  fHistChargeTriggerPerChannel_TVX = CreateHist2D("fHistChargeTriggerPerChannel_TVX","Trigger charge per channel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  fListHist->Add(fHistChargeTriggerPerChannel_TVX);
  
  fHistChargeTriggerPerChannel_PF_TVX = CreateHist2D("fHistChargeTriggerPerChannel_PF_TVX","Trigger charge per channel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  fListHist->Add(fHistChargeTriggerPerChannel_PF_TVX);
  
  fHistChargeTailPerChannel = CreateHist2D("fHistChargeTailPerChannel","Tail charge per channel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  fListHist->Add(fHistChargeTailPerChannel);
  
  fHistChargeTailPerChannel_PF = CreateHist2D("fHistChargeTailPerChannel_PF","Tail charge per channel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  fListHist->Add(fHistChargeTailPerChannel_PF); 
  
  fHistChargeTailPerChannel_TVX = CreateHist2D("fHistChargeTailPerChannel_TVX","Tail charge per channel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  fListHist->Add(fHistChargeTailPerChannel_TVX);
  
  fHistChargeTailPerChannel_PF_TVX = CreateHist2D("fHistChargeTailPerChannel_PF_TVX","Tail charge per channel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  fListHist->Add(fHistChargeTailPerChannel_PF_TVX);
  
  fHistChargePerPM_BB = CreateHist2D("fHistChargePerPM_BB","Charge per PM BB events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  fListHist->Add(fHistChargePerPM_BB);
  
  fHistChargePerPM_Time = CreateHist2D("fHistChargePerPM_Time","Charge per PM Time events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  fListHist->Add(fHistChargePerPM_Time);
  
  fHistChargePerPM_BB_PF = CreateHist2D("fHistChargePerPM_BB_PF","Charge per PM BB events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  fListHist->Add(fHistChargePerPM_BB_PF);
  
  fHistChargePerPM_Time_PF = CreateHist2D("fHistChargePerPM_Time_PF","Charge per PM Time events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  fListHist->Add(fHistChargePerPM_Time_PF);
  
  fHistChargePerPM_BB_TVX = CreateHist2D("fHistChargePerPM_BB_TVX","Charge per PM BB events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  fListHist->Add(fHistChargePerPM_BB_TVX);
  
  fHistChargePerPM_Time_TVX = CreateHist2D("fHistChargePerPM_Time_TVX","Charge per PM Time events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  fListHist->Add(fHistChargePerPM_Time_TVX);
  
  fHistChargePerPM_BB_PF_TVX = CreateHist2D("fHistChargePerPM_BB_PF_TVX","Charge per PM BB events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  fListHist->Add(fHistChargePerPM_BB_PF_TVX);
  
  fHistChargePerPM_Time_PF_TVX = CreateHist2D("fHistChargePerPM_Time_PF_TVX","Charge per PM Time events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  fListHist->Add(fHistChargePerPM_Time_PF_TVX); 

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
if (!fHistChargePerPM_BG) {
    fHistChargePerPM_BG = CreateHist2D("fHistChargePerPM_BG","Charge per PM BG events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
    fListHist->Add(fHistChargePerPM_BG);
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
    fHistDecision = CreateHist2D("fHistDecision","Offline decision in ADA-ADC",4,0 ,4,4,0,4,"ADA","ADC");
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
					kNChannelBins, kChannelMin, kChannelMax,"Leading time [TDC counts]","ADC counts","Channel");
    fListHist->Add(fHistTimeVsChargePerPM_UnCorr);
  } 
/*  
if(!fHistTimeVsChargePerPM_UnCorr_Cut) {
    fHistTimeVsChargePerPM_UnCorr_Cut = CreateHist3D("fHistTimeVsChargePerPM_UnCorr_Cut","Raw Time vs Charge per PM",2600, 400, 3000, 
    					200,-4,0,
					kNChannelBins, kChannelMin, kChannelMax,"Leading time [ns]","ADC counts","Channel");
    fListHist->Add(fHistTimeVsChargePerPM_UnCorr_Cut);
  }
  
if(!fHistTimePairCorrelation) {
    fHistTimePairCorrelation = CreateHist3D("fHistTimePairCorrelation"," per Pair",1200, 3400, 4600, 
    					500,-250,250,
					8, -0.5, 7.5,"t1+t2[TDC counts]","t1-t2[TDC counts]","Channel");
    fListHist->Add(fHistTimePairCorrelation);
  }
if(!fHistTimePairCorrelation_Cut) {
    fHistTimePairCorrelation_Cut = CreateHist3D("fHistTimePairCorrelation_Cut"," per Pair",1200, 3400, 4600, 
    					500,-250,250,
					8, -0.5, 7.5,"t1+t2[TDC counts]","t1-t2[TDC counts]","Channel");
    fListHist->Add(fHistTimePairCorrelation_Cut);
  }        */

  // Post output data.
  PostData(1, fListHist);
  PostData(2, fTimeSlewingTree);
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
  if(fInputEvent->IsIncompleteDAQ()) return;
  
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
  AliESDADfriend* esdADfriend = 0x0;
  if(fESDfriend) esdADfriend = fESDfriend->GetADfriend();
  
  if(fAnalysisUtils->IsPileUpEvent(fESD)) return;
  if(fAnalysisUtils->IsSPDClusterVsTrackletBG(fESD)) return;
  
  fRun=fEvent->GetRunNumber();
  
  /*/
  if (fRun!=fOldRun){
    SetCalibData();
    fOldRun=fRun;
  }
  /*/
  
  //Triggers from VZERO for reference
  AliESDVZERO* esdVZERO = fESD->GetVZEROData();
  if(esdVZERO){
  	UShort_t fTriggerVZERO = esdVZERO->GetTriggerBits();
	if(fTriggerVZERO & (1 << 0) ? kTRUE : kFALSE) fHistTriggerOthers->Fill(0);
	if(fTriggerVZERO & (1 << 1) ? kTRUE : kFALSE) fHistTriggerOthers->Fill(1);
  	}
  Int_t nVZEROBBflags = 0;
  for(Int_t i=0; i<64; i++) nVZEROBBflags += esdVZERO->GetBBFlag(i);
  
  Float_t totChargeADA = 0;
  Float_t totChargeADC = 0;
  Int_t nBBflagsADA = 0;
  Int_t nBBflagsADC = 0;
  Int_t nBGflagsADA = 0;
  Int_t nBGflagsADC = 0;
  
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
	if(esdAD->GetTime(i)> -1024 + 1e-6)fHistChargePerPM_Time->Fill(i,esdAD->GetAdc(i));
	if(fESD->GetHeader()->IsTriggerInputFired("0TVX")){
		if(esdAD->GetBBFlag(i))fHistChargePerPM_BB_TVX->Fill(i,esdAD->GetAdc(i));
		if(esdAD->GetTime(i)> -1024 + 1e-6)fHistChargePerPM_Time_TVX->Fill(i,esdAD->GetAdc(i));
		if(esdAD->GetBBFlag(i)) fHistChargeBBPerPMPerV0Flag->Fill(i,esdAD->GetAdc(i),nVZEROBBflags);
		}
	
	fHistTimePerPM_Corr->Fill(i,esdAD->GetTime(i));
	
	if(i<8){
		fHistTimeVsChargeADC_Corr->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
		if(esdAD->BBTriggerADC(i)) fHistTimeVsChargeADC_Cut->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
		}
        if(i>7){
		fHistTimeVsChargeADA_Corr->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
		if(esdAD->BBTriggerADA(i-8)) fHistTimeVsChargeADA_Cut->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
		}
	fHistWidthPerPM->Fill(i,esdAD->GetWidth(i));
	fHistWidthVsCharge->Fill(esdAD->GetWidth(i),esdAD->GetAdc(i));
	
	
	if(esdAD->GetTime(i) < -1024+1e-6) fHistChargeNoTime->Fill(i,esdAD->GetAdc(i));
	if(esdAD->GetTime(i) < -1024+1e-6 && (esdAD->GetBBFlag(i)||esdAD->GetBGFlag(i)))fHistFlagNoTime->Fill(i);
	
	if(esdADfriend){
		if(esdAD->GetAdc(i)>0)fHistTimeVsChargePerPM_UnCorr->Fill(TMath::Nint(esdADfriend->GetTime(i)/(25./256.)),TMath::Log10(1.0/esdAD->GetAdc(i)),i);
		fAdc[i] = esdAD->GetAdc(i);
		fTdc[i] = TMath::Nint(esdADfriend->GetTime(i)/(25./256.));
		fAdcTail[i] = esdAD->GetAdcTail(i);
		
		if(!esdAD->GetBBFlag(i)) {
			fHistChargeNoFlag->Fill(esdAD->GetAdc(i));
			fHistTimeNoFlag->Fill(i,esdADfriend->GetTime(i));
			}
		fHistTimePerPM_UnCorr->Fill(i,esdADfriend->GetTime(i));
 		if(i<8) fHistTimeVsChargeADC_UnCorr->Fill(esdADfriend->GetTime(i),esdAD->GetAdc(i));
 		if(i>7) fHistTimeVsChargeADA_UnCorr->Fill(esdADfriend->GetTime(i),esdAD->GetAdc(i));
		Int_t nbbFlag = 0;
      		Int_t nbgFlag = 0;
		Int_t charge[21];
		for(Int_t iClock=0; iClock<21; iClock++){
			charge[iClock] = esdADfriend->GetPedestal(i,iClock);
			Bool_t intgr = esdADfriend->GetIntegratorFlag(i,iClock);
			Bool_t bbFlag = esdADfriend->GetBBFlag(i, iClock);
			Bool_t bgFlag = esdADfriend->GetBGFlag(i,iClock);
			if(bbFlag) nbbFlag++;
			if(bgFlag) nbgFlag++;
	
			if(!intgr)fHistChargeVsClockInt0->Fill(i,iClock-10,charge[iClock]);
			if(intgr)fHistChargeVsClockInt1->Fill(i,iClock-10,charge[iClock]);
			fHistBBFlagVsClock->Fill(i,iClock-10, bbFlag);
			fHistBGFlagVsClock->Fill(i,iClock-10, bgFlag);
			}
		fHistBBFlagPerChannel->Fill(i,nbbFlag);
		fHistBGFlagPerChannel->Fill(i,nbgFlag);
		
		Bool_t localPF = kTRUE;
		for(Int_t iClock=0; iClock<10; iClock++) if(esdADfriend->GetBBFlag(i,iClock) || esdADfriend->GetBGFlag(i,iClock))localPF = kFALSE;
		for(Int_t iClock=11; iClock<21; iClock++) if(esdADfriend->GetBBFlag(i,iClock) || esdADfriend->GetBGFlag(i,iClock))localPF = kFALSE;
		
		Float_t triggerCharge = esdAD->GetAdcTrigger(i);
		Float_t tailCharge = esdAD->GetAdcTail(i);
		
		if(esdAD->GetBBFlag(i)){
			fHistChargeTriggerPerChannel->Fill(i,triggerCharge);
			fHistChargeTailPerChannel->Fill(i,tailCharge);
			if(localPF){
				fHistChargeTriggerPerChannel_PF->Fill(i,triggerCharge);
				fHistChargePerPM_BB_PF->Fill(i,esdAD->GetAdc(i));
				fHistChargePerPM_Time_PF->Fill(i,esdAD->GetAdc(i));
				fHistChargeTailPerChannel_PF->Fill(i,tailCharge);
				}
			if(fESD->GetHeader()->IsTriggerInputFired("0TVX")){ 
				fHistChargeTriggerPerChannel_TVX->Fill(i,triggerCharge);
				fHistChargeTailPerChannel_TVX->Fill(i,tailCharge);
				if(esdAD->GetBBFlag(i)) fHistChargeTriggerPerPMPerV0Flag->Fill(i,triggerCharge,nVZEROBBflags);
				if(esdAD->GetBBFlag(i)) fHistChargeTailPerPMPerV0Flag->Fill(i,tailCharge,nVZEROBBflags);
				}
			if(fESD->GetHeader()->IsTriggerInputFired("0TVX") && localPF){
				fHistChargeTriggerPerChannel_PF_TVX->Fill(i,triggerCharge);
				fHistChargePerPM_BB_PF_TVX->Fill(i,esdAD->GetAdc(i));
				fHistChargePerPM_Time_PF_TVX->Fill(i,esdAD->GetAdc(i));
				fHistChargeTailPerChannel_PF_TVX->Fill(i,tailCharge);
				}
			} 
		
		Int_t maxClock = TMath::LocMax(21,charge);
		fHistMaxChargeClock->Fill(i,maxClock-10);
		if(!esdADfriend->GetIntegratorFlag(i,maxClock))fHistMaxChargeValueInt0->Fill(i,charge[maxClock]);
		if( esdADfriend->GetIntegratorFlag(i,maxClock))fHistMaxChargeValueInt1->Fill(i,charge[maxClock]);
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
  fHistMeanTimeADA->Fill(esdAD->GetADATime());
  fHistMeanTimeADC->Fill(esdAD->GetADCTime());
  if(esdAD->GetADATime()!= -1024 && esdAD->GetADCTime()!= -1024){
  	fHistMeanTimeDifference->Fill(esdAD->GetADATime()-esdAD->GetADCTime());
	fHistMeanTimeCorrelation->Fill(esdAD->GetADATime(),esdAD->GetADCTime());
	fHistMeanTimeSumDiff->Fill(esdAD->GetADATime()-esdAD->GetADCTime(),esdAD->GetADATime()+esdAD->GetADCTime());
	}
  fHistDecision->Fill(esdAD->GetADADecision(),esdAD->GetADCDecision());
    
  //Triggers
  UShort_t fTriggerBC = esdAD->GetTriggerBits();
  for(Int_t i = 0; i<6; i++) if(fTriggerBC & (1 << i) ? kTRUE : kFALSE) fHistTriggerMasked->Fill(i);
  for(Int_t i = 12; i<16; i++) if(fTriggerBC & (1 << i) ? kTRUE : kFALSE) fHistTriggerMasked->Fill(i-6);
  
  if(esdADfriend){
  	UShort_t fTriggerUnBC = esdADfriend->GetTriggerInputs();
  	for(Int_t i = 0; i<6; i++) if(fTriggerUnBC & (1 << i) ? kTRUE : kFALSE) fHistTriggerUnMasked->Fill(i);
	for(Int_t i = 12; i<16; i++) if(fTriggerUnBC & (1 << i) ? kTRUE : kFALSE) fHistTriggerUnMasked->Fill(i-6);	
  }
  fTimeSlewingTree->Fill();
  // Post output data.
  PostData(1, fListHist);
  PostData(2, fTimeSlewingTree);
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

