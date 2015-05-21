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
//                 AliAnalysisTaskADQA class
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

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"

#include "AliAnalysisTaskADQA.h"

ClassImp(AliAnalysisTaskADQA)

//________________________________________________________________________
AliAnalysisTaskADQA::AliAnalysisTaskADQA() 
  : AliAnalysisTaskSE(),fListHist(0),fHistTotalChargePerEventADA(0),fHistTotalChargePerEventADC(0),
    fHistChargePerPM_All(0), fHistChargePerPM_BB(0), fHistChargePerPM_BG(0), fHistTimePerPM_Corr(0), fHistTimePerPM_UnCorr(0),
    fHistWidthPerPM(0),fHistTimeVsCharge_Corr(0),fHistTimeVsCharge_UnCorr(0),fHistWidthVsCharge(0),
    fHistNBBflagsADA(0),fHistNBBflagsADC(0),fHistNBBflagsADAVsADC(0),
    fHistNBBCoincidencesADA(0),fHistNBBCoincidencesADC(0),fHistNBBCoincidencesADAVsADC(0),
    fHistNBGflagsADA(0),fHistNBGflagsADC(0),fHistNBGflagsADAVsADC(0),
    fHistNBGCoincidencesADA(0),fHistNBGCoincidencesADC(0),fHistNBGCoincidencesADAVsADC(0),
    fHistChargeNoFlag(0),fHistTimeNoFlag(0), fHistChargeNoTime(0),
    fHistMeanTimeADA(0),fHistMeanTimeADC(0),fHistMeanTimeDifference(0),fHistMeanTimeCorrelation(0),fHistMeanTimeSumDiff(0),fHistDecision(0),
    fHistTrigger(0),
    fHistChargeVsClockInt0(0),fHistChargeVsClockInt1(0),fHistBBFlagVsClock(0),fHistBGFlagVsClock(0),fHistBBFlagPerChannel(0),fHistBGFlagPerChannel(0),fHistMaxChargeClock(0),
    fHistTimeVsChargePerPM_UnCorr(0)
{
  // Dummy constructor
}
//________________________________________________________________________
AliAnalysisTaskADQA::AliAnalysisTaskADQA(const char *name) 
  : AliAnalysisTaskSE(name),fListHist(0),fHistTotalChargePerEventADA(0),fHistTotalChargePerEventADC(0),
    fHistChargePerPM_All(0), fHistChargePerPM_BB(0), fHistChargePerPM_BG(0), fHistTimePerPM_Corr(0), fHistTimePerPM_UnCorr(0),
    fHistWidthPerPM(0),fHistTimeVsCharge_Corr(0),fHistTimeVsCharge_UnCorr(0),fHistWidthVsCharge(0),
    fHistNBBflagsADA(0),fHistNBBflagsADC(0),fHistNBBflagsADAVsADC(0),
    fHistNBBCoincidencesADA(0),fHistNBBCoincidencesADC(0),fHistNBBCoincidencesADAVsADC(0),
    fHistNBGflagsADA(0),fHistNBGflagsADC(0),fHistNBGflagsADAVsADC(0),
    fHistNBGCoincidencesADA(0),fHistNBGCoincidencesADC(0),fHistNBGCoincidencesADAVsADC(0),
    fHistChargeNoFlag(0),fHistTimeNoFlag(0), fHistChargeNoTime(0),
    fHistMeanTimeADA(0),fHistMeanTimeADC(0),fHistMeanTimeDifference(0),fHistMeanTimeCorrelation(0),fHistMeanTimeSumDiff(0),fHistDecision(0),
    fHistTrigger(0),
    fHistChargeVsClockInt0(0),fHistChargeVsClockInt1(0),fHistBBFlagVsClock(0),fHistBGFlagVsClock(0),fHistBBFlagPerChannel(0),fHistBGFlagPerChannel(0),fHistMaxChargeClock(0),
    fHistTimeVsChargePerPM_UnCorr(0)
{
  // Constructor
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskADQA::~AliAnalysisTaskADQA(){
  // Destructor
  if (fListHist) { delete fListHist; fListHist = 0x0; }
}
//________________________________________________________________________
void AliAnalysisTaskADQA::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  const Int_t kNTdcTimeBins = 3062;
  const Float_t kTdcTimeMin = 0.976562;
  const Float_t kTdcTimeMax = 300;
  const Int_t kNTdcWidthBins =  153;
  const Float_t kTdcWidthMin = 2.343750;
  const Float_t kTdcWidthMax = 121.875000;
  
  const Int_t kNChargeChannelBins = 5000;
  const Int_t kChargeChannelMin = 0;
  const Int_t kChargeChannelMax = 5000;
  
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
if (!fHistTimePerPM_Corr) {
    fHistTimePerPM_Corr = CreateHist2D("fHistTimePerPM_Corr","Corrected Time per PM",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax,"PM number","Leading time [ns]");
    fListHist->Add(fHistTimePerPM_Corr);
  } 
if (!fHistTimePerPM_UnCorr) {
    fHistTimePerPM_UnCorr = CreateHist2D("fHistTimePerPM_UnCorr","Raw Time per PM",kNChannelBins, kChannelMin, kChannelMax, kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax,"PM number","Leading time [ns]");
    fListHist->Add(fHistTimePerPM_UnCorr);
  }   
if (!fHistWidthPerPM) {
    fHistWidthPerPM = CreateHist2D("fHistWidthPerPM","Width per PM",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax,"PM number","Time width [ns]");
    fListHist->Add(fHistWidthPerPM);
  }  
if (!fHistTimeVsCharge_Corr) {
    fHistTimeVsCharge_Corr = CreateHist2D("fHistTimeVsCharge_Corr","Corrected Time vs Charge",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsCharge_Corr);
  }
if (!fHistTimeVsCharge_UnCorr) {
    fHistTimeVsCharge_UnCorr = CreateHist2D("fHistTimeVsCharge_UnCorr","Raw Time vs Charge",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
    fListHist->Add(fHistTimeVsCharge_UnCorr);
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
    fHistTimeNoFlag = CreateHist1D("fHistTimeNoFlag","Time in PM without BB flag",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax,"Time","Entries");
    fListHist->Add(fHistTimeNoFlag);
  }   
if (!fHistChargeNoTime) {
    fHistChargeNoTime = CreateHist1D("fHistChargeNoTime","Charge in PM without time measurement",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"Charge","Entries");
    fListHist->Add(fHistChargeNoTime);
  }

if (!fHistMeanTimeADA) {
    fHistMeanTimeADA = CreateHist1D("fHistMeanTimeADA","Mean Time in ADA",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax,"Time","Entries");
    fListHist->Add(fHistMeanTimeADA);
  }   

if (!fHistMeanTimeADC) {
    fHistMeanTimeADC = CreateHist1D("fHistMeanTimeADC","Mean Time in ADC",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax,"Time","Entries");
    fListHist->Add(fHistMeanTimeADC);
  }   

if (!fHistMeanTimeDifference) {
    fHistMeanTimeDifference = CreateHist1D("fHistMeanTimeDifference","Mean Time difference",1024,-150,150,"AD Mean time t_{A} - t_{C} [ns]","Entries");
    fListHist->Add(fHistMeanTimeDifference);
  }  
   
if (!fHistMeanTimeCorrelation) {
    fHistMeanTimeCorrelation = CreateHist2D("fHistMeanTimeCorrelation","Mean Time in ADA-ADC",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax,kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax,"Time ADA","Time ADC");
    fListHist->Add(fHistMeanTimeCorrelation);
  }   

if (!fHistMeanTimeSumDiff) {
    fHistMeanTimeSumDiff = CreateHist2D("fHistMeanTimeSumDiff","Mean Time in ADA-ADC",99, -150.0, 149.707031, 100, 0.0, 400.390625,"AD Mean time t_{A} - t_{C} [ns]","AD Mean time t_{A} + t_{C} [ns]");
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
    fHistDecision->GetXaxis()->SetBinLabel(2, "Fake");
    fHistDecision->GetXaxis()->SetBinLabel(3, "BB");
    fHistDecision->GetXaxis()->SetBinLabel(4, "BG");
    fHistDecision->GetYaxis()->SetBinLabel(1, "Empty");
    fHistDecision->GetYaxis()->SetBinLabel(2, "Fake");
    fHistDecision->GetYaxis()->SetBinLabel(3, "BB");
    fHistDecision->GetYaxis()->SetBinLabel(4, "BG");
    fListHist->Add(fHistDecision);
 }

    
if(!fHistTrigger) {
    fHistTrigger = CreateHist1D("fHistTrigger","Trigger inputs, from online Bits",11,0 ,11,"AD0 Trigger Type","Counts");
    fHistTrigger->SetFillColor(kAzure-8);
    fHistTrigger->SetLineWidth(2);
    fHistTrigger->GetXaxis()->SetLabelSize(0.04);
    fHistTrigger->GetXaxis()->SetNdivisions(808,kFALSE);
    fHistTrigger->GetXaxis()->SetBinLabel(1, "UBA");
    fHistTrigger->GetXaxis()->SetBinLabel(2, "UBC");
    fHistTrigger->GetXaxis()->SetBinLabel(3, "UGA");
    fHistTrigger->GetXaxis()->SetBinLabel(4, "UGC");
    fHistTrigger->GetXaxis()->SetBinLabel(5, "UBA & UBC");
    fHistTrigger->GetXaxis()->SetBinLabel(6, "UBA || UBC");
    fHistTrigger->GetXaxis()->SetBinLabel(7, "(UBA || UBC) & !(UGA || UGC)");
    fHistTrigger->GetXaxis()->SetBinLabel(8, "UGA & UBC");
    fHistTrigger->GetXaxis()->SetBinLabel(9, "UGC & UBA");
    fHistTrigger->GetXaxis()->SetBinLabel(10, "UGA || UGC");
    fHistTrigger->GetXaxis()->SetBinLabel(11, "(UGA & UBC) || (UGC & UBA)");
    fListHist->Add(fHistTrigger);
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
if(!fHistTimeVsChargePerPM_UnCorr) {
    fHistTimeVsChargePerPM_UnCorr = CreateHist3D("fHistTimeVsChargePerPM_UnCorr","Raw Time vs Charge per PM",kNTdcTimeBins, kTdcTimeMin, kTdcTimeMax, 
    					kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,
					kNChannelBins, kChannelMin, kChannelMax,"Leading time [ns]","ADC counts","Channel");
    fListHist->Add(fHistTimeVsChargePerPM_UnCorr);
  } 
   
  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskADQA::UserExec(Option_t *) 
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
	
	fHistTimePerPM_Corr->Fill(i,esdAD->GetTime(i));
	fHistTimeVsCharge_Corr->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
	fHistWidthPerPM->Fill(i,esdAD->GetWidth(i));
	fHistWidthVsCharge->Fill(esdAD->GetWidth(i),esdAD->GetAdc(i));
	if(!esdAD->GetBBFlag(i)) {
		fHistChargeNoFlag->Fill(esdAD->GetAdc(i));
		fHistTimeNoFlag->Fill(esdAD->GetTime(i));
		}
	if(esdAD->GetTime(i) < 1e-6) fHistChargeNoTime->Fill(esdAD->GetAdc(i));
	
	if(esdADfriend){
		fHistTimePerPM_UnCorr->Fill(i,esdADfriend->GetTime(i));
		fHistTimeVsCharge_UnCorr->Fill(esdADfriend->GetTime(i),esdAD->GetAdc(i));
		fHistTimeVsChargePerPM_UnCorr->Fill(esdADfriend->GetTime(i),esdAD->GetAdc(i),i);
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
		fHistMaxChargeClock->Fill(i,TMath::LocMax(21,charge)-10);
		
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
  Bool_t UBA = kFALSE;
  Bool_t UBC = kFALSE;
  Bool_t UGA = kFALSE;
  Bool_t UGC = kFALSE;
  
  if(nBBCoincidencesADA>=1) UBA = kTRUE;
  if(nBBCoincidencesADC>=1) UBC = kTRUE;
  if(nBGCoincidencesADA>=1) UGA = kTRUE;
  if(nBGCoincidencesADC>=1) UGC = kTRUE;

  if(UBA) fHistTrigger->Fill(0);
  if(UBC) fHistTrigger->Fill(1);
  if(UGA) fHistTrigger->Fill(2);
  if(UGC) fHistTrigger->Fill(3);
  if(UBA && UBC) fHistTrigger->Fill(4);
  if(UBA || UBC) fHistTrigger->Fill(5);
  if((UBA || UBC) && !(UGA || UGC)) fHistTrigger->Fill(6);
  if(UGA && UBC) fHistTrigger->Fill(7);
  if(UGC && UBA) fHistTrigger->Fill(8);
  if(UGA || UGC) fHistTrigger->Fill(9);
  if((UGA && UBC) || (UGC && UBA)) fHistTrigger->Fill(10);
  

  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskADQA::Terminate(Option_t *) 
{
 // Draw result to the screen
  // Called once at the end of the query

}
//________________________________________________________________________
TH1F * AliAnalysisTaskADQA::CreateHist1D(const char* name, const char* title,Int_t nBins, 
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
TH2F * AliAnalysisTaskADQA::CreateHist2D(const char* name, const char* title,Int_t nBinsX, 
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
TH3F * AliAnalysisTaskADQA::CreateHist3D(const char* name, const char* title,Int_t nBinsX, 
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

