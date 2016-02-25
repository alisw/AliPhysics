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
//                 AliAnalysisTaskADVVQA class
//            This task is for QAing the AD data from ESD/AOD
//		with VZERO veto events
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
#include "AliADCalibData.h"

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"
#include "AliESDVZERO.h"
#include "AliAnalysisUtils.h"

#include "AliAnalysisTaskADVVQA.h"

ClassImp(AliAnalysisTaskADVVQA)

//________________________________________________________________________
AliAnalysisTaskADVVQA::AliAnalysisTaskADVVQA() 
  : AliAnalysisTaskSE(),fListHist(0),fList_NoVBA_NoVBC(0),fList_VBA_NoVBC(0),fList_NoVBA_VBC(0),fList_VBA_VBC(0),fList_TVX(0),fList_UBA_UBC(0),
  fRun(0),fOldRun(0),fCalibData(0),fAnalysisUtils(0)
{
  // Dummy constructor
}
//________________________________________________________________________
AliAnalysisTaskADVVQA::AliAnalysisTaskADVVQA(const char *name) 
  : AliAnalysisTaskSE(name),fListHist(0),fList_NoVBA_NoVBC(0),fList_VBA_NoVBC(0),fList_NoVBA_VBC(0),fList_VBA_VBC(0),fList_TVX(0),fList_UBA_UBC(0),
  fRun(0),fOldRun(0),fCalibData(0),fAnalysisUtils(new AliAnalysisUtils)
{
  // Constructor
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskADVVQA::~AliAnalysisTaskADVVQA(){
  // Destructor
  if (fListHist) { delete fListHist; fListHist = 0x0; }
  if (fList_NoVBA_NoVBC) { delete fList_NoVBA_NoVBC; fList_NoVBA_NoVBC = 0x0; }
  if (fList_VBA_NoVBC) { delete fList_VBA_NoVBC; fList_VBA_NoVBC = 0x0; }
  if (fList_NoVBA_VBC) { delete fList_NoVBA_VBC; fList_NoVBA_VBC = 0x0; }
  if (fList_VBA_VBC) { delete fList_VBA_VBC; fList_VBA_VBC = 0x0; }
  if (fList_TVX) { delete fList_TVX; fList_TVX = 0x0; }
  if (fList_UBA_UBC) { delete fList_UBA_UBC; fList_UBA_UBC = 0x0; }
}
//________________________________________________________________________
void AliAnalysisTaskADVVQA::UserCreateOutputObjects()
{

  fListHist = new TList();
  fListHist->SetOwner();
  fListHist->SetName("fListHist");
  
  fList_NoVBA_NoVBC = new TList();
  fList_NoVBA_NoVBC->SetOwner();
  fList_NoVBA_NoVBC->SetName("Hists_NoVBA_NoVBC");
  fListHist->Add(fList_NoVBA_NoVBC);
  
  InitHistos(fList_NoVBA_NoVBC,"NoVBA_NoVBC");
  
  fList_VBA_NoVBC = new TList();
  fList_VBA_NoVBC->SetOwner();
  fList_VBA_NoVBC->SetName("Hists_VBA_NoVBC");
  fListHist->Add(fList_VBA_NoVBC);
  
  InitHistos(fList_VBA_NoVBC,"VBA_NoVBC");
  
  fList_NoVBA_VBC = new TList();
  fList_NoVBA_VBC->SetOwner();
  fList_NoVBA_VBC->SetName("Hists_NoVBA_VBC");
  fListHist->Add(fList_NoVBA_VBC);
  
  InitHistos(fList_NoVBA_VBC,"NoVBA_VBC");
  
  fList_VBA_VBC = new TList();
  fList_VBA_VBC->SetOwner();
  fList_VBA_VBC->SetName("Hists_VBA_VBC");
  fListHist->Add(fList_VBA_VBC);
  
  InitHistos(fList_VBA_VBC,"VBA_VBC");
  
  fList_TVX = new TList();
  fList_TVX->SetOwner();
  fList_TVX->SetName("Hists_TVX");
  fListHist->Add(fList_TVX);
  
  InitHistos(fList_TVX,"TVX");
  
  fList_UBA_UBC = new TList();
  fList_UBA_UBC->SetOwner();
  fList_UBA_UBC->SetName("Hists_UBA_UBC");
  fListHist->Add(fList_UBA_UBC);
  
  InitHistos(fList_UBA_UBC,"UBA_UBC");

  // Post output data.
  PostData(1, fListHist);
}
//________________________________________________________________________
void AliAnalysisTaskADVVQA::InitHistos(TList *list,const char* TriggerName)
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
  
  
  TString fHistTotalChargePerEventADA_Name = "fHistTotalChargePerEventADA_";
  fHistTotalChargePerEventADA_Name += TriggerName;
  TH1F *fHistTotalChargePerEventADA = CreateHist1D(fHistTotalChargePerEventADA_Name.Data(),"Total Charge in ADA per event",kNChargeSideBins,kChargeSideMin,kChargeSideMax,"ADC counts","Entries");
  list->Add(fHistTotalChargePerEventADA);
  
  TString fHistTotalChargePerEventADC_Name = "fHistTotalChargePerEventADC_";
  fHistTotalChargePerEventADC_Name += TriggerName;
  TH1F *fHistTotalChargePerEventADC = CreateHist1D(fHistTotalChargePerEventADC_Name.Data(),"Total Charge in ADC per event",kNChargeSideBins,kChargeSideMin,kChargeSideMax,"ADC counts","Entries");
  list->Add(fHistTotalChargePerEventADC);
  
  TString fHistChargePerPM_All_Name = "fHistChargePerPM_All_";
  fHistChargePerPM_All_Name += TriggerName;     
  TH2F *fHistChargePerPM_All = CreateHist2D(fHistChargePerPM_All_Name.Data(),"Charge per PM all events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  list->Add(fHistChargePerPM_All);
 
  TString fHistChargePerPM_BB_Name = "fHistChargePerPM_BB_";
  fHistChargePerPM_BB_Name += TriggerName;     
  TH2F *fHistChargePerPM_BB = CreateHist2D(fHistChargePerPM_BB_Name.Data(),"Charge per PM BB events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  list->Add(fHistChargePerPM_BB);
  
  TString fHistChargePerPM_BG_Name = "fHistChargePerPM_BG_";
  fHistChargePerPM_BG_Name += TriggerName;
  TH2F *fHistChargePerPM_BG = CreateHist2D(fHistChargePerPM_BG_Name.Data(),"Charge per PM BG events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  list->Add(fHistChargePerPM_BG);
  
  TString fHistChargePerPM_Time_Name = "fHistChargePerPM_Time_";
  fHistChargePerPM_Time_Name += TriggerName;
  TH2F *fHistChargePerPM_Time = CreateHist2D(fHistChargePerPM_Time_Name.Data(),"Charge per PM Time events",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"PM number","ADC counts");
  list->Add(fHistChargePerPM_Time);
  
  TString fHistTimePerPM_Corr_Name = "fHistTimePerPM_Corr_";
  fHistTimePerPM_Corr_Name += TriggerName;  
  TH2F *fHistTimePerPM_Corr = CreateHist2D(fHistTimePerPM_Corr_Name.Data(),"Corrected Time per PM",kNChannelBins, kChannelMin, kChannelMax, kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"PM number","Leading time [ns]");
  list->Add(fHistTimePerPM_Corr);
  
  TString fHistTimePerPM_UnCorr_Name = "fHistTimePerPM_UnCorr_";
  fHistTimePerPM_UnCorr_Name += TriggerName;
  TH2F *fHistTimePerPM_UnCorr = CreateHist2D(fHistTimePerPM_UnCorr_Name.Data(),"Raw Time per PM",kNChannelBins, kChannelMin, kChannelMax, kNRawTimeBins, kRawTimeMin, kRawTimeMax,"PM number","Leading time [ns]");
  list->Add(fHistTimePerPM_UnCorr);
  
  TString fHistWidthPerPM_Name = "fHistWidthPerPM_";
  fHistWidthPerPM_Name += TriggerName;
  TH2F *fHistWidthPerPM = CreateHist2D(fHistWidthPerPM_Name.Data(),"Width per PM",kNChannelBins, kChannelMin, kChannelMax, kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax,"PM number","Time width [ns]");
  list->Add(fHistWidthPerPM);
  
  TString fHistTimeVsChargeADA_Cut_Name = "fHistTimeVsChargeADA_Cut_";
  fHistTimeVsChargeADA_Cut_Name += TriggerName;
  TH2F *fHistTimeVsChargeADA_Cut = CreateHist2D(fHistTimeVsChargeADA_Cut_Name.Data(),"Cutted Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
  list->Add(fHistTimeVsChargeADA_Cut);
  
  TString fHistTimeVsChargeADA_Corr_Name = "fHistTimeVsChargeADA_Corr_";
  fHistTimeVsChargeADA_Corr_Name += TriggerName; 
  TH2F *fHistTimeVsChargeADA_Corr = CreateHist2D(fHistTimeVsChargeADA_Corr_Name.Data(),"Corrected Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
  list->Add(fHistTimeVsChargeADA_Corr);
  
  TString fHistTimeVsChargeADA_UnCorr_Name = "fHistTimeVsChargeADA_UnCorr_";
  fHistTimeVsChargeADA_UnCorr_Name += TriggerName; 
  TH2F *fHistTimeVsChargeADA_UnCorr = CreateHist2D(fHistTimeVsChargeADA_UnCorr_Name.Data(),"Raw Time vs Charge",kNRawTimeBins, kRawTimeMin, kRawTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
  list->Add(fHistTimeVsChargeADA_UnCorr);

  TString fHistTimeVsChargeADC_Cut_Name = "fHistTimeVsChargeADC_Cut_";
  fHistTimeVsChargeADC_Cut_Name += TriggerName;
  TH2F *fHistTimeVsChargeADC_Cut = CreateHist2D(fHistTimeVsChargeADC_Cut_Name.Data(),"Cutted Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
  list->Add(fHistTimeVsChargeADC_Cut);

  TString fHistTimeVsChargeADC_Corr_Name = "fHistTimeVsChargeADC_Corr_";
  fHistTimeVsChargeADC_Corr_Name += TriggerName;
  TH2F *fHistTimeVsChargeADC_Corr = CreateHist2D(fHistTimeVsChargeADC_Corr_Name.Data(),"Corrected Time vs Charge",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
  list->Add(fHistTimeVsChargeADC_Corr);

  TString fHistTimeVsChargeADC_UnCorr_Name = "fHistTimeVsChargeADC_UnCorr_";
  fHistTimeVsChargeADC_UnCorr_Name += TriggerName;
  TH2F *fHistTimeVsChargeADC_UnCorr = CreateHist2D(fHistTimeVsChargeADC_UnCorr_Name.Data(),"Raw Time vs Charge",kNRawTimeBins, kRawTimeMin, kRawTimeMax, kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Leading time [ns]","ADC counts");
  list->Add(fHistTimeVsChargeADC_UnCorr);

  TString fHistWidthVsCharge_Name = "fHistWidthVsCharge_";
  fHistWidthVsCharge_Name += TriggerName;
  TH2F *fHistWidthVsCharge = CreateHist2D(fHistWidthVsCharge_Name.Data(),"Width vs Charge",kNTdcWidthBins, kTdcWidthMin, kTdcWidthMax,kNChargeChannelBins/10,kChargeChannelMin,kChargeChannelMax,"Time width [ns]","ADC counts");
  list->Add(fHistWidthVsCharge);

  TString fHistNBBflagsADAVsADC_Name = "fHistNBBflagsADAVsADC_";
  fHistNBBflagsADAVsADC_Name += TriggerName;
  TH2F *fHistNBBflagsADAVsADC = CreateHist2D(fHistNBBflagsADAVsADC_Name.Data(),"Number of BB flags",9,-0.5,8.5,9,-0.5,8.5,"Number of BB flags in ADA","Number of BB flags in ADC");
  list->Add(fHistNBBflagsADAVsADC);

  TString fHistNBGflagsADAVsADC_Name = "fHistNBGflagsADAVsADC_";
  fHistNBGflagsADAVsADC_Name += TriggerName;
  TH2F *fHistNBGflagsADAVsADC = CreateHist2D(fHistNBGflagsADAVsADC_Name.Data(),"Number of BG flags",9,-0.5,8.5,9,-0.5,8.5,"Number of BG flags in ADA","Number of BG flags in ADC");
  list->Add(fHistNBGflagsADAVsADC);

  TString fHistNBBCoincidencesADAVsADC_Name = "fHistNBBCoincidencesADAVsADC_";
  fHistNBBCoincidencesADAVsADC_Name += TriggerName;
  TH2F *fHistNBBCoincidencesADAVsADC = CreateHist2D(fHistNBBCoincidencesADAVsADC_Name.Data(),"Number of BB coincidences",5,-0.5,4.5, 5,-0.5,4.5,"Number of BB coincidences in ADA","Number of BB coincidences in ADC");
  list->Add(fHistNBBCoincidencesADAVsADC);

  TString fHistNBGCoincidencesADAVsADC_Name = "fHistNBGCoincidencesADAVsADC_";
  fHistNBGCoincidencesADAVsADC_Name += TriggerName;
  TH2F *fHistNBGCoincidencesADAVsADC = CreateHist2D(fHistNBGCoincidencesADAVsADC_Name.Data(),"Number of BG coincidences",5,-0.5,4.5, 5,-0.5,4.5,"Number of BG coincidences in ADA","Number of BG coincidences in ADC");
  list->Add(fHistNBGCoincidencesADAVsADC);

  TString fHistChargeNoFlag_Name = "fHistChargeNoFlag_";
  fHistChargeNoFlag_Name += TriggerName;
  TH1F *fHistChargeNoFlag = CreateHist1D(fHistChargeNoFlag_Name.Data(),"Charge in PM without BB flag",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"Charge","Entries");
  list->Add(fHistChargeNoFlag);

  TString fHistTimeNoFlag_Name = "fHistTimeNoFlag_";
  fHistTimeNoFlag_Name += TriggerName;
  TH2F *fHistTimeNoFlag = CreateHist2D(fHistTimeNoFlag_Name.Data(),"Time in PM without BB flag",kNChannelBins, kChannelMin, kChannelMax,kNRawTimeBins, kRawTimeMin, kRawTimeMax,"Channel","Time","Entries");
  list->Add(fHistTimeNoFlag);

  TString fHistChargeNoTime_Name = "fHistChargeNoTime_";
  fHistChargeNoTime_Name += TriggerName;
  TH2F *fHistChargeNoTime = CreateHist2D(fHistChargeNoTime_Name.Data(),"Charge in PM without time measurement",kNChannelBins, kChannelMin, kChannelMax,kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"Channel","Charge","Entries");
  list->Add(fHistChargeNoTime);

  TString fHistFlagNoTime_Name = "fHistFlagNoTime_";
  fHistFlagNoTime_Name += TriggerName;
  TH1F *fHistFlagNoTime = CreateHist1D(fHistFlagNoTime_Name.Data(),"events with BB/BG flag but no time",kNChannelBins, kChannelMin, kChannelMax,"Channel","Entries");
  list->Add(fHistFlagNoTime);


  TString fHistMeanTimeADA_Name = "fHistMeanTimeADA_";
  fHistMeanTimeADA_Name += TriggerName;
  TH1F *fHistMeanTimeADA = CreateHist1D(fHistMeanTimeADA_Name.Data(),"Mean Time in ADA",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time","Entries");
  list->Add(fHistMeanTimeADA);

  TString fHistMeanTimeADC_Name = "fHistMeanTimeADC_";
  fHistMeanTimeADC_Name += TriggerName;
  TH1F *fHistMeanTimeADC = CreateHist1D(fHistMeanTimeADC_Name.Data(),"Mean Time in ADC",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time","Entries");
  list->Add(fHistMeanTimeADC);

  TString fHistMeanTimeDifference_Name = "fHistMeanTimeDifference_";
  fHistMeanTimeDifference_Name += TriggerName;
  TH1F *fHistMeanTimeDifference = CreateHist1D(fHistMeanTimeDifference_Name.Data(),"Mean Time difference",1024,-150,150,"AD Mean time t_{A} - t_{C} [ns]","Entries");
  list->Add(fHistMeanTimeDifference);

  TString fHistMeanTimeCorrelation_Name = "fHistMeanTimeCorrelation_";
  fHistMeanTimeCorrelation_Name += TriggerName;
  TH2F *fHistMeanTimeCorrelation = CreateHist2D(fHistMeanTimeCorrelation_Name.Data(),"Mean Time in ADA-ADC",kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,kNCorrTimeBins, kCorrTimeMin, kCorrTimeMax,"Time ADA","Time ADC");
  list->Add(fHistMeanTimeCorrelation);

  TString fHistMeanTimeSumDiff_Name = "fHistMeanTimeSumDiff_";
  fHistMeanTimeSumDiff_Name += TriggerName;
  TH2F *fHistMeanTimeSumDiff = CreateHist2D(fHistMeanTimeSumDiff_Name.Data(),"Mean Time in ADA-ADC",307, -150.000000, 149.804688, 410, 0.000000, 400.390625,"AD Mean time t_{A} - t_{C} [ns]","AD Mean time t_{A} + t_{C} [ns]");
  list->Add(fHistMeanTimeSumDiff);

  TString fHistDecision_Name = "fHistDecision_";
  fHistDecision_Name += TriggerName;
  TH2F *fHistDecision = CreateHist2D(fHistDecision_Name.Data(),"Offline decision in ADA-ADC",4,0 ,4,4,0,4,"ADA","ADC");
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
  list->Add(fHistDecision);

  TString fHistTriggerMasked_Name = "fHistTriggerMasked_";
  fHistTriggerMasked_Name += TriggerName;
  TH1F *fHistTriggerMasked = CreateHist1D(fHistTriggerMasked_Name.Data(),"Trigger inputs, from FEE (BC masked)",10,0 ,10,"AD0 Trigger Type","Counts");
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
  list->Add(fHistTriggerMasked);

  TString fHistTriggerUnMasked_Name = "fHistTriggerUnMasked_";
  fHistTriggerUnMasked_Name += TriggerName;
  TH1F *fHistTriggerUnMasked = CreateHist1D(fHistTriggerUnMasked_Name.Data(),"Trigger inputs, from FEE (BC masked)",10,0 ,10,"AD0 Trigger Type","Counts");
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
  list->Add(fHistTriggerUnMasked);//31
  
  TString fHistChargeTriggerPerChannel_Name = "fHistChargeTriggerPerChannel_";
  fHistChargeTriggerPerChannel_Name += TriggerName;
  TH2F *fHistChargeTriggerPerChannel = CreateHist2D(fHistChargeTriggerPerChannel_Name.Data(),"Trigger charge per chanel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  list->Add(fHistChargeTriggerPerChannel);//32

  TString fHistChargeTriggerPerChannelPF_Name = "fHistChargeTriggerPerChannel_PF_";
  fHistChargeTriggerPerChannelPF_Name += TriggerName;
  TH2F *fHistChargeTriggerPerChannel_PF = CreateHist2D(fHistChargeTriggerPerChannelPF_Name.Data(),"Trigger charge per chanel",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"ADC counts");
  list->Add(fHistChargeTriggerPerChannel_PF);//33

  TString fHistChargeTriggerADA_Name = "fHistChargeTriggerADA_";
  fHistChargeTriggerADA_Name += TriggerName;
  TH1F *fHistChargeTriggerADA = CreateHist1D(fHistChargeTriggerADA_Name.Data(),"Trigger charge ADA",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
  list->Add(fHistChargeTriggerADA);//34

  TString fHistChargeTriggerADAPF_Name = "fHistChargeTriggerADA_PF_";
  fHistChargeTriggerADAPF_Name += TriggerName;
  TH1F *fHistChargeTriggerADA_PF = CreateHist1D(fHistChargeTriggerADAPF_Name.Data(),"Trigger charge ADA, PF protection",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
  list->Add(fHistChargeTriggerADA_PF);//35
  
  TString fHistChargeTriggerADC_Name = "fHistChargeTriggerADC_";
  fHistChargeTriggerADC_Name += TriggerName;
  TH1F *fHistChargeTriggerADC = CreateHist1D(fHistChargeTriggerADC_Name.Data(),"Trigger charge ADC",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
  list->Add(fHistChargeTriggerADC);//36

  TString fHistChargeTriggerADCPF_Name = "fHistChargeTriggerADC_PF_";
  fHistChargeTriggerADCPF_Name += TriggerName;
  TH1F *fHistChargeTriggerADC_PF = CreateHist1D(fHistChargeTriggerADCPF_Name.Data(),"Trigger charge ADC, PF protection",kNChargeChannelBins,kChargeChannelMin,kChargeChannelMax,"ADC counts");
  list->Add(fHistChargeTriggerADC_PF);//37
  
  TString fHistMaxChargeValueInt0_Name = "fHistMaxChargeValueInt0_";
  fHistMaxChargeValueInt0_Name += TriggerName;
  TH2F *fHistMaxChargeValueInt0 = CreateHist2D(fHistMaxChargeValueInt0_Name.Data(),"Maximum charge value per PM Int0",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"PM number","ADC counts");
  list->Add(fHistMaxChargeValueInt0);//38
  
  TString fHistMaxChargeValueInt1_Name = "fHistMaxChargeValueInt1_";
  fHistMaxChargeValueInt1_Name += TriggerName;
  TH2F *fHistMaxChargeValueInt1 = CreateHist2D(fHistMaxChargeValueInt1_Name.Data(),"Maximum charge value per PM Int1",kNChannelBins, kChannelMin, kChannelMax,1024,0,1024,"PM number","ADC counts");
  list->Add(fHistMaxChargeValueInt1);//39
  
  TString fHistMeanChargeADAPerFlag_Name = "fHistMeanChargeADAPerFlag_";
  fHistMeanChargeADAPerFlag_Name += TriggerName;
  TH2F *fHistMeanChargeADAPerFlag = CreateHist2D(fHistMeanChargeADAPerFlag_Name.Data(),"Mean charge ADA per flag",9,-0.5,8.5,kNChargeSideBins,kChargeSideMin,kChargeSideMax,"Number of BB coincidences","Charge","Entries");
  list->Add(fHistMeanChargeADAPerFlag);//40
  
  TString fHistMeanChargeADCPerFlag_Name = "fHistMeanChargeADCPerFlag_";
  fHistMeanChargeADCPerFlag_Name += TriggerName;
  TH2F *fHistMeanChargeADCPerFlag = CreateHist2D(fHistMeanChargeADCPerFlag_Name.Data(),"Mean charge ADC per flag",9,-0.5,8.5,kNChargeSideBins,kChargeSideMin,kChargeSideMax,"Number of BB coincidences","Charge","Entries");
  list->Add(fHistMeanChargeADCPerFlag);//41

}

//________________________________________________________________________
void AliAnalysisTaskADVVQA::SetCalibData()
{
 
    AliCDBManager *man = AliCDBManager::Instance();
    man->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB");
    man->SetRun(fRun);

    AliCDBEntry *ent = man->Get("AD/Calib/Data");
    fCalibData = (AliADCalibData*)ent->GetObject();
    
}

//________________________________________________________________________
void AliAnalysisTaskADVVQA::FillHistos(TList *list) 
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
  AliESDADfriend* esdADfriend = 0x0;
  if(fESDfriend) esdADfriend = fESDfriend->GetADfriend();
  
  fRun=fEvent->GetRunNumber();
  
  if (fRun!=fOldRun){
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
	
	((TH2F*)(list->At(2)))->Fill(i,esdAD->GetAdc(i));
	if(esdAD->GetBBFlag(i))((TH2F*)(list->At(3)))->Fill(i,esdAD->GetAdc(i));
	if(esdAD->GetBGFlag(i))((TH2F*)(list->At(4)))->Fill(i,esdAD->GetAdc(i));
	if(esdAD->GetTime(i)> -1024 + 1e-6)((TH2F*)(list->At(5)))->Fill(i,esdAD->GetAdc(i));
	
	((TH2F*)(list->At(6)))->Fill(i,esdAD->GetTime(i));
	
	if(i<8){
		((TH2F*)(list->At(13)))->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
		if(esdAD->BBTriggerADC(i)) ((TH2F*)(list->At(12)))->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
		}
        if(i>7){
		((TH2F*)(list->At(10)))->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
		if(esdAD->BBTriggerADA(i-8)) ((TH2F*)(list->At(9)))->Fill(esdAD->GetTime(i),esdAD->GetAdc(i));
		}
	((TH2F*)(list->At(8)))->Fill(i,esdAD->GetWidth(i));
	((TH2F*)(list->At(15)))->Fill(esdAD->GetWidth(i),esdAD->GetAdc(i));
	
	
	if(esdAD->GetTime(i) < -1024+1e-6) ((TH2F*)(list->At(22)))->Fill(i,esdAD->GetAdc(i));
	if(esdAD->GetTime(i) < -1024+1e-6 && (esdAD->GetBBFlag(i)||esdAD->GetBGFlag(i)))((TH1F*)(list->At(23)))->Fill(i);
	
	
	if(esdADfriend){
		if(!esdAD->GetBBFlag(i)) {
			((TH1F*)(list->At(20)))->Fill(esdAD->GetAdc(i));
			((TH2F*)(list->At(21)))->Fill(i,esdADfriend->GetTime(i));
			}
		((TH2F*)(list->At(7)))->Fill(i,esdADfriend->GetTime(i));
 		if(i<8) ((TH2F*)(list->At(14)))->Fill(esdADfriend->GetTime(i),esdAD->GetAdc(i));
 		if(i>7) ((TH2F*)(list->At(11)))->Fill(esdADfriend->GetTime(i),esdAD->GetAdc(i));
		
		Bool_t localPF = kTRUE;
		for(Int_t iClock=0; iClock<10; iClock++) if(esdADfriend->GetBBFlag(i,iClock) || esdADfriend->GetBGFlag(i,iClock)){globalPF = kFALSE; localPF = kFALSE;}
		for(Int_t iClock=11; iClock<21; iClock++) if(esdADfriend->GetBBFlag(i,iClock) || esdADfriend->GetBGFlag(i,iClock)){globalPF = kFALSE; localPF = kFALSE;}
		Int_t k = i + 16*esdADfriend->GetIntegratorFlag(i,10);
		if(esdAD->GetBBFlag(i)) fCharges[i] = esdADfriend->GetPedestal(i,10) - fCalibData->GetPedestal(k);
		((TH2F*)(list->At(32)))->Fill(i,fCharges[i]);
		if(localPF)((TH2F*)(list->At(33)))->Fill(i,fCharges[i]);
		
		Int_t chargeArray[21];
		for(Int_t iClock=0; iClock<21; iClock++) chargeArray[iClock] = esdADfriend->GetPedestal(i,iClock);
		Int_t maxClock = TMath::LocMax(21,chargeArray);
		if(!esdADfriend->GetIntegratorFlag(i,maxClock))((TH2F*)(list->At(38)))->Fill(i,chargeArray[maxClock]);
		if( esdADfriend->GetIntegratorFlag(i,maxClock))((TH2F*)(list->At(39)))->Fill(i,chargeArray[maxClock]);

		}
  }
  for(Int_t i = 0; i<8; i++)chargeADC += fCharges[i];
  for(Int_t i = 8; i<16; i++)chargeADA += fCharges[i];
  ((TH1F*)(list->At(34)))->Fill(chargeADA);
  ((TH1F*)(list->At(36)))->Fill(chargeADC);
  if(globalPF){
  	  ((TH1F*)(list->At(35)))->Fill(chargeADA);
  	  ((TH1F*)(list->At(37)))->Fill(chargeADC);
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
	
  ((TH1F*)(list->At(0)))->Fill(totChargeADA);
  ((TH1F*)(list->At(1)))->Fill(totChargeADC);
  
  ((TH1F*)(list->At(40)))->Fill(nBBflagsADA,totChargeADA);
  ((TH1F*)(list->At(41)))->Fill(nBBflagsADC,totChargeADC);
  
  ((TH2F*)(list->At(16)))->Fill(nBBflagsADA,nBBflagsADC);
  ((TH2F*)(list->At(17)))->Fill(nBGflagsADA,nBGflagsADC);
  ((TH2F*)(list->At(18)))->Fill(nBBCoincidencesADA,nBBCoincidencesADC);
  ((TH2F*)(list->At(19)))->Fill(nBGCoincidencesADA,nBGCoincidencesADC);
  ((TH1F*)(list->At(24)))->Fill(esdAD->GetADATime());
  ((TH1F*)(list->At(25)))->Fill(esdAD->GetADCTime());
  if(esdAD->GetADATime()!= -1024 && esdAD->GetADCTime()!= -1024){
  	((TH1F*)(list->At(26)))->Fill(esdAD->GetADATime()-esdAD->GetADCTime());
	((TH2F*)(list->At(27)))->Fill(esdAD->GetADATime(),esdAD->GetADCTime());
	((TH2F*)(list->At(28)))->Fill(esdAD->GetADATime()-esdAD->GetADCTime(),esdAD->GetADATime()+esdAD->GetADCTime());
	}
  ((TH2F*)(list->At(29)))->Fill(esdAD->GetADADecision(),esdAD->GetADCDecision());
    
  //Triggers
  UShort_t fTriggerBC = esdAD->GetTriggerBits();
  for(Int_t i = 0; i<6; i++) if(fTriggerBC & (1 << i) ? kTRUE : kFALSE) ((TH1F*)(list->At(30)))->Fill(i);
  for(Int_t i = 12; i<16; i++) if(fTriggerBC & (1 << i) ? kTRUE : kFALSE) ((TH1F*)(list->At(30)))->Fill(i-6);
  
  if(esdADfriend){
  	UShort_t fTriggerUnBC = esdADfriend->GetTriggerInputs();
  	for(Int_t i = 0; i<6; i++) if(fTriggerUnBC & (1 << i) ? kTRUE : kFALSE) ((TH1F*)(list->At(31)))->Fill(i);
	for(Int_t i = 12; i<16; i++) if(fTriggerUnBC & (1 << i) ? kTRUE : kFALSE) ((TH1F*)(list->At(31)))->Fill(i-6);	
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskADVVQA::UserExec(Option_t *) 
{
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
  
  if(fAnalysisUtils->IsPileUpEvent(fESD)) return;
  if(fAnalysisUtils->IsSPDClusterVsTrackletBG(fESD)) return;
  
  if(!(fESD->GetHeader()->IsTriggerInputFired("0VBA")) && !(fESD->GetHeader()->IsTriggerInputFired("0VBC"))) FillHistos(fList_NoVBA_NoVBC);
  if( (fESD->GetHeader()->IsTriggerInputFired("0VBA")) && !(fESD->GetHeader()->IsTriggerInputFired("0VBC"))) FillHistos(fList_VBA_NoVBC);
  if(!(fESD->GetHeader()->IsTriggerInputFired("0VBA")) &&  (fESD->GetHeader()->IsTriggerInputFired("0VBC"))) FillHistos(fList_NoVBA_VBC);
  if( (fESD->GetHeader()->IsTriggerInputFired("0VBA")) &&  (fESD->GetHeader()->IsTriggerInputFired("0VBC"))) FillHistos(fList_VBA_VBC);
  if(  fESD->GetHeader()->IsTriggerInputFired("0TVX")) FillHistos(fList_TVX);
  if( (fESD->GetHeader()->IsTriggerInputFired("0UBA")) &&  (fESD->GetHeader()->IsTriggerInputFired("0UBC"))) FillHistos(fList_UBA_UBC);  

  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskADVVQA::Terminate(Option_t *) 
{
 // Draw result to the screen
  // Called once at the end of the query

}
//________________________________________________________________________
TH1F * AliAnalysisTaskADVVQA::CreateHist1D(const char* name, const char* title,Int_t nBins, 
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
TH2F * AliAnalysisTaskADVVQA::CreateHist2D(const char* name, const char* title,Int_t nBinsX, 
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
TH3F * AliAnalysisTaskADVVQA::CreateHist3D(const char* name, const char* title,Int_t nBinsX, 
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

