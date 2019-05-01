#include<TList.h>
#include "TSystem.h"
class AliAnalysisTaskCMEV0;

 void AddTaskCMEV0(Int_t gFilterBit = 768, Int_t gClusterTPC = 70, Double_t fpTLow = 0.2, Double_t fpTHigh = 10.0, 
 Double_t fEtaLow = -0.8, Double_t fEtaHigh = 0.8, TString sAnalysisFile = "AOD", TString sDataSet = "2015", TString sAnalysisType = "AUTOMATIC", 
 TString sEventTrigger = "MB", Bool_t bEventCutsQA = kFALSE, Bool_t bTrackCutsQA = kFALSE,Double_t dVertexLow = -10.,Double_t dVertexHigh = 10., 
 Bool_t bPileUp = kFALSE, Bool_t bPileUpTight = kFALSE, Float_t fPileUpSlope = 3.43, Float_t fPileUpConst = 43.0, TString sCentEstimator = "V0",
 Bool_t bFBeffi = kFALSE,TString sEfficiencyFB = "alien:///alice/cern.ch/user/m/mhaque/gain/FB96_Hijing_LHC15o_HI_CorSec.root",
 TString sFBEffiDimension = "1D", 
 Bool_t bApplyNUA = kFALSE, TString sNUAFile="alien:///alice/cern.ch/user/m/mhaque/gain/Run2015o_Pass1_FB768_pT0p2_5GeV_NUA_Wgt_PosNeg_Run.root", 
 Bool_t bZDCGainEq= kFALSE, TString sZDCFile="alien:///alice/cern.ch/user/m/mhaque/gain/Run2015o_pass1_ZDNP_WgtTotEn_VsCentRun.root", 
 Bool_t bV0MgainCorr= kFALSE, TString sV0MFile="alien:///alice/cern.ch/user/m/mhaque/gain/Run2015_V0GainEq_RbyR_pPb_FAST.root", 
 Bool_t bFillTPCQn= kFALSE, Bool_t bFillNUAhist= kFALSE,Bool_t bFillZDCHist= kFALSE, Bool_t bSkipNestedLoop= kFALSE, 
 Int_t fSetHarmN = 1, Int_t fSetHarmM = 1, Int_t fSetPsiHarm = 2, Bool_t bUseNUAinEP = kFALSE, TString sNUAtype="NewR", Float_t hbtCut = 0.0,
 const char *suffix = "")
{

  //gSystem->Load("libPWGflowBase.so");
  //gSystem->Load("libPWGflowTasks.so");

  gSystem->AddIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -g ");
  gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/OCDB -I$ALICE_ROOT/STEER/macros -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/TRD -I$ALICE_ROOT/ZDC -I$ALICE_ROOT/macros -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/OADB $ALICE_PHYSICS/OADB/macros -I$ALICE_PHYSICS/PWGGA -I$ALICE_PHYSICS/PWGCF -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/TENDER -I$ALICE_PHYSICS/TENDER/Tender -I$ALICE_PHYSICS/TENDER/TenderSupplies -I$ALICE_PHYSICS/PARfiles -I$ALICE_PHYSICS/PWGCF/FLOW/macros I$ALICE_PHYSICS/PWGPP/ZDC -g ");

  //Fixed Track cuts: only vary for systematic check
  Float_t dDCAxy = 2.4;
  Float_t dDCAz  = 3.2;

  Float_t dcentrMin= 0.;
  Float_t dcentrMax=90.;


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();


  TString taskFEname;
  taskFEname.Form("ZDCFlowEventTask%s", suffix);

  AliAnalysisTaskFlowEvent *taskFE = new AliAnalysisTaskFlowEvent(taskFEname,"",bEventCutsQA);

  taskFE->SetQAOn(bEventCutsQA);
  taskFE->SetAnalysisType(sAnalysisType); //sanalysisType = AUTOMATIC see the initializers!!

  if(sDataSet=="2015"||sDataSet=="2015LI"||sDataSet=="2015pPb"||sDataSet=="pPb"){
    taskFE->SelectCollisionCandidates(AliVEvent::kINT7);
  }
  else{
    taskFE->SelectCollisionCandidates(AliVEvent::kMB);
  }

 
  taskFEname.Form("EventCuts%s", suffix);

  AliFlowEventCuts *cutsEvent = new AliFlowEventCuts(taskFEname);
  cutsEvent->SetCheckPileup(kFALSE);
  cutsEvent->SetPrimaryVertexZrange(dVertexLow, dVertexHigh);      // vertex-z cut
  cutsEvent->SetQA(bEventCutsQA);                                  // enable the qa plots

  if(sDataSet=="2015pPb"||sDataSet=="pPb"){
    cutsEvent->SetCutTPCmultiplicityOutliersAOD(kFALSE); 	   // multiplicity outlier cut
  }
  else {
    cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE); 	 
  }

  if(sDataSet=="2015"||sDataSet=="2015LI"||sDataSet=="2015pPb")
  {
   cutsEvent->SetCentralityPercentileRange(dcentrMin, dcentrMax, kTRUE);
  }
  else
  {
   cutsEvent->SetCentralityPercentileRange(dcentrMin, dcentrMax);
   cutsEvent->SetCheckPileup(kTRUE);
  }

  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0); //default is V0

  if(sCentEstimator=="TPC"){
    cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kTPConly);
  }

  if(sDataSet == "2011"){
    cutsEvent->SetLHC11h(kTRUE);
  }
  else if(sDataSet == "2010"){
    cutsEvent->SetLHC10h(kTRUE);
  }


  taskFEname.Form("RefMultCuts%s", suffix);

  AliFlowTrackCuts* RefMultCuts = new AliFlowTrackCuts(taskFEname);
  RefMultCuts->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  RefMultCuts->SetAODfilterBit(gFilterBit);
  RefMultCuts->SetMinimalTPCdedx(-99999);
  RefMultCuts->SetMaxDCAToVertexXY(dDCAxy);
  RefMultCuts->SetMaxDCAToVertexZ(dDCAz);
  RefMultCuts->SetMinNClustersTPC(gClusterTPC);
  RefMultCuts->SetMinChi2PerClusterTPC(0.1);
  RefMultCuts->SetMaxChi2PerClusterTPC(4.);
  RefMultCuts->SetPtRange(fpTLow,fpTHigh);
  RefMultCuts->SetEtaRange(fEtaLow,fEtaHigh);
  RefMultCuts->SetAcceptKinkDaughters(kFALSE);
  


  cutsEvent->SetRefMultCuts(RefMultCuts);
  cutsEvent->SetRefMultMethod(AliFlowEventCuts::kTPConly);


  taskFE->SetCutsEvent(cutsEvent);	//pass these cuts to your flow event task

  taskFEname.Form("RP_cuts_%s", suffix);

  AliFlowTrackCuts *cutsRP  =  new AliFlowTrackCuts(taskFEname);
  cutsRP->SetParamType(AliFlowTrackCuts::kAODFilterBit);  //sets how we want to select the tracks.
  cutsRP->SetAODfilterBit(gFilterBit);                    
  cutsRP->SetPhiMin(0.);
  cutsRP->SetPhiMax(TMath::TwoPi());
  cutsRP->SetDivSigma(kTRUE);
  cutsRP->SetQA(bTrackCutsQA);
  cutsRP->SetMinimalTPCdedx(-99999);
  cutsRP->SetMinNClustersTPC(gClusterTPC);
  cutsRP->SetPtRange(fpTLow,fpTHigh);
  cutsRP->SetEtaGap(-2.0,2.0);  
  cutsRP->SetEtaRange(-5.0,5.0);
  cutsRP->SetRequireCharge(kTRUE);
  cutsRP->SetMinChi2PerClusterTPC(0.1);
  cutsRP->SetMaxChi2PerClusterTPC(4.0);
  cutsRP->SetAcceptKinkDaughters(kFALSE);

  taskFEname.Form("POI_cuts_%s", suffix);
  AliFlowTrackCuts *cutsPOI   = new AliFlowTrackCuts(taskFEname);
  cutsPOI->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  cutsPOI->SetAODfilterBit(gFilterBit);
  cutsPOI->SetMinimalTPCdedx(-99999);
  cutsPOI->SetMinNClustersTPC(gClusterTPC);
  cutsPOI->SetPhiMin(0.);
  cutsPOI->SetPhiMax(TMath::TwoPi());
  cutsPOI->SetPtRange(fpTLow,fpTHigh);
  cutsPOI->SetEtaRange(fEtaLow,fEtaHigh);
  cutsPOI->SetRequireCharge(kTRUE);
  cutsPOI->SetQA(bTrackCutsQA);
  cutsPOI->SetMinChi2PerClusterTPC(0.1);
  cutsPOI->SetMaxChi2PerClusterTPC(4.0);
  cutsPOI->SetAcceptKinkDaughters(kFALSE);

  taskFE->SetCutsRP(cutsRP);
  taskFE->SetCutsPOI(cutsPOI);
  taskFE->SetSubeventEtaRange(-5.5,-4.0, 4.0,5.5);


  // Finally add the task to the manager
  mgr->AddTask(taskFE);


  TString file = AliAnalysisManager::GetCommonFileName();

  TString ContainerFE;
  ContainerFE.Form("FlowEventCont_%s", suffix);


  AliAnalysisDataContainer *coutputFE = mgr->CreateContainer(ContainerFE,AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer    *cinput = mgr->GetCommonInputContainer();  //AOD event

  mgr->ConnectInput( taskFE, 0, cinput); 	//connect the input data (AOD) to the flow event task
  mgr->ConnectOutput(taskFE, 1, coutputFE); 	//get the output of taskFE to a exchange container.





  if(bEventCutsQA || bTrackCutsQA){
    TString taskFEQA = file;      // file is the common outfile filename
    taskFEQA   += ":QAcharge";
    //taskFEQA   += suffix;    

    TString ContainerFEQA;
    ContainerFEQA.Form("FEContQA_%s", suffix);

    AliAnalysisDataContainer *coutputFEQA = mgr->CreateContainer(ContainerFEQA,TList::Class(),AliAnalysisManager::kOutputContainer,taskFEQA.Data());
    mgr->ConnectOutput(taskFE, 2, coutputFEQA);          // kOutputContainer: written to the output file
  }


  //==========================================================================================

  TString TaskZDCflow;
  TaskZDCflow.Form("TaskZDCFlow_%s", suffix);

  AliAnalysisTaskCMEV0 *taskQC_prot = new AliAnalysisTaskCMEV0(TaskZDCflow);

  if(sDataSet == "2015"||sDataSet == "2015LI"||sDataSet == "2015pPb")
  {
    taskQC_prot->SelectCollisionCandidates(AliVEvent::kINT7);
  }
  else{
    taskQC_prot->SelectCollisionCandidates(AliVEvent::kMB);
  }

  taskQC_prot->SetHarmonicN(fSetHarmN);  
  taskQC_prot->SetHarmonicM(fSetHarmM);  
  taskQC_prot->SetPsiHarmonic(fSetPsiHarm);
  taskQC_prot->SetCentEstimator(sCentEstimator);
  taskQC_prot->SetRejectPileUp(bPileUp);  
  taskQC_prot->SetRejectPileUpTight(bPileUpTight); //kTRUE:700,kFALSE:15000
  taskQC_prot->SetDataSet(sDataSet);   
//taskQC_prot->SetAnalysisSet(sAnalysisDef); 
  taskQC_prot->SetApplyNUACorr(bApplyNUA);
  taskQC_prot->SetStoreTPCQnAvg(bFillTPCQn);
  taskQC_prot->SetFillNUAHist(bFillNUAhist);
  taskQC_prot->SetApplyZDCCorr(bZDCGainEq);
  taskQC_prot->SetApplyNUAinEP(bUseNUAinEP);
  taskQC_prot->SetSourceFileNUA(sNUAtype);
  //taskQC_prot->SetDebugLevel(0);
  taskQC_prot->SetFillZDNCalHist(bFillZDCHist);
  taskQC_prot->SetSkipNestedLoop(bSkipNestedLoop);
  taskQC_prot->SetMCEffiDimension(sFBEffiDimension);
  taskQC_prot->SetRemoveNegTrkRndm(kFALSE);
  taskQC_prot->SetApplyV0MCorr(bV0MgainCorr);
  taskQC_prot->SetPileUpCutParam(fPileUpSlope,fPileUpConst);
  taskQC_prot->SetTrackFilterBit(gFilterBit);
  taskQC_prot->SetHBTcutParameter(hbtCut);



  
  if(bFBeffi){
    TFile* FileFBeffi  = TFile::Open(sEfficiencyFB,"READ");
    TList* FBEffiListUse = dynamic_cast<TList*>(FileFBeffi->FindObjectAny("fMcEffiHij"));

    if(FBEffiListUse) {
      //std::cout<<"\n\n Efficiency Histograms found\n:"<<FBEffiListUse->ls()<<"\n\n"<<std::endl;
       taskQC_prot->SetFBEfficiencyList(FBEffiListUse);
    }
    else{
       std::cout<<"\n\n !!!!**** ERROR: FB Efficiency Histograms not found  *****\n\n"<<std::endl;
      //exit(1);
    }   
  }
  if(bApplyNUA && sNUAtype!="OldJ") {
    TFile* fNUAFile = TFile::Open(sNUAFile,"READ");
    if(!fNUAFile) {
      printf("\n\n *** ERROR: NUA wgt file not found! **EXIT** \n\n");
      //exit(1);
    }
    else{ 
      TList* fListNUA = dynamic_cast<TList*>(fNUAFile->FindObjectAny("fNUA_ChPosChNeg"));
      if(fListNUA){
        taskQC_prot->SetInputListNUA(fListNUA);
      }
      else{
        printf("\n\n *** ERROR: NUA file Exist, But fList name is wrong !! **EXIT** \n\n");
      }
    }
  }

  if(bZDCGainEq){
     TFile* fZDCGainFile = TFile::Open(sZDCFile,"READ");
     if(!fZDCGainFile) {
       printf("\n\n *** ERROR: ZDC wgt file not found! **EXIT** \n\n");
       exit(1);
     }

     TList* fZDCWgtUse = dynamic_cast<TList*>(fZDCGainFile->FindObjectAny("fZDN_ZDP_Wgts"));
  
     if(fZDCWgtUse) {
       taskQC_prot->SetGainCorrZDNP(fZDCWgtUse);
     }
     else{
       printf("\n\n !!!!**** ERROR: ZDC Channel wgt List not found **EXIT**!!!\n\n");
       exit(1);
     }
  }

  TList* fListV0MUse=NULL;

  if(bV0MgainCorr){

    TFile* fV0MFile = TFile::Open(sV0MFile,"READ");

     if(!fV0MFile) {
        printf("\n\n *** ERROR: VOM Gain correction file not found! **EXIT** \n\n");
        exit(1);
     } 
     else{
       //TList* fListV0MUse = dynamic_cast<TList*>(fV0MFile->FindObjectAny("fV0MChWgts"));
        fListV0MUse = dynamic_cast<TList*>(fV0MFile->FindObjectAny("fV0MChWgts"));
     }

     if(fListV0MUse) {
        taskQC_prot->SetInputListforV0M(fListV0MUse);
     }
     else{
        printf("\n\n !!!!**** ERROR: VOM Gain List not found **EXIT**!!!\n\n");
        taskQC_prot->SetInputListforV0M(NULL);
       //exit(1);
     }
  }





  mgr->AddTask(taskQC_prot);                      // connect the task to the analysis manager
  mgr->ConnectInput(taskQC_prot, 0, cinput);      // give AOD event to my Task..!!
  mgr->ConnectInput(taskQC_prot, 1, coutputFE);   // give FlowEvent object to my Task..!!

 //mgr->ConnectInput(taskQC_prot, 2,(AliAnalysisDataContainer*) mgr->GetContainers()->FindObject("ZDCEPExchangeContainer"));




  TString outputSP = file;      // file is the common outfile filename
  outputSP += ":ZDCgains";
  //outputSP   += suffix;    should I do this or not?

  TString fZDCCont1;
  fZDCCont1.Form("ResultsZDCCME_%s", suffix);
  AliAnalysisDataContainer *coutputSP1 = mgr->CreateContainer(fZDCCont1,TList::Class(),AliAnalysisManager::kOutputContainer,outputSP.Data());
  mgr->ConnectOutput(taskQC_prot, 1, coutputSP1);

  TString fZDCCont2;
  fZDCCont2.Form("QAandCalibHist_%s", suffix);
  AliAnalysisDataContainer *coutputSP2 = mgr->CreateContainer(fZDCCont2,TList::Class(),AliAnalysisManager::kOutputContainer,outputSP.Data());
  mgr->ConnectOutput(taskQC_prot, 2, coutputSP2);


  TString fZDCCont3;
  fZDCCont3.Form("NUACorrectionHist_%s", suffix);
  AliAnalysisDataContainer *coutputSP3 = mgr->CreateContainer(fZDCCont3,TList::Class(),AliAnalysisManager::kOutputContainer,outputSP.Data());
  mgr->ConnectOutput(taskQC_prot, 3, coutputSP3);

  //if(!mgr->InitAnalysis())  // check if we can initialize the manager
  //{
  //  return;
  //}


  //return taskQC_prot;

  printf("\n\n=======================  Info: AliAnalysisTaskCMEV0 Configured properly =============================\n\n");

}//main ends
