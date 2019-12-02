#if defined(__CLING__)

  R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
  #include <PWGHF/hfe/macros/configs/pp/ConfigHaHFECorrel.C>


#elif defined(__CINT__)
  // ROOT5-specific code here ...
#endif


  AliAnalysisTaskHaHFECorrel *AddTaskHaHFECorrel(Double_t period, Int_t MinNTr, Int_t MaxNTr, Bool_t TRDQA, Bool_t TagEff, Bool_t RecEff, Bool_t OneTimeCheck, Bool_t CorrHadron, Bool_t CorrLP, Bool_t MCTruth,  Bool_t IsMC, Bool_t IsAOD, Bool_t IsHFE, Bool_t UseTender, Bool_t UseEventWeights, Double_t EtaMax, Int_t ITSnCut, Float_t ITSSharedCluster,  Int_t TPCnCut, Int_t TPCnCutdEdx,   Double_t PhotElecPtCut, Int_t PhotElecTPCnCut,Bool_t PhotElecITSrefitCut, Int_t PhotCorrCase, Double_t InvmassCut, Int_t HTPCnCut,   Bool_t HITSrefitCut, Bool_t HTPCrefitCut, Bool_t UseITSsa, Double_t SigmaITScut, Double_t SigmaTOFcut, Double_t SigmaTPCcut, Int_t VarOptE, Int_t VarOptH, Int_t VarOptPhot,  const char * ID="")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHaHFECorrel", "No analysis manager found.");
    return 0;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHaHFECorrel", "This task requires an input event handler");
    return 0x0;
  }

 

  /*
  AliMCEventHandler* mcHand = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHand);
  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if (!mcH) {
    MCthere=kFALSE;
  }
  */

  // VarOption for variations affecting online efficiencies, 0 is default configuration
  




  //TString VarString = Form("_%i", VarOptE);
  //ID+=VarString.Data();

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/ConfigHaHFECorrel.C");
  AliAnalysisTaskHaHFECorrel *taskMB = 
    ConfigHaHFECorrel(period, MinNTr, MaxNTr, TRDQA, TagEff, RecEff, OneTimeCheck,  CorrHadron, CorrLP, MCTruth, IsMC, IsAOD, IsHFE, UseTender, UseEventWeights, EtaMax, ITSnCut, ITSSharedCluster, TPCnCut, TPCnCutdEdx, PhotElecPtCut,PhotElecTPCnCut, PhotElecITSrefitCut, PhotCorrCase, InvmassCut,  HTPCnCut,  HITSrefitCut, HTPCrefitCut, UseITSsa, SigmaITScut, SigmaTOFcut, SigmaTPCcut, VarOptE, VarOptH, VarOptPhot, ID);
  if (!taskMB) {
    Error("AddTaskHaHFECorrel", "No task found.");
  }
  //taskMB->SelectCollisionCandidates(AliVEvent::kINT7);

  
  //Load correction weights for pi0, eta
    if (IsMC) {
      TH1::AddDirectory(kFALSE);
      printf("Loading Pi0EtaCorrectionFiles\n");
      TString CorrectPi0EtaFile;
      if (!IsHFE) CorrectPi0EtaFile="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/Pi0EtaWeightsMinBias.root";
      else if (IsHFE) CorrectPi0EtaFile="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/Pi0EtaWeightsHFE.root";
      TFile *CorrectPi0Eta = TFile::Open(CorrectPi0EtaFile.Data());
      CorrectPi0Eta->ls();
      if (!CorrectPi0Eta->IsZombie()) {   
	TH1F * Pi0W = 0;
	Pi0W =(TH1F*)CorrectPi0Eta->Get("Pi0WeightsIncl");
	TH1F * EtaW=0;
	EtaW = (TH1F*)CorrectPi0Eta->Get("EtaWeightsIncl");
	if (Pi0W) taskMB->SetPi0WeightToData(*Pi0W);
        else printf("Could not load Pi0Weights\n");
	if (EtaW)  taskMB->SetEtaWeightToData(*EtaW);
        else printf("Could not load EtaWeights\n");
      }
      else printf("Could not open Pi0Eta correction file \n");
      TH1::AddDirectory(kTRUE);
    }


    
   // Load correction weights for MC background
   if (IsMC) {
     TH1::AddDirectory(kFALSE);
     printf("Loading MCBGCorrectionFiles\n");
     TString CorrectBGFileName ="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/BGWeights.root";
     TFile *CorrectBGFile = TFile::Open(CorrectBGFileName.Data());
     if (!CorrectBGFile->IsZombie()) {    
       TH2F * BGHist =0;
       BGHist=(TH2F*)CorrectBGFile->Get("BGWeightHist");
       if (BGHist) taskMB->SetBGWeight(*BGHist);
       else printf("Could not load BGWeights\n");
     }
     else printf("Could not open BG correction file \n");
     TH1::AddDirectory(kTRUE);
   }
  




   TH1::AddDirectory(kFALSE);
   printf("Loading SPDnTr files\n");
   TString SPDnTrFileName;
   if (IsMC) SPDnTrFileName="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/SPDProfile_MC.root"; //SPDnTrAvg_MC.root";
   else SPDnTrFileName="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/SPDProfile_Data.root"; //SPDnTrAvg_Data.root";
   TFile *SPDnTrFile  = TFile::Open(SPDnTrFileName.Data());
   if (SPDnTrFile) {    
     TH3F* SPDConfigProfiles= (TH3F*)SPDnTrFile->Get("SPDConfigs_Hist");
     TH1I* SPDConfigHist = (TH1I*) SPDnTrFile->Get("SPDConfigHist");
     if (SPDConfigHist) taskMB->SetSPDConfigHist(*SPDConfigHist);
     else printf("Could not load SPDConfigHist\n");
     if (SPDConfigProfiles) taskMB->SetSPDConfigProfiles(*SPDConfigProfiles);
     else printf("Could not load SPDConfigProfiles\n");
   }
   else printf("Could not open SPDnTrAvg correction file \n");
   TH1::AddDirectory(kTRUE);

   TH1::AddDirectory(kFALSE);
   printf("Loading RecEffFiles: %i for hadrons, %i for electrons \n", VarOptH, VarOptE);
   TString RecEffFileName="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/RecEff_Periods.root";
   TFile *RecEffFile = TFile::Open(RecEffFileName.Data());
   //  RecEffFile->ls();
   //RecEffFile->ls();
   if (RecEffFile) {    
     TH3F * HadRecEff = (TH3F*)RecEffFile->Get(Form("HadRecEff_%i", VarOptH));
     TH3F * EleRecEff = (TH3F*)RecEffFile->Get(Form("EleRecEff_%i", VarOptE));
     if (HadRecEff) taskMB->SetHadRecEff(*HadRecEff);
     else {
       HadRecEff = (TH3F*)RecEffFile->Get(Form("HadRecEff_0"));
       if (HadRecEff) {
   	taskMB->SetHadRecEff(*HadRecEff);
   	printf("WARNING: Could not load HadRecEff - using default values\n");
       }
       else {
   	printf("WARNING: Could not load HadRecEff\n");
   	return NULL;
       }
     }
     if (EleRecEff) taskMB->SetEleRecEff(*EleRecEff);
     else {
       EleRecEff = (TH3F*)RecEffFile->Get(Form("EleRecEff_0"));
       if (EleRecEff) {
   	taskMB->SetEleRecEff(*EleRecEff);
   	printf("WARNING: Could not load EleRecEff - using default values\n");
       }
       else {
   	printf("WARNING: Could not load EleRecEff\n");
   	return NULL;
       }
     }
   }
   else printf("Could not open RecEff correction file \n");
   TH1::AddDirectory(kTRUE);


   TH1::AddDirectory(kFALSE);
   printf("Loading NonTagCorr\n");
   TString NonTagFileName="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/NonTagCorrHist.root";
   TFile *NonTagFile = TFile::Open(NonTagFileName.Data());
   if (!NonTagFile->IsZombie()) {
     NonTagFile->ls();
     TH1F * NonTagCorr = (TH1F*)NonTagFile->Get("NonTagCorr");
     if (NonTagCorr) taskMB->SetNonTagCorr(*NonTagCorr);
     else printf("Could not load NonTagCorr\n");
   }
   else printf("Could not open NonTag correction file \n");
   TH1::AddDirectory(kTRUE);

  
   TH1::AddDirectory(kFALSE);
   printf("Loading EventWeightFile\n");
   TString EventWeightFileName="alien:///alice/cern.ch/user/f/flherrma/HaHFECorrel/TrigVtxEff_Periods.root";
   TFile *EventWeightFile = TFile::Open(EventWeightFileName.Data());
   EventWeightFile->ls();
   if (EventWeightFile) {    
     TH3F * TriggerWeight  = (TH3F*)EventWeightFile->Get("TrigEff");
     if (TriggerWeight) taskMB->SetTriggerWeight(*TriggerWeight); 
     TH2F * VtxWeight;
     VtxWeight = (TH2F*)EventWeightFile->Get("VtxEff");
  //   //    if (!IsMC || !VtxWeight) VtxWeight = (TH1F*)EventWeightFile->Get("VtxWeight");
     if (VtxWeight) taskMB->SetVtxWeight(*VtxWeight);
     if (!VtxWeight || !TriggerWeight) printf("Could no open EventWeight hists");
   }
   else  printf("Could not open EventWeight file \n"); 
   TH1::AddDirectory(kTRUE);

   printf("Setting VarOptions: Electron - %i, Hadron - %i, Phot - %i", VarOptE, VarOptH, VarOptPhot);
   taskMB->SetEleVarOpt(VarOptE);
   taskMB->SetHadVarOpt(VarOptH);
   taskMB->SetPhotVarOpt(VarOptPhot);
 






  mgr->AddTask(taskMB);

  TString containerName1 = mgr->GetCommonFileName();
  containerName1 += ":PWGHF_HaHFECorrel_kINT7_";
  containerName1 += ID;
        
  TString name1 = "histMB_";
  name1 += ID;
        
  TString name2 = "Main_";
  name2 += ID;
  TString name3 = "LP_";
  name3 += ID;
  TString name4 = "Hadron_";
  name4 += ID;
  TString name5 = "QA_";
  name5 += ID;


  AliAnalysisDataContainer *cinput      = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1    = mgr->CreateContainer(name1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  AliAnalysisDataContainer *coutputMain = mgr->CreateContainer(name2.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  AliAnalysisDataContainer *coutputLP   = mgr->CreateContainer(name3.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  AliAnalysisDataContainer *coutputHa   = mgr->CreateContainer(name4.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  AliAnalysisDataContainer *coutputQA   = mgr->CreateContainer(name5.Data(), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());

  mgr->ConnectInput(taskMB, 0, cinput);
  mgr->ConnectOutput(taskMB, 1, coutput1);
  mgr->ConnectOutput(taskMB, 2, coutputMain);
  mgr->ConnectOutput(taskMB, 3, coutputLP);
  mgr->ConnectOutput(taskMB, 4, coutputHa);
  mgr->ConnectOutput(taskMB, 5, coutputQA);

  printf("return Task");


  return taskMB;
}
