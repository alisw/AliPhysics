void runAODProof(Int_t c=2, const char * proofMode = "full")
{ //1 data AOD049
  //2 MC AOD048
  //3 data AOD086
  //4 MC AOD090
  
  //  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libOADB.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  AliAnalysisAlien * handler = new AliAnalysisAlien("test");
  handler->SetOverwriteMode();
  handler->SetRunMode(proofMode);
  handler->SetProofReset(0);
  //handler->SetROOTVersion("v5-33-02a");
  //handler->SetAliROOTVersion("v5-03-11-AN");
  handler->SetAliROOTVersion("v5-04-25-AN");
   
  //handler->SetNproofWorkers(1);
  //handler->SetNproofWorkersPerSlave(4);
  handler->SetProofCluster(Form("%s@alice-caf.cern.ch", gSystem->Getenv("CAFUSER")));
  //handler->SetProofCluster(Form("%s@skaf.saske.sk",gSystem->Getenv("CAFUSER")));
  // Set handler for Real DATA:
  if (c == 1){
    //handler->SetProofDataSet("/alice/data/LHC10h_000138653_p2_AOD049#aodTree");
    //handler->SetProofDataSet("/alice/data/LHC10h_000138653_p2_AOD049#aodTree|/alice/data/LHC10h_000138662_p2_AOD049#aodTree|/alice/data/LHC10h_000138666_p2_AOD049#aodTree|/alice/data/LHC10h_000138795_p2_AOD049#aodTree|/alice/data/LHC10h_000139107_p2_AOD049#aodTree|/alice/data/LHC10h_000139110_p2_AOD049#aodTree");
    handler->SetProofDataSet("/default/lmilano/LHC10h_000138653_AOD049_p2#aodTree|/default/lmilano/LHC10h_000138662_AOD049_p2#aodTree|/default/lmilano/LHC10h_000138666_AOD049_p2#aodTree|/default/lmilano/LHC10h_000138730_AOD049_p2#aodTree|/default/lmilano/LHC10h_000138732_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139507_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139465_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139437_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139107_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139510_AOD049_p2#aodTree");
  }
  if (c == 2){
    //handler->SetProofDataSet("/alice/sim/LHC11a10a_000138653_AOD048#aodTree|/alice/sim/LHC11a10a_000138662_AOD048#aodTree|/alice/sim/LHC11a10a_000138666_AOD048#aodTree|/alice/sim/LHC11a10a_000138795_AOD048#aodTree|/alice/sim/LHC11a10a_000139107_AOD048#aodTree|/alice/sim/LHC11a10a_000139110_AOD048#aodTree");      
  handler->SetProofDataSet("/default/lmilano/LHC11a10a_138653_AOD048#aodTree|/default/lmilano/LHC11a10a_138662_AOD048#aodTree|/default/lmilano/LHC11a10a_138666_AOD048#aodTree|/default/lmilano/LHC11a10a_138730_AOD048#aodTree|/default/lmilano/LHC11a10a_138732_AOD048#aodTree|/default/lmilano/LHC11a10a_139507_AOD048#aodTree|/default/lmilano/LHC11a10a_139465_AOD048#aodTree|/default/lmilano/LHC11a10a_139437_AOD048#aodTree|/default/lmilano/LHC11a10a_139107_AOD048#aodTree|/default/lmilano/LHC11a10a_139510_AOD048#aodTree");      
  }
  if (c == 3){
    handler->SetProofDataSet("/default/lmilano/LHC10h_000138653_AOD086_p2#aodTree|/default/lmilano/LHC10h_000138662_AOD086_p2#aodTree|/default/lmilano/LHC10h_000138666_AOD086_p2#aodTree|/default/lmilano/LHC10h_000139107_AOD086_p2#aodTree|/default/lmilano/LHC10h_000138275_AOD086_p2#aodTree");      
    //handler->SetProofDataSet("/default/lmilano/LHC10h_000138275_AOD086_p2#aodTree");      
  }
  if (c == 4){
    handler->SetProofDataSet("/default/lmilano/LHC11a10a_bis_138653_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_138662_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_138666_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_139107_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_138653_AOD090#aodTree");      
    //handler->SetProofDataSet("/default/lmilano/LHC11a10a_bis_138653_AOD090#aodTree");      
  }
   
  gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
  gROOT->LoadMacro("AliSpectraAODPID.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");

  handler->SetAliRootMode("default");
  handler->SetAdditionalLibs("AliSpectraAODHistoManager.cxx AliSpectraAODHistoManager.h AliSpectraAODPID.cxx AliSpectraAODPID.h AliSpectraAODTrackCuts.cxx AliSpectraAODTrackCuts.h AliSpectraAODEventCuts.cxx AliSpectraAODEventCuts.h AliAnalysisTaskSpectraAOD.cxx AliAnalysisTaskSpectraAOD.h");
  handler->SetAnalysisSource("Histograms.h HistogramNames.h AliSpectraAODHistoManager.cxx+ AliSpectraAODTrackCuts.cxx+ AliSpectraAODEventCuts.cxx+ AliSpectraAODPID.cxx+ AliAnalysisTaskSpectraAOD.cxx+");
  //   handler->SetFileForTestMode("filelist.txt"); // list of local files for testing
  //  handler->SetAliRootMode("");
  handler->SetClearPackages();
   
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  mgr->SetGridHandler(handler);
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
   
  // Add PID task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  Bool_t isMC = kFALSE;
  if (c == 2 || c == 4) isMC = kTRUE;   
  Printf("-------------------------------adding in runAOD AddTaskPIDResponse");
  AliAnalysisTask * taskPID = AddTaskPIDResponse(isMC);
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  //AddTaskPIDqa();
   
  //LOOP OVER SELECTION
  Double_t CentCutMin[4]={0,0,20,20};
  Double_t CentCutMax[4]={100,5,40,40};
  Double_t QvecCutMin[4]={0,0,0,3};
  Double_t QvecCutMax[4]={100,100,2,100};
  // Double_t CentCutMin[4]={5,10,20,30};
  // Double_t CentCutMax[4]={10,20,30,40};
  // Double_t QvecCutMin[4]={0,0,0,0};
  // Double_t QvecCutMax[4]={100,100,100,100};
  using namespace AliSpectraNameSpace;
  AliSpectraAODPID *pid = new AliSpectraAODPID(kNSigmaTPCTOF); 
  pid->SetNSigmaCut(3.);
   
  for(Int_t iCut=1;iCut<2;iCut++){
    AliAnalysisTaskSpectraAOD *task = new AliAnalysisTaskSpectraAOD("TaskAODExercise");
    mgr->AddTask(task);
    //physics selection
    task->SelectCollisionCandidates();     
    // set pid object
    task->SetPID(pid);
    // Set the cuts
    AliSpectraAODEventCuts * vcuts = new AliSpectraAODEventCuts("Event Cuts");
    AliSpectraAODTrackCuts  * tcuts = new AliSpectraAODTrackCuts("Track Cuts");
    if(c==1 || c==2)tcuts->SetTrackType(5); //AOD 046 & 047. Standard Cut with loose DCA
    //if(c==1 || c==2)tcuts->SetTrackType(6); //AOD 046 & 047. Standard Cut with tight DCA
    //if(c==3 || c==4)tcuts->SetTrackType(10); //AOD 086 & 090. Standard Raa cut
    //if(c==3 || c==4)tcuts->SetTrackType(4); //AOD 086 & 090. Jet analysis
    if(c==3 || c==4)tcuts->SetTrackType(7); //AOD 086 & 090. TPC Only
    tcuts->SetEta(.8);
    //tcuts->SetDCA(.1);
    tcuts->SetPt(5);
    tcuts->SetY(.5);
    tcuts->SetPtTOFMatching(.6);   
    tcuts->SetQvecMin(QvecCutMin[iCut]);   
    tcuts->SetQvecMax(QvecCutMax[iCut]);    
    vcuts->SetCentralityCutMax(CentCutMax[iCut]);  
    vcuts->SetCentralityCutMin(CentCutMin[iCut]);
    task->SetEventCuts(vcuts);
    task->SetTrackCuts(tcuts);
    vcuts->PrintCuts();
    tcuts->PrintCuts();
     
    // check for MC or real data
    if (c == 2 || c == 4)
      {
	task->SetIsMC(kTRUE);
	vcuts->SetIsMC(kTRUE);
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("chistpt%d",iCut), AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, 
								    Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer(Form("cvcutpt%d",iCut), AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, 
								    Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer(Form("ctcutpt%d",iCut), AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, 
								    Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	AliAnalysisDataContainer *coutputpt4 = mgr->CreateContainer(Form("cpidpt%d",iCut),  AliSpectraAODPID::Class(),     AliAnalysisManager::kOutputContainer, 
								    Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
      }
    if (c == 1 || c==3)
      {
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("chistpt%d",iCut), AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, 
								    Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer(Form("cvcutpt%d",iCut), AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, 
								    Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer(Form("ctcutpt%d",iCut), AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, 
								    Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	AliAnalysisDataContainer *coutputpt4 = mgr->CreateContainer(Form("cpidpt%d",iCut),  AliSpectraAODPID::Class(),     AliAnalysisManager::kOutputContainer, 
								    Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	 
      }
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutputpt1);
    mgr->ConnectOutput(task, 2, coutputpt2);
    mgr->ConnectOutput(task, 3, coutputpt3);
    mgr->ConnectOutput(task, 4, coutputpt4);
  }
  mgr->SetDebugLevel(2);
   
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof");
}
