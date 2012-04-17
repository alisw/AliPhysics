void runAODProof(Int_t c=2, const char * proofMode = "full")
{

   gEnv->SetValue("XSec.GSI.DelegProxy", "2");

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
   handler->SetAliROOTVersion("v5-03-13-AN");
   
   handler->SetProofCluster(Form("%s@alice-caf.cern.ch", gSystem->Getenv("CAFUSER")));
   //handler->SetProofCluster(Form("%s@skaf.saske.sk",gSystem->Getenv("CAFUSER")));
   
   // Set handler for Real DATA:
   if (c == 1)
     {
       handler->SetProofDataSet("/alice/data/LHC10h_000138653_p2_AOD049#aodTree");
       //handler->SetProofDataSet("/alice/data/LHC10h_000138653_p2_AOD049#aodTree|/alice/data/LHC10h_000138662_p2_AOD049#aodTree|/alice/data/LHC10h_000138666_p2_AOD049#aodTree|/alice/data/LHC10h_000138795_p2_AOD049#aodTree|/alice/data/LHC10h_000139107_p2_AOD049#aodTree|/alice/data/LHC10h_000139110_p2_AOD049#aodTree");
     }
   if (c == 2)
     {
       handler->SetProofDataSet("/alice/sim/LHC11a10a_000138653_AOD048#aodTree|/alice/sim/LHC11a10a_000138662_AOD048#aodTree|/alice/sim/LHC11a10a_000138666_AOD048#aodTree|/alice/sim/LHC11a10a_000138795_AOD048#aodTree|/alice/sim/LHC11a10a_000139107_AOD048#aodTree|/alice/sim/LHC11a10a_000139110_AOD048#aodTree");      
     }
   handler->SetNproofWorkersPerSlave(1);
   handler->SetAliRootMode("default");
   handler->SetAdditionalLibs("AliSpectraAODHistoManager.cxx AliSpectraAODHistoManager.h AliSpectraAODEventCuts.cxx AliSpectraAODEventCuts.h AliSpectraAODTrackCuts.cxx AliSpectraAODTrackCuts.h AliAnalysisTaskSpectraAOD.cxx AliAnalysisTaskSpectraAOD.h");
   handler->SetAnalysisSource("AliSpectraAODHistoManager.cxx+ AliSpectraAODEventCuts.cxx+ AliSpectraAODTrackCuts.cxx+ AliAnalysisTaskSpectraAOD.cxx+");
   handler->SetFileForTestMode("filelist.txt"); // list of local files for testing
   //  handler->SetAliRootMode("");
   handler->SetClearPackages();


   AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
   mgr->SetGridHandler(handler);
   AliAODInputHandler* aodH = new AliAODInputHandler();
   mgr->SetInputEventHandler(aodH);

   gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
   gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
   gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
   gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");


   // Add PID task
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   Bool_t isMC = kFALSE;
   if (c == 2 || c == 3) isMC = kTRUE;   
   AliAnalysisTask * taskPID = AddTaskPIDResponse(isMC);
   mgr->AddTask(taskPID);

   //LOOP OVER SELECTION
   Double_t CentCutMin[4]={0,0,20,20};
   Double_t CentCutMax[4]={100,5,40,40};
   Double_t QvecCutMin[4]={0,0,0,3};
   Double_t QvecCutMax[4]={100,100,2,100};
   for(Int_t iCut=1;iCut<2;iCut++){
     AliAnalysisTaskSpectraAOD *task = new AliAnalysisTaskSpectraAOD("TaskAODExercise");
     mgr->AddTask(task);
     //physics selection
     task->SelectCollisionCandidates();     
     // Set the cuts
     AliSpectraAODEventCuts * vcuts = new AliSpectraAODEventCuts("Event Cuts");
     AliSpectraAODTrackCuts  * tcuts = new AliSpectraAODTrackCuts("Track Cuts");
     tcuts->SetTrackType(5);
     tcuts->SetEta(.8);
     //tcuts->SetDCA(.1);
     tcuts->SetPt(2.5);
     tcuts->SetPtTOFMatching(0.6);   
     tcuts->SetQvecMin(QvecCutMin[iCut]);   
     tcuts->SetQvecMax(QvecCutMax[iCut]);    
     vcuts->SetCentralityCutMax(CentCutMax[iCut]);  
     vcuts->SetCentralityCutMin(CentCutMin[iCut]);
     task->SetEventCuts(vcuts);
     task->SetTrackCuts(tcuts);
     task->SetNSigmaForIdentification(5.); // FIXME
     task->SetYCut(.5);
     vcuts->PrintCuts();
     tcuts->PrintCuts();
     
     // check for MC or real data
     if (c == 2 || c == 3)
       {
	 task->SetIsMC(kTRUE);
	 AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	 AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("chistpt%d",iCut), AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, 
								     Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	 AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer(Form("cvcutpt%d",iCut), AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, 
								     Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	 AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer(Form("ctcutpt%d",iCut), AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, 
								     Form("Pt.AOD.1._MC_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
       }
     if (c == 1)
       {
	 AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	 AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("chistpt%d",iCut), AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, 
								     Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	 AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer(Form("cvcutpt%d",iCut), AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, 
								     Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	 AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer(Form("ctcutpt%d",iCut), AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, 
								     Form("Pt.AOD.1._data_ptcut_Cent%.0fto%.0f_QVec%.1fto%.1f.root",CentCutMin[iCut],CentCutMax[iCut],QvecCutMin[iCut],QvecCutMax[iCut]));
	 
       }
     mgr->ConnectInput(task, 0, cinput);
     mgr->ConnectOutput(task, 1, coutputpt1);
     mgr->ConnectOutput(task, 2, coutputpt2);
     mgr->ConnectOutput(task, 3, coutputpt3);
   }
   mgr->SetDebugLevel(2);
   
   if (!mgr->InitAnalysis()) return;
   mgr->PrintStatus();
   if (c == 3)
   {
      mgr->StartAnalysis("local");
   }
   if (c != 3)
   {
      mgr->StartAnalysis("proof");
   }
}
