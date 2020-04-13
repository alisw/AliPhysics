AliAnalysisTaskQA* AddTask_QA(TString   outputname = "MB_tree",
			      TString   photonCutNumberV0Reader       = "",//"060000084001001500000000",
			      TString   TaskEventCutnumber            = "",//00000003",
			      TString   TaskPhotonCutnumber           = "",//"090000092663743800000000",
			      Bool_t    isMC                          = kFALSE,
			      Int_t     IsHeavyIon                    = 2,
			      Bool_t    hasHistogram                  = kTRUE,
			      Bool_t    enableWriteVariableTree       = kFALSE,// enables V0finding efficiency histograms
			      Bool_t    doEtaShiftV0Reader            = kFALSE,
			      Bool_t    enableV0findingEffi           = kFALSE,              // enables V0finding efficiency histograms
			      TString   NameInactiveBranch            = ""// "X;Y;Z"
			      ){


  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return NULL;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton     = photonCutNumberV0Reader.Data();
  TString cutnumberEvent      = "00000003";
  if(IsHeavyIon==2) cutnumberEvent = "80000003";

  //========= Check V0 Reader in  ANALYSIS manager  =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return NULL;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }
  
  AliConvEventCuts *analysisEventCuts = new AliConvEventCuts();
  analysisEventCuts->SetV0ReaderName(V0ReaderName);
  analysisEventCuts->InitializeCutsFromCutString(TaskEventCutnumber.Data());
  analysisEventCuts->SetFillCutHistograms("",kFALSE);

  AliConversionPhotonCuts *analysisCuts = new AliConversionPhotonCuts();
  analysisCuts->SetV0ReaderName(V0ReaderName);
  analysisCuts->InitializeCutsFromCutString(TaskPhotonCutnumber.Data());
  analysisCuts->SetFillCutHistograms("",kFALSE);

  AliAnalysisTaskQA *task = new AliAnalysisTaskQA(Form("%s_%s_QA",TaskEventCutnumber.Data(),TaskPhotonCutnumber.Data()));
  task->SetEventCuts(analysisEventCuts,IsHeavyIon);
  task->SetConversionCuts(analysisCuts,IsHeavyIon);
  task->FillType(hasHistogram);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetWriteVariableTree(enableWriteVariableTree);
  //activete or deactivate branches depending on the puropose, output file size etc. 
  task->SetTreeInactiveBranch(NameInactiveBranch.Data());
  // task->SetTreeInactiveBranch("eta");
  // task->SetTreeInactiveBranch("phi");
  // task->SetTreeInactiveBranch("theta");
  // task->SetTreeInactiveBranch("chi2");
  // task->SetTreeInactiveBranch("qt");
  // task->SetTreeInactiveBranch("psipair");
  // task->SetTreeInactiveBranch("cosPA");
  // task->SetTreeInactiveBranch("InvMass");
  // task->SetTreeInactiveBranch("X");
  // task->SetTreeInactiveBranch("Y");
  // task->SetTreeInactiveBranch("Z");
  // task->SetTreeInactiveBranch("R");
  // task->SetTreeInactiveBranch("Qual");
  // task->SetTreeInactiveBranch("DCAz");
  // task->SetTreeInactiveBranch("DCAr");
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("%s_%s_%s",outputname.Data(),TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer(Form("tree_%s_%s_%s",outputname.Data(),TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()),
			 TTree::Class(),
			 AliAnalysisManager::kOutputContainer, "AnalysisResults.root");  
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectInput(task,0,cinput);


  //connect containers
  return task;
}

