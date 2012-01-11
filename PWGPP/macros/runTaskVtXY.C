//
// This is the macro to start the task measuring the beam spot size.
//
void runTaskVtXY(Char_t *dataset="/COMMON/COMMON/LHC09a4_run8100X#esdTree") {

  if (!gProof) {
     cerr<<"Proof session has not been open !"<<endl;
     return;    
  }

 // Enable the needed packages
    gProof->UploadPackage("STEERBase");
    gProof->EnablePackage("STEERBase");
    gProof->UploadPackage("ESD");
    gProof->EnablePackage("ESD");
    gProof->UploadPackage("AOD");
    gProof->EnablePackage("AOD");
    gProof->UploadPackage("ANALYSIS");
    gProof->EnablePackage("ANALYSIS");
    gProof->UploadPackage("ANALYSISalice");
    gProof->EnablePackage("ANALYSISalice"); 

 // Create the analysis manager
 mgr = new AliAnalysisManager("testAnalysis");
 // Create, add task
 gProof->Load("AliAnalysisTaskVtXY.cxx++g");
 task = new AliAnalysisTaskVtXY;
 mgr->AddTask(task);
 // Add ESD handler
 AliESDInputHandler* esdHandler = new AliESDInputHandler;
 mgr->SetInputEventHandler(esdHandler);  
 
  // Attach input
 cInput = mgr->CreateContainer("cInput", TChain::Class(), 
                             AliAnalysisManager::kInputContainer);
 mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

 // Attach output
 cOutput= mgr->CreateContainer("cOutput", TList::Class(), 
                            AliAnalysisManager::kOutputContainer, "VtXY.root");
 mgr->ConnectOutput(task, 0, cOutput);

 // Enable debug printouts
 mgr->SetDebugLevel(2);
 // Run analysis
 mgr->InitAnalysis();
 mgr->PrintStatus();

 mgr->StartAnalysis("proof",dataset);

}
