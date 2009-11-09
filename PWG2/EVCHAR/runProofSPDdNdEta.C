/*************************************************************************
* Macro runProofSPDdNdEta                                                *
* To run dN/dEta reconstruction analysis                                 *
*                                                                        *
* Author:  M. Nicassio (INFN Bari)                                       *
* Contact: Maria.Nicassio@ba.infn.it, Domenico.Elia@ba.infn.it           *
**************************************************************************/

void runProofSPDdNdEta (Bool_t kpythia, Int_t MBTrigg, Bool_t kreadmc, Bool_t kppAna, 
                        Bool_t kallStat, Int_t nEntries, Int_t firstEntry, Char_t* dataSet) {

 gEnv->SetValue("XSec.GSI.DelegProxy","2");
 // Connecting to the PROOF cluster
 TProof::Open("alicecaf");

 gProof->ClearPackages();
 // Enable the needed packages
 gProof->UploadPackage("AF-v4-17");
 gProof->EnablePackage("AF-v4-17");

// gProof->ShowEnabledPackages(); 

 // Create the analysis manager
 mgr = new AliAnalysisManager("SPD dN/dEta Analysis");

 // Create, add task
 gProof->Load("AliAnalysisTaskSPDdNdEta.cxx++g");
 task = new AliAnalysisTaskSPDdNdEta();
 mgr->AddTask(task);

 // Add ESD handler
 AliESDInputHandler* esdH = new AliESDInputHandler;
 mgr->SetInputEventHandler(esdH);
 // Add MC handler
 if (kreadmc) {
   AliMCEventHandler*   mcH = new AliMCEventHandler();
   mgr->SetMCtruthEventHandler(mcH);
 } 
 task->SetReadMC(kreadmc);

 task->SetTrigger(MBTrigg);

 task->SetEvtGen(kpythia);
 
 // Attach input
 cInput = mgr->GetCommonInputContainer();
 mgr->ConnectInput(task, 0, cInput);

 // Attach output
 if (kreadmc) cOutput= mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer, "SPDdNdEtaMC.root");

 else cOutput= mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer, "SPDdNdEtaData.root");

 mgr->ConnectOutput(task, 0, cOutput);

 // Enable debug printouts
 mgr->SetDebugLevel(2);

 // Run analysis
 mgr->InitAnalysis();
 mgr->PrintStatus();
 
 if (kallStat)
   mgr->StartAnalysis("proof", dataSet);
 else 
   mgr->StartAnalysis("proof", dataSet, nEntries, firstEntry);

}
