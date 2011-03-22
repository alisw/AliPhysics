AliAnalysisTaskSESelectHF4Prong *AddTaskSelectHF4Prong()
{
  //
  // Test macro for the AliAnalysisTaskSE for heavy-flavour selection
  // and creation of a stand-alone AOD
  // A.Dainese, andrea.dainese@lnl.infn.it
  // F.Colamaria, fabio.colamaria@ba.infn.it
  //

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSelectHF4Prong", "No analysis manager to connect to.");
    return NULL;
  }   

  // Output 
  AliAODHandler *aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAOD.VertexingHF.sa.root");
  aodHandler->SetCreateNonStandardAOD();
  mgr->SetOutputEventHandler(aodHandler);
 
  //Cuts loading
  TFile* filecuts=new TFile("Charm4ProngCutsDef.root");
  if(!filecuts->IsOpen()){
    cout<<"Input file not found: exit"<<endl;
    return;
  }

  AliRDHFCutsD0toKpipipi* RDHFCharm4Prong=new AliRDHFCutsD0toKpipipi();
  RDHFCharm4Prong = (AliRDHFCutsD0toKpipipi*)filecuts->Get("Charm4ProngCuts");
  RDHFCharm4Prong->SetName(Form("Charm4ProngCuts%d",1));

  if(!RDHFCharm4Prong){
    cout<<"Specific AliRDHFCuts not found"<<endl;
    return;
  }

  // Analysis task    
  AliAnalysisTaskSESelectHF4Prong *hfTask = new AliAnalysisTaskSESelectHF4Prong("SelectHFAnalysis",RDHFCharm4Prong);
  hfTask->SetDebugLevel(2);
  mgr->AddTask(hfTask);  

  // Create containers for input/output  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *contHist = mgr->CreateContainer("histos_bin1",TList::Class(),AliAnalysisManager::kOutputContainer,"HistMassInvAndCuts.root");
  AliAnalysisDataContainer *contHist2 = mgr->CreateContainer("histos_bin2",TList::Class(),AliAnalysisManager::kOutputContainer,"HistMassInvAndCuts.root");
  AliAnalysisDataContainer *contHist3 = mgr->CreateContainer("histos_bin3",TList::Class(),AliAnalysisManager::kOutputContainer,"HistMassInvAndCuts.root");
  AliAnalysisDataContainer *contHist4 = mgr->CreateContainer("histos_bin4",TList::Class(),AliAnalysisManager::kOutputContainer,"HistMassInvAndCuts.root");
  AliAnalysisDataContainer *contHist5 = mgr->CreateContainer("histos_bin5",TList::Class(),AliAnalysisManager::kOutputContainer,"HistMassInvAndCuts.root");
  AliAnalysisDataContainer *contHistCuts = mgr->CreateContainer("histoscuts",TList::Class(),AliAnalysisManager::kOutputContainer,"HistMassInvAndCuts.root");

  mgr->ConnectOutput(hfTask,0,mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(hfTask,1,contHist);
  mgr->ConnectOutput(hfTask,2,contHist2);
  mgr->ConnectOutput(hfTask,3,contHist3);
  mgr->ConnectOutput(hfTask,4,contHist4);
  mgr->ConnectOutput(hfTask,5,contHist5);
  mgr->ConnectOutput(hfTask,6,contHistCuts);

  return hfTask;
}
