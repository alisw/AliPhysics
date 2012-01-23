//DEFINITION OF A FEW CONSTANTS
const Int_t    chargeFrCorr = 20;
//----------------------------------------------------

AliAnalysisTaskSEDStarJets *AddTaskDStarJets(Bool_t theMCon=kTRUE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDStarJets2", "No analysis manager to connect to.");
    return NULL;
  } 

  TFile* filecuts=new TFile("DStartoKpipiCuts.root");
  if(!filecuts->IsOpen()){
    cout<<"Input file not found: exit"<<endl;
    return;
  }

  AliRDHFCutsDStartoKpipi* RDHFDStartoKpipi=new AliRDHFCutsDStartoKpipi();
  RDHFDStartoKpipi = (AliRDHFCutsDStartoKpipi*)filecuts->Get("DStartoKpipiCuts");
  RDHFDStartoKpipi->SetName("DStartoKpipiCuts");

  // mm let's see if everything is ok
  if(!RDHFDStartoKpipi){
    cout<<"Specific AliRDHFCuts not found"<<endl;
    return;
  } 

  //CREATE THE TASK
  printf("CREATE TASK\n");
  // create the task
  AliAnalysisTaskSEDStarJets *task = new AliAnalysisTaskSEDStarJets("AliAnalysisTaskSEDStarJets",RDHFDStartoKpipi);
  task->SetMC(theMCon);
  task->SetChargeFractionCorrection(chargeFrCorr);

  // Create and connect containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DStarJet";

  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  
  // ----- output data -----
  
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("charmJetCorr", TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cuts",AliRDHFCutsDStartoKpipi::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //cuts
  mgr->AddTask(task);
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);

  return task ;
}

