void runTask(Float_t etamax=0.5,const char * incollection = 0, const char * outfile = "dndeta.root")
{
  // for running with root only
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 

  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");


  TChain * chain = new TChain ("TE");
  if (incollection == 0) {
    chain->Add("galice.root");
  }
  else if (TString(incollection).Contains("xml")){
    TGrid::Connect("alien://");
    TGridCollection * coll = gGrid->OpenCollection(incollection);
    while(coll->Next()){
      chain->Add(TString("alien://")+coll->GetLFN());
    }
  } else {
    ifstream file_collect(incollection);
    TString line;
    while (line.ReadLine(file_collect) ) {
      chain->Add(line.Data());
    }
  }
  chain->GetListOfFiles()->Print();

  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$(ALICE_ROOT)/include
  // or in each macro
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("dNdeta");

  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);

  // Create tasks
  gROOT->LoadMacro("AliAnalysisTaskHMTFMC.cxx+g");

  AliAnalysisTask *task1 = new AliAnalysisTaskHMTFMC("TaskdNdeta");
  // ((AliAnalysisTaskHMTFMC*)task1)->SetEtaMax(etamax);
  // Enable MC event handler

  AliMCEventHandler* handler = new AliMCEventHandler;
  handler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(handler);

  // Add tasks
  mgr->AddTask(task1);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutput", TList::Class(),    AliAnalysisManager::kOutputContainer, outfile);

  // Connect input/output
  mgr->ConnectInput(task1, 0, cinput);
  mgr->ConnectOutput(task1, 1, coutput1);

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();

  mgr->StartAnalysis("local", chain);
}
