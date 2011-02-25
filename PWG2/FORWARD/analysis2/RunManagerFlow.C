void RunManagerFlow(TString path = "", TString file ="AliAODs.root", Int_t nevents = 100, TString type="", Int_t etabins)
{

  // --- Load libs -------------------------------------------------
  
  gSystem->Load("libVMC");
  gSystem->Load("libGeom");
  gSystem->Load("libXMLIO");
  gSystem->Load("libTree");
  
  gSystem->Load("libSTEERBase");
  
  gSystem->Load("libESD") ;
  gSystem->Load("libAOD") ;
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  gSystem->Load("libPhysics");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG2forward");
  gSystem->Load("libPWG2forward2");

  // --- Add to chain either ESD or AOD ----------------------------

  TChain* chain = 0;

  TChain* chainAOD = new TChain("aodTree");
  if (path.Sizeof() == 1) { // Is it a single AOD file?
    chainAOD->Add(Form("%s", file.Data()));
    chain = chainAOD;
  }
  if (path.Sizeof() > 1) {  // More AOD files, assumes AliAODs1.root numbering
    for (Int_t i=1;i<=9;i++) {
      if (Form("%s/%s%d", path.Data(), file.Data(), i))
        chainAOD->Add(Form("%s/AliAODs%d.root", path.Data(), i));
    }
    chain = chainAOD;
  }

  // --- Initiate the event handlers -------------------------------

  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "A test setup for the analysis train");
  mgr->SetUseProgressBar(kTRUE);

  // AOD input handler
  AliAODInputHandler* aodHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodHandler);

  // Create output handler
  AliAODHandler* aodOut = new AliAODHandler();
  mgr->SetOutputEventHandler(aodOut);
  aodOut->SetOutputFileName("AliAODs.Flow.root");

  // --- Add the tasks ---------------------------------------------
  
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskForwardFlow.C");
  AddTaskForwardFlow(type.Data(), etabins);

  // --- Run the analysis ------------------------------------------

  TStopwatch t;
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    t.Start();
    if (nevents == 0) mgr->StartAnalysis("local", chain);
    if (nevents != 0) mgr->StartAnalysis("local", chain, nevents);
    t.Stop();
    t.Print();
  }  

  // ---------------------------------------------------------------
}
//
// EOF
//
