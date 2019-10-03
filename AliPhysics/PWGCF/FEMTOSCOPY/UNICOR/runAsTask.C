//=============================================================================
void runAsTask(int grid=0) {
  gSystem->Load("libPhysics");
  gSystem->Load("libEG");
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG2unicor");
  
  if (grid) {
    if (!TGrid::Connect("alien://")) return;
    TChain *tr = CreateChainFromTags("wn.xml", "ESD");
  } else {
    gROOT->LoadMacro("makechain.C");
    tr = makechain("esdTree","filelist.txt");
  }

  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/UNICOR/AddTaskUnicor.C");
  AliAnalysisTaskUnicor *mytask = AddTaskUnicor();

  mgr->InitAnalysis();
  mgr->PrintStatus(); 
  mgr->StartAnalysis("local",tr);

  TFile::Open("AnalysisResults.root","read");
  gDirectory->Cd("PWG2UNICOR");
  TList *list = (TList *) gDirectory->Get("unilis");
  char *outfil = "unicor-result-as-anal.root";
  for (int i=0; i<list->GetEntries(); i++) {
    AliUnicorAnal *an = (AliUnicorAnal *) list->At(i);
    if (i==0) an->Save(outfil,"recreate");
    else an->Save(outfil);
    delete an;
  }
}
//=============================================================================
TChain* CreateChainFromTags(const char *xmlfile, const char *type="ESD")
{
  // Create a chain using tags from the xml file.
  TGridCollection* coll = gGrid->OpenCollection(xmlfile);
  if (!coll) {
    ::Error("CreateChainFromTags", "Cannot create an AliEn collection from %s", xmlfile);
    return NULL;
  }
  TGridResult* tagResult = coll->GetGridResult("",kFALSE,kFALSE);
  AliTagAnalysis *tagAna = new AliTagAnalysis(type);
  tagAna->ChainGridTags(tagResult);
  
  AliRunTagCuts      *runCuts = new AliRunTagCuts();
  AliLHCTagCuts      *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts    *evCuts  = new AliEventTagCuts();
  // Check if the cuts configuration file was provided
  if (!gSystem->AccessPathName("ConfigureCuts.C")) {
    gROOT->LoadMacro("ConfigureCuts.C");
    ConfigureCuts(runCuts, lhcCuts, detCuts, evCuts);
  }
  TChain *chain = tagAna->QueryTags(runCuts, lhcCuts, detCuts, evCuts);
  if (!chain || !chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}
//=============================================================================
