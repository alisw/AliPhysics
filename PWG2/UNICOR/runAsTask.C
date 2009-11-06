{
gSystem->Load("libPhysics.so");
gSystem->Load("libEG.so");
gSystem->Load("libTree.so");
gSystem->Load("libVMC.so"); 
gSystem->Load("libSTEERBase.so");
gSystem->Load("libESD.so");
gSystem->Load("libAOD.so");
gSystem->Load("libANALYSIS");
gSystem->Load("libANALYSISalice");
gSystem->Load("libPWG2unicor");

gROOT->LoadMacro("makechain.C");
tr = makechain("esdTree","filelist.txt");

AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
AliVEventHandler* esdH = new AliESDInputHandler;
mgr->SetInputEventHandler(esdH);  

gROOT->LoadMacro("AddTaskUnicor.C");

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
