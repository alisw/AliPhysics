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
tr = makechain("esdTree","/data.local1/alice/data/silvia/0000*AliESDs.root");

AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
AliVEventHandler* esdH = new AliESDInputHandler;
mgr->SetInputEventHandler(esdH);  

gROOT->LoadMacro("AddTaskUnicor.C");

AliAnalysisTaskUnicor *mytask = AddTaskUnicor();

mgr->InitAnalysis();
mgr->PrintStatus(); 
mgr->StartAnalysis("local",tr);
}
