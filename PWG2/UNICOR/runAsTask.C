{
gSystem->Load("libVMC.so"); 
gSystem->Load("libSTEERBase.so");
gSystem->Load("libESD.so");
gSystem->Load("libAOD.so");
gSystem->Load("libANALYSIS");
gSystem->Load("libANALYSISalice");
gSystem->Load("libPWG2unicor");

gROOT->LoadMacro("makechain.C");
tr = makechain("esdTree","/data.local1/alice/data/silvia/0000*AliESDs.root");    // local emergency copy

AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
AliVEventHandler* esdH = new AliESDInputHandler;
mgr->SetInputEventHandler(esdH);  
AliAnalysisTaskUnicor *mytask = new AliAnalysisTaskUnicor();
mgr->AddTask(mytask);
AliAnalysisDataContainer *cinput = mgr->CreateContainer("cinput",TChain::Class(),AliAnalysisManager::kInputContainer);
AliAnalysisDataContainer *coutpt = mgr->CreateContainer("unilis", TList::Class(),AliAnalysisManager::kOutputContainer,"unicor-result-as-list.root");
mgr->ConnectInput (mytask,0,cinput);
mgr->ConnectOutput(mytask,1,coutpt);
mgr->InitAnalysis();
mgr->PrintStatus(); 
mgr->StartAnalysis("local",tr);
}
