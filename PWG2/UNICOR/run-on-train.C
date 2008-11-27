{
gROOT->LoadMacro("makechain.C"); 
tr = makechain("esdTree","/u/sma/data/mc/v4-13-Rev-01/fulltrd_pdc_0.txt",10);

gSystem->Load("libEG.so");
gSystem->Load("libPhysics.so");
gSystem->Load("libTree.so");
gSystem->Load("libVMC.so"); 
gSystem->Load("libSTEERBase.so");
gSystem->Load("libESD.so");
gSystem->Load("libANALYSIS");
gSystem->Load("libUNICOR.so");

AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
AliVEventHandler* esdH = new AliESDInputHandler;
mgr->SetInputEventHandler(esdH);  
AliDAnalysisTask *mytask = new AliDAnalysisTask();
mgr->AddTask(mytask);
AliAnalysisDataContainer *cinput = mgr->CreateContainer("cinput",TChain::Class(),AliAnalysisManager::kInputContainer);
AliAnalysisDataContainer *coutpt = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,"kuku.root");
mgr->ConnectInput (mytask,0,cinput);
mgr->ConnectOutput(mytask,0,coutpt);
mgr->InitAnalysis();
mgr->PrintStatus(); 
mgr->StartAnalysis("local",tr);
}
