
/*_____________________________________________________________
 
 Addtask for simulations taks : HF Correlations
 Jitendra Kumar (jitendra.kumar@cern.ch)
 Andrea Rossi   (andrea.rossi@cern.ch)
 _____________________________________________________________*/


AliAnalysisHFCorrOnFlySim *AddTaskHFCorrOnFlySim(TString suffixName ="", Bool_t useweights=kFALSE, TString namefile="", TString namehistowgt="") {

   AliAnalysisHFCorrOnFlySim* clus = new  AliAnalysisHFCorrOnFlySim("");
   clus->SetEtaRange(-20.0, 20.0);
   clus->SetPtRange(0.3, 1000.0);
   clus->SetYRange(-20., 20.);
   clus->SetMultRange(0,5000);
   clus->SetEventProperties(kTRUE);
   clus->SetPartProperties(kTRUE);
   clus->SetHFCorrelations(kTRUE);
   clus->SetHHCorrelations(kFALSE);
   clus->SetQQbarCorrBetween("c", 1, "c", -1);
   clus->SetQQbarCorrelations(kTRUE);

   if (useweights) {
     TFile *f = TFile::Open(namefile.Data());
     TH1D* hwgt = (TH1D*)f->Get(namehistowgt.Data());
     clus->SetUsePtWeights(kTRUE);
     clus->SetHistoWeights(hwgt);
   }
    
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
     Printf("AliAnalysisHFCorrOnFlySim, No analysis manager to connect to.");
     return NULL;
   }
  
   if(!mgr->GetMCtruthEventHandler()){
     Printf("AliAnalysisHFCorrOnFlySim; This task requires an input MC event handler");
     return NULL;
   }
   
   mgr->AddTask(clus);
  
  // Create containers for input/output
  TString finDirname   = "";
  TString inname       = "cinput";
  TString outBasic     = "BasicPlots";
  TString Specific     = "Specific";
  
  finDirname += suffixName.Data();
  inname           +=   finDirname.Data();
  outBasic         +=   finDirname.Data();
  Specific         +=   finDirname.Data();
   

  //Input and Output Slots:
  AliAnalysisDataContainer *cinputSim = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":KineSimulations";
  //TString outputfile = "AnaKineResults.root";

  AliAnalysisDataContainer *coutputSim1 = mgr->CreateContainer(outBasic,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputSim2 = mgr->CreateContainer(Specific,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  mgr->ConnectInput(clus,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(clus,1,coutputSim1);
  mgr->ConnectOutput(clus,2,coutputSim2);

  return clus;
}
