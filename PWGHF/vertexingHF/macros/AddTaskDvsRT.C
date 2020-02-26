AliAnalysisTaskSEDvsRT *AddTaskDvsRT(Int_t system=0,
                                     Bool_t readMC = kFALSE,
                                     Int_t MCOption = 0,
                                     Int_t pdgSpecies = 421,
                                     TString finDirname = "",
                                     TString filename="",
                                     TString finAnObjname = "D0toKpiCuts",
                                     Bool_t useNsparse = kFALSE,
                                     Bool_t isLcV0 = kFALSE
                                    )
{ 
       //=========================================================
       // Macro for RT-dependence measurement of charmed hadrons
       // Output THn of inv. mass vs pT vs RT in toward and away regions
       // ========================================================
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
     ::Error("AddTaskDvsRT", "No analysis manager to connect to");
  }
  
  Bool_t stdcuts = kFALSE;
  TFile* filecuts;
  if (filename.EqualTo("") ) {
    stdcuts = kTRUE;
    } else { 
    filecuts = TFile::Open(filename.Data());
    if(!filecuts || (filecuts && !filecuts->IsOpen())){
       Printf("FATAL: Input file not found, check your cut object!");
       return 0x0;
       }
    }
  
  
  //Load cuts for chosen meson species  
  AliRDHFCuts *analysisCuts = 0x0;
  TString channel = "";
  
  if(pdgSpecies == 411) { //Dplus
     if (stdcuts) {
        analysisCuts = new AliRDHFCutsDplustoKpipi();
        if (system == 0) analysisCuts->SetStandardCutsPP2010();
        else analysisCuts->SetStandardCutsPbPb2011();
     }
     else analysisCuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(finAnObjname);
     channel = "Dplus";
  } else if(pdgSpecies == 421) { //D0
     if (stdcuts) {
        analysisCuts = new AliRDHFCutsD0toKpi();
        if (system == 0) analysisCuts->SetStandardCutsPP2010();
        else analysisCuts->SetStandardCutsPbPb2011();
     }
     else analysisCuts = (AliRDHFCutsD0toKpi*)filecuts->Get(finAnObjname);
     channel = "D0";
  } else if (pdgSpecies == 413) { //DStar
     if (stdcuts) {
        analysisCuts = new AliRDHFCutsDStartoKpipi();
        if (system == 0) analysisCuts->SetStandardCutsPP2010();
        else analysisCuts->SetStandardCutsPbPb2011();
     }
     else analysisCuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(finAnObjname);
     channel = "DStar";
  } else if (pdgSpecies == 431) { //Ds
     if (stdcuts) {
        analysisCuts = new AliRDHFCutsDstoKKpi();
        if (system == 0) analysisCuts->SetStandardCutsPP2010();
        else analysisCuts->SetStandardCutsPbPb2011();
     }
     else analysisCuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(finAnObjname);
     channel = "Ds";
  } else if (pdgSpecies == 4122) { //Lc
     if (stdcuts) {
        analysisCuts = new AliRDHFCutsLctoV0();
        if (system == 0) analysisCuts->SetStandardCutsPP2010();
        else analysisCuts->SetStandardCutsPbPb2011();
     }
     else analysisCuts = (AliRDHFCutsLctoV0*)filecuts->Get(finAnObjname);
     channel="Lc2pK0S";
  } else {
     Printf("No valid analysis species selected, please check macro parameters");
     return 0x0; 
     }
  
  //Initialise analysis task
  AliAnalysisTaskSEDvsRT *dRTTask = new AliAnalysisTaskSEDvsRT("dRTAnalysis",pdgSpecies,analysisCuts);
  dRTTask->SetReadMC(readMC);
  dRTTask->SetDebugLevel(0);
  dRTTask->SetUseBit(kTRUE);
//  dRTTask->SetDoImpactParameterHistos(kFALSE);
  dRTTask->SetMCOption(MCOption);
  dRTTask->SetLctoV0(isLcV0);
  dRTTask->SetUseNsparse(useNsparse);
  
  //RT-specific settings
  dRTTask->SetEtaCut(0.8);
  dRTTask->SetPtLeadMin(5.0);
  dRTTask->SetAveMultiInTrans(4.939);
  
  //channel-specific inv mass settings
  if (pdgSpecies == 421) { 
     dRTTask->SetMassLimits(1.5648,2.1648);
     dRTTask->SetNMassBins(200);
  } else if (pdgSpecies == 4122) {
     dRTTask->SetMassLimits(pdgSpecies,0.250);
     dRTTask->SetNMassBins(1000);
  } else if (pdgSpecies == 411) dRTTask->SetMassLimits(pdgSpecies,0.2);
  
  //this is where multiplicity estimator settings would go
  //!---
  
  //---!
  
  
  // Create container for input/output
  TString inname = "cinput";
  TString outname = "coutput";
  TString cutsname = "coutputCuts";
  TString normname = "coutputNorm";
  
  inname   += channel.Data();
  outname  += channel.Data();
  cutsname += channel.Data();
  normname += channel.Data();
  inname   += finDirname.Data();
  outname  += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname, TChain::Class(), AliAnalysisManager::kInputContainer);
  
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DvsRT_";
  outputfile += channel.Data();
  outputfile += finDirname.Data();
  AliAnalysisDataContainer *coutputCuts = mgr->CreateContainer(cutsname, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(outname, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutputNorm = mgr->CreateContainer(normname, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  
  mgr->ConnectInput(dRTTask,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(dRTTask,1,coutput);
  mgr->ConnectOutput(dRTTask,2,coutputCuts);
  mgr->ConnectOutput(dRTTask,3,coutputNorm);
  
  return dRTTask;
  
}
