AliAnalysisTaskSEDvsMultiplicity *AddTaskDvsMultiplicity(Int_t system=0,
							 Bool_t readMC=kFALSE,
							 Int_t MCOption=0,
							 Int_t pdgMeson="411",
							 TString finDirname="Loose",
							 TString filename="",
							 TString finAnObjname="AnalysisCuts", 
							 TString estimatorFilename="")
{
  //                                                                                                                                    
  // Test macro for the AliAnalysisTaskSE for D+ candidates 

  //Invariant mass histogram and                                                 
  // association with MC truth (using MC info in AOD)                                                                                   
  //  R. Bala, bala@to.infn.it                                                                                                                                  
  // Get the pointer to the existing analysis manager via the static access method.                                                     
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDvsMultiplicity", "No analysis manager to connect to.");
  }

  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( filename.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
      filecuts=TFile::Open(filename.Data());
      if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
	AliFatal("Input file not found : check your cut object");
      }
  }
  
  

                                      

  
  //Analysis Task

 AliRDHFCuts *analysiscuts=0x0;
  
 TString Name;
 if(pdgMeson=411){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDplustoKpipi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(finAnObjname);
    Name="Dplus";
 }else if(pdgMeson=421){
  
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(finAnObjname);
     Name="D0";
 }else if(pdgMeson=413){
if(stdcuts) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(finAnObjname);
Name="DStar";
 }


  
 AliAnalysisTaskSEDvsMultiplicity *dMultTask = new AliAnalysisTaskSEDvsMultiplicity("dMultAnalysis",pdgMeson,analysiscuts);
  dMultTask->SetReadMC(readMC);
 
  dMultTask->SetDebugLevel(0);
 
  dMultTask->SetUseBit(kTRUE);

  dMultTask->SetDoImpactParameterHistos(kFALSE);

  if(estimatorFilename.EqualTo("") ) {
    printf("Estimator file not provided, multiplcity corrected histograms will not be filled\n");
  }else{
     const Char_t* periodNames[4] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e"};
     TProfile* multEstimatorAvg[4];                       
     TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
     if(!fileEstimator)  {
       AliFatal("File with multiplicity estimator not found\n"); 
       return;
     }
     for(Int_t ip=0; ip<4; ip++) {
       multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("SPDmult10_%s",periodNames[ip]))->Clone(Form("SPDmult10_%s_clone",periodNames[ip])));  
    

     }
     dMultTask->SetMultiplVsZProfileLHC10b(multEstimatorAvg[0]);
     dMultTask->SetMultiplVsZProfileLHC10c(multEstimatorAvg[1]);
     dMultTask->SetMultiplVsZProfileLHC10d(multEstimatorAvg[2]);
     dMultTask->SetMultiplVsZProfileLHC10e(multEstimatorAvg[3]);
  }
  mgr->AddTask(dMultTask);
  
  // Create containers for input/output 

  TString inname = "cinput";
  TString outname = "coutput";
  TString cutsname = "coutputCuts";
  TString normname = "coutputNorm";
 
  inname += finDirname.Data();
  outname += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
 
  
    


TString contname=Form("cinput%s",Name.Data());
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(contname.Data(),TChain::Class(),
							       AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DMult";
  
 contname=Form("coutputCuts%s",Name.Data());
  AliAnalysisDataContainer *coutputCuts = mgr->CreateContainer(contname.Data(),TList::Class(),
								    AliAnalysisManager::kOutputContainer,
								    outputfile.Data());
  
 contname=Form("coutput%s",Name.Data());
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(contname.Data(),TList::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
 contname=Form("coutputNorm%s",Name.Data());  
AliAnalysisDataContainer *coutputNorm = mgr->CreateContainer(contname.Data(),TList::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  
  
  mgr->ConnectInput(dMultTask,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(dMultTask,1,coutput);
  
  mgr->ConnectOutput(dMultTask,2,coutputCuts);

  mgr->ConnectOutput(dMultTask,3,coutputNorm);  
 
  return dMultTask;
}
