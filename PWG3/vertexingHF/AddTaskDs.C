AliAnalysisTaskSEDs *AddTaskDs(Bool_t readMC=kTRUE)
{
  //                                                                           
  // Test macro for the AliAnalysisTaskSE for Ds candidates 


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDs", "No analysis manager to connect to.");
    return NULL;
  }



  //ptbins
  const Int_t nptbins=4;
  Float_t ptbins[nptbins+1]={0.,2.,3.,5.,99999.};

  //setting cut values
  Int_t nvars=14;
  Float_t** ancutsval;
  ancutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){ancutsval[ic]=new Float_t[nptbins];}  
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    ancutsval[0][ipt]=0.2;
    ancutsval[1][ipt]=0.4;
    ancutsval[2][ipt]=0.4;
    ancutsval[3][ipt]=0.;
    ancutsval[4][ipt]=0.;
    ancutsval[5][ipt]=0.005;
    ancutsval[6][ipt]=0.038;
    ancutsval[7][ipt]=0.;
    ancutsval[8][ipt]=0.;
    ancutsval[9][ipt]=0.95;
    ancutsval[10][ipt]=0.;
    ancutsval[11][ipt]=0.1;
    ancutsval[12][ipt]=0.004;
    ancutsval[13][ipt]=0.035;
  }
  Float_t** prcutsval;
  prcutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){prcutsval[ic]=new Float_t[nptbins];}  
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    prcutsval[0][ipt]=0.2;
    prcutsval[1][ipt]=0.4;
    prcutsval[2][ipt]=0.4;
    prcutsval[3][ipt]=0.;
    prcutsval[4][ipt]=0.;
    prcutsval[5][ipt]=0.005;
    prcutsval[6][ipt]=0.06;
    prcutsval[7][ipt]=0.;
    prcutsval[8][ipt]=0.;
    prcutsval[9][ipt]=0.85;
    prcutsval[10][ipt]=0.;
    prcutsval[11][ipt]=0.1;
    prcutsval[12][ipt]=0.1;
    prcutsval[13][ipt]=0.1;
  }

  //Analysis cuts
  AliRDHFCutsDstoKKpi *analysiscuts = new AliRDHFCutsDstoKKpi();
  analysiscuts->SetPtBins(nptbins+1,ptbins);
  analysiscuts->SetCuts(nvars,nptbins,ancutsval);
  AliRDHFCutsDstoKKpi *prodcuts = new AliRDHFCutsDstoKKpi();
  prodcuts->SetPtBins(nptbins+1,ptbins);
  prodcuts->SetCuts(nvars,nptbins,prcutsval);

  // Analysis task                                                                                                                     
  AliAnalysisTaskSEDs *dsTask = new AliAnalysisTaskSEDs("DsAnalysis",prodcuts,analysiscuts);
  dsTask->SetReadMC(readMC);
  dsTask->SetDebugLevel(0);
  mgr->AddTask(dsTask);

  //                                                                                                                                    
  // Create containers for input/output                                                                                                 
  AliAnalysisDataContainer *cinputDs = mgr->CreateContainer("cinputDs",TChain::Class(),
							    AliAnalysisManager::kInputContainer);

 TString outputfile = AliAnalysisManager::GetCommonFileName(); 
 outputfile += ":PWG3_D2H_InvMassDs";
 
 AliAnalysisDataContainer *coutputDsCuts = mgr->CreateContainer("coutputDsCuts",TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
 								outputfile.Data());

 
  AliAnalysisDataContainer *coutputDs = mgr->CreateContainer("coutputDs",TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputfile.Data());


  mgr->ConnectInput(dsTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(dsTask,1,coutputDs);

  mgr->ConnectOutput(dsTask,2,coutputDsCuts);
  return dsTask;
}
