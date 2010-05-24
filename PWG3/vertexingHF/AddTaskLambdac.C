AliAnalysisTaskSELambdac *AddTaskLambdac(const char *name,Bool_t storeNtuple,Bool_t readMC,Bool_t MCPid,Bool_t realPid,Bool_t resPid,Bool_t useKF)
{
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLambdac", "No analysis manager to connect to.");
    return NULL;
  }

  //pt bins
  const Int_t nptbins=4;
    Float_t* ptbins;
    ptbins=new Float_t[nptbins+1];
    ptbins[0]=0.;
    ptbins[1]=2.;
    ptbins[2]=3.;
    ptbins[3]=5.;
    ptbins[4]=99999.;
    const Int_t nvars=12;

  Float_t** prodcutsval;
  prodcutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){prodcutsval[ic]=new Float_t[nptbins];}
  for(Int_t ipt=0;ipt<nptbins;ipt++){
   prodcutsval[0][ipt]=0.2;
   prodcutsval[1][ipt]=0.4;
   prodcutsval[2][ipt]=0.4;
   prodcutsval[3][ipt]=0.;
   prodcutsval[4][ipt]=0.;
   prodcutsval[5][ipt]=0.01;
   prodcutsval[6][ipt]=0.06;
   prodcutsval[7][ipt]=0.02;
   prodcutsval[8][ipt]=0.;
   prodcutsval[9][ipt]=0.85;
   prodcutsval[10][ipt]=0.;
   prodcutsval[11][ipt]=0.1;
  }

  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  Int_t ic=0;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
       anacutsval[ic][ipt]=0.2;
          }
  Int_t ic=1;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.8;
	     }
  Int_t ic=2;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.4;
	     }
  Int_t ic=3;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.;
	     }
  Int_t ic=4;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.;
	     }
  Int_t ic=5;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.;
	     }
  Int_t ic=6;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.1;
	     }
  Int_t ic=7;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.;
	     }
  Int_t ic=8;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.;
	     }
  Int_t ic=9;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.;
	     }
  Int_t ic=10;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.;
	     }
  Int_t ic=11;
     for(Int_t ipt=0;ipt<nptbins;ipt++){
          anacutsval[ic][ipt]=0.1;
	     }

  anacutsval[1][0]=0.6;

  anacutsval[2][2]=0.5;
  anacutsval[2][3]=0.5;

  anacutsval[9][2]=0.9;
  anacutsval[9][3]=0.9;


  AliRDHFCutsLctopKpi *prodcuts = new AliRDHFCutsLctopKpi();
    prodcuts->SetPtBins(nptbins+1,ptbins);
      prodcuts->SetCuts(nvars,nptbins,prodcutsval);
  //Analysis cuts
  AliRDHFCutsLctopKpi *analysiscuts = new AliRDHFCutsLctopKpi();
  ((AliRDHFCuts*)analysiscuts)->SetPtBins(nptbins+1,ptbins);
  ((AliRDHFCuts*)analysiscuts)->SetCuts(nvars,nptbins,anacutsval);
  // Aanalysis task                                                                                                                     
  AliAnalysisTaskSELambdac *lambdacTask = new AliAnalysisTaskSELambdac("LambdacAnalysis",storeNtuple,analysiscuts,prodcuts);
  lambdacTask->SetReadMC(readMC);
  if(MCPid) lambdacTask->SetMCPid();
  if(resPid) lambdacTask->SetResonantPid();
  if(realPid) lambdacTask->SetRealPid();
  lambdacTask->SetDebugLevel(0);
  if(useKF) {
    lambdacTask->SetUseKF();
     Float_t cuts[10]={0.1,0.1,1.5,0.5,0.1,1.5,0.5,0.1,1.5,0.5};
      lambdacTask->SetCutsKF(cuts);
       }
  mgr->AddTask(lambdacTask);

  //                                                                                                                                    
  // Create containers for input/output                                                                                                 
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassLambdac";
  AliAnalysisDataContainer *cinputLambdac = mgr->CreateContainer("cinputLambdac",TChain::Class(),
                                                          AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutputLambdacCuts = mgr->CreateContainer("coutputLambdacCuts",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  AliAnalysisDataContainer *coutputLambdac = mgr->CreateContainer("coutputLambdac",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  if(storeNtuple){
    AliAnalysisDataContainer *coutputLambdac2 = mgr->CreateContainer("coutputLambdac2",TNtuple::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								 "InvMassLambdac_nt1.root");

    coutputLambdac2->SetSpecialOutput();
  }

  mgr->ConnectInput(lambdacTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(lambdacTask,1,coutputLambdac);
  mgr->ConnectOutput(lambdacTask,2,coutputLambdacCuts);
  
  if(storeNtuple){
    mgr->ConnectOutput(lambdacTask,3,coutputLambdac2);
  }
  return lambdacTask;
}
