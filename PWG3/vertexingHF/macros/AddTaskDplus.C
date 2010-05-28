AliAnalysisTaskSEDplus *AddTaskDplus(Bool_t storeNtuple=kTRUE,
				     Bool_t readMC=kTRUE)
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
    ::Error("AddTaskDplus", "No analysis manager to connect to.");
    return NULL;
  }

  //ptbins
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
    prodcutsval[11][ipt]=10000000.0;
    
  }



 Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
 //Double_t cutsDplus[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};
   Int_t ic=0;
   for(Int_t ipt=0;ipt<nptbins;ipt++){
     anacutsval[ic][ipt]=0.2;
   }
   Int_t ic=1;
   for(Int_t ipt=0;ipt<nptbins;ipt++){
     anacutsval[ic][ipt]=0.4;
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
     anacutsval[ic][ipt]=0.01;
   }
   Int_t ic=11;
   for(Int_t ipt=0;ipt<nptbins;ipt++){
     anacutsval[ic][ipt]=10000000000.;
   }
   anacutsval[6][0]=0.022100;
   anacutsval[6][1]=0.034;
   anacutsval[6][2]=0.020667;
   anacutsval[6][3]=0.023333;

   anacutsval[7][0]=0.08;
   anacutsval[7][1]=0.09;
   anacutsval[7][2]=0.095;
   anacutsval[7][3]=0.115;

   anacutsval[8][0]=0.5;
   anacutsval[8][1]=1.0;
   anacutsval[8][2]=0.5;
   anacutsval[8][3]=0.5;

   anacutsval[9][0]=0.979;
   anacutsval[9][1]=0.9975;
   anacutsval[9][2]=0.995;
   anacutsval[9][3]=0.9975;

   anacutsval[10][0]=0.0055;
   anacutsval[10][1]=0.0028;
   anacutsval[10][2]=0.000883;
   anacutsval[10][3]=0.000883;



 
//Production cuts

  AliRDHFCutsDplustoKpipi *prodcuts = new AliRDHFCutsDplustoKpipi();
  prodcuts->SetPtBins(nptbins+1,ptbins);
  prodcuts->SetCuts(nvars,nptbins,prodcutsval);

  //Analysis cuts
  AliRDHFCutsDplustoKpipi *analysiscuts = new AliRDHFCutsDplustoKpipi();
  analysiscuts->SetPtBins(nptbins+1,ptbins);
  analysiscuts->SetCuts(nvars,nptbins,anacutsval);


  // Aanalysis task
  AliAnalysisTaskSEDplus *dplusTask = new AliAnalysisTaskSEDplus("DplusAnalysis",analysiscuts,prodcuts,storeNtuple);
  dplusTask->SetReadMC(readMC);
  dplusTask->SetDoLikeSign(kTRUE);
  dplusTask->SetDebugLevel(0);

  mgr->AddTask(dplusTask);

 // Create containers for input/output 

  AliAnalysisDataContainer *cinputDplus = mgr->CreateContainer("cinputDplus",TChain::Class(),
                                                          AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassDplus";

 AliAnalysisDataContainer *coutputDplusCuts = mgr->CreateContainer("coutputDplusCuts",TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
 								outputfile.Data());
 
  AliAnalysisDataContainer *coutputDplus = mgr->CreateContainer("coutputDplus",TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								outputfile.Data());

  if(storeNtuple){
    AliAnalysisDataContainer *coutputDplus2 = mgr->CreateContainer("coutputDplus2",TNtuple::Class(),
                                                           AliAnalysisManager::kOutputContainer,
								 "InvMassDplus_nt1.root");

    coutputDplus2->SetSpecialOutput();
  }
  mgr->ConnectInput(dplusTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(dplusTask,1,coutputDplus);

  mgr->ConnectOutput(dplusTask,2,coutputDplusCuts);
  
  if(storeNtuple){
    mgr->ConnectOutput(dplusTask,3,coutputDplus2);
  }
  return dplusTask;
}
