AliAnalysisTaskSEDplus *AddTaskDplus(Bool_t storeNtuple=kFALSE,
				     Bool_t readMC=kFALSE)
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
  const Int_t nptbins=14;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=1.;
  ptbins[2]=2.;
  ptbins[3]=3.;
  ptbins[4]=4.;
  ptbins[5]=5.;
  ptbins[6]=6.;
  ptbins[7]=7.;
  ptbins[8]=8.;
  ptbins[9]=9.;
  ptbins[10]=10.;
  ptbins[11]=12.;
  ptbins[12]=15.;
  ptbins[13]=20.;
  ptbins[14]=99999.;
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
   anacutsval[6][1]=0.022100;
   anacutsval[6][2]=0.034;
   anacutsval[6][3]=0.020667;
   anacutsval[6][4]=0.020667;
   anacutsval[6][5]=0.023333;
   anacutsval[6][6]=0.023333;
   anacutsval[6][7]=0.023333;
   anacutsval[6][8]=0.023333;
   anacutsval[6][9]=0.023333;
   anacutsval[6][10]=0.023333;
   anacutsval[6][11]=0.023333;
   anacutsval[6][12]=0.023333;
   anacutsval[6][13]=0.023333;
   anacutsval[6][14]=0.023333;
   
   anacutsval[7][0]=0.08;
   anacutsval[7][1]=0.08;
   anacutsval[7][2]=0.09;
   anacutsval[7][3]=0.095;
   anacutsval[7][4]=0.095;
   anacutsval[7][5]=0.115;
   anacutsval[7][6]=0.115;
   anacutsval[7][7]=0.115;
   anacutsval[7][8]=0.115;
   anacutsval[7][9]=0.115;
   anacutsval[7][10]=0.115;
   anacutsval[7][11]=0.115;
   anacutsval[7][12]=0.115;
   anacutsval[7][13]=0.115;
   anacutsval[7][14]=0.115;


   anacutsval[8][0]=0.5;
   anacutsval[8][1]=0.5;
   anacutsval[8][2]=1.0;
   anacutsval[8][3]=0.5;
   anacutsval[8][4]=0.5;
   anacutsval[8][5]=0.5;
   anacutsval[8][6]=0.5;
   anacutsval[8][7]=0.5;
   anacutsval[8][8]=0.5;
   anacutsval[8][9]=0.5;
   anacutsval[8][10]=0.5;
   anacutsval[8][11]=0.5;
   anacutsval[8][12]=0.5;
   anacutsval[8][13]=0.5;
   anacutsval[8][14]=0.5;


   anacutsval[9][0]=0.979;
   anacutsval[9][1]=0.979;
   anacutsval[9][2]= 0.9975;//0.99 ; 0.9975;
   anacutsval[9][3]= 0.995;   //0.99;  //0.995;
   anacutsval[9][4]= 0.995;   //0.99;  //0.995;
   anacutsval[9][5]=0.9975; //0.99;
   anacutsval[9][6]=0.9975; 
   anacutsval[9][7]=0.9975;
   anacutsval[9][8]=0.9975;
   anacutsval[9][9]=0.9975;
   anacutsval[9][10]=0.9975;
   anacutsval[9][11]=0.9975;
   anacutsval[9][12]=0.9975;
   anacutsval[9][13]=0.9975;
   anacutsval[9][14]=0.9975;

   anacutsval[10][0]=0.0055;
   anacutsval[10][1]=0.0055;
   anacutsval[10][2]=0.0028;//0.00400
   anacutsval[10][3]=0.000883;
   anacutsval[10][4]=0.000883;
   anacutsval[10][5]=0.000883;
   anacutsval[10][6]=0.000883;
   anacutsval[10][7]=0.000883;
   anacutsval[10][8]=0.000883;
   anacutsval[10][9]=0.000883;
   anacutsval[10][10]=0.000883;
   anacutsval[10][11]=0.000883;
   anacutsval[10][12]=0.000883;
   anacutsval[10][13]=0.000883;
   anacutsval[10][14]=0.000883;
   

 
//Production cuts

  AliRDHFCutsDplustoKpipi *prodcuts = new AliRDHFCutsDplustoKpipi();
  prodcuts->SetPtBins(nptbins+1,ptbins);
  prodcuts->SetCuts(nvars,nptbins,prodcutsval);

  //Analysis cuts
  AliRDHFCutsDplustoKpipi *analysiscuts = new AliRDHFCutsDplustoKpipi();
  analysiscuts->SetPtBins(nptbins+1,ptbins);
  analysiscuts->SetCuts(nvars,nptbins,anacutsval);
  analysiscuts->SetUsePID(kTRUE);
  
  //  analysiscuts->SetTPCPID(kTRUE);
  //analysiscuts->SetTOFPID(kTRUE);

  // Aanalysis task
  AliAnalysisTaskSEDplus *dplusTask = new AliAnalysisTaskSEDplus("DplusAnalysis",analysiscuts,prodcuts,storeNtuple);
  dplusTask->SetReadMC(readMC);
  dplusTask->SetDoLikeSign(kTRUE);
  //  dplusTask->SetUseTPCpid(kTRUE);
  //dplusTask->SetUseTOFpid(kTRUE);
  dplusTask->SetDebugLevel(0);
   dplusTask->SetMassLimits(0.2);
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
