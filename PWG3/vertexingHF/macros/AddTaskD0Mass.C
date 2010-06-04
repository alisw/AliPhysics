AliAnalysisTaskSED0Mass *AddTaskD0Mass(TString finname="D0toKpiCuts.root",Int_t flag=0/*0 = D0,1 = LS*/,Bool_t readMC=kFALSE,Bool_t cutOnDistr=kFALSE)
{
  //
  // AddTask for the AliAnalysisTaskSE for D0 candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD) and cut variables distributions
  // C.Bianchin  chiara.bianchin@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
    return NULL;
  }   

  TString filename="",out1name="",out2name="",out3name="",out4name="",out5name="",inname="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":PWG3_D2H_";
  if(flag==0){
    filename+="D0InvMass";
    if(cutOnDistr) filename+="C"; 
    //list mass
    out1name="coutputmassD0Mass";
    if(cutOnDistr) out1name+="C"; 
    //list distr
    out2name="coutputmassD0distr";
    if(cutOnDistr) out2name+="C"; 
    //hist entries
    out3name="nEntriesD0";
    if(cutOnDistr) out3name+="C"; 
    //list checks
    out4name="checksD0";
    if(cutOnDistr) out4name+="C"; 
    //cuts object
    out5name="cutsD0";
    if(cutOnDistr) out5name+="C"; 

    inname="cinputmassD0_0";
    if(cutOnDistr) inname+="C"; 

  } else {
    filename+="D0InvMassLikeSign";
    if(cutOnDistr) filename+="C"; 
    //list mass
    out1name="coutputmassLSMass";
    if(cutOnDistr) out1name+="C"; 
    //list distr
    out2name="coutputmassLSdistr";
    if(cutOnDistr) out2name+="C"; 
    //hist entries
    out3name="nEntriesLS";
    if(cutOnDistr) out3name+="C"; 
    //list checks
    out4name="checksLS";
    if(cutOnDistr) out4name+="C"; 
    //cuts object
    out5name="cutsLS";
    if(cutOnDistr) out5name+="C"; 

    inname="cinputmassD0_1";
    if(cutOnDistr) inname+="C"; 
  }

   //setting my cut values

    //cuts order
    //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
    //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
    //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
    //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);


  TFile* filecuts=new TFile(finname.Data());
  if(!filecuts->IsOpen()){
    cout<<"Input file not found: exit"<<endl;
    return;
  }

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi = (AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiCuts");
  RDHFD0toKpi->SetName(Form("D0toKpiCuts%d",flag));

  if(!RDHFD0toKpi){
    cout<<"Specific AliRDHFCuts not found"<<endl;
    return;
  }

  // Aanalysis task    
  TString taskname="MassAndDistrAnalysis";
  if (flag==0)taskname.Prepend("D0");
  else taskname.Prepend("LS");
  AliAnalysisTaskSED0Mass *massD0Task = new AliAnalysisTaskSED0Mass(taskname.Data(),RDHFD0toKpi);
  massD0Task->SetDebugLevel(2);
  massD0Task->SetArray(flag);
  massD0Task->SetReadMC(readMC);
  massD0Task->SetCutOnDistr(cutOnDistr);
  mgr->AddTask(massD0Task);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputmassD0 = mgr->CreateContainer(inname,TChain::Class(), 
							  AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputmassD01 = mgr->CreateContainer(out1name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //mass
  AliAnalysisDataContainer *coutputmassD02 = mgr->CreateContainer(out2name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //distr
  AliAnalysisDataContainer *coutputmassD03 = mgr->CreateContainer(out3name,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //nev
  AliAnalysisDataContainer *coutputmassD04 = mgr->CreateContainer(out4name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //check
  AliAnalysisDataContainer *coutputmassD05 = mgr->CreateContainer(out5name,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  
  mgr->ConnectInput(massD0Task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massD0Task,1,coutputmassD01);
  mgr->ConnectOutput(massD0Task,2,coutputmassD02);
  mgr->ConnectOutput(massD0Task,3,coutputmassD03);
  mgr->ConnectOutput(massD0Task,4,coutputmassD04);
  mgr->ConnectOutput(massD0Task,5,coutputmassD05);


  return massD0Task;
}
