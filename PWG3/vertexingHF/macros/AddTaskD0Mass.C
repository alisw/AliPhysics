AliAnalysisTaskSED0Mass *AddTaskD0Mass(Int_t flag=0/*0 = D0,1 = LS*/,Bool_t readMC=kFALSE,Bool_t cutOnDistr=kFALSE)
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
  TString cutobjname="mycuts";
  cutobjname+=flag;

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(4); // default is 5
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
 // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);

  //cout<<"nvars = "<<RDHFD0toKpi->GetNVars()<<endl;
  const Int_t nvars=9;
  //nvars=RDHFD0toKpi->GetNVars();  
  //cout<<"Nvars = "<<nvars<<"\t"<<RDHFD0toKpi->GetNVars()<<endl;
  RDHFD0toKpi->SetName(cutobjname);
  RDHFD0toKpi->SetTitle(cutobjname);

  const Int_t nptbins=5;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=1.;
  ptbins[2]=2.;
  ptbins[3]=3.;
  ptbins[4]=5.;
  ptbins[5]=10.;
  
  RDHFD0toKpi->SetPtBins(nptbins+1,ptbins);
  

  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }
  //cout<<"\tnptbins = "<<nptbins<<endl;
  /*
  //setting PPR cut values
  rdcutsvalPPR[0][0]=0.7;
  rdcutsvalPPR[1][0]=0.04;
  rdcutsvalPPR[2][0]=0.8;
  rdcutsvalPPR[3][0]=0.5;
  rdcutsvalPPR[4][0]=0.5;
  rdcutsvalPPR[5][0]=0.05;
  rdcutsvalPPR[6][0]=0.05;
  rdcutsvalPPR[7][0]=-0.0002;
  rdcutsvalPPR[8][0]=0.5;

  rdcutsvalPPR[0][1]=rdcutsvalPPR[0][2]=0.7;
  rdcutsvalPPR[1][1]=rdcutsvalPPR[1][2]=0.02;
  rdcutsvalPPR[2][1]=rdcutsvalPPR[2][2]=0.8;
  rdcutsvalPPR[3][1]=rdcutsvalPPR[3][2]=0.7;
  rdcutsvalPPR[4][1]=rdcutsvalPPR[4][2]=0.7;
  rdcutsvalPPR[5][1]=rdcutsvalPPR[5][2]=0.05;
  rdcutsvalPPR[6][1]=rdcutsvalPPR[6][2]=0.05;
  rdcutsvalPPR[7][1]=rdcutsvalPPR[7][2]=-0.0002;
  rdcutsvalPPR[8][1]=rdcutsvalPPR[8][2]=0.6;

  rdcutsvalPPR[0][3]=0.7;
  rdcutsvalPPR[1][3]=0.02;
  rdcutsvalPPR[2][3]=0.8;
  rdcutsvalPPR[3][3]=0.7;
  rdcutsvalPPR[4][3]=0.7;
  rdcutsvalPPR[5][3]=0.05;
  rdcutsvalPPR[6][3]=0.05;
  rdcutsvalPPR[7][3]=-0.0001;
  rdcutsvalPPR[8][3]=0.8;

  rdcutsvalPPR[0][4]=0.7;
  rdcutsvalPPR[1][4]=0.02;
  rdcutsvalPPR[2][4]=0.8;
  rdcutsvalPPR[3][4]=0.7;
  rdcutsvalPPR[4][4]=0.7;
  rdcutsvalPPR[5][4]=0.05;
  rdcutsvalPPR[6][4]=0.05;
  rdcutsvalPPR[7][4]=-0.00005;
  rdcutsvalPPR[8][4]=0.8;
  */
  //setting my cut values

  rdcutsvalmine[0][0]=0.7;
  rdcutsvalmine[1][0]=0.04;
  rdcutsvalmine[2][0]=0.8;
  rdcutsvalmine[3][0]=0.5;
  rdcutsvalmine[4][0]=0.5;
  rdcutsvalmine[5][0]=0.05;
  rdcutsvalmine[6][0]=0.05;
  rdcutsvalmine[7][0]=-0.00025;
  rdcutsvalmine[8][0]=0.7;

  rdcutsvalmine[0][1]=rdcutsvalmine[0][2]=0.7;
  rdcutsvalmine[1][1]=rdcutsvalmine[1][2]=0.02;
  rdcutsvalmine[2][1]=rdcutsvalmine[2][2]=0.8;
  rdcutsvalmine[3][1]=rdcutsvalmine[3][2]=0.7;
  rdcutsvalmine[4][1]=rdcutsvalmine[4][2]=0.7;
  rdcutsvalmine[5][1]=rdcutsvalmine[5][2]=1.;
  rdcutsvalmine[6][1]=rdcutsvalmine[6][2]=1.;
  rdcutsvalmine[7][1]=rdcutsvalmine[7][2]=-0.00025;
  rdcutsvalmine[8][1]=rdcutsvalmine[8][2]=0.8;

  rdcutsvalmine[0][3]=0.7;
  rdcutsvalmine[1][3]=0.02;
  rdcutsvalmine[2][3]=0.8;
  rdcutsvalmine[3][3]=0.7;
  rdcutsvalmine[4][3]=0.7;
  rdcutsvalmine[5][3]=0.05;
  rdcutsvalmine[6][3]=0.05;
  rdcutsvalmine[7][3]=-0.00015;
  rdcutsvalmine[8][3]=0.8;

  rdcutsvalmine[0][4]=0.7;
  rdcutsvalmine[1][4]=0.02;
  rdcutsvalmine[2][4]=0.8;
  rdcutsvalmine[3][4]=0.7;
  rdcutsvalmine[4][4]=0.7;
  rdcutsvalmine[5][4]=0.05;
  rdcutsvalmine[6][4]=0.05;
  rdcutsvalmine[7][4]=-0.00015;
  rdcutsvalmine[8][4]=0.9;

  cout<<"Filled array ("<<nvars<<","<<nptbins<<")"<<endl;
  /*
  for(Int_t j=0;j<nvars;j++){
    for(Int_t k=0;k<nptbins;k++){
      cout<<rdcutsvalmine[j][k]<<"\t";
    }
    cout<<endl;
  }
  */

  //cout<<"\tbefore SetCuts : npt = "<<RDHFD0toKpi->GetNPtBins()<<endl;
  RDHFD0toKpi->SetCuts(nvars,nptbins,rdcutsvalmine);
  //  RDHFD0toKpi->PrintAll();

  // Aanalysis task    
  TString taskname="MassAndDistrAnalysis";
  if (flag==0)taskname.Prepend("D0");
  else taskname.Prepend("LS");
  AliAnalysisTaskSED0Mass *massD0Task = new AliAnalysisTaskSED0Mass(taskname.Data(),RDHFD0toKpi);
  massD0Task->SetDebugLevel(0);
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
  //AliAnalysisDataContainer *coutputmassD05 = mgr->CreateContainer(out5name,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kParamContainer, filename.Data()); //cuts
AliAnalysisDataContainer *coutputmassD05 = mgr->CreateContainer(out5name,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  //AliAnalysisDataContainer *coutputmassD05 = mgr->CreateContainer(out5name,TList::Class(),AliAnalysisManager::kParamContainer, filename.Data()); //cuts
  //AliAnalysisDataContainer *coutputmassD05 = mgr->CreateContainer(out5name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
   

  mgr->ConnectInput(massD0Task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massD0Task,1,coutputmassD01);
  mgr->ConnectOutput(massD0Task,2,coutputmassD02);
  mgr->ConnectOutput(massD0Task,3,coutputmassD03);
  mgr->ConnectOutput(massD0Task,4,coutputmassD04);
  mgr->ConnectOutput(massD0Task,5,coutputmassD05);


  return massD0Task;
}
