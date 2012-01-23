AliAnalysisTaskSECharmFraction* AddTaskSECharmFraction(TString fileout="d0D0.root",Int_t switchMC[5],Bool_t readmc=kFALSE,Bool_t usepid=kTRUE,Bool_t likesign=kFALSE,TString cutfile="D0toKpiCharmFractCuts.root",TString containerprefix="c",Int_t ppPbPb=0,Int_t analysLevel=2, Float_t minC=0., Float_t maxC=20.,Float_t minCloose=40., Float_t maxCloose=80.)
{  
  //
  // Configuration macro for the task to analyze the fraction of prompt charm
  // using the D0 impact parameter
  // andrea.rossi@ts.infn.it
  //
  //==========================================================================

  //######## !!! THE SWITCH FOR MC ANALYSIS IS NOT IMPLEMENTED YET!!! ##########à
  switchMC[0]=1;
  switchMC[1]=1;
  switchMC[2]=1;
  switchMC[3]=1;
  switchMC[4]=1;
  Int_t last=0;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCharmFraction", "No analysis manager to connect to.");
    return NULL;
  }   
  
  TString str,containername;
  if(fileout=="standard"){
    fileout=AliAnalysisManager::GetCommonFileName();
    fileout+=":PWG3_D2H_";
    fileout+="d0D0";
    if(containerprefix!="c")fileout+=containerprefix;
    str="d0D0";
  }
  else if(fileout=="standardUp"){
    fileout=AliAnalysisManager::GetCommonFileName();
    fileout+=":PWG3_D2H_Up_";
    fileout+="d0D0";
    if(containerprefix!="c")fileout+=containerprefix;
    str="d0D0";    
  }
  else {
    str=fileout;
    str.ReplaceAll(".root","");
  }
  str.Prepend("_");

  AliAnalysisTaskSECharmFraction *hfTask;
  if(!gSystem->AccessPathName(cutfile.Data(),kFileExists)){
    TFile *f=TFile::Open(cutfile.Data());
    AliRDHFCutsD0toKpi *cutTight= (AliRDHFCutsD0toKpi*)f->Get("D0toKpiCutsStandard");
    cutTight->PrintAll();
    AliRDHFCutsD0toKpi *cutLoose= (AliRDHFCutsD0toKpi*)f->Get("D0toKpiCutsLoose");
    cutLoose->PrintAll();
    hfTask = new AliAnalysisTaskSECharmFraction("AliAnalysisTaskSECharmFraction",cutTight,cutLoose);
  }
  else {
    //hfTask = new AliAnalysisTaskSECharmFraction("AliAnalysisTaskSECharmFraction");
    AliRDHFCutsD0toKpi *cutTight=new AliRDHFCutsD0toKpi("D0toKpiCutsStandard");
    AliRDHFCutsD0toKpi *cutLoose=new AliRDHFCutsD0toKpi("D0toKpiCutsLoose");
    if(ppPbPb==1){
      cutTight->SetStandardCutsPbPb2010();
      cutTight->SetMinCentrality(minC);
      cutTight->SetMaxCentrality(maxC);
      cutLoose->SetStandardCutsPbPb2010();
      cutLoose->SetMinCentrality(minCloose);
      cutLoose->SetMaxCentrality(maxCloose);
    }
    else {
      cutTight->SetStandardCutsPP2010();
      cutLoose->SetStandardCutsPP2010();
    }
    hfTask = new AliAnalysisTaskSECharmFraction("AliAnalysisTaskSECharmFraction",cutTight,cutLoose);  
    cutLoose->PrintAll();
  }
  
  if(ppPbPb==1){// Switch Off recalctulation of primary vertex w/o candidate's daughters
    // a protection that must be kept here to be sure 
    //that this is done also if the cut objects are provided by outside 
    Printf("AddTaskSECharmFraction: Switch Off recalculation of primary vertex w/o candidate's daughters (PbPb analysis) \n");
    AliRDHFCutsD0toKpi *cloose=hfTask->GetLooseCut();
    AliRDHFCutsD0toKpi *ctight=hfTask->GetTightCut();
    cloose->SetRemoveDaughtersFromPrim(kFALSE);
    ctight->SetRemoveDaughtersFromPrim(kFALSE);
    if(analysLevel<2){
      printf("Cannot activate the filling of all the histograms for PbPb analysis \n changing analysis level to 2 \n");
      analysLevel=2;
    }
    // Activate Default PID for proton rejection (TEMPORARY)
    // cloose->SetUseDefaultPID(kTRUE);
    // ctight->SetUseDefaultPID(kTRUE);
  }

  hfTask->SetReadMC(readmc);
  hfTask->SetNMaxTrForVtx(2);
  hfTask->SetAnalyzeLikeSign(likesign);
  hfTask->SetUsePID(usepid);
  hfTask->SetStandardMassSelection();
  hfTask->SetAnalysisLevel(analysLevel);

  // hfTask->SignalInvMassCut(0.27);

  /*  ############### HERE THE POSSIBILITY TO SWITCH ON/OFF THE TLISTS AND MC SELECTION WILL BE SET #########à

  hfTask->SetUseCuts(setD0usecuts);
  hfTask->SetCheckMC(setcheckMC);
  hfTask->SetCheckMC_D0(setcheckMC_D0);
  hfTask->SetCheckMC_2prongs(setcheckMC_2prongs);
  hfTask->SetCheckMC_prompt(setcheckMC_prompt);
  hfTask->SetCheckMC_fromB(setcheckMC_fromB);
  hfTask->SetCheckMC_fromDstar(setSkipD0star);
  hfTask->SetStudyPureBackground(setStudyPureBack);*/
  //  hfTask->SetSideBands(0);
  //  hfTask->SetDebugLevel(2);
  mgr->AddTask(hfTask);
 
 
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput =   mgr->GetCommonInputContainer();
  //mgr->CreateContainer("cinput",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(hfTask,0,cinput);
  

  //Now container for general properties histograms
  containername="outputNentries";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputNentries = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,1,coutputNentries);

  containername="outputSignalType";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSignalType = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,2,coutputSignalType);


  containername="outputSignalType_LsCuts";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSignalType_LsCuts = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,3,coutputSignalType_LsCuts);


  containername="outputSignalType_TghCuts";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSignalType_TghCuts = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,4,coutputSignalType_TghCuts);


  containername="outputNormalizationCounter";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputNormCounter = mgr ->CreateContainer(containername.Data(), AliNormalizationCounter::Class(), 
								       AliAnalysisManager::kOutputContainer, 
								       fileout.Data());
  mgr->ConnectOutput(hfTask, 5, coutputNormCounter);
  
  //Now Container for MC TList
  containername="listMCproperties";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistMCprop = mgr->CreateContainer(containername.Data(),TList::Class(),
							       AliAnalysisManager::kOutputContainer, 
							       fileout.Data());
  mgr->ConnectOutput(hfTask,6,clistMCprop);
  
  // Now container for TLists 
  last=7;
  //##########  NO CUTS TLISTS CONTAINER ##############à
  containername="listNCsign";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCsign = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistNCsign);
  last++;


  containername="listNCback";  
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCback = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistNCback);
  last++;

  containername="listNCfromB";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCfromB = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistNCfromB);
  last++;


  containername="listNCfromDstar";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCfromDstar = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistNCfromDstar);
  last++;


  containername="listNCother";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistNCother = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistNCother);
  last++;


  //######### LOOSE CUTS TLISTS CONTAINER #############
  containername="listLSCsign";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCsign = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCsign);
  last++;


  containername="listLSCback";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCback = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCback);
  last++;

  containername="listLSCfromB";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCfromB = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCfromB);
  last++;


  containername="listLSCfromDstar";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCfromDstar = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCfromDstar);
  last++;


  containername="listLSCother";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistLSCother = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistLSCother);
  last++;



  //######### TIGHT CUTS TLISTS CONTAINER #############
    containername="listTGHCsign";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCsign = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCsign);
  last++;


  containername="listTGHCback";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCback = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCback);
  last++;

  containername="listTGHCfromB";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCfromB = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCfromB);
  last++;


  containername="listTGHCfromDstar";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCfromDstar = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCfromDstar);
  last++;


  containername="listTGHCother";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTGHCother = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  mgr->ConnectOutput(hfTask,last,clistTGHCother);
  last++;
  
  // Container for Cuts Objects
  containername="cutsObjectTight";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *cCutsObjectTight = mgr->CreateContainer(containername,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer,fileout.Data()); //cuts
  mgr->ConnectOutput(hfTask,last,cCutsObjectTight);
  last++;
  
  containername="cutsObjectLoose";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *cCutsObjectLoose = mgr->CreateContainer(containername,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer,fileout.Data()); //cuts
  mgr->ConnectOutput(hfTask,last,cCutsObjectLoose);

  return hfTask;
}
