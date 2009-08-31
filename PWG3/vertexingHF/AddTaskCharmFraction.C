AliAnalysisTaskCharmFraction* AddTaskCharmFraction(
         const char* fileout="d0D0.root",
	 Bool_t sideband=kFALSE,
	 Bool_t setD0usecuts=kTRUE,
	 Bool_t setcheckMC=kTRUE,
	 Bool_t setcheckMC_prompt=kTRUE,
	 Bool_t setcheckMC_fromB=kFALSE,
	 Bool_t setcheckMC_D0=kTRUE,
	 Bool_t setcheckMC_2prongs=kTRUE,
	 Bool_t setSkipD0star=kFALSE,
	 Bool_t setStudyPureBack=kFALSE)
{  
  //
  // Configuration macro for the task to analyze the fraction of prompt charm
  // using the D0 impact parameter
  // andrea.rossi@ts.infn.it
  //
  //==========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCharmFraction", "No analysis manager to connect to.");
    return NULL;
  }   

  TString str=fileout,containername;
  str.ReplaceAll(".root","");
  str.Prepend("_");

  AliAnalysisTaskCharmFraction *hfTask;
  if(!sideband) {
    hfTask = new AliAnalysisTaskCharmFraction("CharmFraction",10);
  } else {
    hfTask= new AliAnalysisTaskCharmFraction("CharmFractionSideB",10);
    hfTask->SetSideBands(-2.);
  }

  hfTask->SetUseCuts(setD0usecuts);
  hfTask->SetCheckMC(setcheckMC);
  hfTask->SetCheckMC_D0(setcheckMC_D0);
  hfTask->SetCheckMC_2prongs(setcheckMC_2prongs);
  hfTask->SetCheckMC_prompt(setcheckMC_prompt);
  hfTask->SetCheckMC_fromB(setcheckMC_fromB);
  hfTask->SetCheckMC_fromDstar(setSkipD0star);
  hfTask->SetStudyPureBackground(setStudyPureBack);
  //  hfTask->SetSideBands(0);
  //  hfTask->SetDebugLevel(2);
  mgr->AddTask(hfTask);
 
  //Now the same for sidebands
  /*AliAnalysisTaskCharmFraction *hfTaskSideB 

  
  hfTaskSideB->SetUseCuts(fD0usecuts);
  hfTaskSideB->SetCheckMC(fcheckMC);
  hfTaskSideB->SetCheckMC_D0(fcheckMC_D0);
  hfTaskSideB->SetCheckMC_2prongs(fcheckMC_2prongs);
  hfTaskSideB->SetCheckMC_prompt(fcheckMC_prompt);
  hfTaskSideB->SetCheckMC_fromB(fcheckMC_fromB);
  hfTaskSideB->SetCheckMC_fromDstar(fSkipD0star);
  hfTaskSideB->SetStudyPureBackground(fStudyPureBack);

  //  hfTaskSideB->SetDebugLevel(2);
  mgr->AddTask(hfTaskSideB);
  */
 

  // Create containers for input/output
  AliAnalysisDataContainer *cinput =   mgr->GetCommonInputContainer();//mgr->CreateContainer("cinput",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(hfTask,0,cinput);
  //  mgr->ConnectInput(hfTaskSideB,0,cinput);

  //Now container for general properties histograms
  containername="coutputCptd0d0";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputCptd0d0 = mgr->CreateContainer(containername.Data(),TH2::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,0,coutputCptd0d0);

  containername="coutputSecVtxXY";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSecVtxXY = mgr->CreateContainer(containername.Data(),TH2::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,1,coutputSecVtxXY);


  containername="coutputd0d0";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputd0d0 = mgr->CreateContainer(containername.Data(),TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,2,coutputd0d0);

  containername="coutputCpt";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputCpt = mgr->CreateContainer(containername.Data(),TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,3,coutputCpt);

  containername="coutputSecVtxZ";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSecVtxZ = mgr->CreateContainer(containername.Data(),TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,4,coutputSecVtxZ);

  containername="coutputSecVtxX";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSecVtxX = mgr->CreateContainer(containername.Data(),TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,5,coutputSecVtxX);

  containername="coutputSecVtxY";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputSecVtxY = mgr->CreateContainer(containername.Data(),TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,6,coutputSecVtxY);

  containername="coutputSecVtxPhi";
  containername.Append(str.Data());
 AliAnalysisDataContainer *coutputSecVtxPhi = mgr->CreateContainer(containername.Data(),TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);
  mgr->ConnectOutput(hfTask,7,coutputSecVtxPhi);


  //Now container for d0D0  
  AliAnalysisDataContainer **coutput=new AliAnalysisDataContainer*[10];
  containername="coutputAll";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputAll = mgr->CreateContainer(containername.Data(),TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);


  TString name="coutput";
  for(Int_t j=0;j<10;j++){
    containername=name;
    containername+=j;
    containername.Append(str.Data());
    coutput[j] = mgr->CreateContainer(containername.Data(),TH1::Class(),
				      AliAnalysisManager::kOutputContainer, 
				      fileout);
    
    mgr->ConnectOutput(hfTask,j+8,coutput[j]);
  }
  mgr->ConnectOutput(hfTask,18,coutputAll);
  //Now container for MC d0D0  
  AliAnalysisDataContainer **coutputMC=new AliAnalysisDataContainer*[10];
  containername="coutputAllMC";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputAllMC = mgr->CreateContainer(containername.Data(),TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);


  name="coutputMC";
  for(Int_t j=0;j<10;j++){
    containername=name;
    containername+=j;
    containername.Append(str.Data());
    coutputMC[j] = mgr->CreateContainer(containername.Data(),TH1::Class(),
				      AliAnalysisManager::kOutputContainer, 
				      fileout);
    
    mgr->ConnectOutput(hfTask,j+19,coutputMC[j]);
  }
  mgr->ConnectOutput(hfTask,29,coutputAllMC);

  //Now container for histo with d0 with respect to True Vtx
  AliAnalysisDataContainer **coutputd0VtxTrue=new AliAnalysisDataContainer*[10];
  containername="coutputd0VtxTrueAll";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputd0VtxTrueAll = mgr->CreateContainer(containername.Data(),TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout);


  name="coutputd0VtxTrue";
  for(Int_t j=0;j<10;j++){
    containername=name;
    containername+=j;
    containername.Append(str.Data());
    coutputd0VtxTrue[j] = mgr->CreateContainer(containername.Data(),TH1::Class(),
				      AliAnalysisManager::kOutputContainer, 
				      fileout);
    
    mgr->ConnectOutput(hfTask,j+30,coutputd0VtxTrue[j]);
  }
  mgr->ConnectOutput(hfTask,40,coutputd0VtxTrueAll);
  //INV MASS
  containername="coutputD0InvMassAll";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputD0InvMassAll = mgr->CreateContainer(containername.Data(),TH1::Class(),
								       AliAnalysisManager::kOutputContainer, 
								       fileout);
  mgr->ConnectOutput(hfTask,41,coutputD0InvMassAll);
  containername="coutputD0MCInvMassAll";
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputD0MCInvMassAll = mgr->CreateContainer(containername.Data(),TH1::Class(),
									 AliAnalysisManager::kOutputContainer, 
									 fileout);
  mgr->ConnectOutput(hfTask,42,coutputD0MCInvMassAll);
  
 ////////
 //NOW THE SAME FOR SIDE BANDS
 /*
 //Now container for general properties histograms

 AliAnalysisDataContainer *coutputSBCptd0d0 = mgr->CreateContainer("coutputSBCptd0d0",TH2::Class(),
								   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
 mgr->ConnectOutput(hfTaskSideB,0,coutputSBCptd0d0);

  AliAnalysisDataContainer *coutputSBSecVtxXY = mgr->CreateContainer("coutputSBSecVtxXY",TH2::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
  mgr->ConnectOutput(hfTaskSideB,1,coutputSBSecVtxXY);


  AliAnalysisDataContainer *coutputSBd0d0 = mgr->CreateContainer("coutputSBd0d0",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
  mgr->ConnectOutput(hfTaskSideB,2,coutputSBd0d0);

  AliAnalysisDataContainer *coutputSBCpt = mgr->CreateContainer("coutputSBCpt",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
  mgr->ConnectOutput(hfTaskSideB,3,coutputSBCpt);

  AliAnalysisDataContainer *coutputSBSecVtxZ = mgr->CreateContainer("coutputSBSecVtxZ",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
  mgr->ConnectOutput(hfTaskSideB,4,coutputSBSecVtxZ);

  AliAnalysisDataContainer *coutputSBSecVtxX = mgr->CreateContainer("coutputSBSecVtxX",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
  mgr->ConnectOutput(hfTaskSideB,5,coutputSBSecVtxX);

  AliAnalysisDataContainer *coutputSBSecVtxY = mgr->CreateContainer("coutputSBSecVtxY",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
  mgr->ConnectOutput(hfTaskSideB,6,coutputSBSecVtxY);

 AliAnalysisDataContainer *coutputSBSecVtxPhi = mgr->CreateContainer("coutputSBSecVtxPhi",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
  mgr->ConnectOutput(hfTaskSideB,7,coutputSBSecVtxPhi);


  //Now container for d0D0SideB  
  AliAnalysisDataContainer **coutputSB=new AliAnalysisDataContainer*[10];
  AliAnalysisDataContainer *coutputSBAll = mgr->CreateContainer("coutputSBAll",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");


  TString name="coutputSB",strname;
  for(Int_t j=0;j<10;j++){
    strname=name;
    strname+=j;
    coutputSB[j] = mgr->CreateContainer(strname.Data(),TH1::Class(),
				      AliAnalysisManager::kOutputContainer, 
				      "d0D0SideB.root");
    
    mgr->ConnectOutput(hfTaskSideB,j+8,coutputSB[j]);
  }
  mgr->ConnectOutput(hfTaskSideB,18,coutputSBAll);
  //Now container for MC d0D0SideB  
  AliAnalysisDataContainer **coutputSBMC=new AliAnalysisDataContainer*[10];
  AliAnalysisDataContainer *coutputSBAllMC = mgr->CreateContainer("coutputSBAllMC",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");


  name="coutputSBMC";
  for(Int_t j=0;j<10;j++){
    strname=name;
    strname+=j;
    coutputSBMC[j] = mgr->CreateContainer(strname.Data(),TH1::Class(),
				      AliAnalysisManager::kOutputContainer, 
				      "d0D0SideB.root");
    
    mgr->ConnectOutput(hfTaskSideB,j+19,coutputSBMC[j]);
  }
  mgr->ConnectOutput(hfTaskSideB,29,coutputSBAllMC);

  //Now container for histo with d0 with respect to True Vtx
  AliAnalysisDataContainer **coutputSBd0VtxTrue=new AliAnalysisDataContainer*[10];
  AliAnalysisDataContainer *coutputSBd0VtxTrueAll = mgr->CreateContainer("coutputSBd0VtxTrueAll",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");


  name="coutputSBd0VtxTrue";
  for(Int_t j=0;j<10;j++){
    strname=name;
    strname+=j;
    coutputSBd0VtxTrue[j] = mgr->CreateContainer(strname.Data(),TH1::Class(),
				      AliAnalysisManager::kOutputContainer, 
				      "d0D0SideB.root");
    
    mgr->ConnectOutput(hfTaskSideB,j+30,coutputSBd0VtxTrue[j]);
  }
  mgr->ConnectOutput(hfTaskSideB,40,coutputSBd0VtxTrueAll);

//INV MASS
 AliAnalysisDataContainer *coutputSBD0InvMassAll = mgr->CreateContainer("coutputSBD0InvMassAll",TH1::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
 mgr->ConnectOutput(hfTaskSideB,41,coutputSBD0InvMassAll);

 AliAnalysisDataContainer *coutputSBD0MCInvMassAll = mgr->CreateContainer("coutputSBD0MCInvMassAll",TH1::Class(),
								      AliAnalysisManager::kOutputContainer, 
							   "d0D0SideB.root");
 mgr->ConnectOutput(hfTaskSideB,42,coutputSBD0MCInvMassAll);

 */

  return hfTask;
}
