AliAnalysisTaskSE *AddTaskDFilterAndCorrelations(
  AliAnalysisTaskSEDmesonsFilterCJ::ECandidateType cand = AliAnalysisTaskSEDmesonsFilterCJ::kDstartoKpipi,
  TString filename = "DStartoKpipiCuts.root",
  Bool_t theMCon = kFALSE,
  Bool_t reco = kTRUE /*must be true if theMCon is false*/,
  TString suffix = "",
  TString jetArrname = "",
  TString trackArrname = "PicoTracks",
  TString rhoname="",
  Int_t leadingHadType = 0 /*0 = charged, 1 = neutral, 2 = both*/,
  Float_t R = 0.4,
  Float_t jptcut = 10.,
  const char *cutType = "TPC",
  Int_t thnsparse=1, /*-1 = no thnsparse, 0 = heavy, 1 = light*/
  Double_t percjetareacut = -1.,
  AliAnalysisTaskEmcal::TriggerType trType=AliAnalysisTaskEmcal::kND,
  Int_t typeDjet=2,
  TString subwagons=""
)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSEDmesonsFilterCJ", "No analysis manager to connect to.");
    return NULL;
  } 

  Bool_t useStdC = kFALSE;
  TFile* filecuts=TFile::Open(filename);
  if(!filecuts || (filecuts && !filecuts->IsOpen())) {
    cout<<"Input file not found: use std cuts"<<endl;
    useStdC = kTRUE;
  }

  AliRDHFCuts *analysiscuts=0x0;
  switch (cand) {
  case 0 :
    if(useStdC) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPP2010();
    } else
      analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get("D0toKpiCuts");
    break;
  case 1 :
    if(useStdC) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      analysiscuts->SetStandardCutsPP2010();
    } else
      analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get("DStartoKpipiCuts");
    analysiscuts->SetName("DStartoKpipiCuts");
    break;
  }
  
  if (!analysiscuts) { // mm let's see if everything is ok
    AliFatal("Specific AliRDHFCuts not found");
    return;
  } 

  printf("CREATE TASK\n"); //CREATE THE TASK

  TString candname="DStar"; 
  if(cand==0)  candname="D0";
  TString sR = Form("R%.0f",R*10);
  
  TString taskFiltername="DmesonsFilterCJ";
  taskFiltername+=candname;
  taskFiltername+=suffix;
  if(theMCon) taskFiltername+="MC";
  if(!reco)   taskFiltername+="gen";
  
  AliAnalysisTaskSEDmesonsFilterCJ* taskFilter = mgr->GetTask(taskFiltername.Data());
  Bool_t bTaskFilter=kTRUE;
  if (!taskFilter){
     bTaskFilter=kFALSE;
     taskFilter = new AliAnalysisTaskSEDmesonsFilterCJ(taskFiltername.Data(),analysiscuts,cand);
     if(!theMCon) reco=kTRUE;
     taskFilter->SetMC(theMCon); //D meson settings
     taskFilter->SetUseReco(reco);
     mgr->AddTask(taskFilter);
  } else Printf("Task %s already exist, continue",taskFiltername.Data());

  // create the task
  TString taskCorrName="TaskFlavourJetCorrelations";
  taskCorrName+=candname;
  taskCorrName+=suffix;
  if(theMCon) taskCorrName+="MC";
  if(!reco)   taskCorrName+="gen";
  taskCorrName+=cutType;
  taskCorrName+=Form("PTj%.0f",jptcut);
  taskCorrName+=sR;
  taskCorrName+=Form("Typ%d",typeDjet);

  AliAnalysisTaskFlavourJetCorrelations *taskCorr = new AliAnalysisTaskFlavourJetCorrelations(taskCorrName.Data(), 
     analysiscuts, cand);
  
  taskCorr->SetJetArrayName(jetArrname);
  taskCorr->SetTrackArrayName(trackArrname);
  //taskCorr->SetRadius(R);
  AliParticleContainer *trackCont  = taskCorr->AddParticleContainer(trackArrname);
  trackCont->SetClassName("AliVTrack");
  
  AliJetContainer *jetCont = taskCorr->AddJetContainer(jetArrname,cutType,R);
  if(jetCont) {
     jetCont->ConnectParticleContainer(trackCont);
     //jetCont->SetJetAcceptanceType(cutType);
     jetCont->SetJetPtCut(jptcut);
     jetCont->SetPercAreaCut(percjetareacut);
     jetCont->SetRhoName(rhoname); 
  }
  taskCorr->SetMC(theMCon);
  taskCorr->SetUseReco(reco);
  taskCorr->SetTypeDJetSelection(typeDjet);
  if(theMCon && trType!=AliAnalysisTaskEmcal::kND){
     taskCorr->SetCaloTriggerPatchInfoName("EmcalTriggers");
     taskCorr->SetTriggerTypeSel(trType);
  }
  if(thnsparse==-1)taskCorr->TurnOffTHnSparse();
  if(thnsparse==0) taskCorr->HeavyTHnSparse();
  if(thnsparse==1) taskCorr->LightTHnSparse();   
    
  mgr->AddTask(taskCorr);

  if(theMCon) {
     suffix+="MC";
     if(reco) suffix+="rec";  
  }
  
  
  // Create and connect containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  TString outputfileF = outputfile, outputfileC = outputfile;
  outputfileF += ":PWG3_D2H_DmesonsForJetCorrelations";
  outputfileC += ":PWG3_D2H_DEmcalJet";
  
  outputfileF += suffix;
  outputfileC += suffix;

  TString nameContainerF0="histograms";
  TString nameContainerF1="cuts";
  
  TString nameContainerC0="hCor";
  TString nameContainerC1="cutsJ";

  TString nameContainerFC2="Dcandidates";
  TString nameContainerFC3="DSBcandidates";
    
  TString nameContainerFC4 = "DcandidatesAndTracks";
  TString nameContainerFC5 = "DSBcandidatesAndTracks";

  nameContainerF0  += candname;
  nameContainerF1  += candname;
  nameContainerC0  += candname;
  nameContainerC1  += candname;
  nameContainerFC2 += candname;
  nameContainerFC3 += candname;
  nameContainerFC4 += candname;
  nameContainerFC5 += candname;
  
  nameContainerF0  += suffix;
  nameContainerF1  += suffix;
  nameContainerC0  += suffix;
  nameContainerC1  += suffix;
  nameContainerFC2 += suffix;
  nameContainerFC3 += suffix;
  nameContainerFC4 += suffix;
  nameContainerFC5 += suffix;
  
  nameContainerC0+=sR;
  nameContainerC1+=sR;
  if(typeDjet!=2) { //no particular name for default
     nameContainerC0+=Form("Typ%d",typeDjet);
     nameContainerC1+=Form("Typ%d",typeDjet);
  }
  
  nameContainerC0 += subwagons;
  nameContainerC1 += subwagons;
  
  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  //cinput0->SetName(Form("in%s%s",candname.Data(),suffix.Data()));
  
  // ----- output data -----
  AliAnalysisDataContainer *coutputF0;
  AliAnalysisDataContainer *coutputF1;
  AliAnalysisDataContainer *coutputFC2;
  AliAnalysisDataContainer *coutputFC3;
  AliAnalysisDataContainer *coutputFC4;
  AliAnalysisDataContainer *coutputFC5;
  
  if(!bTaskFilter){
  coutputF0 = mgr->CreateContainer(nameContainerF0, TList::Class(),AliAnalysisManager::kOutputContainer,outputfileF.Data());
  
  coutputF1 = mgr->CreateContainer(nameContainerF1, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfileF.Data());
  
  coutputFC2 = mgr->CreateContainer(nameContainerFC2, TClonesArray::Class(),AliAnalysisManager::kExchangeContainer, outputfileF.Data()); // exchange
  
  coutputFC3 = mgr->CreateContainer(nameContainerFC3, TClonesArray::Class(),AliAnalysisManager::kExchangeContainer, outputfileF.Data()); // exchange

  coutputFC4 = mgr->CreateContainer(nameContainerFC4, TClonesArray::Class(), AliAnalysisManager::kExchangeContainer, outputfile.Data()); // exchange
      
  coutputFC5 = mgr->CreateContainer(nameContainerFC5, TClonesArray::Class(), AliAnalysisManager::kExchangeContainer, outputfile.Data()); // exchange

  } else {
     TObjArray * cnt = mgr->GetContainers();
     coutputF0 = (AliAnalysisDataContainer*)cnt->FindObject(nameContainerF0);
     coutputF1 = (AliAnalysisDataContainer*)cnt->FindObject(nameContainerF1);
     coutputFC2= (AliAnalysisDataContainer*)cnt->FindObject(nameContainerFC2);
     coutputFC3= (AliAnalysisDataContainer*)cnt->FindObject(nameContainerFC3);
     coutputFC4= (AliAnalysisDataContainer*)cnt->FindObject(nameContainerFC4);
     coutputFC5= (AliAnalysisDataContainer*)cnt->FindObject(nameContainerFC5);
  }
  
  AliAnalysisDataContainer *coutputC0 = mgr->CreateContainer(nameContainerC0, TList::Class(),AliAnalysisManager::kOutputContainer,outputfileC.Data());

  AliAnalysisDataContainer *coutputC1 = mgr->CreateContainer(nameContainerC1, AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer, outputfileC.Data());
  
    
  mgr->ConnectInput(taskFilter,0,cinput0);
  mgr->ConnectInput(taskCorr,0,cinput0);
  
  mgr->ConnectOutput(taskFilter,1,coutputF0);
  mgr->ConnectOutput(taskFilter,2,coutputF1);
  mgr->ConnectOutput(taskFilter,3,coutputFC2);
  mgr->ConnectOutput(taskFilter,4,coutputFC3);
  mgr->ConnectOutput(taskFilter,5,coutputFC4);
  mgr->ConnectOutput(taskFilter,6,coutputFC5);
  
  mgr->ConnectInput(taskCorr,1,coutputFC2);
  mgr->ConnectInput(taskCorr,2,coutputFC3);
  mgr->ConnectInput(taskCorr,3,coutputFC4);
  mgr->ConnectInput(taskCorr,4,coutputFC5);
  mgr->ConnectOutput(taskCorr,1,coutputC0);
  mgr->ConnectOutput(taskCorr,2,coutputC1);
  //if(cand==1) mgr->ConnectOutput(task,4,coutput4);

  Printf("Input and Output connected to the manager");
  return taskCorr; 
}

