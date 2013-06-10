AliAnalysisTaskEmcalJetTriggerQA* AddTaskEmcalJetTriggerQA(TString kTracksName = "PicoTracks", 
							   TString kClusName = "caloClusterCorr",
							   Double_t R = 0.4, 
							   Double_t ptminTrack = 0.15, 
							   Double_t etminClus = 0.3, 
							   Int_t rhoType = 0,
							   UInt_t type = AliAnalysisTaskEmcal::kEMCAL,
							   TString trigClass = "",
							   TString kEmcalCellsName = "") {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalJetTriggerQA","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalJetTriggerQA", "This task requires an input event handler");
      return NULL;
    }


  
 TString strJetsFull = "";
 TString strJetsCh = "";
 TString wagonName = Form("DiJetR%03d%spT%04d%sET%04dRhoType%d",(int)(100*R),kTracksName.Data(),(int)(1000.*ptminTrack),kClusName.Data(),(int)(1000.*etminClus),rhoType);

 if(kClusName.IsNull()) {  //particle level jets
   strJetsFull = Form("Jet_AKTFullR%03d_%s_pT%04d",(int)(100*R),kTracksName.Data(),(int)(1000.*ptminTrack));
   strJetsCh = Form("Jet_AKTChargedR%03d_%s_pT%04d",(int)(100*R),kTracksName.Data(),(int)(1000.*ptminTrack));
   wagonName = Form("JetTriggerQAR%03d%spT%04dRhoType%dTrigClass%s",(int)(100*R),kTracksName.Data(),(int)(1000.*ptminTrack),rhoType,trigClass.Data());
 }
 else if(kTracksName.IsNull()) { //neutral/calo jets
   strJetsFull = Form("Jet_AKTNeutralR%03d_%s_pT%04d_%s_ET%04d",(int)(100*R),"PicoTracks",(int)(1000.*ptminTrack),kClusName.Data(),(int)(1000.*etminClus));
   strJetsCh = Form("Jet_AKTFullR%03d_%s_pT%04d_%s_ET%04d",(int)(100*R),"PicoTracks",(int)(1000.*ptminTrack),"caloClustersCorr",(int)(1000.*etminClus));
   wagonName = Form("JetTriggerQANeutralR%03d%sET%04dRhoType%dTrigClass%s",(int)(100*R),kClusName.Data(),(int)(1000.*etminClus),rhoType,trigClass.Data());
 }
 else { //full jets
   strJetsFull = Form("Jet_AKTFullR%03d_%s_pT%04d_%s_ET%04d",(int)(100*R),kTracksName.Data(),(int)(1000.*ptminTrack),kClusName.Data(),(int)(1000.*etminClus));
   strJetsCh = Form("Jet_AKTChargedR%03d_%s_pT%04d_%s_ET%04d",(int)(100*R),kTracksName.Data(),(int)(1000.*ptminTrack),"caloClustersCorr",(int)(1000.*etminClus));
   wagonName = Form("JetTriggerQAR%03d%spT%04d%sET%04dRhoType%dTrigClass%s",(int)(100*R),kTracksName.Data(),(int)(1000.*ptminTrack),kClusName.Data(),(int)(1000.*etminClus),rhoType,trigClass.Data());
 }

  AliAnalysisTaskEmcalJetTriggerQA *task = new AliAnalysisTaskEmcalJetTriggerQA(wagonName);
  task->SetTracksName(kTracksName.Data());
  task->SetClusName(kClusName.Data());
  task->SetJetsName(strJetsFull.Data());
  task->SetJetsChName(strJetsCh.Data());
  task->SetJetRadius(R);
  task->SetJetPtCut(0.15);
  task->SetPercAreaCut(0.557);
  task->SetTrackPtCut(ptminTrack);
  task->SetClusPtCut(etminClus);
  task->SetAnaType(type);
  task->SetTriggerClass(trigClass.Data());
  task->SetCaloCellsName(kEmcalCellsName.Data());

  if(rhoType==1) {
    task->SetRhoName("RhoSparse_Scaled");
    task->SetRhoChName("RhoSparse");
  }

  task->SelectCollisionCandidates();

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  AliAnalysisDataContainer *coutput1 = 0x0;

  TString containerName1 = Form("%s",wagonName.Data());

  TString outputfile = Form("%s:%s",AliAnalysisManager::GetCommonFileName(),wagonName.Data());

  coutput1 = mgr->CreateContainer(containerName1, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);

  return task;  

}
