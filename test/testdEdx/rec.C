const char * recoStorage="local:///cvmfs/alice.gsi.de/alice/data/2010/OCDB";
Int_t run=0;


void ModifyRecoParam(TObjArray* recoArray, Bool_t useIonTail, Double_t crossTalkCorrection){
  //
  // Modify reco param - and store it in the OCDB in local directory
  //
  AliCDBManager * man  =  AliCDBManager::Instance();
  for (Int_t i=0; i<4; i++){
    AliTPCRecoParam* p = ( AliTPCRecoParam*)recoArray->At(i);
    p->SetUseIonTailCorrection(useIonTail);
    p->SetCrosstalkCorrection(crossTalkCorrection);
  }
  TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDBrec";
  AliCDBStorage*pocdbStorage = AliCDBManager::Instance()->GetStorage(localStorage.Data());  
  AliCDBMetaData *metaData= new AliCDBMetaData();
  AliCDBId*   id1=new AliCDBId("TPC/Calib/RecoParam/", man->GetRun(), man->GetRun());
  pocdbStorage->Put(recoArray, (*id1), metaData);
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/RecoParam/",localStorage.Data());
}


void rec(Bool_t useIonTail, Double_t crossTalkCorrection) {
  //
  // run reconstruction
  // Parameters:
  //   useIonTail - switch for using ion tail correction - OCDB entry will be overritten in  
  //
  // stard reco setting
  //
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("TPC/Calib/RecoParam/",recoStorage);
  man->SetRun(run);
  AliReconstruction reco;
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));
  reco.SetRunPlaneEff(kTRUE);
  reco.SetRunQA(":"); 
  reco.SetRunGlobalQA(kFALSE);
  reco.ResetCheckRecoCDBvsSimuCDB();

  //
  //Switch Iontail in RecoParam. Creation of a new OCDB entry
  //
  AliCDBEntry* e = man->Get("TPC/Calib/RecoParam/",run); // get default
  // modify content
  TObjArray* recoArray = (TObjArray*)e->GetObject();
  ModifyRecoParam(recoArray, useIonTail, crossTalkCorrection);
  //
  //
  // Run reconstruction
  //
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
  //
  // Print the OCDB setup which we used
  //
  AliCDBEntry* ocdbEntry = man->Get("TPC/Calib/RecoParam/",run);
  TObjArray* recoArray = (TObjArray*)ocdbEntry->GetObject();
  for (Int_t i=0; i<4; i++){
    AliTPCRecoParam* recoParam = ( AliTPCRecoParam*)recoArray->At(i);
    printf("ipar=%d\t%d\t%f\n",i,recoParam->GetUseIonTailCorrection(), recoParam->GetCrosstalkCorrection());
  } 
}
