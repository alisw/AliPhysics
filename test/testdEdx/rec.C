void rec() {
  AliReconstruction reco;

  // Switch Iontail in RecoParam. Ceration of a new OCDB entry
  //AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("local:///cvmfs/alice.gsi.de/alice/data/2010/OCDB");
  //AliCDBEntry* e = man->Get("TPC/Calib/RecoParam/",136844);
  //TObjArray* a = (TObjArray*)e->GetObject();
  //for (Int_t i=0; i<4; i++){
  // AliTPCRecoParam* p = ( AliTPCRecoParam*)a->At(i);
  // p->SetUseIonTailCorrection(kTRUE);
  //}


  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));
  reco.SetRunPlaneEff(kTRUE);

  
  reco.SetRunQA(":"); 
  reco.SetRunGlobalQA(kFALSE);


  TStopwatch timer;
  reco.ResetCheckRecoCDBvsSimuCDB();
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
