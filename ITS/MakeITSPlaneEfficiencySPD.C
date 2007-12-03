void MakeITSPlaneEfficiencySPD(Int_t firstRun=0,Int_t lastRun=AliCDBRunRange::Infinity()){
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  }
  
  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSPlaneEff");
  md1->SetResponsible("Giuseppe Bruno");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("head 19/11/07"); //root version

  AliCDBId idplaneeffSPD("ITS/PlaneEff/PlaneEffSPD",firstRun, lastRun);
  
  AliITSPlaneEffSPD* planeeffSPD = new AliITSPlaneEffSPD();

//  planeeffSPD->SetOwner(kFALSE);
  

  // loop over SPD chip
  Bool_t BFound=kFALSE;
  //for(Int_t key=0;key<1200;key++){
  for(UInt_t mod=0;mod<240;mod++){
  for(UInt_t chip=0;chip<5;chip++){
  // suppose to have 1000 tracks in each chip and an average efficiency of 99%
    for(Int_t j=0; j<1000; j++) {
      BFound=kFALSE;
      if (gRandom->Uniform(0,1000)>10) BFound=kTRUE;
      planeeffSPD->UpDatePlaneEff(BFound,mod,chip);
    }
  }}
  if(AliCDBManager::Instance()->GetDefaultStorage()->Put(planeeffSPD, idplaneeffSPD, md1))
  {cout << "Riuscito " << endl;}
}