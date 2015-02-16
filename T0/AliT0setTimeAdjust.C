void AliT0setTimeAdjust(Int_t firstRun, Int_t lastRun)
{
  // Set T0 alignment OCDB entry calib/TimeAdjust
  AliT0CalibSeasonTimeShift *clb = new AliT0CalibSeasonTimeShift();
  Float_t par[4] = {0,0,0,0};
  Float_t sigmas[4] = {0,0,0,0};

  clb->SetT0Par(par, sigmas);
  clb->Print();

  AliCDBMetaData md;
  md.SetResponsible("Alla");
  TString fPath="T0/Calib/TimeAdjust";
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
    AliCDBId id(fPath.Data(),firstRun, lastRun );
    storage->Put(clb, id, &md);
  }
  
}
