void ConfigOCDB(Int_t run=0) {
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(run);
  // magnetic field
  if ( !TGeoGlobalMagField::Instance()->GetField() ) {
    printf("Loading field map...\n");
    AliGRPManager grpMan;
    if( !grpMan.ReadGRPEntry() ) {
      printf("Cannot get GRP entry\n");
    }
    if( !grpMan.SetMagField() ) {
      printf("Problem with magnetic field setup\n");
    }
  }
  if ( !TGeoGlobalMagField::Instance()->GetField()){
    AliMagF::BMap_t smag = AliMagF::k5kG;
    Double_t bzfac = 1;
    AliMagF* magF= new AliMagF("Maps","Maps", bzfac, 1., smag);
    TGeoGlobalMagField::Instance()->SetField(magF);
  }
  //TPC calib
  AliTPCcalibDB *db=AliTPCcalibDB::Instance();
  db->SetRun(run);
  //geometry
  AliGeomManager::LoadGeometry();
  // init geometry in parameters
  if (db->GetParameters()) db->GetParameters()->ReadGeoMatrices();
}


