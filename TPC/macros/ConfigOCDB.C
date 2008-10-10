//
// Macro to be invoked before analysis 
// Setup TPC OCDB entries
//
// To be used on the proof
//


void ConfigOCDB(){
  //
  // 
  //
  // import geometry
  //
  printf("SETUP OCBD for PROOF\n");
  TGeoManager::Import("/u/miranov/proof/geometry.root");
  AliGeomManager::LoadGeometry("/u/miranov/proof/geometry.root");
  //
  //
  // Setup magnetic field
  //
  AliMagF* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field,0);
  //
  //
  //
  AliCDBManager::Instance()->SetRun(1);
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Parameters","local://$ALICE_ROOT/TPC/Calib/Parameters");

  AliTPCClusterParam * param = AliTPCcalibDB::Instance()->GetClusterParam();
  AliTPCClusterParam::SetInstance(param);
  AliTPCcalibDB::Instance()->SetExBField(0);
  //
  //
  //
  printf("END of SETUP OCBD for PROOF\n");


}
