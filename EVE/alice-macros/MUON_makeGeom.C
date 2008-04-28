{
  AliMpCDB::LoadMpSegmentation2(); 
  gAlice->Init("$ALICE_ROOT/MUON/Config.C");
  //gAlice->Init("Config.C");

  gGeoManager->Export("geometry.root");

}
