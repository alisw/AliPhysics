void merge(Int_t run=117048)
{
  AliCDBManager *cdbm = AliCDBManager::Instance();
  cdbm->SetRun(run);
  cdbm->SetDefaultStorage("local://$ALICE_ROOT/OCDB");     
  cdbm->SetSpecificStorage("GRP/GRP/Data",Form("local://%s","$ALICE_ROOT/test/rawmerge"));  
  //cdbm->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));  

  AliGeomManager::LoadGeometry("geometry.root");
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));

  AliSimulation sim;
  sim.SetSeed(12345);
  sim.SetRunGeneration(0);
  sim.SetRunSimulation(0);
  sim.SetRunHLT("");
  sim.SetNumberOfEvents(1);
  sim.SetMakeSDigits("");
  sim.SetMakeDigits("ITS TPC");
  sim.SetMakeDigitsFromHits("");
  sim.SetGAliceFile("galice.root");
  sim.SetGeometryFile("geometry.root");
  sim.MergeWith("$ALICE_ROOT/test/rawmerge/bkg/galice.root",1);
  sim.Run();
  return;
}
