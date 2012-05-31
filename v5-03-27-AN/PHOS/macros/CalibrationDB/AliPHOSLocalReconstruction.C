AliPHOSLocalReconstruction()
{
  //Run PHOS clusterization using information from calibration database.
  // Author: Boris Polishchuk (Boris.Polichtchouk at cern.ch)

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("PHOS","local://CalibDB");
 
  // Make clusterization, calibration parameters will be taken from CDB

  AliPHOSClusterizerv1 *clu = new AliPHOSClusterizerv1("galice.root");
  clu->SetEventRange(0,-1);
  clu->Exec("deb tim");
}
