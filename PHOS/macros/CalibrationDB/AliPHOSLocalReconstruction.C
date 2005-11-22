AliPHOSLocalReconstruction()
{
  //Run PHOS clusterization using information from calibration database.
  // Author: Boris Polishchuk (Boris.Polichtchouk at cern.ch)

  // Open local calibration data base
  AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");
 
  // Make clusterization, calibration parameters will be taken from CDB
  AliPHOSClusterizerv1 *clu = new AliPHOSClusterizerv1("galice.root");
  clu->SetEventRange(0,-1);
  clu->Exec("deb tim");
}
