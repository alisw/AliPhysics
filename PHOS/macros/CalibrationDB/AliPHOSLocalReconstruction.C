AliPHOSLocalReconstruction()
{
  //Run PHOS clusterization using information from calibration database.

  // Open local calibration data base
  AliCDBLocal *loc = new AliCDBLocal("CalibDB");
 
  // Make clusterization, calibration parameters will be taken from CDB
  AliPHOSClusterizerv1 *clu = new AliPHOSClusterizerv1("galice.root");
  clu->SetEventRange(0,-1);
  clu->Exec("deb tim");
}
