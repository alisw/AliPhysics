

//Run EMCAL clusterization using information from calibration database.
// Author: Gustavo Conesa

AliEMCALLocalReconstruction()
{

  // Open local calibration data base
  AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");
  
  // Make clusterization, calibration parameters will be taken from CDB

  AliReconstruction rec ; 
  rec.SetRunLocalReconstruction("EMCAL") ;
  rec.Run() ;
}
