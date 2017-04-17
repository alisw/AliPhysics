///
/// \file AliEMCALLocalReconstruction.C
/// \ingroup EMCAL_CalibDB
/// \brief Simple macro to test EMCAL Reconstruction reading local OCDB
///
/// Run EMCAL clusterization using information from calibration database.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

AliEMCALLocalReconstruction()
{
  // Open local calibration data base
  AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");
  
  // Make clusterization, calibration parameters will be taken from CDB
  AliReconstruction rec ; 
  rec.SetRunLocalReconstruction("EMCAL") ;
  rec.Run() ;
}
