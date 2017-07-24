///
/// \file AliEMCALDecalibrate.C
/// \ingroup EMCAL_CalibDB
/// \brief Decalibrate exiting digits
///
/// Digitize and decalibrate events assuming that SDigits 
/// have been already produced.
/// Decalibration coefficients are located in the local file
/// DeCalibDB/EMCAL/Calib/Data/Run_xxx.root
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
//// based on PHOS macro by  Boris Polishchuk (Boris.Polichtchouk at cern.ch)
///

/// Main method
/// \param nevents number of events to process
void AliEMCALDecalibrate(Int_t nevents=3)
{
  // Use "decalibration" database to simulate decalibrated EMCAL.
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("EMCAL","local://DeCalibDB");

  // Make digitization, calibration parameters will be taken from CDB
  AliSimulation sim ; 
  sim.SetRunGeneration(kFALSE) ;
  sim.SetMakeSDigits("") ;
  sim.SetMakeDigits("EMCAL") ;
  sim.Run(nevents) ;  
}
