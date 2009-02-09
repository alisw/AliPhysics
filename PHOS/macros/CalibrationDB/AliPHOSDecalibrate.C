/* $Id$ */

// Digitize and decalibrate events assuming that SDigits 
// have been already produced.
// Author: Boris Polichtchouk (Boris.Polichtchouk@cern.ch)

void AliPHOSDecalibrate(Int_t nevents=1)
{

  //Use "decalibration" database to simulate decalibrated PHOS.

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("PHOS","local://DeCalibDB");

  // Make digitization, calibration parameters will be taken from CDB

  AliSimulation sim ; 
  sim.SetRunGeneration(kFALSE) ;
  sim.SetMakeSDigits("") ;
  sim.SetMakeDigits("PHOS") ;
  sim.Run(nevents) ;  

}
