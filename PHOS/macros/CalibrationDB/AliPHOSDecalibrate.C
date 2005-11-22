/* $Id$ */

// Digitize and decalibrate events assuming that SDigits 
// have been already produced.
// Decalibration coefficients are located in the local file
// deCalibDB/PHOS/Calib/GainFactors_and_Pedestals/Run_xxx.root
// Author: Boris Polichtchouk (Boris.Polichtchouk@cern.ch)

void AliPHOSDecalibrate(Int_t nevents=1)
{

  //Load (de)calibration database into aliroot session
  //and set it to AliPHOSGetter.
  
  AliPHOSCalibData* deCal  = (AliPHOSCalibData*)(AliCDBManager::Instance()
    ->GetStorage("local://deCalibDB")->Get("PHOS/Calib/GainFactors_and_Pedestals",1)
    ->GetObject());
  
  AliPHOSGetter* gime = AliPHOSGetter::Instance("galice.root");
  gime->SetCalibData(deCal);


  AliSimulation sim ; 
  sim.SetRunGeneration(kFALSE) ;
  sim.SetMakeSDigits("") ;
  sim.SetMakeDigits("PHOS") ;
  sim.Run(nevents) ;  
}
