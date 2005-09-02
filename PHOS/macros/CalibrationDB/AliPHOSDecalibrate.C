/* $Id$ */

// Digitize and decalibrate events assuming that SDigits 
// have been already produced.
// Decalibration coefficients are located in the local file
// DeCalibDB/PHOS/Calib/GainFactors_and_Pedestals/Run0-10_v0.root

void AliPHOSDecalibrate(Int_t nevents=1)
{
  //Load calibration database into aliroot session
  //and set it to AliPHOSGetter.
  AliCDBLocal *loc = new AliCDBLocal("DeCalibDB");
 
  AliPHOSCalibData* clb = (AliPHOSCalibData*)AliCDBStorage::Instance()
    ->Get("PHOS/Calib/GainFactors_and_Pedestals",gAlice->GetRunNumber());
  
  AliPHOSGetter* gime = AliPHOSGetter::Instance("galice.root");
  gime->SetCalibData(clb);


  AliSimulation sim ; 
  sim.SetRunGeneration(kFALSE) ;
  sim.SetMakeSDigits("") ;
  sim.SetMakeDigits("PHOS") ;
  sim.Run(nevents) ;  
}
