
// Digitize and decalibrate events assuming that SDigits 
// have been already produced.
// Decalibration coefficients are located in the local file
// DeCalibDB/EMCAL/Calib/GainFactors_and_Pedestals/Run_xxx.root
// Author: Boris Polichtchouk (Boris.Polichtchouk@cern.ch)
// Modified from PHOS script for EMCAL by Gustavo Conesa

void AliEMCALDecalibrate(Int_t nevents=3)
{

  //Load (de)calibration database into aliroot session
  //and set it to AliEMCALGetter.
  
  AliEMCALCalibData* deCal  = (AliEMCALCalibData*)(AliCDBManager::Instance()
    ->GetStorage("local://DeCalibDB")->Get("EMCAL/Calib/GainFactors_and_Pedestals",1)
    ->GetObject());
  
  //Loader  
  AliRunLoader* rl=0x0;
  
  rl = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "update");
  if (rl == 0x0)
    {
      gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
      return;
    }

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));
  //  rl->LoadSDigits("EMCAL");
  emcalLoader->SetCalibData(deCal);

  AliSimulation sim ; 
  sim.SetRunGeneration(kFALSE) ;
  sim.SetMakeSDigits("") ;
  sim.SetMakeDigits("EMCAL") ;
  sim.Run(nevents) ;  
}
