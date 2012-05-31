void MakePHOSCalibration(){

   //Fills PHOS re-calibration parameters into OADB
   //Each run-dependent object contains list of 3 objects:
   //calibration parameters for pass1, pass2 and pass3 reconstruction.
   //"$ALICE_ROOT/OADB/PHOS/PHOSRecalibration.root"
   
  AliOADBContainer calibContainer("phosRecalibration");

  AliCDBManager * man = AliCDBManager::Instance();
  man->SetRun(140000) ;
  man->SetDefaultStorage("local://OCDB");
  AliPHOSCalibData* phosCalibData = new AliPHOSCalibData(-1);

  // -- LHC10h --
  TObjArray * lhc10aAll = new TObjArray(3); 
  lhc10aAll->SetName("PHOSRecalibration_LHC10b");
  lhc10aAll->AddAt(phosCalibData,2) ; //pass 3 reconstruction
  calibContainer.AppendObject(lhc10aAll,114737,117223) ;

  calibContainer.WriteToFile("PHOSCalibrations.root");


}