void AddPHOSCalibration(){

   //Fills PHOS re-calibration parameters into OADB
   //Each run-dependent object contains list of 3 objects:
   //calibration parameters for pass1, pass2 and pass3 reconstruction.
   //"$ALICE_ROOT/OADB/PHOS/PHOSRecalibration.root"

  //Init existing calibrations
  AliOADBContainer calibContainer("phosRecalibration");
  calibContainer.InitFromFile("$ALICE_ROOT/OADB/PHOS/PHOSCalibrations.root","phosRecalibration");

  printf("So far we stored the following list of calibrations\n") ;
  Int_t n=calibContainer.GetNumberOfEntries() ;
  printf("-------------------------------------------------------\n") ;
  for(Int_t i=0;i<n;i++){
    TObjArray* a= (TObjArray*)calibContainer.GetObjectByIndex(i);
    printf("Entry(%d): %s,  runs %d-%d \n",i,a->GetName(),calibContainer.LowerLimit(i),calibContainer.UpperLimit(i)) ;
    printf("      pass 1: ("); 
    if(a->At(0)) printf("OK"); else printf("no") ;
    printf(")     pass 2: ("); 
    if(a->At(1)) printf("OK"); else printf("no") ;
    printf(")     pass 3: ("); 
    if(a->At(2)) printf("OK"); else printf("no") ;
    printf(")\n"); 
  }
  printf("--------------------------------------------------------\n") ;

  
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetRun(140000) ;
  man->SetDefaultStorage("local://OCDB");
  AliPHOSCalibData* phosCalibData = new AliPHOSCalibData(-1);

/*  
  // -- LHC10c --
  TObjArray * lhc10cAll = new TObjArray(3); 
  lhc10cAll->SetName("PHOSRecalibration_LHC10c");
  lhc10cAll->AddAt(phosCalibData,1) ; //pass 2 reconstruction
  lhc10cAll->AddAt(phosCalibData,2) ; //pass 2 reconstruction
  calibContainer.AppendObject(lhc10cAll,118503,121040) ;

  // -- LHC10d --
  TObjArray * lhc10dAll = new TObjArray(3); 
  lhc10dAll->SetName("PHOSRecalibration_LHC10d");
  lhc10dAll->AddAt(phosCalibData,1) ; //pass 2 reconstruction
  lhc10dAll->AddAt(phosCalibData,2) ; //pass 2 reconstruction
  calibContainer.AppendObject(lhc10dAll,122195,126437) ;

  // -- LHC10e --
  TObjArray * lhc10eAll = new TObjArray(3); 
  lhc10eAll->SetName("PHOSRecalibration_LHC10e");
  lhc10eAll->AddAt(phosCalibData,0) ; //pass 2 reconstruction
  lhc10eAll->AddAt(phosCalibData,1) ; //pass 2 reconstruction
  calibContainer.AppendObject(lhc10eAll,127712,130850) ;

  // -- LHC10f --
  TObjArray * lhc10fAll = new TObjArray(3); 
  lhc10fAll->SetName("PHOSRecalibration_LHC10f");
  lhc10fAll->AddAt(phosCalibData,0) ; //pass 2 reconstruction
  lhc10fAll->AddAt(phosCalibData,1) ; //pass 2 reconstruction
  calibContainer.AppendObject(lhc10eAll,133004,135029) ;
 
  // -- LHC10g --
  TObjArray * lhc10gAll = new TObjArray(3); 
  lhc10gAll->SetName("PHOSRecalibration_LHC10g");
  lhc10gAll->AddAt(phosCalibData,0) ; //pass 2 reconstruction
  lhc10gAll->AddAt(phosCalibData,1) ; //pass 2 reconstruction
  calibContainer.AppendObject(lhc10gAll,135654,136377) ;
*/
  // -- LHC10g --
  TObjArray * lhc10hAll = new TObjArray(3); 
  lhc10hAll->SetName("PHOSRecalibration_LHC10h");
  lhc10hAll->AddAt(phosCalibData,1) ; //pass 2 reconstruction
  calibContainer.AppendObject(lhc10hAll,136851,139517) ;
  
  calibContainer.WriteToFile("PHOSCalibrations.root");


}