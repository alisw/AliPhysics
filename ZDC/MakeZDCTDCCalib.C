MakeZDCTDCCalib(){
   // Create OCDB object for ZDC TDC calibration
   const char* macroname = "MakeZDCTDCCalib.C";
 
   AliZDCTDCCalib *tdcCalib = new AliZDCTDCCalib();
   
   for(Int_t j=0; j<6; j++){  
      tdcCalib->SetMeanTDC(j, 0.);
      tdcCalib->SetWidthTDC(j, 2.5);
   }
  
   if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
     // save in file
     const char* filename = "ZDCTDCCalib.root";
     TFile f(filename, "RECREATE");
     if(!f){
       Error(macroname,"cannot open file for output\n");
       return;
     }
     Info(macroname,"Saving alignment objects to the file %s", filename);
     f.cd();
     f.WriteObject(tdcCalib,"ZDCTDCCalib","kSingleKey");
     f.Close();
   }
   else{
     // save in CDB storage
     TString Storage = gSystem->Getenv("STORAGE");
     if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
       Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
       return;
     }
     Info(macroname,"Saving alignment objects in CDB storage %s", Storage.Data());
     AliCDBManager* cdb = AliCDBManager::Instance();
     AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
     if(!storage){
       Error(macroname,"Unable to open storage %s\n",Storage.Data());
       return;
     }
     AliCDBMetaData* md = new AliCDBMetaData();
     md->SetResponsible("Chiara Oppedisano");
     md->SetComment("Calibration objects for ZDC saturation calibration");
     md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
     md->SetObjectClassName("AliZDCTDCCalib");
     AliCDBId id("ZDC/Calib/TDCCalib",0,AliCDBRunRange::Infinity());
     storage->Put(tdcCalib,id,md);
   }
   
   tdcCalib->Delete(); 

}
