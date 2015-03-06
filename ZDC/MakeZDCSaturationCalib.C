MakeZDCSaturationCalib(){
   // Create OCDB object for ZDC saturation calibration
   const char* macroname = "MakeZDCSaturationCalib.C";
 
   AliZDCSaturationCalib *satCalib = new AliZDCSaturationCalib();
   Float_t satcalibZNA[4] = {0.,0.,0.,0.};
   Float_t satcalibZNC[4] = {0.,0.,0.,0.};
   satCalib->SetZNASatCalib(satcalibZNA);
   satCalib->SetZNCSatCalib(satcalibZNC);
  
   if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
     // save in file
     const char* filename = "ZDCSaturationCalib.root";
     TFile f(filename, "RECREATE");
     if(!f){
       Error(macroname,"cannot open file for output\n");
       return;
     }
     Info(macroname,"Saving alignment objects to the file %s", filename);
     f.cd();
     f.WriteObject(satCalib,"ZDCSaturation","kSingleKey");
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
     md->SetObjectClassName("AliZDCSaturationCalib");
     AliCDBId id("ZDC/Calib/SaturationCalib",0,AliCDBRunRange::Infinity());
     storage->Put(satCalib,id,md);
   }
   
   satCalib->Delete(); 

}
