MakeZDCEnergyCalib(){
   // Create OCDB object for ZDC pedestal
   const char* macroname = "MakeZDCEnergyCalib.C";
 
   AliZDCEnCalib *enCalib = new AliZDCEnCalib();
   enCalib->SetEnCalib(0, 1.);
   enCalib->SetEnCalib(1, 1.);
   enCalib->SetEnCalib(2, 1.);
   enCalib->SetEnCalib(3, 1.);
   enCalib->SetEnCalib(4, 1.);
   enCalib->SetEnCalib(5, 1.);
  
   if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
     // save in file
     const char* filename = "ZDCEnergyCalib.root";
     TFile f(filename, "RECREATE");
     if(!f){
       Error(macroname,"cannot open file for output\n");
       return;
     }
     Info(macroname,"Saving alignment objects to the file %s", filename);
     f.cd();
     f.WriteObject(enCalib,"ZDCPedestals","kSingleKey");
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
     md->SetComment("Calibration objects for ZDC mapping");
     md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
     md->SetObjectClassName("AliZDCPedestals");
     AliCDBId id("ZDC/Calib/EnergyCalib",0,AliCDBRunRange::Infinity());
     storage->Put(enCalib,id,md);
   }
   
   enCalib->Delete(); 

}
