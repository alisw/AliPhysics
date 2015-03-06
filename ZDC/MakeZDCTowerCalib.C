MakeZDCTowerCalib(){
   // Create OCDB object for ZDC tower calibration
   const char* macroname = "MakeZDCTowerCalib.C";
 
   AliZDCTowerCalib *towerCalib = new AliZDCTowerCalib();
   
   for(Int_t j=0; j<5; j++){  
      towerCalib->SetZN1EqualCoeff(j, 1.);
      towerCalib->SetZP1EqualCoeff(j, 1.);
      towerCalib->SetZN2EqualCoeff(j, 1.);
      towerCalib->SetZP2EqualCoeff(j, 1.);  
   }
  
   if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
     // save in file
     const char* filename = "ZDCTowerCalib.root";
     TFile f(filename, "RECREATE");
     if(!f){
       Error(macroname,"cannot open file for output\n");
       return;
     }
     Info(macroname,"Saving alignment objects to the file %s", filename);
     f.cd();
     f.WriteObject(towerCalib,"ZDCTower","kSingleKey");
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
     md->SetObjectClassName("AliZDCTowerCalib");
     AliCDBId id("ZDC/Calib/TowerCalib",0,AliCDBRunRange::Infinity());
     storage->Put(towerCalib,id,md);
   }
   
   towerCalib->Delete(); 

}
