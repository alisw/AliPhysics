MakeZDCChMapCalib(){
   // Create OCDB object for ZDC mapping
   const char* macroname = "MakeZDCChMapCalib.C";
 
   const int kNch=48; //const int kNScch=32;
   Int_t readMap[kNch][6]; 

   FILE *file;
   if((file = fopen("./ShuttleInput/ZDCChMapping.dat","r")) != NULL){
     for(Int_t j=0; j<kNch; j++){       
       for(Int_t k=0; k<6; k++){
   	 int read = fscanf(file,"%d",&readMap[j][k]);
    	 if(read == 0) printf(" Failing in reading data from mapping file\n");
       }
     }
     fclose(file);
   }

   AliZDCChMap *mapCalib = new AliZDCChMap();
   for(Int_t k=0; k<kNch; k++){
      mapCalib->SetADCModule(k,readMap[k][1]);
      mapCalib->SetADCChannel(k,readMap[k][2]);
      mapCalib->SetDetector(k,readMap[k][4]);
      mapCalib->SetSector(k,readMap[k][5]);
    }
  
   if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
     // save in file
     const char* filename = "ZDCChMap.root";
     TFile f(filename, "RECREATE");
     if(!f){
       Error(macroname,"cannot open file for output\n");
       return;
     }
     Info(macroname,"Saving alignment objects to the file %s", filename);
     f.cd();
     f.WriteObject(mapCalib,"ZDCChMap","kSingleKey");
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
     md->SetObjectClassName("AliZDCChMap");
     AliCDBId id("ZDC/Calib/ChMap",0,AliCDBRunRange::Infinity());
     storage->Put(mapCalib,id,md);
   }
   
   mapCalib->Delete(); 

}
