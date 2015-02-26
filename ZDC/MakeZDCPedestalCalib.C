MakeZDCPedestalCalib(){
   // Create OCDB object for ZDC pedestal
   const char* macroname = "MakeZDCPedestalCalib.C";
 
   const Int_t nch = 24;
   Float_t MeanPed[2*nch], MeanPedWidth[2*nch], 
           MeanPedOOT[2*nch], MeanPedWidthOOT[2*nch],
           CorrCoeff0[2*nch], CorrCoeff1[2*nch];
   for(int j=0; j<2*nch; j++){
     if(j<nch){
       MeanPed[j] = MeanPedOOT[j] = 50.;
       MeanPedWidth[j] = MeanPedWidthOOT[j] = 4.;
       CorrCoeff0[j] = 0.;
       CorrCoeff1[j] = 1.;
     }
     else{
       MeanPed[j] = MeanPedOOT[j] = 8.*50.;
       MeanPedWidth[j] = MeanPedWidthOOT[j] = 8.*4..;
       CorrCoeff0[j] = 0.;
       CorrCoeff1[j] = 1.;
     }
   }
  
   AliZDCPedestals *calibPed = new AliZDCPedestals("ZDC");
   calibPed->SetMeanPed(MeanPed);
   calibPed->SetMeanPedWidth(MeanPedWidth);
   calibPed->SetOOTPed(MeanPedOOT);
   calibPed->SetOOTPedWidth(MeanPedWidthOOT);
   calibPed->SetPedCorrCoeff(CorrCoeff0, CorrCoeff1); 

   if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
     // save in file
     const char* filename = "ZDCPedestals.root";
     TFile f(filename, "RECREATE");
     if(!f){
       Error(macroname,"cannot open file for output\n");
       return;
     }
     Info(macroname,"Saving alignment objects to the file %s", filename);
     f.cd();
     f.WriteObject(calibPed,"ZDCPedestals","kSingleKey");
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
     md->SetComment("Calibration objects for ZDC pedestals");
     md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
     md->SetObjectClassName("AliZDCPedestals");
     AliCDBId id("ZDC/Calib/Pedestals",0,AliCDBRunRange::Infinity());
     storage->Put(calibPed,id,md);
   }
   
   calibPed->Delete(); 

}
