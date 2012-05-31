void MakeEMCALPF(){
  const char* macroname = "MakeEMCALPF.C";
  TFile *f2 = new TFile("peakfindervectors2.root",  "read" ); 

  //  AliCaloPeakFinderVectors *pfv =  (AliCaloPeakFinderVectors* )f2->GetKey( "AliCaloPeakFinderVectors"); 

  AliCaloPeakFinderVectors pfv =  *((AliCaloPeakFinderVectors* )f2->GetKey( "AliCaloPeakFinderVectors")); 

  f2->Close();
  {
    //    TString Storage = "local://home/perthi/aliroot-current/OCDB/";
    TString Storage = "local://OCDB/";  

    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname ,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    
    Info(macroname,"Saving PF objects in CDB storage %s",
	 Storage.Data());
    AliCDBManager* cdb = AliCDBManager::Instance();


    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }

    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Per Thomas Hille");
    md->SetComment("Peak-Finder vectors for EMCAL");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("EMCAL/Config/PeakFinder",0,AliCDBRunRange::Infinity());
    
    // if(pfv == 0)
    //   {
    // 	cout << "  ERROR !!!!!!!!!" << endl;
    //   }
 
    // else
    
      {
	storage->Put( &pfv,id,md);
      }

    //  delete md;
  }

  
  //
}

