enum SurveyDataType_t { kSurvey = 0, kDummy = 1};

void MakeEMCALFullMisAlignment(TString geoname = "EMCAL_FIRSTYEAR",TString surveyFilename = "emcal_survey_FIRSTYEAR.txt",SurveyDataType_t type = kSurvey){
  // Create TClonesArray of full misalignment objects for EMCAL
  //
  const char* macroname = "MakeEMCALFullMisAlignment.C";
  if(geoname=="")geoname=AliEMCALGeometry::GetDefaultGeometryName();
  const AliEMCALGeometry *emcalGeom = AliEMCALGeometry::GetInstance(geoname,"");
  if(!emcalGeom) {
    Error("MakeEMCALFullMisAlignment","Cannot obtain AliEMCALGeometry singleton\n");
    return;
  }
  TClonesArray *array = new TClonesArray("AliAlignObjParams",emcalGeom->GetNumberOfSuperModules());
  TClonesArray &alobj = *array;

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  
  AliCDBStorage* storage;
  
  if( TString(gSystem->Getenv("TOCDB")) == TString("kTRUE") ){
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBPath path("GRP","Geometry","Data");
    AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    AliGeomManager::SetGeometry(geom);
  }else{
    AliGeomManager::LoadGeometry(); //load geom from default CDB
				    //storage
  }    


  AliEMCALSurvey emcalSurvey(surveyFilename,type);
  emcalSurvey.CreateAliAlignObjParams(alobj);

  // *************************    2nd step    ***************

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "EMCALfullMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"EMCALAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Jennifer Klay");
    md->SetComment("Full misalignment for EMCAL_FIRSTYEAR based on survey information");
    md->AddDateToComment();
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("EMCAL/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

