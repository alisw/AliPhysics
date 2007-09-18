void MakeSTRUCTResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for all STRUCTures
  // (presently this includes only FRAME)
  // 
  const char* macroname = "MakeSTRUCTResMisAlignment.C";

  TClonesArray *array = new TClonesArray("AliAlignObjParams",20);

  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t dvoluid = AliGeomManager::LayerToVolUID(iLayer,iIndex); //dummy vol id

  const char* basepath ="ALIC_1/B077_1/BSEGMO";
  TString segmpath;

  for(Int_t sm=0; sm<18; sm++){
    segmpath=basepath;
    segmpath+=sm;
    segmpath+="_1";
    cout<<segmpath.Data()<<endl;
    new((*array)[sm]) AliAlignObjParams(segmpath.Data(),dvoluid,0.,0.,0.,0.,0.,0.,kTRUE);
  }


  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "STRUCTresMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects in %s", filename);
    f.cd();
    f.WriteObject(array,"STRUCTAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    Info(macroname,"Saving alignment objects in CDB storage %s",Storage.Data());
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Grosso Raffaele");
    md->SetComment("Residual misalignment for STRUCT. Presently includes objects for FRAME.");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("GRP/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

