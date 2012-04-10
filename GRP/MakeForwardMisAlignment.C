void MakeForwardMisAlignment(Int_t inputRun){
  // Creates misalignment of the beam pipe, frame and
  // FMD2&3.
  // The misalignment is based on the ITS global alignment
  // object with the exception of the beam pipe where
  // only the translations are used ignoring the rotation
  // matrix.
  const char* macroname = "MakeForwardMisAlignment.C";

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(inputRun);
  AliCDBEntry *inputEntry = cdb->Get("ITS/Align/Data");
  TClonesArray *inputArray = (TClonesArray*)inputEntry->GetObject();
  AliAlignObjParams *itsAlign = (AliAlignObjParams*)inputArray->At(0);
  itsAlign->Print();
  printf("\n\n\n");
  Double_t inputTrans[3];
  itsAlign->GetTranslation(inputTrans);
  Double_t inputAngles[3];
  itsAlign->GetAngles(inputAngles);

  TClonesArray *array = new TClonesArray("AliAlignObjParams",20);

  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t dvoluid = AliGeomManager::LayerToVolUID(iLayer,iIndex); //dummy vol id

  const char* baseSymName ="FRAME/Sector"; //base of symbolic name corresponding to base of path "ALIC_1/B077_1/BSEGMO";
  TString symname;

  for(Int_t sm=0; sm<18; sm++){
    symname = baseSymName;
    symname += sm;
    cout<<symname.Data()<<endl;
    new((*array)[sm]) AliAlignObjParams(symname.Data(),dvoluid,0.,0.,0.,0.,0.,0.,kTRUE);
  }

  new((*array)[18]) AliAlignObjParams("CP1",dvoluid,inputTrans[0],inputTrans[1],0.,0.,0.,0.,kTRUE);
  array->At(18)->Print();
  printf("\n\n\n");

  // FMD
  TClonesArray *arrayFMD = new TClonesArray("AliAlignObjParams",4);
  new((*arrayFMD)[0]) AliAlignObjParams("FMD/FMD2_T",0,inputTrans[0],inputTrans[1],inputTrans[2],inputAngles[0],inputAngles[1],inputAngles[2],kTRUE);
  new((*arrayFMD)[1]) AliAlignObjParams("FMD/FMD2_B",0,inputTrans[0],inputTrans[1],inputTrans[2],inputAngles[0],inputAngles[1],inputAngles[2],kTRUE);
  new((*arrayFMD)[2]) AliAlignObjParams("FMD/FMD3_T",0,inputTrans[0],inputTrans[1],inputTrans[2],inputAngles[0],inputAngles[1],inputAngles[2],kTRUE);
  new((*arrayFMD)[3]) AliAlignObjParams("FMD/FMD3_B",0,inputTrans[0],inputTrans[1],inputTrans[2],inputAngles[0],inputAngles[1],inputAngles[2],kTRUE);
  arrayFMD->Print();
  printf("\n\n\n");

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "ForwardMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects in %s", filename);
    f.cd();
    f.WriteObject(array,"STRUCTAlignObjs","kSingleKey");
    f.WriteObject(arrayFMD,"FMDAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    Info(macroname,"Saving alignment objects in CDB storage %s",Storage.Data());
    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Grosso Raffaele");
    md->SetComment("Misalignment for FRAME and beam pipe");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("GRP/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
    AliCDBMetaData* mdFMD = new AliCDBMetaData();
    mdFMD->SetResponsible("FMD Experts");
    mdFMD->SetComment("Misalignment for FMD using ITS object");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId idFMD("FMD/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(arrayFMD,idFMD,mdFMD);
  }

  array->Delete();

}

