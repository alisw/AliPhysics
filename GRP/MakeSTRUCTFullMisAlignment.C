void MakeSTRUCTFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for all STRUCTures
  // (presently this includes only FRAME)
  // Full alignment objects for FRAME segments are built by
  // averaging in each supermodule the values produced by the
  // ANSYS finite-elements simulation
  //
  const char* macroname = "MakeSTRUCTFullMisAlignment.C";
  TClonesArray *array = new TClonesArray("AliAlignObjParams",20);

  // the following hardcoded values have been obtained by taking the output of the ANSYS simulation
  // (finite elements simulation of deformations for given loads) and recalculating the avarege
  // displacement of the center of each space-frame sector
  ((*array)[0]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO0_1",0,0.13375,-0.0275,0,0,0,0,kTRUE);
  ((*array)[1]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO1_1",0,0.25125,0.00125,0,0,0,0,kTRUE);
  ((*array)[2]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO2_1",0,0.2325,-0.03625,0,0,0,0,kTRUE);
  ((*array)[3]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO3_1",0,0.17,-0.14375,0,0,0,0,kTRUE);
  ((*array)[4]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO4_1",0,0.1475,-0.2,0,0,0,0,kTRUE);
  ((*array)[5]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO5_1",0,0.12625,-0.13375,0,0,0,0,kTRUE);
  ((*array)[6]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO6_1",0,0.06375,-0.02125,0,0,0,0,kTRUE);
  ((*array)[7]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO7_1",0,0.0475,0.015,0,0,0,0,kTRUE);
  ((*array)[8]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO8_1",0,0.1775,-0.0175,0,0,0,0,kTRUE);
  ((*array)[9]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO9_1",0,0.32375,-0.03125,0,0,0,0,kTRUE);
  ((*array)[10]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO10_1",0,0.35125,-0.0325,0,0,0,0,kTRUE);
  ((*array)[11]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO11_1",0,0.285,-0.105,0,0,0,0,kTRUE);
  ((*array)[12]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO12_1",0,0.205,-0.24375,0,0,0,0,kTRUE);
  ((*array)[13]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO13_1",0,0.1775,-0.32,0,0,0,0,kTRUE);
  ((*array)[14]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO14_1",0,0.1525,-0.25375,0,0,0,0,kTRUE);
  ((*array)[15]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO15_1",0,0.07625,-0.1225,0,0,0,0,kTRUE);
  ((*array)[16]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO16_1",0,0.0075,-0.04875,0,0,0,0,kTRUE);
  ((*array)[17]) = new AliAlignObjParams("ALIC_1/B077_1/BSEGMO17_1",0,0.01375,-0.04125,0,0,0,0,kTRUE);

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "STRUCTfullMisalignment.root";
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
    md->SetComment("Full misalignment for FRAME, including dead weight and full weight deformations derived from ANSYS finiteelements simulation: dispmap corrected excel worksheet from Werner Riegler");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("GRP/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

