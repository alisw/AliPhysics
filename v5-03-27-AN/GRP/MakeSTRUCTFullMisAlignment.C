void MakeSTRUCTFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for all STRUCTures
  // (presently this includes only FRAME)
  // Full alignment objects for FRAME segments are built by
  // averaging in each supermodule the values produced by the
  // ANSYS finite-elements simulation
  //
  const char* macroname = "MakeSTRUCTFullMisAlignment.C";
  TClonesArray *array = new TClonesArray("AliAlignObjParams",20);

  const char* baseSymName ="FRAME/Sector"; //base of symbolic name corresponding to base of path "ALIC_1/B077_1/BSEGMO";
  TString symname;

  // the following hardcoded values have been obtained by taking the output of the ANSYS simulation
  // (finite elements simulation of deformations for given loads) and recalculating the avarege
  // displacement of the center of each space-frame sector
  Double_t dx[18]={0.13375,0.25125,0.2325,0.17,0.1475,0.12625,0.06375,0.0475,0.1775,0.32375,0.35125,0.285,0.205,0.1775,0.1525,0.07625,0.0075,0.01375};
  Double_t dy[18]={-0.0275,0.00125,-0.03625,-0.14375,-0.2,-0.13375,-0.02125,0.015,-0.0175,-0.03125,-0.0325,-0.105,-0.24375,-0.32,-0.25375,-0.1225,-0.04875,-0.04125};

  for(Int_t sm=0; sm<18; sm++){
    symname = baseSymName;
    symname += sm;
    cout<<symname.Data()<<endl;
    ((*array)[sm]) = new AliAlignObjParams(symname.Data(),0,dx[sm],dy[sm],0.,0.,0.,0.,kTRUE);
  }

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

