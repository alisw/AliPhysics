const char* GetARversion();

void MakeAllDETsFullMisAlignment(Char_t* CDBstorage = "local://$HOME/Full"){
   // Make full misalignment objects for all detectors
  // Pass different "CDBstorage" argument if needed (e.g. to fill
  // conditions' data base on alien) or set it to null string to have
  // the objects saved locally on file 
  // This macro defines the default name and place for the detector-macros
  // in charge of producing the full misalignment objects as 
  // $ALICE_ROOT/DET/MakeDETFullMisAlignment.C
  //
  const char* macroname="MakeAllDETsFullMisAlignment.C";

  TString strStorage(CDBstorage);
  if(strStorage.IsNull()){
    gSystem->Setenv("TOCDB","kFALSE");
  }else{  
    gSystem->Setenv("TOCDB","kTRUE");
    gSystem->Setenv("STORAGE",strStorage.Data());
    gSystem->Setenv("ARVERSION",GetARversion());
  }

  // Load geometry from CDB updating it if we are producing the
  // alignment objects for the CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  
  if(strStorage.IsNull()){ //if we produce the objects into a file
    AliGeomManager::LoadGeometry(); //load geom from default CDB storage
  }else{ // if we produce the objects in a CDB storage
    // update geometry in it
    Info(macroname,"Updating geometry in CDB storage %s",strStorage.Data());
    gROOT->ProcessLine(".L $ALICE_ROOT/GRP/UpdateCDBIdealGeom.C");
    UpdateCDBIdealGeom(strStorage.Data());
    // load the same geometry from given CDB storage
    AliCDBPath path("GRP","Geometry","Data");
    AliCDBStorage* storage = cdb->GetStorage(strStorage.Data());
    AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroname,"Couldn't load geometry data from CDB!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    if (!geom) Fatal(macroname,"Couldn't find TGeoManager in the specified CDB entry!");
    AliGeomManager::SetGeometry(geom);
  }
   
  // run macro for non-sensitive modules
  // (presently generates only FRAME alignment objects)
  gSystem->Exec("aliroot -b -q $ALICE_ROOT/GRP/MakeSTRUCTFullMisAlignment.C");

  // run macros for sensitive modules
  TString sModules="EMCAL,FMD,HMPID,ITS,MUON,PMD,PHOS,T0,TRD,TPC,TOF,VZERO,ZDC";
  TObjArray *detArray = sModules.Tokenize(',');
  TIter iter(detArray);
  TObjString *ostr;
  TString exec_det_macro;

  while((ostr = (TObjString*) iter.Next())){
    TString str(ostr->String());
    exec_det_macro="aliroot -b -q $ALICE_ROOT/";
    exec_det_macro+=str;
    exec_det_macro+="/Make";
    exec_det_macro+=str;
    exec_det_macro+="FullMisAlignment.C";
    
    gSystem->Exec(exec_det_macro.Data());
  }

  return;
}

const char* GetARversion(){
  // Get AliRoot version from $ALICE_ROOT/CVS/Repository file
  // It's the best we can do without a GetVersion() method
  TFile *fv= TFile::Open("$ALICE_ROOT/CVS/Repository?filetype=raw","READ");
  Int_t size = fv->GetSize();
  char *buf = new Char_t[size];
  memset(buf, '\0', size);
  fv->Seek(0);
  const char* alirootv;
  if ( fv->ReadBuffer(buf, size) ) {
    Printf("Error reading AliRoot version from file to buffer!");
    alirootv = "";
  }
  if(buf=="AliRoot"){
    alirootv="HEAD";
  }else{
    alirootv = buf;
  }
  return alirootv;
}
