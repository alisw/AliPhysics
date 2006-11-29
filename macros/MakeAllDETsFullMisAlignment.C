void MakeAllDETsFullMisAlignment(Char_t* CDBstorage = "local://$HOME/Full"){
  // Make full misalignment objects for all detectors
  // Pass different "CDBstorage" argument if needed (e.g. to fill
  // conditions' data base on alien) or set it to null string to have
  // the objects saved locally on file 
  // This macro defines the default name and place for the detector-macros
  // in charge of producing the full misalignment objects as 
  // $ALICE_ROOT/DET/MakeDETFullMisAlignment.C
  //
  TString strStorage(CDBstorage);
  if(strStorage.IsNull()){
    gSystem->Setenv("$TOCDB","kFALSE");
  }else{  
    gSystem->Setenv("$TOCDB","kTRUE");
    gSystem->Setenv("$STORAGE",strStorage.Data());
    gSystem->Setenv("$ARVERSION","v4-04-Release");
  }

  // if not already present, create geometry file needed by those detectors
  // producing their objects in the local RS
  if(gSystem->AccessPathName("./geometry.root")){
    gAlice->Init();
    gGeoManager->Export("geometry.root");
  }else{
    TGeoManager::Import("geometry.root");
  }

  TString dets = "EMCAL,FMD,ITS,MUON,PHOS,PMD,HMPID,T0,TOF,TPC,TRD,VZERO,ZDC";
//   TString dets = "ABSO,DIPO,FMD,FRAME,HALL,ITS,MAG,MUON,PHOS,PIPE,PMD,HMPID,SHIL,T0,TOF,TPC,TRD,ZDC,EMCAL,CRT,VZERO";
  TObjArray *detArray = dets.Tokenize(',');
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
