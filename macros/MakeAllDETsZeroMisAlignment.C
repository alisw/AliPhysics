void MakeAllDETsZeroMisAlignment(Char_t* CDBstorage = "local://$HOME/Zero"){
  // Make zero misalignment objects for all detectors
  // Pass different "CDBstorage" argument if needed (e.g. to fill
  // conditions' data base on alien) or set it to null string to have
  // the objects saved locally on file 
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

  TString dets="EMCAL,FMD,HMPID,ITS,MUON,PMD,PHOS,T0,TRD,TPC,TOF,VZERO,ZDC";
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
    exec_det_macro+="ZeroMisAlignment.C";
    
    gSystem->Exec(exec_det_macro.Data());
  }

  return;
}
