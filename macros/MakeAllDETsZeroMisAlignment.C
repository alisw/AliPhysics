void MakeAllDETsZeroMisAlignment(Char_t* CDBstorage = "local://$HOME/Zero"){
  // Make zero misalignment objects for all detectors
  // Pass different "CDBstorage" argument if needed (e.g. to fill
  // conditions' data base on alien) or set it to null string to have
  // the objects saved locally on file 
  //
  TString strStorage(CDBstorage);
  if(strStorage.IsNull()){
    gSystem->Setenv("TOCDB","kFALSE");
  }else{  
    gSystem->Setenv("TOCDB","kTRUE");
    gSystem->Setenv("STORAGE",strStorage.Data());
    gSystem->Setenv("ARVERSION","v4-05-08");
  }

  if(!AliGeomManager::GetGeometry()){
    if(!(AliCDBManager::Instance())->IsDefaultStorageSet())
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
      AliCDBManager::Instance()->SetRun(0);
    AliGeomManager::LoadGeometry();
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
