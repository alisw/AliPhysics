/*
  Create refernce map

*/

void CreateRefMap(){
  // aliroot -b -q $ALICE_ROOT/TPC/CalibMacros/CreateRefMap.C
  //
  //cdb storage - output stored in the working directory
  //

  TString storage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  //set reference run numbers for the various ocdb entries
  Int_t pedestalRun=121642;
  Int_t noiseRun=121642;
  Int_t pulserRun=121645;
  Int_t ceRun=120818;
  Int_t altroRun=120503;
  Int_t qaRun=123537;
  Int_t rawRun=123537;
  //comment why the update was done
  TString comment("");
  //find first run for which the data are valid
  Int_t first=0;
  first=TMath::Max(first,pedestalRun);
  first=TMath::Max(first,noiseRun);
  first=TMath::Max(first,pulserRun);
  first=TMath::Max(first,ceRun);
  first=TMath::Max(first,altroRun);
  first=TMath::Max(first,qaRun);
  first=TMath::Max(first,rawRun);
//   first=0;
  //create the map
  TMap map;
  map.Add(new TObjString("TPC/Calib/Pedestals"),new TObjString(Form("%d",pedestalRun)));
  map.Add(new TObjString("TPC/Calib/PadNoise"),new TObjString(Form("%d",noiseRun)));
  map.Add(new TObjString("TPC/Calib/Pulser"),new TObjString(Form("%d",pulserRun)));
  map.Add(new TObjString("TPC/Calib/CE"),new TObjString(Form("%d",ceRun)));
  map.Add(new TObjString("TPC/Calib/AltroConfig"),new TObjString(Form("%d",altroRun)));
  map.Add(new TObjString("TPC/Calib/QA"),new TObjString(Form("%d",qaRun)));
  map.Add(new TObjString("TPC/Calib/Raw"),new TObjString(Form("%d",rawRun)));
  map.Add(new TObjString("Comment"), new TObjString(comment));
  //create meta data
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TMap");
  metaData->SetResponsible("Jens Wiechula (Jens.Wiechula@cern.ch)");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("5-26-00"); //root version
  metaData->SetComment("Map for reference run numbers");
  //store object
  AliCDBId id1("TPC/Calib/Ref", first, AliCDBRunRange::Infinity());
  //
  gStorage = AliCDBManager::Instance()->GetStorage(storage.Data());
  gStorage->Put(&map, id1, metaData);
}

