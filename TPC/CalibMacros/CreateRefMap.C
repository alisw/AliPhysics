CreateRefMap(){
  //cdb storage
  TString storage="local:///data/Work/data/2009/myOCDB";
  //set reference run numbers for the various ocdb entries
  Int_t pedestalRun=86876;
  Int_t noiseRun=86876;
  Int_t pulserRun=83680;
  Int_t ceRun=83680;
  Int_t altroRun=83680;
  Int_t qaRun=83680;
  Int_t rawRun=83680;
  //comment why the update was done
  TString comment("Update of Pedestal Referenc: 3FECs exchanged.");
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
  metaData->SetAliRootVersion("5-24-00"); //root version
  metaData->SetComment("Map for reference run numbers");
  //store object
  AliCDBId id1("TPC/Calib/Ref", first, AliCDBRunRange::Infinity());
  //
  gStorage = AliCDBManager::Instance()->GetStorage(storage.Data());
  gStorage->Put(&map, id1, metaData);
}