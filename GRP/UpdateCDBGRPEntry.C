void UpdateCDBGRPEntry() {
  // produce the GRP default entry in CDB
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  AliCDBId id("GRP/GRP/Data",0,AliCDBRunRange::Infinity());
  AliCDBMetaData *metadata= new AliCDBMetaData();

  // Get root version
  const char* rootv = gROOT->GetVersion();

  // Get AliRoot version from file to set it in the metadata of the entry
  TFile *fv= TFile::Open("$ALICE_ROOT/CVS/Repository?filetype=raw","READ");
  Int_t size = fv->GetSize();
  char *buf = new Char_t[size];
  memset(buf, '\0', size);
  fv->Seek(0);
  char* alirootv;
  if ( fv->ReadBuffer(buf, size) ) {
    Printf("Error reading AliRoot version from file to buffer!");
    alirootv = "";
  }
  if(buf=="AliRoot"){
    alirootv="HEAD";
  }else{
    alirootv = buf;
    metadata->SetResponsible("Panos.Christakoglou@cern.ch");
    metadata->SetComment("Default values for the GRP monitored parameters");
    metadata->SetAliRootVersion(alirootv);
  }
  
  TList *list = GetGRPList();
  Printf("Storing in CDB the default values for the GRP %d parameters produced with root %s and AliRoot version %s",list->GetEntries(),rootv,alirootv);

  man->Put(list,id,metadata);
}

//_______________________________________//
TList *GetGRPList() {
  TList *list = new TList();

  TMap *mapDAQ1 = new TMap();
  mapDAQ1->Add(new TObjString("fAliceStartTime"),new TObjString("0"));
  list->Add(mapDAQ1);
  TMap *mapDAQ2 = new TMap();
  mapDAQ2->Add(new TObjString("fAliceStopTime"),new TObjString("9999"));
  list->Add(mapDAQ2);
  TMap *mapDAQ3 = new TMap();
  mapDAQ3->Add(new TObjString("fAliceBeamEnergy"),new TObjString("14"));
  list->Add(mapDAQ3);
  TMap *mapDAQ4 = new TMap();
  mapDAQ4->Add(new TObjString("fAliceBeamType"),new TObjString("pp"));
  list->Add(mapDAQ4);
  TMap *mapDAQ5 = new TMap();
  mapDAQ5->Add(new TObjString("fNumberOfDetectors"),new TObjString("15"));
  list->Add(mapDAQ5);
  TMap *mapDAQ6 = new TMap();
  mapDAQ6->Add(new TObjString("fDetectorMask"),new TObjString("1048575"));
  list->Add(mapDAQ6);
  TMap *mapDAQ7 = new TMap();
  mapDAQ7->Add(new TObjString("fLHCPeriod"),new TObjString("LHC07a"));
  list->Add(mapDAQ7);

  TMap *mapDCS1 = new TMap();
  mapDCS1->Add(new TObjString("fLHCState"),new TObjString("test"));
  list->Add(mapDCS1);
  TMap *mapDCS2 = new TMap();
  mapDCS2->Add(new TObjString("fLHCCondition"),new TObjString("test"));
  list->Add(mapDCS2);
  TMap *mapDCS3 = new TMap();
  mapDCS3->Add(new TObjString("fLHCLuminosity"),new TObjString("0"));
  list->Add(mapDCS3);
  TMap *mapDCS4 = new TMap();
  mapDCS4->Add(new TObjString("fBeamIntensity"),new TObjString("0"));
  list->Add(mapDCS4);
   TMap *mapDCS5 = new TMap();
  mapDCS5->Add(new TObjString("fL3Current"),new TObjString("0"));
  list->Add(mapDCS5);
  TMap *mapDCS6 = new TMap();
  mapDCS6->Add(new TObjString("fL3Polarity"),new TObjString("0"));
  list->Add(mapDCS6);
  TMap *mapDCS7 = new TMap();
  mapDCS7->Add(new TObjString("fDipoleCurrent"),new TObjString("0"));
  list->Add(mapDCS7);
  TMap *mapDCS8 = new TMap();
  mapDCS8->Add(new TObjString("fDipolePolarity"),new TObjString("0"));
  list->Add(mapDCS8);
  TMap *mapDCS9 = new TMap();
  mapDCS9->Add(new TObjString("fCavernTemperature"),new TObjString("0"));
  list->Add(mapDCS9);
  TMap *mapDCS10 = new TMap();
  mapDCS10->Add(new TObjString("fCavernPressure"),new TObjString("0"));
  list->Add(mapDCS10);

  return list;
}
