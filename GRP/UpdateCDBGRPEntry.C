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
  
  TList *list = new TList();
  TMap *mappp = GetGRPList("pp");
  list->Add(mappp);
  TMap *mappbpb = GetGRPList("PbPb");
  list->Add(mappbpb);

  Printf("Storing in CDB the default values for the GRP %d parameters produced with root %s and AliRoot version %s",list->GetEntries(),rootv,alirootv);

  man->Put(list,id,metadata);
}

//_______________________________________//
TMap *GetGRPList(const char* system) {
  TString fSystem = system;
  TMap *map = new TMap();
  map->SetName(system);

  //DAQ
  map->Add(new TObjString("fAliceStartTime"),new TObjString("0"));
  map->Add(new TObjString("fAliceStopTime"),new TObjString("9999"));
  if(fSystem == "pp")
    map->Add(new TObjString("fAliceBeamEnergy"),new TObjString("14000"));
  else map->Add(new TObjString("fAliceBeamEnergy"),new TObjString("5500"));
  map->Add(new TObjString("fAliceBeamType"),new TObjString(system));
  map->Add(new TObjString("fNumberOfDetectors"),new TObjString("15"));
  map->Add(new TObjString("fDetectorMask"),new TObjString("1048575"));
  map->Add(new TObjString("fLHCPeriod"),new TObjString("LHC07a"));

  //DCS
  map->Add(new TObjString("fLHCState"),new TObjString("test"));
  map->Add(new TObjString("fLHCCondition"),new TObjString("test"));
  map->Add(new TObjString("fLHCLuminosity"),new TObjString("0"));
  map->Add(new TObjString("fBeamIntensity"),new TObjString("0"));
  map->Add(new TObjString("fL3Current"),new TObjString("0"));
  map->Add(new TObjString("fL3Polarity"),new TObjString("0"));
  map->Add(new TObjString("fDipoleCurrent"),new TObjString("0"));
  map->Add(new TObjString("fDipolePolarity"),new TObjString("0"));
  map->Add(new TObjString("fCavernTemperature"),new TObjString("0"));
  map->Add(new TObjString("fCavernPressure"),new TObjString("0"));

  return map;
}
