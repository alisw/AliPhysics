
void AliTRDanalyzeBlackEventsLHC(const char *filename) {


  TString strFile = filename;
  strFile += "?EventType=7";
  AliRawReaderRoot *reader = new AliRawReaderRoot(strFile.Data(), 0);
  reader->SelectEquipment(0, 1024, 1041);
  reader->Select("TRD");

  //AliTRDRawStreamTB::SupressWarnings(kTRUE);
  //AliTRDrawStreamTB::SetForceCleanDataOnly();
  AliTRDrawStreamTB::AllowCorruptedData();

  AliTRDrawStreamTB::DisableStackNumberChecker();
  AliTRDrawStreamTB::DisableStackLinkNumberChecker();
  AliTRDrawStreamTB::DisableSkipData();

  AliTRDrawStreamTB *raw = new AliTRDrawStreamTB(reader); 
  //raw->Init();
  //raw->SetRawVersion(3);

  AliTRDqaBlackEvents *qa = new AliTRDqaBlackEvents();
  qa->Init();

  int counter = 0;
  while (reader->NextEvent()) {
    cout << "next event " << counter++ <<  endl;
    cout << qa->AddEvent(raw) << endl;
  }

  cout << "Processing" << endl;
  qa->Process("qaTRD_black.root");
}

// drawing 
//
// AliTRDqaBlackEvents qa;
// qa.DrawSm("qaTRD_black.root", 0);
//
// qa.DrawChamber("qaTRD_black.root", 2);
//

