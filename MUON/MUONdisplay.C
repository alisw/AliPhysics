MUONdisplay (Int_t nevent=0, TString fileName="galice.root") {
 
  // Getting runloader 
  AliRunLoader * RunLoader = AliRunLoader::Open(fileName,"MUONFolder","READ");
  if (RunLoader == 0x0) {
    Error("MUONdisplay","Inut file %s error!",fileName);
    return;   
  }
  RunLoader->LoadHeader();
  RunLoader->LoadKinematics("READ");

  // Getting MUONloader 
  AliLoader * MUONLoader  = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadHits("READ");
  MUONLoader->LoadDigits("READ");
  MUONLoader->LoadRecPoints("READ");

  if (RunLoader->GetAliRun() == 0x0) RunLoader->LoadgAlice();
  gAlice = RunLoader->GetAliRun();

  // Getting Module MUON  
  AliMUON * MUON  = (AliMUON *) gAlice->GetDetector("MUON");
  if (!MUON) {
    Error("MUONdisplay","Module MUON not found in the input file");
    return;
  }
  // Getting Muon data
  AliMUONData * muondata = MUON->GetMUONData(); 
  muondata->SetLoader(MUONLoader);
  muondata->SetTreeAddress("H,D,RC,GLT");

// Create Event Display object
   AliMUONDisplay *muondisplay = new AliMUONDisplay(750);

// Display first event
   RunLoader->GetEvent(nevent);
   muondisplay->ShowNextEvent(0);
}
