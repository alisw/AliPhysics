MUONdisplay (Int_t nevent=0, TString fileName="galice.root") {
 
  // Getting runloader 
  AliRunLoader * RunLoader = AliRunLoader::Open(fileName.Data(),"MUONFolder","READ");
  if (RunLoader == 0x0) {
    Error("MUONdisplay","Inut file %s error!",fileName.Data());
    return;   
  }
  RunLoader->LoadHeader();
  RunLoader->LoadKinematics("READ");

  // Getting MUONloader 
  AliLoader * MUONLoader  = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadHits("READ");
  MUONLoader->LoadDigits("READ");
  MUONLoader->LoadRecPoints("READ");

  //  if (RunLoader->GetAliRun() == 0x0) 
  RunLoader->LoadgAlice();
  gAlice = RunLoader->GetAliRun();


// Create Event Display object
   AliMUONDisplay *muondisplay = new AliMUONDisplay(750, MUONLoader);
   RunLoader->GetEvent(nevent);
 
   muondisplay->ShowNextEvent(0);
}
