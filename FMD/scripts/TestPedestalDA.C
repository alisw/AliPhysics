void TestPedestalDA(Char_t* fileName, Int_t runNumber){
  
  //This script runs the pedestal DA using the class AliFMDPedestalDA
  
  
  gSystem->Load("libFMDutil");
  Bool_t old = kTRUE;
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(runNumber);
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  AliLog::SetModuleDebugLevel("FMD", 1);
  AliFMDParameters::Instance()->Init();
  AliFMDParameters::Instance()->UseRcuTrailer(!old);
  AliFMDParameters::Instance()->UseCompleteHeader(!old);
  
  AliRawReader *reader = new AliRawReaderDate(fileName,-1);
  TStopwatch timer;
  timer.Start();
  AliFMDPedestalDA pedestalDA;
  //pedestalDA.SetSaveDiagnostics(kTRUE);
  pedestalDA.Run(reader);
  timer.Stop();
  timer.Print();

  
}
