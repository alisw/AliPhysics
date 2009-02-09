void TestGainDA(Char_t* fileName, Int_t runNumber){
  
  //This script runs the gain DA using the class AliFMDGainDA
  
  gSystem->Load("libFMDutil");
  Bool_t old = kTRUE;
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(runNumber);
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliLog::SetModuleDebugLevel("FMD", 1);
  AliFMDParameters::Instance()->Init();
  AliFMDParameters::Instance()->UseRcuTrailer(!old);
  AliFMDParameters::Instance()->UseCompleteHeader(!old);
  
  AliRawReader *reader = new AliRawReaderDate(fileName);
  TStopwatch timer;
  timer.Start();
  AliFMDGainDA gainDA;
  gainDA.SetSaveDiagnostics(kTRUE);
  gainDA.Run(reader);
  
  timer.Stop();
  timer.Print();

  
}
