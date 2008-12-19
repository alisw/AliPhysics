void AliITSSD2D(TString inFile, TString outFile);

void AliITSSDigits2Digits(TString inFile= "galice.root", TString outFile = "")
 {
    // This macro takes SDigits and produces Digits. No merging is done
    // and only one galice.root file is used. 
    // Dynamically link some shared libs 
    TStopwatch timer;
    if(gAlice)
     {
       delete AliRunLoader::GetRunLoader();
       delete gAlice;
       gAlice = 0x0;
    } // end if gAlice
    cout << "Creating digits from summable digits for the ITS..." << endl;
    AliITSSD2D(inFile,outFile);
    timer.Stop(); 
    timer.Print();
}

void AliITSSD2D(TString inFile, TString outFile){
  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  manager->SetInputStream(0,inFile);
  if(outFile != "")manager->SetOutputFile(outFile);
  AliITSDigitizer *dITS  = new AliITSDigitizer(manager);
  manager->Exec("");
  delete manager;
}






