void AliITSsDigits2Digits(const char *filename = "galice.root"){
    // This macro takes SDigits and produces Digits. No merging is done
    // and only one galice.root file is used. 
    // Dynamically link some shared libs 
    TStopwatch timer;

    if(gAlice){
	delete gAlice;
	gAlice = 0;
    } // end if gAlice
    AliRunDigitizer * manager = new AliRunDigitizer(1,1);
    manager->SetInputStream(0,filename);
    AliITSDigitizer *dITS  = new AliITSDigitizer(manager);
    timer.Start();
    manager->Exec("");
    timer.Stop(); timer.Print();
}
