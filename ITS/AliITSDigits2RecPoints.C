void AliITSDigits2RecPoints(Int_t evNumber1=0,Int_t evNumber2=0, const char *filename="galice.root"){
    TStopwatch timer;

    if(gAlice){
	delete gAlice;
	gAlice = 0;
    } // end if gAlice
    cout << "Creating reconstructed points from digits for the ITS..." << endl;
    AliITSreconstruction *itsr = new AliITSreconstruction(filename);
    timer.Start();
    itsr->Init();
    itsr->Exec(); 
    timer.Stop(); timer.Print();    
    delete itsr;
}
