Int_t AliITSHits2SDigitsDefault(const char *inFile = "galice.root"){

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    } // end if

    // Connect the Root Galice file containing Geometry, Kine and Hits
  
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
    if (file) {file->Close(); delete file;}
    cout << "AliITSHits2SDigitsDefault" << endl;
    file = new TFile(inFile,"UPDATE");
    if (!file->IsOpen()) {
	cerr<<"Can't open "<<inFile<<" !" << endl;
	return 1;
    } // end if !file
    file->ls();

    // Get AliRun object from file or return if not on file
    if (gAlice) delete gAlice;
    gAlice = (AliRun*)file->Get("gAlice");
    if (!gAlice) {
	cerr << "AliITSITSHits2Digits.C : AliRun object not found on file"
	    << endl;
	return 2;
    } // end if !gAlice

    gAlice->GetEvent(0);
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");      
    if (!ITS) {
	cerr<<"AliITSHits2DigitsDefault.C : AliITS object not found on file"
	    << endl;
	return 3;
    }  // end if !ITS
    if(!(ITS->GetITSgeom())){
	cerr << " AliITSgeom not found. Can't digitize with out it." << endl;
	return 4;
    } // end if

    if(!gAlice->TreeS()){ 
	cout << "Having to create the SDigits Tree." << endl;
	gAlice->MakeTree("S");
    } // end if !gAlice->TreeS()
    //make branch
    ITS->MakeBranch("S");
    ITS->SetTreeAddress();
    cout << "Digitizing ITS..." << endl;

    TStopwatch timer;
    Long_t size0 = file->GetSize();
    timer.Start();
    ITS->Hits2SDigits();
    timer.Stop(); timer.Print();

    delete gAlice;   gAlice=0;
    file->Close();
    Long_t size1 = file->GetSize();
    cout << "File size before = " << size0 << " file size after = " << size1;
    cout << "Increase in file size is " << size1-size0 << " Bytes" << endl;
    delete file;
    return 0;
};

