Int_t AliITSHits2DigitsDefault(const char *inFile = "galice.root"){
    ////////////////////////////////////////////////////////////////////
    //      This macro will take hits from a galice.root file and 
    // produce digits for the ITS using the standard detector 
    // simulations. It will measure the time required to do so and the
    // increase in the galice.root file. There is only one input, that
    // of the name of the root file containing the hits and to which the
    // digits will be written to. This macro will process all of the 
    // events on the root file.
    ////////////////////////////////////////////////////////////////////

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    } // end if

    // Connect the Root Galice file containing Geometry, Kine and Hits
  
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
    if (file) {file->Close(); delete file;}
    cout << "AliITSHits2DigitsDefault" << endl;
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

    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");      
    if (!ITS) {
	cerr<<"ITSHits2Digits.C : AliITS object not found on file" << endl;
	return 3;
    }  // end if !ITS
    if(!(ITS->GetITSgeom())){
	cerr << " AliITSgeom not found. Can't digitize with out it." << endl;
	return 4;
    } // end if
    cout << "Digitizing ITS..." << endl;

    TStopwatch timer;
    Long_t size0 = file->GetSize();
    timer.Start();
    gAlice->Hits2Digits("ITS");
    timer.Stop(); timer.Print();

    file->Close();
    Long_t size1 = file->GetSize();
    cout << "File size before = " << size0 << " file size after = " << size1;
    cout << "Increase in file size is " << size1-size0 << " Bytes" << endl;
    delete file;
    return 0;
};

