Int_t AliITSHitsToDigitsDefault(Int_t evNumber1=0,Int_t evNumber2=0,
				const char *inFile="galice.root"){
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    } // end if

    // Connect the Root Galice file containing Geometry, Kine and Hits

    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
    if (file) {file->Close(); delete file;}
    cout << "AliITSHitsToDigitsDefault" << endl;
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
    AliITS *ITS  = (AliITS*) gAlice->GetDetector("ITS");
    if (!ITS){
	cerr << "ITS not found in gAlice. Exiting." << endl;
	return 3;
    } // end if
    if(!(ITS->GetITSgeom())){
	cerr << " AliITSgeom not found. Can't digitize with out it." << endl;
	return 4;
    } // end if

    if(!gAlice->TreeD()){ 
	cout << "Having to create the Digits Tree." << endl;
	gAlice->MakeTree("D");
    } // end if !gAlice->TreeD()
    //make branch
    ITS->MakeBranch("D");
    ITS->SetTreeAddress();

    cout<<"Digitizing ITS..." << endl;
    TStopwatch timer;
    Long_t size0 = file->GetSize();

    for (Int_t nev=evNumber1; nev<= evNumber2; nev++) {
	cout << "nev         " <<nev<<endl;
	if(nev>0) {
	    gAlice->SetEvent(nev);
	    if(!gAlice->TreeD()) gAlice->MakeTree("D");
	    ITS->MakeBranch("D");
	} // end if nev>0
	if (nev < evNumber1) continue;
	timer.Start();
	ITS->HitsToDigits(nev,0,-1," ","All"," ");
	timer.Stop(); timer.Print();
    } // event loop

    delete gAlice; gAlice=0;
    file->Close();
    Long_t size1 = file->GetSize();
    cout << "File size before = " << size0 << " file size after = " << size1;
    cout << "Increase in file size is " << size1-size0 << " Bytes" << endl;
    delete file;
    return 0;
};
