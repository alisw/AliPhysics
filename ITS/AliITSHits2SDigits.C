Int_t AliITSHits2SDigits(const char *inFile = "galice.root"){

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
    if(!file) file = new TFile(fileNameSDigitsSig.Data());
    TDatime *ct0 = new TDatime(2002,04,26,00,00,00), ct = file->GetCreationDate();

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

    if(ct0->GetDate()>ct.GetDate()){
	// For old files, must change SDD noise.
	AliITS *ITS = (AliITS*) gAlice->GetDetector("ITS");
	AliITSresponseSDD *resp1 = ITS->DetType(1)->GetResponseModel();
	resp1->SetNoiseParam();
	resp1->SetNoiseAfterElectronics();
	Float_t n,b;
	Int_t cPar[8];
	resp1->GetNoiseParam(n,b);
	n = resp1->GetNoiseAfterElectronics();
	cPar[0]=0;
	cPar[1]=0;
	cPar[2]=(Int_t)(b + 2.*n + 0.5);
	cPar[3]=(Int_t)(b + 2.*n + 0.5);
	cPar[4]=0;
	cPar[5]=0;
	cPar[6]=0;
	cPar[7]=0;
	resp1->SetCompressParam(cPar);
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

