Int_t AliITSHitsToSDigitsDefault(const char *inFile="galice.root"){
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
    cout << "input file " << inFile << endl;
    if (file) file->Close(); 
    if (!file) file = new TFile(inFile,"UPDATE");
    file->ls();

    // Get AliRun object from file or create it if not on file

    if (!gAlice) {
	gAlice = (AliRun*)file->Get("gAlice");
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    } // end if !gAlice 

    AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
    if (!ITS) return;

    // Set the simulation models
    AliITSgeom *geom = ITS->GetITSgeom();

    Int_t nbgr_ev=0;

    if(!gAlice->TreeS()){ 
	cout << "Having to create the SDigits Tree." << endl;
	gAlice->MakeTree("S");
    } // end if !gAlice->TreeS()
    //make branch
    ITS->MakeBranch("S");
    ITS->SetTreeAddress();
    gAlice->GetEvent(0);
    cout<<"SDigitizing ITS..." << endl;
    TStopwatch timer;
    Long_t size0 = file->GetSize();

    for (Int_t nev=evNumber1; nev<= evNumber2; nev++) {
	cout << "nev         " <<nev<<endl;
	if(nev>0) {
	    nparticles = gAlice->GetEvent(nev);
	    gAlice->SetEvent(nev);
	    if(!gAlice->TreeD()) gAlice->MakeTree("D");
	    ITS->MakeBranch("D");
	} // end if nev>0
	cout << "nparticles  " <<nparticles<<endl;
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;

	Int_t nbgr_ev=0;
	if(nsignal) nbgr_ev=Int_t(nev/nsignal);
	timer.Start();
	ITS->HitsToDigits(nev,nbgr_ev,size," ","All"," ");
	timer.Stop(); timer.Print();
    } // event loop

    file->Close();
    Long_t size1 = file->GetSize();
    cout << "File size before = " << size0 << " file size after = " << size1;
    cout << "Increase in file size is " << size1-size0 << " Bytes" << endl;
    delete file;
    return 0;
};
