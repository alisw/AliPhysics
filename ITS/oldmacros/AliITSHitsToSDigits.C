Int_t AliITSHitsToSDigits(const char *inFile="galice.root"){
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

    // Connect the Root Galice file containing Geometry, Kine and Hit
  
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
    if (file) {file->Close(); delete file;}
    cout << "AliITSHits2SDigitsDefault" << endl;
    file = new TFile(inFile,"UPDATE");
    if (!file->IsOpen()) {
	cerr<<"Can't open "<<inFile<<" !" << endl;
	return 1;
    } // end if !file
    file->ls();

    // Get AliRun object from file or create it if not on file
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
    } // end ifs

    // SPD
    cout << "Changing from Default SPD simulation, and responce." << endl;
    AliITSDetType *iDetType=ITS->DetType(0);
    AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->
	GetSegmentationModel();
    AliITSresponseSPD *res0 = (AliITSresponseSPD*)iDetType->GetResponseModel();
    AliITSsimulationSPD *sim0=new AliITSsimulationSPD(seg0,res0);
    ITS->SetSimulationModel(0,sim0);
    // test
    cout << "SPD dimensions " << seg0->Dx() << " " << seg0->Dz() << endl;
    cout << "SPD npixels " << seg0->Npz() << " " << seg0->Npx() << endl;
    cout << "SPD pitches " << seg0->Dpz(0) << " " << seg0->Dpx(0) << endl;
    // end test

    // SDD
    cout << "Changing from Default SDD simulation, and responce." << endl;
    //Set response functions
    Float_t baseline = 10.;
    Float_t noise = 1.75;
    // SDD compression param: 2 fDecrease, 2fTmin, 2fTmax or disable,
    // 2 fTolerance
    AliITSDetType *iDetType=ITS->DetType(1);
    AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
    if (!res1) {
	res1=new AliITSresponseSDD();
	ITS->SetResponseModel(1,res1);
    } // end if !res1
    Float_t fCutAmp = baseline + 2.*noise;
    Int_t cp[8]={0,0,fCutAmp,fCutAmp,0,0,0,0}; //1D

    //res1->SetZeroSupp("2D");
    res1->SetZeroSupp("1D");
    res1->SetNoiseParam(noise,baseline);
    res1->SetDo10to8(kTRUE);
    res1->SetCompressParam(cp);
    res1->SetMinVal(4);
    res1->SetDiffCoeff(3.6,40.);
    AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->
	GetSegmentationModel();
    if (!seg1) {
	seg1 = new AliITSsegmentationSDD(ITS->GetITSgeom(),res1);
	ITS->SetSegmentationModel(1,seg1);
    } // end if !seg1
    AliITSsimulationSDD *sim1 = new AliITSsimulationSDD(seg1,res1);
    sim1->SetDoFFT(1);
    sim1->SetCheckNoise(kFALSE);
    ITS->SetSimulationModel(1,sim1);

    // SSD
    cout << "Changing from Default SSD simulation, and responce." << endl;
    AliITSDetType *iDetType = ITS->DetType(2);
    AliITSsegmentationSSD *seg2 = (AliITSsegmentationSSD*)iDetType->
	GetSegmentationModel();
    AliITSresponseSSD *res2 = (AliITSresponseSSD*)iDetType->GetResponseModel();
    res2->SetSigmaSpread(3.,2.);
    AliITSsimulationSSD *sim2 = new AliITSsimulationSSD(seg2,res2);
    ITS->SetSimulationModel(2,sim2);

    if(!gAlice->TreeS()){ 
	cout << "Having to create the SDigits Tree." << endl;
	gAlice->MakeTree("S");
    } // end if !gAlice->TreeS()
    //make branch
    ITS->MakeBranch("S");
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
	ITS->HitsToSDigits(nev,0,-1," ","All"," ");
	timer.Stop(); timer.Print();
    } // event loop

    delete gAlice;   gAlice=0;
    file->Close();
    Long_t size1 = file->GetSize();
    cout << "File size before = " << size0 << " file size after = " << size1;
    cout << "Increase in file size is " << size1-size0 << " Bytes" << endl;
    delete file;
    return 0;
};
