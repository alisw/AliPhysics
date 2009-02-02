// 
// Macro to convert ITS local-coordinate points
// into globa lones
// 

Int_t AliITSL2GConvertPointsV2 
(const char* in_name = "ITS.RecPoints.root", 
 const char* out_name = "ITS.Neural.PointsV2.root", Int_t nev = 0)
{
	TStopwatch timer;

	// Open output file
	TFile *in = new TFile(in_name);
	TFile *out = new TFile(out_name, "recreate");
	
	// Load event files
	if (gAlice) {
		delete AliRunLoader::Instance();
		delete gAlice;
		gAlice=0;
	} 
	AliRunLoader* rl = AliRunLoader::Open("galice.root");
    if (rl == 0x0) {
		cerr << "AliITSL2GConvertPoints.C : Can not open session." << endl;
		return 3;
	}
    Int_t retval = rl->LoadgAlice();
	if (retval) {
		cerr << "AliITSL2GConvertPoints.C : LoadgAlice returned error" << endl;
		return 3;
	}
	gAlice=rl->GetAliRun();
	AliITSLoader* gime = (AliITSLoader*)rl->GetLoader("ITSLoader");
	if (gime == 0x0) {
		cerr << "AliITSL2GConvertPoints.C : can not get ITS loader" << endl;
		return 3;
	}
	AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");
	if (!ITS) {
		cerr << "AliITSL2GConvertPoints.C : AliITS object not found on file" << endl;
		return 3;
	}  // end if !ITS
	AliITSgeom *geom = (AliITSgeom*)ITS->GetITSgeom();
	if(!geom) {
		cerr << "AliITSL2GConvertPoints.C : AliITSgeom not found." << endl;
		return 4;
	} // end if
	gime->LoadRecPoints("read");
	rl->GetEvent(nev);
	TTree *TR = gime->TreeR();
	if (!TR) {
		cerr << "AliITSL2GConvertPoints.C : Can't get the clusters tree." << endl;
		return 4;
	}
		
	// Tree of recpoints
	Int_t nModules = 0;
	TTree *TR = (TTree*)in->Get(Form("Event%d/TreeR", nev));
	nModules = (Int_t)TR->GetEntries();
	if (!nModules) {
		cout << "Empty TreeR!!!" << endl;
		return;
	}

	timer.Start();
	
	// Converts and stores the ITS points into global coordinate format
	Int_t pos = 0;
	AliITSRecPoint *local = 0;
	AliITSNeuralPoint *global = 0;
	TTree *TP = new TTree("TreeP", "Event points in global coords");
	TP->Branch("pos", &pos, "pos/I");
	TP->Branch("Points", "AliITSNeuralPoint", &global);
	
	TObjArray *localArray = 0;
	TR->SetBranchAddress("Clusters", &localArray);
	Int_t module, layer, i, j, count, index;
	Double_t locPos[3], globPos[3], locErr[3][3], globErr[3][3];
	
	cout << geom->GetModuleIndex(1,1,1) << endl;
	cout << geom->GetModuleIndex(2,1,1) << endl;
	cout << geom->GetModuleIndex(3,1,1) << endl;
	cout << geom->GetModuleIndex(4,1,1) << endl;
	cout << geom->GetModuleIndex(5,1,1) << endl;
	cout << geom->GetModuleIndex(6,1,1) << endl;
	
	for(module = 0; module < nModules; module++) {
		TR->GetEvent(module);
		count = (Int_t)localArray->GetEntriesFast();
		for (index = 0; index < count; index++) {
			local = (AliITSRecPoint*)localArray->At(index);
			cout << module << " - " << local->GetDetectorIndex() << endl;
			global = new AliITSNeuralPoint(local, geom, module, index);
			global->SetUser(-1);
			global->ConfMap(0.0, 0.0);
			TP->Fill();
			pos++;
		}
	}

	timer.Stop();
	timer.Print();
	cout << TP->GetEntries() << " points collected" << endl;
	
	out->cd();
	out->mkdir(Form("Event%d", nev));
	out->cd(Form("Event%d", nev));
	TP->Write(Form("TreeP", nev));
	out->Close();
}

