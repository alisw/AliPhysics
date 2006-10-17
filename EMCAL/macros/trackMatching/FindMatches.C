//
// This macro performs the matching operation between ESD tracks 
// and EMCAL clusters (stored as AliESDCaloCluster objects).
// Due to the necessity to know the magnetic field, it is supposed that
// this macro runs in the directory where event files are stored
// (galice.root --> recovery magnetic field), and saves its output in
// a directory specified in the argument (default = work directory)
// 

void FindMatches(const char *fileOut = "matchESD.root")
{
	//
	// Initialize AliRun manager.
	//
	if (gAlice) {
		delete gAlice;
		gAlice = 0;
	}
	
	//
	// Initialize AliRunLoader, in order
	// to read magnetic field from simulation.
	//
	AliRunLoader *rl = AliRunLoader::Open("galice.root");
	if (!rl) return;
	rl->LoadgAlice();
	gAlice = rl->GetAliRun();
	AliMagF *magf = gAlice->Field();
	Bool_t constField = (magf->Type() == 1);
	AliTracker::SetFieldMap(magf, constField);
	
	//
	// Open ESD file and recoveries TTree of ESD objects.
	//
	TFile *esdFile = new TFile("AliESDs.root");
	TTree *esdTree = (TTree*)esdFile->Get("esdTree");
	AliESD *esd = 0;
	esdTree->SetBranchAddress("ESD", &esd);
	Long64_t nEvents = esdTree->GetEntries();
	
	//
	// Set this important flag to the tracker.
	// This example works with and AliESD file already saved
	// after the tracking procedure done during AliReconstruction::Run().
	// All tracks saved in such a file, are already propagated to the vertex,
	// while EMCAL matching needs to be done using track status in the outermost
	// available position.
	// Then, we must use the "outer" parameters in the ESD track, and we must
	// tell to the track to copy the outer parameters from its AliESDtrack seed.
	// 
	AliEMCALTrack::SetUseOuterParams(kTRUE);
	
	//
	// Instantiate match finder, and set some cuts.
	// Distances are in centimeters, angles in degrees.
	//
	AliEMCALTracker *mf = new AliEMCALTracker;
	mf->SetCutX(5.0);
	mf->SetCutY(5.0);
	mf->SetCutZ(5.0);
	mf->SetCutAlpha(-50., 50.);  // --> i.e. exclude this cut
	mf->SetCutAngle(10.);
	mf->SetMaxDistance(5.0);
	
	//
	// Define how to manage energy loss correction.
	// Actually, these corrections are exlcluded.
	//
	// mf->SetCorrection(0.0604557, 34.5437); // at the moment, we exclude energy loss correction
	mf->SetTrackCorrectionMode("NONE");
	mf->SetNumberOfSteps(0);
	//TGeoManager::Import("misaligned_geometry.root");	
	
	//
	// Before starting looping on files, the output file is created.
	// It will be structurally identical to the ESD source tree, with the 
	// only difference that here the "fEMCALindex" datamember of the ESD tracks
	// will have been set to a meaningful value.
	//
	TFile *outFile = TFile::Open(fileOut, "RECREATE");
	TTree *outTree = new TTree("esdTree", "ESD with matched clusters");
	outTree->Branch("ESD", "AliESD", &esd);
	
	//
	// Loop on events.
	// The whole track matching procedure is compiled in the 
	// method "PropagateBack" which accepts and input ESD collection.
	// This method modifies the passed object then, the same object is linked
	// to source tree and target tree branch.
	//
	Int_t nTracks, nStep1, nStep2, nSaved;
	for (Long64_t iev = 0; iev < nEvents; iev++) {
		cout << "Finding matches in event " << iev + 1 << "/" << nEvents << endl;
		esdTree->GetEntry(iev);
		mf->PropagateBack(esd);
		outTree->Fill();
	}
	
	// 
	// Save processed tree and close output file.
	//
	outFile->cd();
	outTree->Write("esdTree", TObject::kOverwrite);
	outFile->Close();
}
