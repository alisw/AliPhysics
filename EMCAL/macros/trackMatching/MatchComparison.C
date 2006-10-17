//
// This macros compares the matches found by algorithm with the true matches
// found with the "SaveTrueMatchesSimple.C" macro.
// It saves 4 histogram, which contain the Pt distribution of:
//   - all found matches
//   - all correctly found matches
//   - all wrong (fake) found matches
//   - all true matches
// which will then be availabel for computing efficiency and contamination.
//

class match_t
{
	public:
	
	Int_t    label;       // GEANT label of particle
	Int_t    indexT;      // index of track in ESD collection
	Int_t    indexC;      // index of cluster in ESD collection
	Double_t p[3];        // track momentum
	Double_t v[3];        // track vertex
};

void MatchComparison()
{
	//
	// Initialize AliRun manager
	//
	if (gAlice) {
		delete gAlice;
		gAlice = 0;
	}
	
	//
	// Initialize run loader and load Kinematics
	//
	AliRunLoader *runLoader = AliRunLoader::Open("galice.root");
	if (!runLoader) return;
	runLoader->LoadgAlice();
	gAlice = runLoader->GetAliRun();
	runLoader->LoadKinematics();
	
	//
	// Initialize histograms with their error computation
	//
	TH1D *hgood = new TH1D("hgood", "Well matched tracks", 20, 0.0,  10.0);
	TH1D *hfake = new TH1D("hfake", "Fake matched tracks", 20, 0.0,  10.0);
	TH1D *htrue = new TH1D("htrue", "True matches", 20, 0.0,  10.0);
	TH1D *hfound = new TH1D("hfound", "Found matches", 20, 0.0,  10.0);
	hgood->Sumw2();
	hfake->Sumw2();
	htrue->Sumw2();
	hfound->Sumw2();
	
	//
	// Open file containing true matches,
	// retrieve the Tree and link to a cursor.
	//
	TFile *fileTrue = TFile::Open("true-matches.root");
	match_t trueMatch;
	
	//
	// Open file of found matches,
	// link the modified ESD container.
	//
	TFile *fileFound = TFile::Open("matchESD.root");
	TTree *treeFound = (TTree*)fileFound->Get("esdTree");
	AliESD *esd = 0;
	treeFound->SetBranchAddress("ESD", &esd);
	Long64_t nEvents = treeFound->GetEntries();
	
	//
	// Loop on all events
	//
	Int_t im, it, ic, nTrueMatches, nTracks;
	Int_t label, trkLabel, cluLabel;
	for (Long64_t iev = 0; iev < nEvents; iev++) {
		
		// get true matches tree of given event
		TTree *treeTrue = (TTree*)fileTrue->Get(Form("tm_%d", iev));
		treeTrue->SetBranchAddress("matches", &trueMatch);
		nTrueMatches = treeTrue->GetEntries();
		
		// set TTree pointers to selected event
		runLoader->GetEvent(iev);
		treeFound->GetEntry(iev);
		AliStack *stack = runLoader->Stack();
		nTracks = esd->GetNumberOfTracks();
		
		// read all true pairs
		for (im = 0; im < nTrueMatches; im++) {
			treeTrue->GetEntry(im);
			AliESDtrack *track = esd->GetTrack(trueMatch.indexT);
			label = TMath::Abs(track->GetLabel());
			TParticle *p = stack->Particle(label);
			htrue->Fill(p->Pt());
		}
		
		// compare found matches
		for (Int_t it = 0; it < nTracks; it++) {
			AliESDtrack *track = esd->GetTrack(it);
			ic = track->GetEMCALcluster();
			if (ic == -999999999) continue;
			ic = TMath::Abs(ic);
			AliESDCaloCluster *cl = esd->GetCaloCluster(ic);
			if (!cl) continue;
			trkLabel = track->GetLabel();
			cluLabel = cl->GetPrimaryIndex();
			if (trkLabel == cluLabel && trkLabel > 0) {
				TParticle *p = stack->Particle(TMath::Abs(trkLabel));
				hgood->Fill(p->Pt());
				hfound->Fill(p->Pt());
			}
			else  {
				TParticle *p = stack->Particle(TMath::Abs(trkLabel));
				hfake->Fill(p->Pt());
				hfound->Fill(p->Pt());
			}
		}
	}
	
	cout << "True matches : " << htrue->GetEntries() << endl;
	cout << "Found matches: " << hfound->GetEntries() << endl;
	cout << "Good matches : " << hgood->GetEntries() << endl;
	cout << "Fake matches : " << hfake->GetEntries() << endl;
	
	TFile *fout = TFile::Open("match-comparison.root", "RECREATE");
	hgood->Write();
	hfake->Write();
	htrue->Write();
	hfound->Write();
	fout->Close();
}
