//
// This macro generates a simple TTree containing
// all true matches collected from one event.
//
// For each match it is stored:
//  - label ID of the track
//  - label ID of the EMCAL cluster
//  - momentum (3D) of the track
//  - vertex (3D) of the track
//
// All tracks which come from kinks are excluded 
// from this computation, and fake tracks are collected.
// If desired, it is possible to change these settings 
// operating on the Bool_t variables listed below.
//  

Bool_t rejectKinks = kTRUE;     // switch to TRUE to include kinks in computation
Bool_t rejectFakes = kFALSE;    // switch to TRUE to include fake tracks in computation

class match_t
{
	public:
	
	Int_t    label;       // GEANT label of particle
	Int_t    indexT;      // index of track in ESD collection
	Int_t    indexC;      // index of cluster in ESD collection
	Double_t p[3];        // track momentum
	Double_t v[3];        // track vertex
};

//
// Read AliESDs.root file and saves all true pairs of trak-cluster.
//
void SaveTrueMatchesSimple_compiled(const char *outFileName)
{
	//
	// open ESD file, retrieve tree and link to branch cursor
	//
	TFile *srcFile = TFile::Open("AliESDs.root");
	if (!srcFile) return;
	TTree *srcTree = (TTree*)srcFile->Get("esdTree");
	if (!srcTree) return;
	AliESD *esd = 0;
	srcTree->SetBranchAddress("ESD", &esd);
	Long64_t nEvents = srcTree->GetEntries();
	
	//
	// Open output file and create output tree
	//
	TFile *outFile = new TFile(outFileName, "RECREATE");
	match_t match;

	//
	// Loop on events and store true matches
	//
	Bool_t isKink;
	Int_t label, count, nTracks, firstCluster, lastCluster;
	for (Long64_t iev = 0; iev < nEvents; iev++) {

		srcTree->GetEntry(iev);
		cout << "Event " << iev + 1 << " of " << nEvents << ": " << endl;
		
		nTracks = esd->GetNumberOfTracks();
		firstCluster = esd->GetFirstEMCALCluster();
		lastCluster = esd->GetFirstEMCALCluster() + esd->GetNumberOfEMCALClusters();
		cout << "Tracks found      : " << nTracks << endl;
		cout << "EMC clusters found: " << lastCluster - firstCluster << endl;
		
		// create new matches tree
		TTree *outTree = new TTree(Form("tm_%d", iev), Form("True matches from event %d", iev));
		outTree->Branch("matches", &match, "label/I:indexT/I:indexC/I:p[3]/D:v[3]/D");
		
		// external loop on tracks
		for (Int_t it = 0; it < nTracks; it++) {
			AliESDtrack *esdTrack = esd->GetTrack(it);
			// start check to reject kinks
			if (rejectKinks) {
				isKink = kFALSE;
				for (Int_t i = 0; i < 3; i++) {
					if (esdTrack->GetKinkIndex(i) > 0) isKink = kTRUE;
				}
				if (isKink) continue;
			}
			// get track GEANT label (to be checked for match)
			label = esdTrack->GetLabel();
			if (rejectFakes) {
				if (label < 0) continue;
			}
			else {
				label = TMath::Abs(label);
			}
			// store track data into candidate match variable
			// anyway, it will be stored in the tree only
			// if it matches a cluster
			match.indexT = it;
			esdTrack->GetPxPyPz(match.p);
			esdTrack->GetXYZ(match.v);
			// internal loop on clusters
			// a counter counts how many true matches are
			// generated for the same track
			count = 0;
			for (Int_t ic = firstCluster; ic < lastCluster; ic++) {
				AliESDCaloCluster *cl = esd->GetCaloCluster(ic);
				// reject pseudo-clusters & unmatched clusters
				if (cl->GetClusterType() != AliESDCaloCluster::kClusterv1) continue;
				if (cl->GetPrimaryIndex() != label) continue;
				// if the method reaches this point, we
				// have found a match to be stored
				match.label = label;
				match.indexC = ic;
				outTree->Fill();
				count++;
			}
			// alert for multiple matches
			if (count > 0) {
				cout << "Found " << count << " clusters which match track " << it << " in ESD";
				if (count > 1) cout << " --> MULTIPLE MATCH";
				cout << endl;
			}
		} // end loop on tracks
		
		outFile->cd();
		outTree->Write();
		delete outTree;
	}
	
	outFile->Close();
}
