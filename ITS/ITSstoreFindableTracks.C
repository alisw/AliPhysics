void ITSstoreFindableTracks(const char *nfile = "galice", Int_t evnum = 0)
{
	TFile *froot = new TFile(Form("%s.root", nfile), "READ");
	
	gAlice = (AliRun*)froot->Get("gAlice");
	if (!gAlice) {
		cout << "gAlice not found in file!!!" << endl;
		return;
	}
	
	AliITS *ITS = (AliITS*)gAlice->GetModule("ITS");
	
	Int_t nparts = gAlice->GetEvent(evnum);
	cout << "Tracks saved in event " << evnum << ": " << nparts << endl;
	
	TClonesArray *recPoints = ITS->RecPoints();
	TTree        *treeRec   = gAlice->TreeR();
   Int_t         ntracks   = gAlice->GetNtrack(); //FCA correction
	Int_t         nmodules  = treeRec->GetEntries();
	Int_t         modmin[6];
	Int_t         modmax[6];
	
	Int_t layer, nlads, ndets;
	AliITSgeom *gm = ITS->GetITSgeom();
	for (layer = 0; layer < 6; layer++) {
		nlads  = gm->GetNladders(layer+1);
		ndets  = gm->GetNdetectors(layer+1);
		modmin[layer] = gm->GetModuleIndex(layer+1, 1, 1);
		modmax[layer] = gm->GetModuleIndex(layer+1, nlads, ndets);
	}
	
	Int_t track, *goodITS[6];
	for (layer = 0; layer < 6; layer++) {
		goodITS[layer] = new Int_t[ntracks];
		for(track = 0; track < ntracks; track++) goodITS[layer][track] = 0;
	}
   
	Int_t irec, index, npoints = 0;
	for (Int_t layer = 1; layer <= 6; layer++) {
		for (Int_t mod = modmin[layer-1]; mod <= modmax[layer-1]; mod++) {
			ITS->ResetRecPoints();
			treeRec->GetEvent(mod);
			cout << "\rLayer " << layer << ", module: " << mod << flush;
			TObjArrayIter *iter = (TObjArrayIter*)recPoints->MakeIterator();
			while ((recp = (AliITSRecPoint*)iter.Next())) {
				for (index = 0; index < 3; index++) {
					track = recp->GetLabel(index);
					if(track < 0) continue;
					if(track > ntracks) {
						cerr << " Error on track number \n";
						continue;
					}
					goodITS[layer-1][track] = 1;
				} //loop over tracks
			} //loop over points
		} //loop over modules
		cout << endl;
	} //loop over layers
	
	// histos filled with neural results
	TFile     *fnew = new TFile(Form("%s_fnd.root", nfile), "RECREATE");
	TParticle *part = 0;
	TTree     *tree = new TTree("tree", "Tree of tracks");
	
	Int_t count = 0, prim = 0;
	Double_t pt = 0.0;
	tree->Branch("count", &count, "count/I");
	tree->Branch("prim", &prim, "prim/I");
	tree->Branch("pt", &pt, "pt/D");
	
	for(track = 0; track < ntracks; track++) {
		prim = 0;
		count = 0;
		part = gAlice->Particle(track);
		if (!(track % 10000)) cout << "\rStoring track " << track << flush;
		if (part->GetFirstMother() == -1) prim = 1;
		for (layer = 0; layer < 6; layer++) count += goodITS[layer][track];
		pt = part->Pt();
		tree->Fill();
	}
	
	fnew->cd();
	tree->Write("tree", TObject::kOverwrite);
	fnew->Close();
}
