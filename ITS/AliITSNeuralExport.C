void AliITSNeuralExport(Int_t min_count_4_good = 5)
{
	TFile *f_in = new TFile("its_neural_new.root");
	TTree *t_in = (TTree*)f_in->Get("TreeT");
	
	AliITSNeuralTrack *in_track = 0;
	t_in->SetBranchAddress("Tracks", &in_track);
	Int_t in_entries = (Int_t)t_in->GetEntries();

	TFile *f_out = new TFile("its_neural_tracks.root", "RECREATE");
	TTree *t_out = new TTree("TreeT0", "Tree of exported tracks from neural tracking into V1 IO track");

	AliITSIOTrack *out_track = 0;
	t_out->Branch("ITStracks", "AliITSIOTrack", &out_track);

	for (Int_t i = 0; i < in_entries; i++) {
		t_in->GetEntry(i);
		out_track = in_track->ExportIOtrack(min_count_4_good);
		t_out->Fill();
	}

	f_out->cd();
	t_out->Write();
	f_out->Close();
}
