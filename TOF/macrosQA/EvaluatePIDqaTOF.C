void EvaluatePIDqaTOF(const char* gridDir){

	//
	// macro to evaluate PID QA for TOF
	//

	TGrid::Connect("alien://");
	TFile* f = TFile::Open(Form("%s/QAresults.root",gridDir));
	TDirectoryFile* d = (TDirectoryFile*)f->Get("PIDqa");
	TList* list = (TList*)d->Get("PIDqa");
	TList* listTOF = (TList*)list->FindObject("TOF");
	TH1F* hT0MakerEff = (TH1F*)listTOF->FindObject("hT0MakerEff");
	TH1F* hnTracksAt_TOF = (TH1F*)listTOF->FindObject("hnTracksAt_TOF");

	new TCanvas();
	hT0MakerEff->Draw();
	new TCanvas();
	hnTracksAt_TOF->Draw();

	TH1F* efficiencyT0Maker = (TH1F*)hT0MakerEff->Clone("efficiencyT0Maker");
	efficiencyT0Maker->Divide(hT0MakerEff, hnTracksAt_TOF, 1, 1, "b");

	new TCanvas();
	
	efficiencyT0Maker->Draw();

}
