Double_t EvaluatePIDqaTOF(const char* gridDir, Int_t irun, Int_t color, TLegend* leg, TCanvas* cmatchAll, Option_t* opt="hist same err"){

	//
	// macro to evaluate PID QA for TOF
	//

	TGrid::Connect("alien://");
	TFile* f = TFile::Open(Form("%s/QAresults.root",gridDir));
	TDirectoryFile* d = (TDirectoryFile*)f->Get("PIDqa");
	TList* list = (TList*)d->Get("PIDqa");
	TList* listTOF = (TList*)list->FindObject("TOF");
	TH1F* hT0MakerEff = (TH1F*)listTOF->FindObject("hT0MakerEff");
	hT0MakerEff->Sumw2();
	TH1F* hnTracksAt_TOF = (TH1F*)listTOF->FindObject("hnTracksAt_TOF");
	hnTracksAt_TOF->Sumw2();
	TH1F* efficiencyT0Maker = (TH1F*)hT0MakerEff->Clone("efficiencyT0Maker");
	efficiencyT0Maker->Divide(hT0MakerEff, hnTracksAt_TOF, 1, 1, "b");
	efficiencyT0Maker->SetLineColor(color);
	efficiencyT0Maker->SetMarkerColor(color);
	efficiencyT0Maker->GetXaxis()->SetTitle("nTracks at TOF");
	efficiencyT0Maker->GetYaxis()->SetTitle("T0TOF efficiency");
	cmatchAll->cd();
	efficiencyT0Maker->DrawCopy(opt);

	leg->AddEntry(efficiencyT0Maker,Form("run %d",irun),"l");

	TCanvas* cmatch = new TCanvas(Form("cmatch_%d",irun),Form("cmatch_%d",irun),50,50,750,550);
	efficiencyT0Maker->DrawCopy();
	cmatch->Print(Form("T0tofEfficiency_run_%d.png",irun));
	cmatch->Print(Form("T0tofEfficiency_run_%d.root",irun));
	efficiencyT0Maker->Fit("pol0");
	return (Double_t)(efficiencyT0Maker->GetFunction("pol0")->GetParameter(1));
}

void MakeTrendT0Tof(Int_t nruns, Int_t* runs, const char* gridDirBase, const char* pass){

	//
	// macro to make the trending of the T0tof efficiency for the given run list
	// e.g.:
	// .L EvaluatePIDqaTOF.C
	// Int_t runs[2] = {123456,234567}
	// MakeTrendT0Tof(2,runs,"alien:///alice/data/2011/LHC11h","ESDs/pass1_HLT")
	//

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	if (nruns > 12) {
		Printf("More colors needed! change the macro...");
		return;
	}
	Int_t colors[12] = {kOrange, kRed, kMagenta, kBlue, kAzure+10, kGreen, kYellow+2, kOrange+2, kPink-9, kViolet+2, kBlue-2, kAzure+1};
	TCanvas* cmatchAll = new TCanvas(Form("cmatchAll"),Form("cmatchAll"),50,50,750,550);
	TLegend * leg = new TLegend(0.6,0.2,0.8,0.2+0.05*nruns);
	leg->SetFillColor(0);
	//leg->SetBorderSize(0);
	for (Int_t irun = 0; irun<nruns; irun++){
      		TString path(gridDirBase);
		path+="/000";
		path+=Form("%d",runs[irun]);
		path+="/";
		path+=pass;
		Printf("path for run %d = %s",runs[irun],path.Data());		
		Double_t eff = 0;
		if (irun == 0) eff = EvaluatePIDqaTOF(path.Data(), runs[irun], colors[irun], leg, cmatchAll, "hist err");
		else eff = EvaluatePIDqaTOF(path.Data(), runs[irun], colors[irun], leg, cmatchAll);
		Printf("the average efficiency for run %d is %f", runs[irun], eff);
	}
	cmatchAll->cd();
	leg->Draw();
	cmatchAll->Print("T0tof_efficiency.png");
	cmatchAll->Print("T0tof_efficiency.root");
	return;
}
		
