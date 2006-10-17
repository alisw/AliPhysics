//
// This macro reads comparison files and makes plots.
//

void DrawComparison
(const char *fileName = "match-comparison.root")
{
	TH1D *hg = new TH1D("hg", "Good matches / true matches", 20, 0.0,  10.0);
	TH1D *hf = new TH1D("hf", "Fake matches / found matches", 20, 0.0,  10.0);
	
	// 
	// Open file
	//
	TFile *file = TFile::Open(fileName);
	if (!file) return;
	
	TH1D *hgood = (TH1D*)file->Get("hgood");
	TH1D *hfake = (TH1D*)file->Get("hfake");
	TH1D *htrue = (TH1D*)file->Get("htrue");
	TH1D *hfound = (TH1D*)file->Get("hfound");

	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	
	hg->Divide(hgood, htrue, 100.0, 1.0, "b");
	hf->Divide(hfake, hfound, 100.0, 1.0, "b");
	
	TCanvas *c = new TCanvas("c", "", 0, 0, 800, 600);
	
	hg->SetMarkerStyle(21);
	hf->SetMarkerStyle(25);
	
	hg->SetXTitle("p_{T} (GeV/c)");

	hg->GetXaxis()->SetRangeUser(0.0, 6.0);
	hf->GetXaxis()->SetRangeUser(0.0, 6.0);
	
	hg->SetMaximum(120.0);
	hg->SetMinimum(0.0);
	hg->Draw("PE1");
	hf->Draw("PE1same");
}	
