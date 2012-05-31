//
// This macro reads comparison files and makes plots.
//

void DrawComparison(const char *fileName = "match-comparison.root")
{
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	// 
	// Open file
	//
	TFile *file = TFile::Open(fileName);
	if (!file) return;
	
	TH1D *hgood  = (TH1D*)file->Get("hgood");
	TH1D *hfake  = (TH1D*)file->Get("hfake");
	TH1D *htrue  = (TH1D*)file->Get("htrue");
	TH1D *hfound = (TH1D*)file->Get("hfound");

        TH1D *hg = (TH1D*)hgood->Clone("hg");
        TH1D *hf = (TH1D*)hfake->Clone("hf");
	
	hg->Divide( htrue );
	hf->Divide( hfound);
	hg->Scale(100);
	hf->Scale(100);
	hg->SetMarkerStyle(21);
	hf->SetMarkerStyle(25);
	
	TCanvas *c = new TCanvas("c", "", 0, 0, 800, 600);

	TLegend* leg = new TLegend(0.6,0.8,0.88,0.88);
	leg->SetFillColor(10);
	leg->AddEntry(hg,"Good/True","p");
	leg->AddEntry(hf,"Fake/Found","p");
			
	hg->SetXTitle("p_{T} (GeV/c)");
	hg->SetYTitle("efficiency (%)");
	hg->SetTitle("Track-EMCAL Cluster Matching");
	hg->SetMaximum(120.0);
	hg->SetMinimum(0.0);
	hg->Draw("PE1");
	hf->Draw("PE1same");
	leg->Draw();
	
}	
