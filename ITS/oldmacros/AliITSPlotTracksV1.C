#include <iostream.h>
#include <fstream.h>

void AliITSPlotTracksV1()
{
	TFile *file = new TFile("AliITSComparisonV1.root");
	TTree *tree = (TTree*)file->Get("Eval"); 
	
	if (!tree) {
		cerr << "Evaluation tree not found!" << endl;
		return;
	}
	
	// histogram definition - I (efficiency)
	TH1D *hFindables = new TH1D("hFindables", "Findable tracks", 10, 0, 2);
	TH1D *hGood = new TH1D("hGood", "Good found tracks", 10, 0, 2);
	TH1D *hFake = new TH1D("hFake", "Fake found tracks", 10, 0, 2);
	TH1D *hRatioG = new TH1D("hRatioG", "", 10, 0, 2);
	TH1D *hRatioF = new TH1D("hRatioF", "", 10, 0, 2);
	hGood->Sumw2(); 
	hFake->Sumw2();  
	hFindables->Sumw2();
	hRatioG->SetLineColor(kBlue); 
	hRatioG->SetLineWidth(2);
	hRatioF->SetLineColor(kRed); 
	hRatioF->SetLineWidth(2);
	
	// histograms definition - II (resolution)
	TH1D *hPhi    = new TH1D("hPhi", "#Phi resolution", 50, -15., 15.); 
	TH1D *hLambda = new TH1D("hLambda", "#lambda resolution", 50, -15., 15.); 
	TH1D *hPt     = new TH1D("hPt", "Relative Pt resolution", 40, -10., 10.);
	TH1D *hDtot   = new TH1D("hDtot", "Total impact parameter distribution", 100, 0., 2000.); 
	TH1D *hDr     = new TH1D("hDr", "Transv. impact parameter distribution", 50, -1000., 1000.);  
	TH1D *hDz     = new TH1D("hDz", "Long. impact parameter distribution", 50, -1000., 1000.);  
  	hPhi->SetFillColor(4);
	hLambda->SetFillColor(4);
	hPt->SetFillColor(2); 
	hDtot->SetFillColor(6);
	hDr->SetFillColor(kGreen);
	hDz->SetFillColor(kBlue);
	
	// Evaluation tree settings
	Int_t labITS, labTPC, signC;
	Int_t i, j, tot = (Int_t)tree->GetEntries();
	Double_t difpt, diflambda, difphi, Dz, Dr, Dtot, ptg;
	tree->SetBranchAddress("labITS"   , &labITS   );
	tree->SetBranchAddress("labTPC"   , &labTPC   );
	tree->SetBranchAddress("difpt"    , &difpt    ); 
	tree->SetBranchAddress("diflambda", &diflambda);
	tree->SetBranchAddress("difphi"   , &difphi   );
	tree->SetBranchAddress("Dz"       , &Dz       );
	tree->SetBranchAddress("Dr"       , &Dr       );
	tree->SetBranchAddress("Dtot"     , &Dtot     );
	tree->SetBranchAddress("ptg"      , &ptg      );
	tree->SetBranchAddress("signC"    , &signC    );
	
	// Filling the histogram of findable tracks (w.r. to momentum)
	for(i = 0; i < tot; i++) {
		tree->GetEntry(i);
		hFindables->Fill(ptg);
	}
	
	// Filling the evaluation histograms
	Int_t neglabs = 0;
	for(i = 0; i < tot; i++) {
		tree->GetEntry(i);
//		if(signC<0) continue;
		if (labITS < 0) neglabs++;
		if (labITS >= 0) {
			hGood->Fill(ptg); 
			hPt->Fill(difpt);
			hLambda->Fill(diflambda);
			hPhi->Fill(difphi);
			hDtot->Fill(Dtot);
			hDr->Fill(Dr);
			hDz->Fill(Dz);
		}
		else {
			hFake->Fill(ptg);
			neglabs++;
		}
	}
  
	// Drawing
	cerr << "Findable tracks   : " << hFindables->GetEntries() << endl;
	cerr << "Good found tracks : " << hGood->GetEntries() << endl;
	cerr << "Fake found tracks : " << hFake->GetEntries() << endl;
	
	gStyle->SetOptStat(111110);
	gStyle->SetOptFit(1);
	
	TCanvas *c1 = new TCanvas("c1","Parameter resolutions",0,0,700,700);
	c1->Divide(2, 2, 0.001, 0.001);
	c1->cd(1); hPhi->SetXTitle("(mrad)"); hPhi->Draw(); hPhi->Fit("gaus", "Q"); 
	c1->cd(2); hLambda->SetXTitle("(mrad)"); hLambda->Draw(); hLambda->Fit("gaus", "Q");
	c1->cd(3); hPt->SetXTitle("(%)"); hPt->Draw(); hPt->Fit("gaus", "Q");
	c1->cd(4); hDtot->SetXTitle("(micron)"); hDtot->Draw();
	c1->Update();
	
	TCanvas *c2 = new TCanvas("c2","Impact parameters resolution",100,100,700,400);
	c2->Divide(2, 1, 0.001, 0.001);
	c2->cd(1); hDr->SetXTitle("(micron)"); hDr->Draw(); hDr->Fit("gaus", "Q"); 
	c2->cd(2); hDz->SetXTitle("(micron)"); hDz->Draw(); hDz->Fit("gaus", "Q"); 
	c2->Update();
	
	TCanvas *c3 = new TCanvas("c3","Momentum distributions",200,200,800,500);
	c3->Divide(3, 1, 0.001, 0.001);
	c3->cd(1); hGood->Draw(); 
	c3->cd(2); hFake->Draw(); 
	c3->cd(3); hFindables->Draw();
	c3->Update();
	
	TCanvas *c4 = new TCanvas("c4","Tracking efficiency",300,300,800,500);
	hRatioG->Divide(hGood, hFindables, 1., 1.);
	hRatioF->Divide(hFake, hFindables, 1., 1.);
	hRatioG->SetMaximum(1.4);
	hRatioG->SetYTitle("Tracking efficiency");
	hRatioG->SetXTitle("Pt (GeV/c)");
	hRatioG->Draw();  // to not plot the erro bars    hg->Draw("histo");
	hRatioF->SetFillColor(1);
	hRatioF->SetFillStyle(3013);
	hRatioF->SetLineColor(2);
	hRatioF->SetLineWidth(2);
	hRatioF->Draw("same");  // to not plot the error bars  hRatioF->Draw("histosame");
	// a line to mark the best efficiency
	TLine *line1 = new TLine(0,1.0,2,1.0); line1->SetLineStyle(4);
	line1->Draw("same");
	TLine *line2 = new TLine(0,0.9,2,0.9); line2->SetLineStyle(4);
	line2->Draw("histosame");
	// a text explanation
	TText *text = new TText(0.461176,0.248448,"Fake tracks");
	text->SetTextSize(0.05);
	text->Draw();
	text = new TText(0.453919,1.11408,"Good tracks");
	text->SetTextSize(0.05);
	text->Draw();
	c4->Update();
	
	cout << "neglabs = " << neglabs << endl;  // provvisoria
}
