//____________________________________________________________________
//
// Script to compare the output of GEANT 3.21 to FLUKA 2. 
//
/** @ingroup FMD_simple_script
 */
void
Compare() 
{
  TFile* fluka_file = TFile::Open("fluka/FMD.Hits.root", "READ");
  TFile* geant_file = TFile::Open("geant321/FMD.Hits.root", "READ");
  if (!fluka_file || !geant_file) {
    std::cerr << "Failed to open one or more of the input files" 
	      << std::endl;
    return;
  }
  

  fluka_file->cd();
  gDirectory->cd("Event0");
  TTree* fluka_tree = static_cast<TTree*>(gDirectory->Get("TreeH"));

  geant_file->cd();
  gDirectory->cd("Event0");
  TTree* geant_tree = static_cast<TTree*>(gDirectory->Get("TreeH"));

  if (!fluka_tree || !geant_tree) {
    std::cerr << "Failed to get one or more of the trees" 
	      << std::endl;
    return;
  }
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelFont(132, "XY");
  gStyle->SetTitleFont(132, "XY");
  gStyle->SetTextFont(132);
  gStyle->SetNdivisions(10, "XY");
  
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetTopMargin(.05);
  c->SetRightMargin(.05);
  
  TH1* fluka_dist = new TH1F("fluka_dist", "FLUKA Energy deposition", 
			     100, 0, 1);
  fluka_dist->SetLineColor(2);
  fluka_dist->SetFillColor(2);
  fluka_dist->SetFillStyle(3001);
  fluka_dist->GetXaxis()->SetTitle("Energy deposited");
  fluka_dist->GetYaxis()->SetTitle("Frequency");
  fluka_tree->Draw("FMD.fEdep>>fluka_dist");

  TH1* geant_dist = new TH1F("geant_dist", "GEANT Energy deposition", 
			     100, 0, 1);
  geant_dist->SetLineColor(3);
  geant_dist->SetFillColor(3);
  geant_dist->SetFillStyle(3002);
  geant_tree->Draw("FMD.fEdep>>geant_dist", "", "SAME");
  
  TLegend* l = new TLegend(.3, .8, .95, .95);
  l->SetFillColor(0);
  l->SetBorderSize(1);
  l->SetTextFont(132);
  l->AddEntry(fluka_dist, Form("%s - %d counts", 
			       fluka_dist->GetTitle(), 
			       Int_t(fluka_dist->Integral())), "LF");
  l->AddEntry(geant_dist, Form("%s - %d counts", 
			       geant_dist->GetTitle(), 
			       Int_t(geant_dist->Integral())), "LF");
  l->Draw();
  
  c->Modified();
  c->cd();
  c->Print("fluka_vs_geant321.C");
}

//____________________________________________________________________
//
// EOF
//
