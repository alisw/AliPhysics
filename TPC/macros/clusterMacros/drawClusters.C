/// \file drawClusters.C
/// \brief Macro to draw clusters TPC clusters out of THnSparse which have been created with readClusters.C
///
/// Used for the HLT-TPC cluster verification
///
/// Usage:
///
/// ~~~
/// aliroot -b -l -q drawClusters.C'("<FileName in results>", "<Simulation Id>", "<Simulation Version>", "<PADGROUP>",<SCALE>)'
/// ~~~
///
/// * PADGROUP:
///   * -1: all pads
///   * 0: inner pad region
///   * 1: middle pad region
///   * 2: outter pad region
///
/// * SCALE:
///   * 0:
///   * 1:
///
/// Example:
/// ~~~
/// aliroot -b -l -q drawClusters.C'("results_recPoints_Pythia_20a.root","Pythia","20a")'
/// ~~~
///
/// Will read:
///  * $CWD/results/results_friends_Pythia_20a.root
///
/// Will write:
///  * $CWD/results/images/xyz_Pythia_20a_all_scaled.png
///  * $CWD/results/images/param_Pythia_20a_all_scaled.png
///  * $CWD/results/rootfiles/xyz_Pythia_20a_all_scaled.root
///  * $CWD/results/rootfiles/param_Pythia_20a_all_scaled.root
///
/// \author Jochen Thaeder <jochen@thaeder.de>

// ==================================================================================
void drawClusters( Char_t *file = "results_recPoints_20000.root",
		   Char_t *type = "Pythia", Char_t *version = "8",
		   Int_t padGroup = -1, Bool_t bScale = kTRUE) {

  /// --------------------------------------------------------
  /// -- Setup
  /// --------------------------------------------------------

  // -- Setup style
  SetupStyle();

  // -- Setup legend
  TH1F* tmp1 =  new TH1F("tmp1","",1,1,2);
  tmp1->SetLineColor(kBlack);
  TH1F* tmp2 =  new TH1F("tmp2","",1,1,2);
  tmp2->SetLineColor(kSpring);
  TH1F* tmp3 =  new TH1F("tmp3","",1,1,2);
  tmp3->SetLineColor(kRed);

  TLegend* leg = new TLegend(0.50, 0.66, 0.84, 0.89, "");
  leg->AddEntry(tmp1, "Offline Cluster","L");
  leg->AddEntry(tmp2, "HLT Cluster","L");
  leg->AddEntry(tmp3, "HLT Redux Cluster","L");
  leg->SetFillColor(kWhite);

  // -- Setup canvas
  Int_t cX = 1380;
  Int_t cY =  800;

  TH1D *hh0 = NULL;
  TH1D *hh1 = NULL;
  TH1D *hh2 = NULL;
  TH2D *hhh = NULL;

  // --------------------------------------------------------
  // -- Open file and get THnSparse
  // --------------------------------------------------------

  TFile* results = TFile::Open(Form("results/%s",file));

  THnSparseF* spo    = static_cast<THnSparseF*>(results->Get("spo"));
  THnSparseF* sphhw  = static_cast<THnSparseF*>(results->Get("sphhw"));
  THnSparseF* sphhwR = static_cast<THnSparseF*>(results->Get("sphhwR"));

  if (padGroup == 0) {
    spo->GetAxis(10)->SetRangeUser(-0.5,62.4);
    sphhw->GetAxis(10)->SetRangeUser(-0.5,62.4);
    sphhwR->GetAxis(10)->SetRangeUser(-0.5,62.4);
  }
  else if (padGroup == 1) {
    spo->GetAxis(10)->SetRangeUser(62.6,126.4);
    sphhw->GetAxis(10)->SetRangeUser(62.6,126.4);
    sphhwR->GetAxis(10)->SetRangeUser(62.6,126.4);
  }
  else if (padGroup == 2) {
    spo->GetAxis(10)->SetRangeUser(126.6,159.4);
    sphhw->GetAxis(10)->SetRangeUser(126.6,159.4);
    sphhwR->GetAxis(10)->SetRangeUser(126.6,159.4);
  }

  // --------------------------------------------------------
  // -- Fill Canvas 0
  // --------------------------------------------------------

  TCanvas* cs0 = new TCanvas("cs0", "all - xyz", 10, 10, cX, cY);
  cs0->Divide(3,2);

  cs0->cd(1);
  PrintHist(spo, sphhw, sphhwR, 0, "Local X", bScale);
  leg->Draw("same");

  cs0->cd(2);
  PrintHist(spo, sphhw, sphhwR, 1, "Local Y", bScale);

  cs0->cd(3);
  PrintHist(spo, sphhw, sphhwR, 2, "Local Z", bScale);

  cs0->cd(4);
  PrintHist(spo, sphhw, sphhwR, 10, "Row", bScale);

  cs0->cd(5);
  PrintHist(spo, sphhw, sphhwR, 11, "Pad", bScale);

  cs0->cd(6);
  PrintHist(spo, sphhw, sphhwR, 9, "Timebin", bScale);

  // --------------------------------------------------------
  // -- Fill Canvas 1
  // --------------------------------------------------------

  TCanvas* cs1 = new TCanvas("cs1", "all - param", 10, 10, cX, cY);
  cs1->Divide(3,2);

  cs1->cd(1);
  PrintHist(spo, sphhw, sphhwR, 12, "Sector", bScale);

  cs1->cd(2);
  gPad->SetLogy();
  PrintHist(spo, sphhw, sphhwR, 13, "Q_{tot}", bScale);

  cs1->cd(3);
  gPad->SetLogy();
  PrintHist(spo, sphhw, sphhwR, 14, "Q_{max}", bScale);
  leg->Draw("same");

  cs1->cd(4);
  gPad->SetLogy();
  PrintHist(spo, sphhw, sphhwR, 6, "#sigma Y^{2}", bScale);

  cs1->cd(5);
  gPad->SetLogy();
  PrintHist(spo, sphhw, sphhwR, 7, "#sigma Z^{2}", bScale);

  cs1->cd(6);
  gPad->SetLogy();
  PrintHist(spo, sphhw, sphhwR, 8, "Local X / Local Y", bScale);

  // --------------------------------------------------------
  // -- Write outFiles
  // --------------------------------------------------------

  gSystem->Exec("if [ ! -d ./results/image ] ; then mkdir -p results/images ; fi");
  gSystem->Exec("if [ ! -d ./results/rootfiles ] ; then mkdir -p results/rootfiles ; fi");

  TString name(Form("%s_%s",type, version));

  if (padGroup == -1)
    name += "_all";
  else
    name += Form("_pad%d", padGroup);

  if (bScale == kTRUE)
    name += "_scaled";

  cs0->SaveAs(Form("results/images/xyz_%s.png",   name.Data()));
  cs1->SaveAs(Form("results/images/param_%s.png", name.Data()));

  cs0->SaveAs(Form("results/rootfiles/xyz_%s.root",   name.Data()));
  cs1->SaveAs(Form("results/rootfiles/param_%s.root", name.Data()));

  return;
}

// -----------------------------------------------
void PrintHist(THnSparseF* spo, THnSparseF* sphhw, THnSparseF* sphhwR,
	       Int_t proj, const Char_t *title, Bool_t bScale) {

  TH1D *hh1 = sphhw->Projection(proj);
  if (hh1) {
    if (bScale) hh1->Scale(1./hh1->Integral());
    hh1->SetTitle(title);
    hh1->SetLineColor(kSpring);
    hh1->DrawCopy(); // drawBase
  }

  TH1D *hh0 = spo->Projection(proj);
  if (hh0) {
    if (bScale) hh0->Scale(1./hh0->Integral());
    hh0->SetTitle(title);
    hh0->SetLineColor(kBlack);
    hh0->DrawCopy("same");
  }

  if (hh1)
    hh1->DrawCopy("same");

  TH1D *hh2 = sphhwR->Projection(proj);
  if (hh2) {
    if (bScale) hh2->Scale(1./hh2->Integral());
    hh2->SetTitle(title);
    hh2->SetLineColor(kRed);
    hh2->DrawCopy("same");
  }
}


/*
  hhh = spo->Projection(4,3);
  hhh->SetTitle("Global X vs Global Z");
  hhh->DrawCopy();
  hhh = sphhw->Projection(4,3);
  hhh->SetLineColor(kSpring);
  hhh->DrawCopy("same");
  hhh = sphhwR->Projection(4,3);
  hhh->SetLineColor(kRed);
  hhh->DrawCopy("same");

  hhh = spo->Projection(5,3);
  hhh->SetTitle("Global X vs Global Y");
  hhh->DrawCopy();
  hhh = sphhw->Projection(5,3);
  hhh->SetLineColor(kSpring);
  hhh->DrawCopy("same");
  hhh = sphhwR->Projection(5,3);
  hhh->SetLineColor(kRed);
  hhh->DrawCopy("same");

  hhh = spo->Projection(5,4);
  hhh->SetTitle("Global Y vs Global Z");
  hhh->DrawCopy();
  hhh = sphhw->Projection(5,4);
  hhh->SetLineColor(kSpring);
  hhh->DrawCopy("same");
  hhh = sphhwR->Projection(5,4);
  hhh->SetLineColor(kRed);
  hhh->DrawCopy("same");
*/


// ==================================================================================
void SetupStyle() {
  /// -- setup style

  gROOT->SetStyle("Plain");

  gStyle->SetHatchesSpacing(0.8);
  gStyle->SetHatchesLineWidth(1);

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(10);

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetFillStyle(1001);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(10);

  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(0);

  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);

  gStyle->SetLegendBorderSize(0);

  Int_t font = 42;

  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);

  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);

  gStyle->SetTickLength(0.02,"xy");
  gStyle->SetEndErrorSize(3);

  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelFont(font,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");

  gStyle->SetTitleFont(font,"xyz");
  gStyle->SetTitleOffset(1.3,"xyz");
  gStyle->SetTitleSize(0.04,"xyz");
  gStyle->SetTitleSize(0.04);

  gStyle->SetMarkerSize(1.2);
  gStyle->SetPalette(1,0);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetLineWidth(1);

  return;
}
