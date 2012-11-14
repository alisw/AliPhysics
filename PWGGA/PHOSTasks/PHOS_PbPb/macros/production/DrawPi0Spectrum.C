// Macro for drawing (raw) pi0 spectrum
// Uses files produced by MakeMmixPi0.C
// author: Yuri Kharlov <Yuri.Kharlov@cern.ch>

const char* fileName = "output/Pi0_FitResult.root";
const char *centrality[] = {"0-10%","10-40%","40-80%"};

DrawPi0Spectrum()
{
  PPRstyle();
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadLeftMargin(0.12);

  //DrawPi0SpectrumCentral();
  //DrawPi0SpectrumSemiCentral(1);
  // DrawPi0SpectrumSemiCentral(2);

   DrawPi0SpectrumCentralTrigger();
  // DrawPi0SpectrumSemiCentralTrigger(2);
}

//-----------------------------------------------------------------------------
DrawPi0SpectrumCentral()
{
  f1=new TFile("output/Pi0_FitResult.root");

  h1 = (TH1D*)f1->Get("MixAll_cent0_kCentral_yr2int");
  h2 = (TH1D*)f1->Get("MixCPV_cent0_kCentral_yr2int");
  h3 = (TH1D*)f1->Get("MixDisp_cent0_kCentral_yr2int");

  Int_t nPt = h1->GetNbinsX();
  for (Int_t iPt=1; iPt<=nPt; iPt++) {
    Double_t binWidth = h1->GetBinWidth(iPt);
    h1->SetBinContent(iPt,h1->GetBinContent(iPt)/binWidth);
    h1->SetBinError  (iPt,h1->GetBinError  (iPt)/binWidth);
    h2->SetBinContent(iPt,h2->GetBinContent(iPt)/binWidth);
    h2->SetBinError  (iPt,h2->GetBinError  (iPt)/binWidth);
    h3->SetBinContent(iPt,h3->GetBinContent(iPt)/binWidth);
    h3->SetBinError  (iPt,h3->GetBinError  (iPt)/binWidth);
  }

  h1->SetXTitle("p_{T} (GeV/c)");
  h2->SetXTitle("p_{T} (GeV/c)");
  h3->SetXTitle("p_{T} (GeV/c)");
  h1->SetYTitle("(1/N_{AA}) dN/p_{T} (GeV/c)^{-1}");
  h2->SetYTitle("(1/N_{AA}) dN/p_{T} (GeV/c)^{-1}");
  h3->SetYTitle("(1/N_{AA}) dN/p_{T} (GeV/c)^{-1}");
  h1->SetMarkerStyle(20);
  h2->SetMarkerStyle(20);
  h3->SetMarkerStyle(20);
  h1->SetMarkerColor(kRed);
  h2->SetMarkerColor(kBlue);
  h3->SetMarkerColor(kGreen+2);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h3->SetLineColor(kGreen+2);

  TH1D *effCPV = (TH1D*)h2->Clone("effCPV");
  effCPV->Divide(h2,h1,1.,1.,"B");
  effCPV->SetYTitle("PID/no PID");
  effCPV->SetMarkerColor(kBlue);
  effCPV->SetLineColor(kBlue);

  TH1D *effDisp = (TH1D*)h3->Clone("effDisp");
  effDisp->Divide(h3,h1,1.,1.,"B");
  effDisp->SetYTitle("PID/no PID");
  effDisp->SetMarkerColor(kGreen+2);
  effDisp->SetLineColor(kGreen+2);

  TCanvas *c11 = new TCanvas("c11","c11");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  l11 = new TLegend(0.6,0.75,0.89,0.89);
  l11->SetFillColor(kWhite);
  l11->SetBorderSize(1);
  l11->AddEntry(h1,"No PID","lp");
  l11->AddEntry(h2,"PID=CPV","lp");
  l11->AddEntry(h3,"PID=Disp","lp");
  l11->Draw();
  c11->Print(Form("LHC11h_pi0Spectrum_PIDspectra.eps"));

  TCanvas *c12 = new TCanvas("c12","c12");
  gPad->SetGridx();
  gPad->SetGridy();
  effCPV ->SetAxisRange(0.,1.09,"Y");
  effCPV ->Draw();
  effDisp->Draw("same");
  l12 = new TLegend(0.75,0.80,0.94,0.94);
  l12->SetFillColor(kWhite);
  l12->SetBorderSize(1);
  l12->AddEntry(effCPV ,"PID=CPV","lp");
  l12->AddEntry(effDisp,"PID=Disp","lp");
  l12->Draw();
  c12->Print(Form("LHC11h_pi0Spectrum_PIDeffi.eps"));
}

//-----------------------------------------------------------------------------
void DrawPi0SpectrumSemiCentral(const Int_t cent=1)
{
  f1=new TFile("output/Pi0_FitResult.root");
  //f3=new TFile("output/LHC11h_Pi0_FitResult.root");

  h1 = (TH1D*)f1->Get(Form("MixDisp_cent%d_kPHOSPb_yr2int",cent));
  h3 = (TH1D*)f1->Get(Form("MixDisp_cent%d_kSemiCentral_yr2",cent));

  TPaveText *txt = new TPaveText(0.6,0.9,0.8,0.99,"NDC");
  txt->SetFillColor(kWhite);
  txt->SetBorderSize(1);
  txt->AddText(Form("Centrality %s",centrality[cent]));

  Int_t nPt = h1->GetNbinsX();
  for (Int_t iPt=1; iPt<=nPt; iPt++) {
    Double_t binWidth = h1->GetBinWidth(iPt);
    h1->SetBinContent(iPt,h1->GetBinContent(iPt)/binWidth);
    h1->SetBinError  (iPt,h1->GetBinError  (iPt)/binWidth);
    h3->SetBinContent(iPt,h3->GetBinContent(iPt)/binWidth);
    h3->SetBinError  (iPt,h3->GetBinError  (iPt)/binWidth);
  }

  h1->SetXTitle("p_{T} (GeV/c)");
  h3->SetXTitle("p_{T} (GeV/c)");
  h1->SetYTitle("(1/N_{AA}) dN/p_{T} (GeV/c)^{-1}");
  h3->SetYTitle("(1/N_{AA}) dN/p_{T} (GeV/c)^{-1}");
  h1->SetMarkerStyle(20);
  h3->SetMarkerStyle(20);
  h1->SetMarkerColor(kRed);
  h3->SetMarkerColor(kBlue);
  h1->SetLineColor(kRed);
  h3->SetLineColor(kBlue);

  TH1D *ratioSemiCent = (TH1D*)h1->Clone(Form("ratioSemiCent%d",cent));
  ratioSemiCent->Divide(h3);
  ratioSemiCent->SetYTitle("kPHOS/kSemiCentral");
  ratioSemiCent->SetMarkerColor(kBlack);
  ratioSemiCent->SetLineColor(kBlack);
  //ratioSemiCent->SetAxisRange(0.,150.,"Y");
  ratioSemiCent->Fit("pol0","Q","",10.,30.);
  Double_t suppr = ratioSemiCent->GetFunction("pol0")->GetParameter(0);
  h1->Scale(1./suppr);

  TCanvas *c21 = new TCanvas("c21","c21");
  gPad->SetLogy();
  h3->SetMinimum(2e-09);
  h3->Draw();
  h1->Draw("same");
  l21 = new TLegend(0.5,0.75,0.89,0.89);
  l21->SetFillColor(kWhite);
  l21->SetBorderSize(1);
  l21->AddEntry(h1,"PHOS trigger","lp");
  l21->AddEntry(h3,"SemiCentral trigger","lp");
  l21->Draw();
  txt->Draw();
  c21->Print(Form("LHC11h_pi0Spectrum_SemiCent%d.eps",cent));

  TCanvas *c22 = new TCanvas("c22","c22");
  ratioSemiCent->DrawClone();
  txt->Draw();
  c22->Print(Form("LHC11h_pi0Ratio_SemiCent%d.eps",cent));
}

DrawPi0SpectrumCentralTrigger()
{
  f1=new TFile("output/Pi0_FitResult.root");
  // f2=new TFile("out/kMB/Pi0_FitResult_0.root");
  // f3=new TFile("out/kPHOSPb/Pi0_FitResult_0.root");

  h1 = (TH1D*)f1->Get("MixCPV_cent0_kCentral_yr2int");
  h2 = (TH1D*)f1->Get("MixCPV_cent0_kMB_yr2int");
  h3 = (TH1D*)f1->Get("MixCPV_cent0_kPHOSPb_yr2int");

  Int_t nPt = h1->GetNbinsX();
  for (Int_t iPt=1; iPt<=nPt; iPt++) {
    Double_t binWidth = h1->GetBinWidth(iPt);
    h1->SetBinContent(iPt,h1->GetBinContent(iPt)/binWidth);
    h1->SetBinError  (iPt,h1->GetBinError  (iPt)/binWidth);
    h2->SetBinContent(iPt,h2->GetBinContent(iPt)/binWidth);
    h2->SetBinError  (iPt,h2->GetBinError  (iPt)/binWidth);
    h3->SetBinContent(iPt,h3->GetBinContent(iPt)/binWidth);
    h3->SetBinError  (iPt,h3->GetBinError  (iPt)/binWidth);
  }

  h1->GetYaxis()->SetRangeUser(10.e**-8, 0.2);
  h1->SetXTitle("p_{T} (GeV/c)");
  h2->SetXTitle("p_{T} (GeV/c)");
  h3->SetXTitle("p_{T} (GeV/c)");
  h1->SetTitle("Raw production, PID=CPV, centrality: 0-10%");
  h1->SetYTitle("(1/N_{AA}) dN/dp_{T} (GeV/c)^{-1}");
  h2->SetYTitle("(1/N_{AA}) dN/dp_{T} (GeV/c)^{-1}");
  h3->SetYTitle("(1/N_{AA}) dN/dp_{T} (GeV/c)^{-1}");
  h1->SetMarkerStyle(20);
  h2->SetMarkerStyle(21);
  h3->SetMarkerStyle(22);
  h1->SetMarkerColor(kRed);
  h2->SetMarkerColor(kBlue);
  h3->SetMarkerColor(kGreen+2);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h3->SetLineColor(kGreen+2);

  TH1D *effCPV = (TH1D*)h2->Clone("effCPV");
  effCPV->Divide(h2,h1,1.,1.,"B");
  effCPV->SetYTitle("PID/no PID");
  effCPV->SetMarkerColor(kBlue);
  effCPV->SetLineColor(kBlue);

  TH1D *effDisp = (TH1D*)h3->Clone("effDisp");
  effDisp->Divide(h3,h1,1.,1.,"B");
  effDisp->SetYTitle("PID/no PID");
  effDisp->SetMarkerColor(kGreen+2);
  effDisp->SetLineColor(kGreen+2);

  TCanvas *c11 = new TCanvas("c11","c11");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();

  // TF1* func = new TF1("func", "exp([0]+[1]*x)/x**4", 2, 30);
  // func->SetLineColor(kRed);
  // //h1->Fit(func, "", "", 2, 10);
  // //h1->Fit(func, "+", "", 8, 20);
  // h1->Fit(func, "+", "", 5, 30);

  // func->SetLineColor(kBlue);
  // h2->Fit(func, "", "", 2, 13);

  // func->SetLineColor(kGreen+2);
  // h3->Fit(func, "", "", 10, 25);


  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  l11 = new TLegend(0.6,0.75,0.89,0.89);
  l11->SetFillColor(kWhite);
  l11->SetBorderSize(1);
  l11->AddEntry(h1,"kCentral","lp");
  l11->AddEntry(h2,"kMB","lp");
  l11->AddEntry(h3,"kPHOSPb","lp");
  l11->Draw();
  c11->Print(Form("LHC11h_pi0Spectrum_PIDspectra.eps"));

  TCanvas *c12 = new TCanvas("c12","c12");
  gPad->SetGridx();
  gPad->SetGridy();
  effCPV ->SetAxisRange(0.,1.09,"Y");
  effCPV ->Draw();
  effDisp->Draw("same");
  l12 = new TLegend(0.75,0.80,0.94,0.94);
  l12->SetFillColor(kWhite);
  l12->SetBorderSize(1);
  l12->AddEntry(effCPV ,"PID=CPV","lp");
  l12->AddEntry(effDisp,"PID=Disp","lp");
  l12->Draw();
  c12->Print(Form("LHC11h_pi0Spectrum_PIDeffi.eps"));
}


void DrawPi0SpectrumSemiCentralTrigger(int cent = 1)
{
  f1=new TFile("output/Pi0_FitResult.root");
  // f2=new TFile(Form("out/kMB/Pi0_FitResult_%d.root", cent));
  // f3=new TFile(Form("out/kPHOSPb/Pi0_FitResult_%d.root", cent));

  h1 = (TH1D*)f1->Get(Form("MixCPV_cent%d_kSemiCentral_yr2int", cent));
  h2 = (TH1D*)f1->Get(Form("MixCPV_cent%d_kMB_yr2int", cent));
  h3 = (TH1D*)f1->Get(Form("MixCPV_cent%d_kPHOSPb_yr2int", cent));

  Int_t nPt = h1->GetNbinsX();
  for (Int_t iPt=1; iPt<=nPt; iPt++) {
    Double_t binWidth = h1->GetBinWidth(iPt);
    h1->SetBinContent(iPt,h1->GetBinContent(iPt)/binWidth);
    h1->SetBinError  (iPt,h1->GetBinError  (iPt)/binWidth);
    h2->SetBinContent(iPt,h2->GetBinContent(iPt)/binWidth);
    h2->SetBinError  (iPt,h2->GetBinError  (iPt)/binWidth);
    h3->SetBinContent(iPt,h3->GetBinContent(iPt)/binWidth);
    h3->SetBinError  (iPt,h3->GetBinError  (iPt)/binWidth);
  }

  h1->GetYaxis()->SetRangeUser(10.e**-8, 0.2);
  h1->SetXTitle("p_{T} (GeV/c)");
  h2->SetXTitle("p_{T} (GeV/c)");
  h3->SetXTitle("p_{T} (GeV/c)");
  h1->SetTitle(Form("Raw production, PID=CPV, centrality: %s", centrality[cent]));
  h1->SetYTitle("(1/N_{AA}) dN/dp_{T} (GeV/c)^{-1}");
  h2->SetYTitle("(1/N_{AA}) dN/dp_{T} (GeV/c)^{-1}");
  h3->SetYTitle("(1/N_{AA}) dN/dp_{T} (GeV/c)^{-1}");
  h1->SetMarkerStyle(20);
  h2->SetMarkerStyle(21);
  h3->SetMarkerStyle(22);
  h1->SetMarkerColor(kRed);
  h2->SetMarkerColor(kBlue);
  h3->SetMarkerColor(kGreen+2);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h3->SetLineColor(kGreen+2);

  TH1D *effCPV = (TH1D*)h2->Clone("effCPV");
  effCPV->Divide(h2,h1,1.,1.,"B");
  effCPV->SetYTitle("PID/no PID");
  effCPV->SetMarkerColor(kBlue);
  effCPV->SetLineColor(kBlue);

  TH1D *effDisp = (TH1D*)h3->Clone("effDisp");
  effDisp->Divide(h3,h1,1.,1.,"B");
  effDisp->SetYTitle("PID/no PID");
  effDisp->SetMarkerColor(kGreen+2);
  effDisp->SetLineColor(kGreen+2);

  TCanvas *c11 = new TCanvas("c11","c11");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();

  // TF1* func = new TF1("func", "exp([0]+[1]*x)/x**4", 2, 30);
  // func->SetLineColor(kRed);
  // //h1->Fit(func, "", "", 2, 10);
  // //h1->Fit(func, "+", "", 8, 20);
  // h1->Fit(func, "+", "", 5, 30);

  // func->SetLineColor(kBlue);
  // h2->Fit(func, "", "", 2, 13);

  // func->SetLineColor(kGreen+2);
  // h3->Fit(func, "", "", 10, 25);


  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  l11 = new TLegend(0.6,0.75,0.89,0.89);
  l11->SetFillColor(kWhite);
  l11->SetBorderSize(1);
  l11->AddEntry(h1,"kSemiCentral","lp");
  l11->AddEntry(h2,"kMB","lp");
  l11->AddEntry(h3,"kPHOSPb","lp");
  l11->Draw();
  c11->Print(Form("LHC11h_pi0Spectrum_PIDspectra.eps"));

  TCanvas *c12 = new TCanvas("c12","c12");
  gPad->SetGridx();
  gPad->SetGridy();
  effCPV ->SetAxisRange(0.,1.09,"Y");
  effCPV ->Draw();
  effDisp->Draw("same");
  l12 = new TLegend(0.75,0.80,0.94,0.94);
  l12->SetFillColor(kWhite);
  l12->SetBorderSize(1);
  l12->AddEntry(effCPV ,"PID=CPV","lp");
  l12->AddEntry(effDisp,"PID=Disp","lp");
  l12->Draw();
  c12->Print(Form("LHC11h_pi0Spectrum_PIDeffi.eps"));
}

//-----------------------------------------------------------------------------
PPRstyle()
{

  gStyle->SetPalette(1);
  // gStyle->SetCanvasBorderMode(-1);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetCanvasColor(10);

  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderSize(1);
  // gStyle->SetFrameBorderMode(-1);
  // gStyle->SetFrameLineWidth(1.2);
  gStyle->SetFrameLineColor(1);

  gStyle->SetHistFillColor(0);
  gStyle->SetHistLineWidth(2);
  gStyle->SetHistLineColor(1);

  gStyle->SetPadColor(10);
  gStyle->SetPadBorderSize(1);
  // gStyle->SetPadBorderMode(-1);

  gStyle->SetStatColor(10);
  gStyle->SetTitleColor(kBlack,"X");
  gStyle->SetTitleColor(kBlack,"Y");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetStatColor(kWhite);

  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetLabelSize(0.04,"Z");
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelFont(42,"Z");
  gStyle->SetStatFont(42);

  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.2,"Y");

  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

}
