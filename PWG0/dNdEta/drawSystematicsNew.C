
Int_t markers[] = {20,21,22,23,28,29};
Int_t colors[]  = {1,2,3,4,6,8,102};

void DrawpiKpAndCombinedZOnly(Float_t upperPtLimit=0.99)
{
  //gROOT->ProcessLine(".L drawPlots.C");
  gSystem->Load("libPWG0base");

  const char* fileNames[] = { "systematics.root", "systematics.root", "systematics.root", "correction_map.root" };
  const char* folderNames[] = { "correction_0", "correction_1", "correction_2", "dndeta_correction" };
  const char* legendNames[] = { "#pi", "K", "p", "standard" };
  Int_t folderCount = 3;



  TH2F* h2DCorrections[4];
  TH1F* h1DCorrections[4];
  for (Int_t i=0; i<4; i++) {
    TFile::Open(fileNames[i]);
    AlidNdEtaCorrection* correctionTmp = new AlidNdEtaCorrection(folderNames[i],folderNames[i]);
    correctionTmp->LoadHistograms();
    
    //    h2DCorrections[i] = correctionTmp->GetTrack2ParticleCorrection()->GetTrackCorrection()->Get2DCorrectionHistogram("yz",-1,1);

    h1DCorrections[i] = correctionTmp->GetTrack2ParticleCorrection()->GetTrackCorrection()->Get1DCorrectionHistogram("z",-10,10,0.8,0.8);
  }

  TH2F* null = new TH2F("","",100,0.1,10,100,0.5,9.99);
  null->SetXTitle("p_{T} (GeV/c)");
  null->SetXTitle("Correction");

  null->Draw();

  //h1DCorrections[0]->SetMaximum(5);
  //h1DCorrections[0]->Draw();
  //  h2DCorrections[0]->Draw("colz");
  for (Int_t i=1; i<4; i++) {
    h1DCorrections[i]->Draw("same");
  }
  
}


void DrawEffectOfChangeInCrossSection() {
  
  //get the data
  TFile* fin = TFile::Open("systematics_xsections.root");

  const Char_t* changes[]  = {"pythia","ddmore","ddless","sdmore","sdless", "dmore", "dless"};
  Int_t colors[] = {1,2,103,102,4,6,1};

  TH1F* hRatios[7];
  for(Int_t i=0; i<7; i++) {
    hRatios[i] = (TH1F*)fin->Get(Form("ratio_vertexReco_triggerBias_%s",changes[i]));
    hRatios[i]->SetLineWidth(2);
    hRatios[i]->SetLineColor(colors[i]);				 
    hRatios[i]->SetMarkerStyle(22);
    hRatios[i]->SetMarkerSize(0.8);
  }

  TPad* p = DrawCanvasAndPad("syst_changeInXsection",700,400);
  p->SetRightMargin(0.2);
  p->SetLeftMargin(0.13);

  TH2F* null = new TH2F("","",100,-1.05,1.05,100,0.93,1.06);
  null->GetXaxis()->SetTitle("#eta");
  null->GetYaxis()->SetTitle("Ratio pythia/modified cross-sections");
  null->Draw();

  TLatex* text[7];

  for(Int_t i=1; i<7; i++) {
    hRatios[i]->Draw("same");
    
    TString str(hRatios[i]->GetTitle());
    str.Remove(0,str.First("DD"));
    str.Remove(str.First(")"),1);
    text[i] = new TLatex(1.08,hRatios[i]->GetBinContent(hRatios[i]->FindBin(0.))-0.002,str.Data());
    text[i]->SetTextColor(colors[i]);
    text[i]->Draw();
  }
}


void DrawEffectOfChangeInComposition() {
  
  //get the data
  TFile* fin = TFile::Open("systematics_composition.root");

  Int_t colors[] = {1,2,103,102,4,6,1};

  TH1F* hRatios[6];

  Float_t etaRange = 0.899;

  for (Int_t i=0; i<6; i++) {
    hRatios[i] = (TH1F*)fin->Get(Form("ratio_%d",i));

    hRatios[i]->SetLineWidth(2);
    hRatios[i]->SetLineColor(colors[i]);
    hRatios[i]->SetMarkerStyle(22);
    hRatios[i]->SetMarkerSize(0.8);

    hRatios[i]->GetXaxis()->SetRangeUser(-etaRange, etaRange);
  }

  TPad* p = DrawCanvasAndPad("syst_changeOfComposition",700,400);
  p->SetRightMargin(0.2);
  p->SetLeftMargin(0.13);

  TH2F* null = new TH2F("","",100,-1.05,1.05,100,0.97,1.03);
  null->GetXaxis()->SetTitle("#eta");
  null->GetYaxis()->SetTitle("Ratio pythia/modified composition");
  null->Draw();

  TLatex* text[6];

  for(Int_t i=0; i<6; i++) {
    hRatios[i]->Draw("same");
    
    TString str(hRatios[i]->GetTitle());
    str.Remove(0,16);
    text[i] = new TLatex(1.08,hRatios[i]->GetBinContent(hRatios[i]->FindBin(0))-0.002,str.Data());
    text[i]->SetTextColor(colors[i]);
    text[i]->SetTextSize(0.053);
    
    text[i]->Draw();
  }


}

TPad* DrawCanvasAndPad(const Char_t* name, Int_t sizeX=600, Int_t sizeY=500) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);

  gStyle->SetTextSize(0.04);
  gStyle->SetTitleSize(0.05,"xyz");
  //gStyle->SetTitleFont(133, "xyz");
  //gStyle->SetLabelFont(133, "xyz");
  //gStyle->SetLabelSize(17, "xyz");
  gStyle->SetLabelOffset(0.01, "xyz");

  gStyle->SetTitleOffset(1.1, "y");
  gStyle->SetTitleOffset(1.1, "x");
  gStyle->SetEndErrorSize(0.0);

  //##############################################

  //making canvas and pads
  TCanvas *c = new TCanvas(name,name,sizeX,sizeY);

  TPad* p1 = new TPad("pad1","", 0, 0.0, 1.0, 1.0, 0, 0, 0);

  p1->SetBottomMargin(0.15);
  p1->SetTopMargin(0.03);
  p1->SetLeftMargin(0.15);
  p1->SetRightMargin(0.03);
  
  p1->SetGridx();
  p1->SetGridy();

  p1->Draw();
  p1->cd();

  return p1;
}

void MisalignmentShowRawTrackPlots()
{
  gSystem->Load("libPWG0base");
  TFile* file = TFile::Open("resC-resC/analysis_esd_raw.root");
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("dndeta");

  TFile* file2 = TFile::Open("resC-fullA/analysis_esd_raw.root");
  dNdEtaAnalysis* fdNdEtaAnalysis2 = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis2->LoadHistograms("dndeta");

  TH3* track1 = fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->Clone("track1");
  TH3* track2 = fdNdEtaAnalysis2->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->Clone("track2");

  // normalize to number of events;
  TH2* event1 = fdNdEtaAnalysis->GetData()->GetEventCorrection()->GetMeasuredHistogram();
  TH2* event2 = fdNdEtaAnalysis2->GetData()->GetEventCorrection()->GetMeasuredHistogram();
  Int_t event1Count = event1->Integral();
  Int_t event2Count = event2->Integral();
  track1->Scale(1.0 / event1Count);
  track2->Scale(1.0 / event2Count);

  const Float_t innerLimit = 0.49;
  const Float_t outerLimit = 0.99;

  track1->GetYaxis()->SetRangeUser(-outerLimit, outerLimit);
  track2->GetYaxis()->SetRangeUser(-outerLimit, outerLimit);
  AliPWG0Helper::CreateDividedProjections(track1, track2, "ze1");
  TH1* fullRange = gROOT->FindObject("track1_ze1_div_track2_ze1");

  track1->GetYaxis()->SetRangeUser(-innerLimit, innerLimit);
  track2->GetYaxis()->SetRangeUser(-innerLimit, innerLimit);
  AliPWG0Helper::CreateDividedProjections(track1, track2, "ze2");
  TH1* central = gROOT->FindObject("track1_ze2_div_track2_ze2");
  central->SetLineColor(1);
  central->SetMarkerStyle(21);

  for (Int_t x=1; x<track1->GetXaxis()->GetNbins(); ++x)
    for (Int_t y=track1->GetYaxis()->FindBin(-innerLimit); y<track1->GetYaxis()->FindBin(innerLimit); ++y)
      for (Int_t z=1; z<track1->GetZaxis()->GetNbins(); ++z)
      {
        track1->SetBinContent(x, y, z, 0);
        track1->SetBinError(x, y, z, 0);
        track2->SetBinContent(x, y, z, 0);
        track2->SetBinError(x, y, z, 0);
      }

  track1->GetYaxis()->SetRangeUser(-outerLimit, outerLimit);
  track2->GetYaxis()->SetRangeUser(-outerLimit, outerLimit);
  AliPWG0Helper::CreateDividedProjections(track1, track2, "ze3");
  TH1* peripheral = gROOT->FindObject("track1_ze3_div_track2_ze3");
  peripheral->SetLineColor(2);
  peripheral->SetMarkerStyle(22);
  peripheral->SetMarkerColor(2);

  TH2* tmp = new TH2F("tmp", ";p_{T} [GeV/c]      ;#frac{tracks residual misalignment}{tracks full misalignment}", 1, 0.1, 10, 1, 0.9, 1.3);

  tmp->SetStats(kFALSE);
  //tmp->GetXaxis()->SetNoExponent();

  Float_t ptStart = 0.1;

  fullRange->GetXaxis()->SetRangeUser(ptStart, 9.9);
  central->GetXaxis()->SetRangeUser(ptStart, 9.9);
  peripheral->GetXaxis()->SetRangeUser(ptStart, 9.9);

  TCanvas* canvas = new TCanvas("MisalignmentShowRawTrackPlots", "MisalignmentShowRawTrackPlots", 700, 400);
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();

  TLegend* legend = new TLegend(0.2, 0.7, 0.4, 0.8);

  legend->AddEntry(central, "|#eta| < 0.5");
  legend->AddEntry(peripheral, "0.5 < |#eta| < 1.0");

  legend->SetFillColor(0);

  tmp->Draw();
  //fullRange->Draw("SAME");
  central->Draw("SAME");
  peripheral->Draw("SAME");

  legend->Draw();

  canvas->SaveAs("syst_mis_ntracks.eps");
}


void drawdNdEtaRatios(const char* canvasName, Int_t n, const char** files, const char** dirs, const char** names, Int_t* histID)
{
  gSystem->Load("libPWG0base");

  TCanvas* canvas = new TCanvas(canvasName, canvasName, 1000, 500);
  canvas->Divide(2, 1);

  TLegend* legend = new TLegend(0.63, 0.73, 0.98, 0.98);
  legend->SetFillColor(0);

  TH1* base = 0;

  for (Int_t i = 0; i < n; ++i)
  {
    TFile::Open(files[i]);

    dNdEtaAnalysis* tmp = new dNdEtaAnalysis(dirs[i], dirs[i]);
    tmp->LoadHistograms();

    TH1* hist = tmp->GetdNdEtaPtCutOffCorrectedHistogram(histID[i]);

    if (i == 0)
      base = hist;

    legend->AddEntry(hist, names[i]);

    hist->SetMarkerColor(colors[i]);
    hist->SetMarkerStyle(markers[i]);

    canvas->cd(1);
    hist->DrawCopy((i == 0) ? "" : "SAME");

    if (i != 0)
    {
      canvas->cd(2);
      hist->Divide(hist, base, 1, 1, "B");
      hist->GetYaxis()->SetRangeUser(0.98, 1.02);
      hist->DrawCopy((i == 1) ? "" : "SAME");
    }
  }

  legend->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void drawdNdEtaRatios(const char* canvasName, Int_t n, const char** files1, const char** files2, const char** dirs1, const char** dirs2, const char** names, Int_t* histID)
{
  gSystem->Load("libPWG0base");

  TCanvas* canvas = new TCanvas(canvasName, canvasName, 700, 400);
  canvas->SetLeftMargin(0.12);

  TLegend* legend = new TLegend(0.35, 0.7, 0.65, 0.85);
  legend->SetFillColor(0);

  for (Int_t i = 0; i < n; ++i)
  {
    TFile::Open(files1[i]);

    dNdEtaAnalysis* tmp1 = new dNdEtaAnalysis(dirs1[i], dirs1[i]);
    tmp1->LoadHistograms();

    TFile::Open(files2[i]);

    dNdEtaAnalysis* tmp2 = new dNdEtaAnalysis(dirs2[i], dirs2[i]);
    tmp2->LoadHistograms();

    TH1* hist1 = tmp1->GetdNdEtaPtCutOffCorrectedHistogram(histID[i]);
    TH1* hist2 = tmp2->GetdNdEtaPtCutOffCorrectedHistogram(histID[i]);

    TH1* division = hist1->Clone();

    division->Divide(hist1, hist2, 1, 1, "B");

    division->SetMarkerColor(colors[i]);
    division->SetMarkerStyle(markers[i]);

    legend->AddEntry(division, names[i]);

    division->SetTitle("");
    division->GetYaxis()->SetTitle("#frac{dN_{ch}/d#eta using MC vtx}{dN_{ch}/d#eta using ESD vtx}");
    division->SetStats(kFALSE);
    division->GetYaxis()->SetTitleOffset(1.3);
    division->GetXaxis()->SetRangeUser(-0.99, 0.99);
    division->GetYaxis()->SetRangeUser(0.981, 1.02);
    division->DrawCopy((i == 0) ? "" : "SAME");
  }

  gPad->SetGridx();
  gPad->SetGridy();

  legend->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

void vertexShiftStudy(Int_t histID)
{
  const char* files[] = { "maps/idealA/mc-vertex/analysis_esd.root", "results/idealC-idealA/analysis_esd.root", "maps/idealA/mc-vertex-shift-0.05/analysis_esd.root", "maps/idealA/mc-vertex-shift-0.1/analysis_esd.root", "maps/idealA/mc-vertex-shift-dep/analysis_esd.root" };
  const char* dirs[] = { "dndeta", "dndeta", "dndeta", "dndeta", "dndeta" };
  const char* names[] = { "mc vtx", "esd vtx", "+ 0.05 cm", "+ 0.1 cm", "old vtx shift" };
  Int_t hist[] = { histID, histID, histID, histID, histID };

  drawdNdEtaRatios("syst_vertex_shift1", 5, files, dirs, names, hist);

  const char* files1[] = { "maps/idealA/mc-vertex/analysis_esd.root", "maps/idealA/mc-vertex/analysis_esd.root", "maps/idealA/mc-vertex/analysis_esd.root" };
  const char* files2[] = { "results/idealC-idealA/analysis_esd.root", "results/idealC-idealA/analysis_esd.root", "results/idealC-idealA/analysis_esd.root" };
  const char* dirs1[] = { "dndeta", "dndeta", "dndeta"};
  const char* names[] = { "|vtx-z| < 10 cm", "-10 cm < vtx-z < 0 cm", "0 cm < vtx-z < 10 cm" };
  Int_t hist[] = { 0, 1, 2 };

  drawdNdEtaRatios("syst_vertex_shift2", 3, files1, files2, dirs, dirs, names, hist);
}

void vertexShift()
{
  TFile::Open("vertex.root");

  TH2* hist = gFile->Get("fVertexCorr");
  TProfile* prof = hist->ProfileX();

  prof->SetStats(kFALSE);
  prof->SetTitle(";MC vtx-z [cm];mean (ESD vtx-z - MC vtx-z) [cm]");
  prof->GetYaxis()->SetTitleOffset(1.2);

  prof->SetLineWidth(2);

  TCanvas* canvas = new TCanvas("syst_vertex_shift", "syst_vertex_shift", 700, 400);

  gPad->SetGridx();
  gPad->SetGridy();

  prof->Draw();

  canvas->SaveAs(Form("%s.eps", canvas->GetName()));
}

