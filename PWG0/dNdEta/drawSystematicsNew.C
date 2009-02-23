
Int_t markers[] = {20,21,22,23,28,29};
Int_t colors[]  = {1,2,3,4,6,8,102};

void loadlibs()
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
}

void InitPad()
{
  if (!gPad)
    return;

  gPad->Range(0, 0, 1, 1);
  gPad->SetLeftMargin(0.15);
  //gPad->SetRightMargin(0.05);
  //gPad->SetTopMargin(0.13);
  //gPad->SetBottomMargin(0.1);

  gPad->SetGridx();
  gPad->SetGridy();
}

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

TPad* DrawChange(Bool_t spd, const char* basename, const char** changes, Int_t nChanges, Int_t nDraw, Int_t* colors, const char** names = 0, Float_t scale = 0.10)
{  
  Float_t etaMax = 1.05;
  if (spd)
    etaMax = 1.79;

  TH1F* hRatios[100];
  for(Int_t i=0; i<nChanges; i++) {
    hRatios[i] = (TH1F*)gFile->Get(Form("%s%s",basename,changes[i]));
    hRatios[i]->SetLineWidth(1);
    hRatios[i]->SetMarkerStyle(22);
    hRatios[i]->SetMarkerSize(0.8);

    Float_t average = hRatios[i]->Integral(hRatios[i]->FindBin(-1), hRatios[i]->FindBin(1)) / (hRatios[i]->FindBin(1) - hRatios[i]->FindBin(-1) + 1);
    Printf("%s: %.2f %%" , hRatios[i]->GetTitle(), (average - 1) * 100);
  }
  
  TPad* p = DrawCanvasAndPad("syst_changeInXsection",700,400);
  p->SetRightMargin(0.2);
  p->SetLeftMargin(0.13);
  
  TH2F* null = new TH2F("","",100,-etaMax, etaMax,100,1. - scale,1. + scale);
  null->GetXaxis()->SetTitle("#eta");
  null->GetYaxis()->SetTitle(hRatios[0]->GetYaxis()->GetTitle());
  null->Draw();
  
  line = new TLine(-etaMax, 1, etaMax, 1);
  line->Draw();

  TLatex* text[100];

  for(Int_t i=1; i<nDraw; i++) {
    hRatios[i]->SetLineColor(colors[i]);
    hRatios[i]->SetMarkerColor(colors[i]);
    hRatios[i]->Draw("HISTPL SAME");
    
    TString str(hRatios[i]->GetTitle());
    
    if (names)
      str = names[i];
    
    text[i] = new TLatex(etaMax + 0.03,hRatios[i]->GetBinContent(hRatios[i]->FindBin(etaMax-0.1))-0.002,str.Data());
    text[i]->SetTextAlign(11);
    text[i]->SetTextColor(colors[i]);
    text[i]->Draw();
  }
  
  return p;
}

void DrawEffectOfChangeInCrossSection(Bool_t spd = kFALSE, const char* fileName = "systematics_vtxtrigger_compositions.root") 
{
  TFile::Open(fileName);

  const Char_t* changes[]  = {"pythia","ddmore","ddless","sdmore","sdless", "dmore", "dless", "sdmoreddless", "sdlessddmore", "ddmore25","ddless25","sdmore25","sdless25", "dmore25", "dless25", "sdmoreddless25", "sdlessddmore25" };
  //const Char_t* changes[]  = { "pythia", "qgsm", "phojet" };
  //const Int_t nChanges = 3;
  Int_t colors[] = {1,1,4,1,2,2,4,2,1};

  c = DrawChange(spd, "ratio_vertexReco_triggerBias_", changes, 17, 9, colors, 0);
  c->SaveAs("cross_sections.eps");
}

void DrawEffectOfChangeInComposition(Bool_t spd = kFALSE, const char* fileName = "new_compositions_analysis.root") 
{
  TFile::Open(fileName);

  const Char_t* changes[]  = { "PythiaRatios", "KBoosted", "KReduced", "pBoosted", "pReduced", "KBoostedpBoosted", "KReducedpReduced", "KBoostedpReduced", "KReducedpBoosted"};
  const char*   names[]    = { "",             "K #times 1.5", "K #times 0.5", "p #times 1.5", "p #times 0.5", "K #times 1.5, p #times 1.5", "K #times 0.5, p #times 0.5", "K #times 1.5, p #times 0.5", "K #times 0.5, p #times 1.5" };
  //const Char_t* changes[]  = { "PythiaRatios", "PiBoosted",      "PiReduced", "KBoosted", "KReduced", "pBoosted", "pReduced", "othersBoosted", "othersReduced" };
  //const char*   names[]    = { "",             "#pi #times 1.5", "#pi #times 0.5", "K #times 1.5", "K #times 0.5", "p #times 1.5", "p #times 0.5",  "others #times 1.5", "others #times 0.5" };
  Int_t colors[] = {1,1,2,2,1,2,1,4,4};

  c = DrawChange(spd, "", changes, 9, 9, colors, names, 0.03);
  c->SaveAs("compositions.eps");
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

void MisalignmentShowRawTrackPlots(const char* dirName = "fdNdEtaAnalysisESD")
{
  loadlibs();

  TFile* file = TFile::Open("fullA-simrec/MB2/analysis_esd_raw.root");
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis(dirName, dirName);
  fdNdEtaAnalysis->LoadHistograms();

  TFile* file2 = TFile::Open("fullA-sim/analysis_esd_raw.root");
  dNdEtaAnalysis* fdNdEtaAnalysis2 = new dNdEtaAnalysis(dirName, dirName);
  fdNdEtaAnalysis2->LoadHistograms();

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

  TH2* tmp = new TH2F("tmp", ";p_{T} [GeV/c]      ;#frac{tracks full misalignment during rec.}{tracks ideal misalignment during rec.}", 1, 0.1, 10, 1, 0.9, 1.3);

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
  loadlibs();
  gROOT->ProcessLine(".L $ALICE_ROOT/PWG0/dNdEta/drawPlots.C");

  TCanvas* canvas = new TCanvas(canvasName, canvasName, 1000, 500);
  canvas->Divide(2, 1);
  canvas->cd(2)->SetGridx();
  canvas->cd(2)->SetGridy();

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
      PrintIntegratedDeviation(hist, base, names[i]);

      canvas->cd(2);
      hist->Divide(hist, base, 1, 1);
      hist->GetYaxis()->SetRangeUser(0.9, 1.1);
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

void MaterialBudgetChange()
{
  const char* files[] =
    { "Material-normal-mcvtx/analysis_esd.root",
      "Material-increased-mcvtx/analysis_esd.root",
      "Material-decreased-mcvtx/analysis_esd.root" };

  const char* dirs[] = { "dndeta", "dndeta", "dndeta"};
  const char* names[] = { "no change", "+ 10 % material", "- 10 % material" };
  Int_t hist[] = { 0, 0, 0 };

  drawdNdEtaRatios("MaterialBudgetChange", 3, files, dirs, names, hist);
}

void MisalignmentChange()
{
  const char* files[] =
    { "maps/v4-09-Release/tpc-only/fullC-simrec--fullA-sim-mcvtx/analysis_esd.root",
      "maps/v4-09-Release/tpc-only/fullC-simrec--fullA-simrec-mcvtx/analysis_esd.root" };

  const char* dirs[] = { "dndeta", "dndeta", "dndeta"};
  const char* names[] = { "no change", "increased in sim/rec", "increased in sim" };
  Int_t hist[] = { 0, 0, 0 };

  drawdNdEtaRatios("MisalignmentChange", 2, files, dirs, names, hist);
}

void dNdEtaVertexRanges()
{
  const char* files[] =
    { "analysis_esd.root",
      "analysis_esd.root",
      "analysis_esd.root" };

  const char* dirs[] = { "dndeta", "dndeta", "dndeta"};
  const char* names[] = { "|vtx-z| < 10 cm", "-10 cm < vtx-z < 0 cm", "0 cm < vtx-z < 10 cm" };
  Int_t hist[] = { 0, 1, 2 };

  drawdNdEtaRatios("dNdEtaVertexRanges", 3, files, dirs, names, hist);
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

void CompareRawTrackPlots(const char* fileName1, const char* fileName2, Float_t ptCut = 0.0, Int_t multCut = 1)
{
  loadlibs();

  const char* dirName = "fdNdEtaAnalysisESD";

  TFile* file = TFile::Open(fileName1);
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis(dirName, dirName);
  fdNdEtaAnalysis->LoadHistograms();

  TFile* file2 = TFile::Open(fileName2);
  dNdEtaAnalysis* fdNdEtaAnalysis2 = new dNdEtaAnalysis(dirName, dirName);
  fdNdEtaAnalysis2->LoadHistograms();

  TH3* track1 = fdNdEtaAnalysis->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->Clone("track1");
  TH3* track2 = fdNdEtaAnalysis2->GetData()->GetTrackCorrection()->GetMeasuredHistogram()->Clone("track2");

  TH2* event1 = fdNdEtaAnalysis->GetData()->GetEventCorrection()->GetMeasuredHistogram();
  TH2* event2 = fdNdEtaAnalysis2->GetData()->GetEventCorrection()->GetMeasuredHistogram();
  Int_t event1Count = event1->Integral(event1->GetXaxis()->FindBin(-9.9), event1->GetXaxis()->FindBin(9.9), multCut, event1->GetNbinsY() + 1);
  Int_t event2Count = event2->Integral(event2->GetXaxis()->FindBin(-9.9), event2->GetXaxis()->FindBin(9.9), multCut, event1->GetNbinsY() + 1);

  Float_t nTrack1 = track1->Integral(track1->GetXaxis()->FindBin(-9.9), track1->GetXaxis()->FindBin(9.9), track1->GetYaxis()->FindBin(-0.79), track1->GetYaxis()->FindBin(0.79), track1->GetZaxis()->FindBin(ptCut), track1->GetZaxis()->GetNbins());
  Float_t nTrack2 = track2->Integral(track2->GetXaxis()->FindBin(-9.9), track2->GetXaxis()->FindBin(9.9), track2->GetYaxis()->FindBin(-0.79), track2->GetYaxis()->FindBin(0.79), track2->GetZaxis()->FindBin(ptCut), track2->GetZaxis()->GetNbins());

  Printf("%d tracks in %d events in first sample; %d tracks in %d events in second sample", (Int_t) nTrack1, event1Count, (Int_t) nTrack2, event2Count);

  // normalize to number of events;
  nTrack1 /= event1Count;
  nTrack2 /= event2Count;

  Printf("There are %.2f tracks/event in the first sample and %.2f tracks/event in the second sample. %.2f %% difference (with pt cut at %.2f GeV/c)", nTrack1, nTrack2, 100.0 * (nTrack1 - nTrack2) / nTrack1, ptCut);

  gROOT->cd();

  AliPWG0Helper::CreateDividedProjections(track1, track2);

  new TCanvas; gROOT->FindObject("track1_yx_div_track2_yx")->Draw("COLZ");
  new TCanvas; gROOT->FindObject("track1_zx_div_track2_zx")->Draw("COLZ");
  new TCanvas; gROOT->FindObject("track1_zy_div_track2_zy")->Draw("COLZ");

  for (Int_t i=0; i<3; i++)
  {
    char c = 'x' + (char) i;

    /*proj1 = track1->Project3D(Form("%ce2", c));
    proj2 = track2->Project3D(Form("%ce2", c));
    AliPWG0Helper::NormalizeToBinWidth(proj1);
    AliPWG0Helper::NormalizeToBinWidth(proj2);

    new TCanvas;
    proj1->DrawCopy();
    proj2->SetLineColor(2);
    proj2->SetMarkerColor(2);
    proj2->DrawCopy("SAME");*/

    AliPWG0Helper::CreateDividedProjections(track1, track2, Form("%ce2", c));
    TH1* pt = gROOT->FindObject(Form("track1_%ce2_div_track2_%ce2", c, c));
    new TCanvas; pt->DrawCopy();
    gPad->SetGridx(); gPad->SetGridy();
  }

  event1_x = event1->ProjectionX("event1_x");
  event1_y = event1->ProjectionY("event1_y");
  event2_x = event2->ProjectionX("event2_x");
  event2_y = event2->ProjectionY("event2_y");

  new TCanvas; event1_x->DrawCopy(); event2_x->SetLineColor(2); event2_x->DrawCopy("SAME"); 

  event1_x->Divide(event2_x);
  event1_y->Divide(event2_y);

  new TCanvas; event1_x->DrawCopy();
  new TCanvas; event1_y->DrawCopy();

  event1->Divide(event2);
  new TCanvas;
  event1->Draw("COLZ");
  event1->SetMinimum(0.5);
  event1->SetMaximum(2);


}

void MagnitudeOfCorrection(const char* fileName, const char* dirName = "dndeta", Float_t ptCut = 0.3)
{
  loadlibs();

  TFile* file = TFile::Open(fileName);
  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis(dirName, dirName);
  fdNdEtaAnalysis->LoadHistograms();
  fdNdEtaAnalysis->GetData()->PrintInfo(ptCut);
}

Double_t ConvSigma1To2D(Double_t sigma)
{
  return TMath::Sqrt( - TMath::Log( 1 - TMath::Erf(sigma / TMath::Sqrt(2)) ) * 2);
}

Double_t ConvDistance1DTo2D(Double_t distance)
{
  return TMath::ErfInverse(1 - TMath::Exp(-distance * distance / 2)) * TMath::Sqrt(2);
}

Double_t Sigma2VertexCount(TH2F* tracks, Double_t nSigma)
{
  Double_t count = 0;

  //nSigma = ConvSigma1To2D(nSigma);

  for (Int_t x=1; x<=tracks->GetNbinsX(); ++x)
    for (Int_t y=1; y<=tracks->GetNbinsY(); ++y)
    {
      Double_t impactX = tracks->GetXaxis()->GetBinCenter(x);
      Double_t impactY = tracks->GetYaxis()->GetBinCenter(y);

      Float_t d = TMath::Sqrt(impactX*impactX + impactY*impactY);

      d = ConvDistance1DTo2D(d);

      if (d < nSigma)
        count += tracks->GetBinContent(x, y);
    }

  return count;
}

TH2F* Sigma2VertexGaussianTracksHist()
{
  TH2F* tracks = new TH2F("Sigma2Vertex_tracks", "Sigma2Vertex_tracks", 200, -5, 5, 200, -5, 5);

  TF2* gaussian2D = new TF2("gaussian2D", "xgausn(0) * ygausn(3)", -5, 5, -5, 5);
  gaussian2D->SetParameters(1, 0, 1, 1, 0, 1);

  for (Int_t x=1; x<=tracks->GetNbinsX(); ++x)
    for (Int_t y=1; y<=tracks->GetNbinsY(); ++y)
      tracks->SetBinContent(x, y, gaussian2D->Eval(tracks->GetXaxis()->GetBinCenter(x), tracks->GetYaxis()->GetBinCenter(y)));

  //normalize
  tracks->Scale(1.0 / tracks->Integral());

  return tracks;
}

TH1F* Sigma2VertexGaussian()
{
  TH2F* tracks = Sigma2VertexGaussianTracksHist();

  TCanvas* canvas = new TCanvas("Sigma2VertexGaussian", "Sigma2VertexGaussian", 1000, 1000);
  canvas->Divide(2, 2);

  canvas->cd(1);
  tracks->Draw("COLZ");

  TH1F* ratio = new TH1F("Sigma2Vertex_ratio", "Sigma2Vertex_ratio;n sigma;included", 50, 0.05, 5.05);
  for (Double_t nSigma = 0.1; nSigma < 5.05; nSigma += 0.1)
    ratio->Fill(nSigma, Sigma2VertexCount(tracks, nSigma));
  ratio->SetMarkerStyle(21);

  canvas->cd(2);
  ratio->DrawCopy("P");

  TH1F* ratio2 = new TH1F("Sigma2Vertex_ratio2", "Sigma2Vertex_ratio2;nSigma;% included 4 sigma / % included n sigma", 50, 0.05, 5.05);
  Double_t sigma3 = Sigma2VertexCount(tracks, 4);
  for (Double_t nSigma = 0.1; nSigma < 5.05; nSigma += 0.1)
    ratio2->Fill(nSigma, sigma3 / ratio->GetBinContent(ratio->FindBin(nSigma)));
  ratio2->SetMarkerStyle(21);

  canvas->cd(3);
  ratio2->DrawCopy("P");

  canvas->SaveAs("Sigma2Vertex.eps");

  return ratio2;
}

TH1F** Sigma2VertexSimulation(const char* fileName = "systematics.root")
{
  TFile* file = TFile::Open(fileName);

  TH1F* sigmavertex = dynamic_cast<TH1F*> (file->Get("fSigmaVertexTracks"));
  TH1F* sigmavertexPrim = dynamic_cast<TH1F*> (file->Get("fSigmaVertexPrim"));
  if (!sigmavertex || !sigmavertexPrim)
  {
    printf("Could not read histogram(s)\n");
    return;
  }

  // calculate ratio
  TH1F* ratio = new TH1F("sigmavertexsimulation_ratio", "sigmavertexsimulation_ratio;N#sigma;% included in 4 #sigma / % included in N#sigma", sigmavertex->GetNbinsX(), sigmavertex->GetXaxis()->GetXmin(), sigmavertex->GetXaxis()->GetXmax());

  // calculate contamination
  TH1F* contamination = ratio->Clone("sigmavertexsimulation_contamination");
  contamination->SetTitle("sigmavertexsimulation_contamination;N#sigma;1 + N_{secondaries} / N_{all}");

  for (Int_t i=1; i<=sigmavertex->GetNbinsX(); ++i)
  {
    ratio->SetBinContent(i, sigmavertex->GetBinContent(sigmavertex->GetXaxis()->FindBin(4)) / sigmavertex->GetBinContent(i));
    contamination->SetBinContent(i, 1 + (sigmavertex->GetBinContent(i) - sigmavertexPrim->GetBinContent(i)) / sigmavertex->GetBinContent(i));
  }

  // print stats
  for (Float_t sigma = 2.0; sigma < 5.25; sigma += 0.5)
  {
    Float_t error1 = 1 - ratio->GetBinContent(sigmavertex->GetXaxis()->FindBin(sigma)) / ratio->GetBinContent(sigmavertex->GetXaxis()->FindBin(sigma - 0.5));
    Float_t error2 = -1 + ratio->GetBinContent(sigmavertex->GetXaxis()->FindBin(sigma)) / ratio->GetBinContent(sigmavertex->GetXaxis()->FindBin(sigma + 0.5));
    Float_t cont = -1 + contamination->GetBinContent(sigmavertex->GetXaxis()->FindBin(sigma));
    Printf("%.2f sigma --> syst. error = - %.2f %% + %.2f %%, cont. = %.2f %%", sigma, error1 * 100, error2 * 100, cont * 100);
  }

  TCanvas* canvas = new TCanvas("Sigma2VertexSimulation", "Sigma2VertexSimulation", 1000, 500);
  canvas->Divide(2, 1);

  canvas->cd(1);
  sigmavertex->SetMarkerStyle(21);
  sigmavertex->Draw("P");

  canvas->cd(2);
  ratio->SetMarkerStyle(21);
  ratio->DrawCopy("P");

  contamination->DrawCopy("SAME");

  TH1F** returnContainer = new TH1F*[2];
  returnContainer[0] = ratio;
  returnContainer[1] = contamination;

  return returnContainer;
}

void Sigma2VertexCompare(const char* fileName = "systematics.root")
{
  TH1F* ratio1 = Sigma2VertexGaussian();

  TH1F** hists = Sigma2VertexSimulation(fileName);
  TH1F* ratio2 = hists[0];
  TH1F* contamination = hists[1];

  ratio1->SetStats(kFALSE);
  ratio2->SetStats(kFALSE);

  ratio1->SetMarkerStyle(0);
  ratio2->SetMarkerStyle(0);

  ratio1->SetLineWidth(2);
  ratio2->SetLineWidth(2);

  TLegend* legend = new TLegend(0.6, 0.8, 0.95, 0.95);
  legend->SetFillColor(0);
  legend->AddEntry(ratio1, "Gaussian");
  legend->AddEntry(ratio2, "Simulation");
  legend->AddEntry(contamination, "1 + Contamination");

  ratio2->SetTitle("");
  ratio2->GetYaxis()->SetTitleOffset(1.5);
  ratio2->GetXaxis()->SetRangeUser(2, 5);

  TCanvas* canvas = new TCanvas("Sigma2VertexCompare", "Sigma2VertexCompare", 500, 500);
  InitPad();

  ratio2->SetMarkerStyle(21);
  ratio1->SetMarkerStyle(22);

  ratio2->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratio2->SetLineColor(kRed);
  ratio2->SetMarkerColor(kRed);
  ratio2->Draw("PL");
  ratio1->Draw("SAMEPL");

  contamination->Draw("SAME");

  legend->Draw();

  canvas->SaveAs("Sigma2VertexCompare.eps");
}
