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

  for (Int_t i=0; i<6; i++) {
    hRatios[i] = (TH1F*)fin->Get(Form("ratio_%d",i));
		
    hRatios[i]->SetLineWidth(2);
    hRatios[i]->SetLineColor(colors[i]);				 
    hRatios[i]->SetMarkerStyle(22);
    hRatios[i]->SetMarkerSize(0.8);
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
    text[i] = new TLatex(1.08,hRatios[i]->GetBinContent(hRatios[i]->FindBin(0.9))-0.002,str.Data());
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
