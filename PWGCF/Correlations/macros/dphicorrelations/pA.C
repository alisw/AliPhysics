// 0       1       2       3       4  5  6 	        7     8             9                            10         11
// nsyield,asyield,nswidth,aswidth,v2,v3,nsasyieldratio,v3/v2,remainingpeak,remainingjetyield/ridgeyield,chi2(v2v3),baseline
const char* graphTitles[] = { "NS Ridge Yield", "AS Ridge Yield", "NS Width", "AS Width", "v2", "v3", "NS Yield / AS Yield", "v3 / v2", "remaining peak in %", "remaining jet / NS yield in %", "chi2/ndf", "baseline" };
const Int_t NGraphs = 6; // pt index
const Int_t NHists = 12; 
TGraphErrors*** graphs = 0;

const char* kCorrFuncTitle = "#frac{1}{#it{N}_{trig}} #frac{d^{2}#it{N}_{assoc}}{d#Delta#etad#Delta#varphi} (rad^{-1})";
// const char* kCorrFuncTitle = "1/N_{trig} dN_{assoc}/d#Delta#etad#Delta#varphi (1/rad)";
const char* kProjYieldTitlePhi = "1/#it{N}_{trig} d#it{N}_{assoc}/d#Delta#varphi per #Delta#eta (rad^{-1})";
const char* kProjYieldTitleEta = "1/#it{N}_{trig} d#it{N}_{assoc}/d#Delta#eta per #Delta#varphi (rad^{-1})";
// const char* kProjYieldTitlePhiOrEta = "1/N_{trig} dN_{assoc}/d#Delta#varphi (1/rad) , dN_{assoc}/d#Delta#eta";
const char* kProjYieldTitlePhiOrEta = "#frac{1}{#it{N}_{trig}} #frac{d#it{N}_{assoc}}{d#Delta#varphi} (rad^{-1}) , d#it{N}_{assoc}/d#Delta#eta";

Float_t etaMax = 1.8;

void CreateGraphStructure()
{
  graphs = new TGraphErrors**[NGraphs];
  for (Int_t i=0; i<NGraphs; i++)
  {
    graphs[i] = new TGraphErrors*[NHists];
    for (Int_t j=0; j<NHists; j++)
    {
      graphs[i][j] = new TGraphErrors;
      graphs[i][j]->GetYaxis()->SetTitle(graphTitles[j]);
    }
  }
}

void WriteGraphs(const char* outputFileName = "graphs.root")
{
  TFile::Open(outputFileName, "RECREATE");
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
    {
      graphs[i][j]->GetYaxis()->SetTitle(graphTitles[j]);
      graphs[i][j]->Write(Form("graph_%d_%d", i, j));
    }

  gFile->Close();
}

void ReadGraphs(const char* fileName = "graphs.root")
{
  CreateGraphStructure();
  TFile* file = TFile::Open(fileName);
  for (Int_t i=0; i<NGraphs; i++)
    for (Int_t j=0; j<NHists; j++)
      graphs[i][j] = (TGraphErrors*) file->Get(Form("graph_%d_%d", i, j));
}

void AddPoint(TGraphErrors* graph, Float_t x, Float_t y, Float_t xe, Float_t ye)
{
	graph->SetPoint(graph->GetN(), x, y);
	graph->SetPointError(graph->GetN() - 1, xe, ye);
}

void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text, Float_t textSize = 0.06)
{
  TLatex* latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextSize(textSize);
  latex->SetTextColor(color);
  latex->Draw();
}

void PadFor2DCorr()
{
  gPad->SetPad(0, 0, 1, 1);
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetTheta(61.62587);
  gPad->SetPhi(45);
}

void DivideGraphs(TGraphErrors* graph1, TGraphErrors* graph2)
{
//   graph1->Print();
//   graph2->Print();

  graph1->Sort();
  graph2->Sort();

  for (Int_t bin1 = 0; bin1 < graph1->GetN(); bin1++)
  {
    Float_t x = graph1->GetX()[bin1];

    Int_t bin2 = 0;
    for (bin2 = 0; bin2<graph2->GetN(); bin2++)
      if (graph2->GetX()[bin2] >= x)
        break;

    if (bin2 == graph2->GetN())
            bin2--;

    if (bin2 > 0)
      if (TMath::Abs(graph2->GetX()[bin2-1] - x) < TMath::Abs(graph2->GetX()[bin2] - x))
        bin2--;

    if (graph1->GetY()[bin1] == 0 || graph2->GetY()[bin2] == 0 || bin2 == graph2->GetN())
    {
      Printf("%d %d removed", bin1, bin2);
      graph1->RemovePoint(bin1--);
      continue;
    }

    Float_t graph2Extrapolated = graph2->GetY()[bin2];
    if (TMath::Abs(x - graph2->GetX()[bin2]) > 0.0001)
    {
      Printf("%f %f %d %d not found", x, graph2->GetX()[bin2], bin1, bin2);
      graph1->RemovePoint(bin1--);
      continue;
    }

    Float_t value = graph1->GetY()[bin1] / graph2Extrapolated;
    Float_t error = value * TMath::Sqrt(TMath::Power(graph1->GetEY()[bin1] / graph1->GetY()[bin1], 2) + TMath::Power(graph2->GetEY()[bin2] / graph2->GetY()[bin2], 2));

    graph1->GetY()[bin1] = value;
    graph1->GetEY()[bin1] = error;

//     Printf("%d %d %f %f %f %f", bin1, bin2, x, graph2Extrapolated, value, error);
  }
}

void GraphShiftX(TGraphErrors* graph, Float_t offset)
{
  for (Int_t i=0; i<graph->GetN(); i++)
    graph->GetX()[i] += offset;
}

TGraphErrors* ReadHepdata(const char* fileName, Bool_t errorsAreAdded = kFALSE, Int_t skipYErrors = 0, Int_t skipXerrors = 1)
{
  // expected format: x [x2] y [ye] [ye2] [xe]
  //
  // skipYErrors:   0 --> ye present
  //                1 --> no errors ye
  //                2 --> y and ye are lower and upper error, i.e. y' = (y + ye) / 2 and ye = (ye - y) / 2
  //                3 --> ye and ye2 are stat and syst error, will be added in quadrature
  // 
  // skipXerrors:   0 --> xe present
  //                1 --> no errors xe
  //                2 --> x2 present, xe not present and is calculated from x2 - x
  
  ifstream fin(fileName);

  graph = new TGraphErrors(0);

  Double_t sum = 0;

  while (fin.good())
  {
    char buffer[2000];
    if (fin.peek() == '#')
    {
      fin.getline(buffer, 2000);
      continue;
    }
  
    Double_t x = -1;
    Double_t x2 = -1;
    Double_t y = -1;
    Double_t ye = 0;
    Double_t xe = 0;

    fin >> x;
    
    if (skipXerrors == 2)
    {
      fin >> x2;
      xe = (x2 - x + 1) / 2;
      x = x + (x2 - x) / 2;
    }
    
    fin >> y;

    if (y == -1)
      continue;

    if (skipYErrors == 0)
    {
      ye = -1;
      fin >> ye;
      if (ye == -1)
        continue;
    }
    else if (skipYErrors == 2)
    {
      ye = -1;
      fin >> ye;
      if (ye == -1)
        continue;
      
      Double_t newy = (y + ye) / 2;
      ye = (ye - y) / 2;
      y = newy;
    }
    else if (skipYErrors == 3)
    {
      ye = -1;
      fin >> ye;
      if (ye == -1)
        continue;
      
      Double_t ye2 = -1;
      fin >> ye2;
      if (ye2 == -1)
        continue;

      ye = TMath::Sqrt(ye*ye + ye2*ye2);
    }

    if (skipXerrors == 0)
    {
      xe = -1;
      fin >> xe;
      if (xe == -1)
        continue;
    }

    //Printf("%f %f %f %f", x, y, xe, ye);

    if (errorsAreAdded)
      ye -= y;

    graph->SetPoint(graph->GetN(), x, y);
    graph->SetPointError(graph->GetN()-1, xe, ye);

    sum += y;
    
    // read rest until end of line...
    fin.getline(buffer, 2000);
  }
  fin.close();

  Printf("%s: %f", fileName, sum);

  return graph;
}

TH2* SubtractEtaGapNS(TH2* hist, Float_t etaLimit, Float_t outerLimit, Bool_t drawEtaGapDist = kFALSE)
{
  TString histName(hist->GetName());
  Int_t etaBins = 0;

  TH1D* etaGap = hist->ProjectionX(histName + "_1", TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)), hist->GetYaxis()->FindBin(-etaLimit - 0.01));
//   Printf("%f", etaGap->GetEntries());
  if (etaGap->GetEntries() > 0)
    etaBins += hist->GetYaxis()->FindBin(-etaLimit - 0.01) - TMath::Max(1, hist->GetYaxis()->FindBin(-outerLimit + 0.01)) + 1;

  TH1D* tracksTmp = hist->ProjectionX(histName + "_2", hist->GetYaxis()->FindBin(etaLimit + 0.01), TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)));
//   Printf("%f", tracksTmp->GetEntries());
  if (tracksTmp->GetEntries() > 0)
    etaBins += TMath::Min(hist->GetYaxis()->GetNbins(), hist->GetYaxis()->FindBin(outerLimit - 0.01)) - hist->GetYaxis()->FindBin(etaLimit + 0.01) + 1;
  
  etaGap->Add(tracksTmp);

  // get per bin result
  if (etaBins > 0)
    etaGap->Scale(1.0 / etaBins);
 
  for (Int_t i=1; i<=etaGap->GetNbinsX()/2; i++)
  {
//     Printf("%d -> %d", i, etaGap->GetNbinsX()+1-i);
    etaGap->SetBinContent(etaGap->GetNbinsX()+1-i, etaGap->GetBinContent(i));
    etaGap->SetBinError(etaGap->GetNbinsX()+1-i, etaGap->GetBinError(i));
  }
  
  if (drawEtaGapDist)
  {
    TH1D* centralRegion = hist->ProjectionX(histName + "_3", hist->GetYaxis()->FindBin(-etaLimit + 0.01), hist->GetYaxis()->FindBin(etaLimit - 0.01));
    
//    centralRegion->Scale(1.0 / (hist->GetYaxis()->FindBin(etaLimit - 0.01) - hist->GetYaxis()->FindBin(-etaLimit + 0.01) + 1));
    centralRegion->Scale(hist->GetXaxis()->GetBinWidth(1));

    TCanvas* c = new TCanvas("SubtractEtaGap", "SubtractEtaGap", 800, 800);
    gPad->SetLeftMargin(0.13);
    centralRegion->SetStats(0);
    TString label(centralRegion->GetTitle());
    label.ReplaceAll(".00", " GeV/#it{c}");
    label.ReplaceAll(".0", " GeV/#it{c}");
    centralRegion->SetTitle(label);
    centralRegion->SetLineColor(3);
    centralRegion->Draw();
//     centralRegion->GetYaxis()->SetTitle(kProjYieldTitlePhi);
    centralRegion->GetYaxis()->SetTitleOffset(1.6);
    TH1* copy = etaGap->DrawCopy("SAME");
    copy->Scale(hist->GetXaxis()->GetBinWidth(1));
    copy->Scale((hist->GetYaxis()->FindBin(etaLimit - 0.01) - hist->GetYaxis()->FindBin(-etaLimit + 0.01) + 1));
    copy->SetLineColor(2);
    TLegend* legend = new TLegend(0.41, 0.73, 0.69, 0.85);
    legend->SetFillColor(0);
    legend->AddEntry(centralRegion, Form("|#Delta#eta| < %.1f", etaLimit), "L");
    legend->AddEntry(copy, Form("%.1f < |#Delta#eta| < %.1f (scaled)", etaLimit, outerLimit), "L");
    legend->Draw();
    
//     DrawLatex(0.705, 0.62, 1, "Pb-Pb 2.76 TeV", 0.025);
//     DrawLatex(0.705, 0.58, 1, "Stat. unc. only", 0.025);
  }
  
//   new TCanvas; etaGap->DrawCopy();
  
  TH2* histTmp2D = (TH2*) hist->Clone("histTmp2D");
  histTmp2D->Reset();
  
  for (Int_t xbin=1; xbin<=histTmp2D->GetNbinsX(); xbin++)
    for (Int_t y=1; y<=histTmp2D->GetNbinsY(); y++)
      histTmp2D->SetBinContent(xbin, y, etaGap->GetBinContent(xbin));
    
  hist->Add(histTmp2D, -1);  
  return histTmp2D;
}

TH1* GetProjections(Int_t i, Int_t j, Int_t centr, char** label, Float_t etaBegin = 1.0, TH1** etaProj = 0)
{
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, centr));
  if (!hist1)
    return 0;
  hist1 = (TH2*) hist1->Clone(Form("%s_%.1f", hist1->GetName(), etaBegin));

  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
  
  hist1->Rebin2D(2, 1); hist1->Scale(0.5);
  
//   new TCanvas; hist1->Draw("surf1");
  
  if (etaBegin > 0)
    SubtractEtaGapNS(hist1, etaBegin, etaMax, kTRUE);
  
  tokens = TString(hist1->GetTitle()).Tokenize("-");
  centralityStr = new TString;
  if (tokens->GetEntries() > 2)
    *centralityStr = tokens->At(2)->GetName();
  if (tokens->GetEntries() > 3)
    *centralityStr = *centralityStr + "-" + tokens->At(3)->GetName();
  *label = centralityStr->Data();
//   Printf("%s", label);
  
  proj1x = ((TH2*) hist1)->ProjectionX(Form("proj1x_%d_%d_%d_%.1f", i, j, centr, etaBegin), hist1->GetYaxis()->FindBin(-etaMax+0.01), hist1->GetYaxis()->FindBin(etaMax-0.01));
  proj1x->GetXaxis()->SetTitleOffset(1);
  proj1x->Scale(hist1->GetYaxis()->GetBinWidth(1));

  proj1y = ((TH2*) hist1)->ProjectionY(Form("proj1y_%d_%d_%d_%.1f", i, j, centr, etaBegin), hist1->GetXaxis()->FindBin(-0.5), hist1->GetXaxis()->FindBin(0.5));
  proj1y->Scale(hist1->GetXaxis()->GetBinWidth(1));
  proj1y->GetXaxis()->SetTitleOffset(1);
  proj1y->GetXaxis()->SetRangeUser(-1.99, 1.99);
  Float_t etaPhiScale = 1.0 * (hist1->GetXaxis()->FindBin(0.5) - hist1->GetXaxis()->FindBin(-0.5) + 1) / (hist1->GetYaxis()->FindBin(etaMax-0.01) - hist1->GetYaxis()->FindBin(-etaMax+0.01) + 1);
  
  if (gStudySystematic == 20)
  {
    Printf(">>>>>>>>>>>> Applying non-closure systematics <<<<<<<<<<<<");
    file2 = TFile::Open("non_closure.root");
    non_closure = (TH1*) file2->Get(Form("non_closure_all_%d_%d_%d", i, j, 0));
    for (Int_t bin=1; bin<=non_closure->GetNbinsX(); bin++)
      non_closure->SetBinError(bin, 0);
    
    proj1x->Multiply(non_closure);  
  }
  
  Float_t zyam = 0;
  if (0)
  {  
    clone = (TH1*) proj1x->Clone();
    clone->Fit("pol0", "0", "", TMath::Pi()/2 - 0.2, TMath::Pi()/2);
    zyam = clone->GetFunction("pol0")->GetParameter(0);
  }
  else
  {
//     zyam = (proj1x->GetBinContent(proj1x->FindBin(TMath::Pi()/2)) + proj1x->GetBinContent(proj1x->FindBin(-TMath::Pi()/2))) / 2;
//     zyam = proj1x->GetBinContent(proj1x->FindBin(TMath::Pi()/2));
    zyam = proj1x->GetBinContent(proj1x->FindBin(1.3));
//     zyam = proj1x->GetMinimum();
  }
    
  proj1x->Add(new TF1("func", "-1", -100, 100), zyam);
  proj1y->Add(new TF1("func", "-1", -100, 100), zyam * etaPhiScale);
  
  proj1x->Scale(1.0 / (2.0 * etaMax));
  
  if (etaProj != 0)
    *etaProj = proj1y;
  return proj1x;
}

TH1* GetProjectionsNew(Int_t i, Int_t j, Int_t centr, char** label, Float_t etaBegin = 1.0, TH1** etaProj = 0)
{
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, centr));
  if (!hist1)
    return 0;
  hist1 = (TH2*) hist1->Clone(Form("%s_%.1f", hist1->GetName(), etaBegin));

  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
  
  hist1->Rebin2D(2, 1); hist1->Scale(0.5);
  
//   new TCanvas; hist1->Draw("surf1");
  
  tokens = TString(hist1->GetTitle()).Tokenize("-");
  centralityStr = new TString;
  if (tokens->GetEntries() > 2)
    *centralityStr = tokens->At(2)->GetName();
  if (tokens->GetEntries() > 3)
    *centralityStr = *centralityStr + "-" + tokens->At(3)->GetName();
  *label = centralityStr->Data();
//   Printf("%s", label);
  
  if (etaBegin > 0)
  {
    proj1x = ((TH2*) hist1)->ProjectionX(Form("proj1x_%d_%d_%d_%.1f", i, j, centr, etaBegin), hist1->GetYaxis()->FindBin(-etaMax+0.01), hist1->GetYaxis()->FindBin(etaMax-0.01));
    proj1x->GetXaxis()->SetTitleOffset(1);
    proj1x->Scale(hist1->GetYaxis()->GetBinWidth(1));

    proj1xR1 = ((TH2*) hist1)->ProjectionX(Form("proj2x_%d_%d_%d_%.1f", i, j, centr, etaBegin), hist1->GetYaxis()->FindBin(-etaMax+0.01), hist1->GetYaxis()->FindBin(-etaBegin-0.01));
    proj1xR2 = ((TH2*) hist1)->ProjectionX(Form("proj3x_%d_%d_%d_%.1f", i, j, centr, etaBegin), hist1->GetYaxis()->FindBin(etaBegin+0.01), hist1->GetYaxis()->FindBin(etaMax-0.01));
    proj1xR1->Add(proj1xR2);
    
    proj1xR1->GetXaxis()->SetTitleOffset(1);
    proj1xR1->Scale(hist1->GetYaxis()->GetBinWidth(1));
    
    proj1xR1->Scale((1.0 * hist1->GetYaxis()->FindBin(etaMax-0.01) - hist1->GetYaxis()->FindBin(-etaMax+0.01) + 1) / (hist1->GetYaxis()->FindBin(-etaBegin-0.01) - hist1->GetYaxis()->FindBin(-etaMax+0.01) + 1 + hist1->GetYaxis()->FindBin(etaMax-0.01) - hist1->GetYaxis()->FindBin(etaBegin+0.01) + 1));
    
    // mirror
    for (Int_t i=1; i<=proj1xR1->GetNbinsX()/2; i++)
    {
  //     Printf("%d -> %d", i, etaGap->GetNbinsX()+1-i);
      proj1xR1->SetBinContent(proj1xR1->GetNbinsX()+1-i, proj1xR1->GetBinContent(i));
      proj1xR1->SetBinError(proj1xR1->GetNbinsX()+1-i, proj1xR1->GetBinError(i));
    }
    
    proj1x->Add(proj1xR1, -1);
    
    proj1x->Scale(1.0 / (2.0 * etaMax));
  
//     new TCanvas; proj1xR1->DrawCopy();
  }
  else
  {
    proj1x = ((TH2*) hist1)->ProjectionX(Form("proj1x_%d_%d_%d_%.1f", i, j, centr, etaBegin), hist1->GetYaxis()->FindBin(-etaMax+0.01), hist1->GetYaxis()->FindBin(etaMax-0.01));
    proj1x->GetXaxis()->SetTitleOffset(1);
    proj1x->Scale(hist1->GetYaxis()->GetBinWidth(1));

    proj1x->Scale(1.0 / (2.0 * etaMax));
  }

  if (gStudySystematic == 20)
  {
    Printf(">>>>>>>>>>>> Applying non-closure systematics <<<<<<<<<<<<");
    file2 = TFile::Open("non_closure.root");
    non_closure = (TH1*) file2->Get(Form("non_closure_all_%d_%d_%d", i, j, 0));
    for (Int_t bin=1; bin<=non_closure->GetNbinsX(); bin++)
      non_closure->SetBinError(bin, 0);
    
    proj1x->Multiply(non_closure);  
  }
  
  Float_t zyam = 0;
  if (0)
  {  
    clone = (TH1*) proj1x->Clone();
    clone->Fit("pol0", "0", "", TMath::Pi()/2 - 0.2, TMath::Pi()/2);
    zyam = clone->GetFunction("pol0")->GetParameter(0);
  }
  else
  {
//     zyam = (proj1x->GetBinContent(proj1x->FindBin(TMath::Pi()/2)) + proj1x->GetBinContent(proj1x->FindBin(-TMath::Pi()/2))) / 2;
//     zyam = proj1x->GetBinContent(proj1x->FindBin(TMath::Pi()/2));
    zyam = proj1x->GetBinContent(proj1x->FindBin(1.3));
//     zyam = proj1x->GetMinimum();
  }
    
  proj1x->Add(new TF1("func", "-1", -100, 100), zyam);
  
  return proj1x;
}

void DrawEtaDep(const char* fileName, Int_t i, Int_t j, Int_t centr)
{
  c = new TCanvas;
  gPad->SetGridx();
  gPad->SetGridy();
  
  TFile::Open(fileName);

  legend = new TLegend(0.8, 0.6, 0.99, 0.99);
  legend->SetFillColor(0);

  for (Int_t n=0; n<4; n++)
  {
    Float_t eta = 0.8 + n * 0.2;
    const char* label = 0;
    TH1* hist = GetProjections(i, j, centr, &label, eta);
    hist->SetStats(0);	
    hist->SetLineColor(n+1);
    hist->Draw((n == 0) ? "" : "SAME");
    legend->AddEntry(hist, Form("|#Delta#eta| > %.1f", eta));
  }
  legend->Draw();
}

void DrawProjectionsTim(const char* fileName, Int_t i, Int_t j, Float_t eta = 1.0, Bool_t etaPhi = kFALSE)
{
    gStyle->SetErrorX(0.0);
    c = new TCanvas;
//    gPad->SetGridx();
//    gPad->SetGridy();
    gPad->SetTopMargin(0.025);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);

    TFile::Open(fileName);
    
//    TLegend *legend = new TLegend(0.45, 0.60, 0.65, 0.90);
    TLegend *legend = new TLegend(0.45, 0.45, 0.65, 0.90);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(62);
    legend->SetTextSize(0.04);
    legend->SetFillStyle(0);
    
//    TLegend *legend2 = new TLegend(0.65, 0.75, 0.85, 0.90);
//    legend2->SetFillColor(0);
//    legend2->SetBorderSize(0);
//    legend2->SetTextFont(62);
//    legend2->SetTextSize(0.04);
//    legend2->SetFillStyle(0);
    
//     Int_t colors[] = { 1, 2, 1, 4, 6, 2 };
    Int_t colors[] = { kRed+1, kOrange+7, kBlack, kGreen+2, kAzure+2, kBlack };
    Int_t markers[] = { 20, 21, 1, 22, 23, 1 };
    
    TH1* first = 0;
    
    Float_t min = 100;
    Float_t max = -100;
    
    Int_t centSeq[] = { 0, 1, 3, 4, 2, 5 };
    
    for (Int_t otcentr=0; otcentr<6; otcentr++)
    {
        Int_t centr = centSeq[otcentr];
        /*    if (centr >= 5)
         continue;*/
        const char* label = 0;
        TH1* etaHist = 0;
        TH1* hist = GetProjectionsNew(i, j, centr, &label, (centr == 2 || centr == 5 || centr == 4) ? -1 : eta, &etaHist);
        if (!hist)
            continue;
        if (etaPhi)
            hist = etaHist;
        hist->SetStats(0);
        hist->GetXaxis()->CenterTitle();
        hist->GetXaxis()->SetLabelSize(0.05);
        hist->GetXaxis()->SetTitleSize(0.05);
        hist->GetXaxis()->SetTitleOffset(1.1);
        hist->GetXaxis()->SetTitle("#Delta#varphi (rad)");
        hist->GetYaxis()->SetNdivisions(506);
        hist->GetYaxis()->CenterTitle();
        hist->GetYaxis()->SetLabelSize(0.05);
        hist->GetYaxis()->SetTitleSize(0.05);
        hist->GetYaxis()->SetTitleOffset(1.35);
	
	if (centr == 2 || centr == 5)
	{
	  hist->SetLineWidth(2);
	  if (centr == 2)
	    hist->SetLineStyle(2);
	}
	  
//         hist->GetYaxis()->SetTitle("#frac{1}{#it{N}_{trig}}#frac{d#it{N}_{assoc}}{d#Delta#varphi} (rad^{-1})");
	hist->GetYaxis()->SetTitle("1/#it{N}_{trig} d#it{N}_{assoc}/d#Delta#varphi per #Delta#eta - const (rad^{-1})");

        hist->SetLineColor(colors[centr]);
        hist->SetMarkerColor(colors[centr]);
        hist->SetMarkerStyle(markers[centr]);
        if (etaPhi)
            hist->GetXaxis()->SetRangeUser(-1.79, 1.79);
        c->cd();
        
        // -----
        tokens = TString(hist->GetTitle()).Tokenize("-");
        sPtTRange = new TString;
        *sPtTRange = tokens->At(0)->GetName();
        *sPtTRange = *sPtTRange + "GeV/#it{c}";
        sPtTRange->ReplaceAll(".0", "");
        sPtARange = new TString;
        *sPtARange = tokens->At(1)->GetName();
        *sPtARange = *sPtARange + "GeV/#it{c}";
        sPtARange->ReplaceAll(".00", "");
        sPtARange->ReplaceAll(".0", "");
        sPtARange->ReplaceAll(" 1", "1");
        sPtTRange->ReplaceAll("p_", "#it{p}_");
        sPtARange->ReplaceAll("p_", "#it{p}_");
        cout << sPtTRange->Data() << endl;
        cout << sPtARange->Data() << endl;
        hist->SetTitle("");
        // -----
        
        if (centr == 0)
	  hist->Draw("");
	else if (centr == 2 || centr == 5)
	  hist->Draw("HISTE SAME");
	else
	  hist->Draw("SAME");
	
        min = TMath::Min(min, hist->GetMinimum());
        max = TMath::Max(max, hist->GetMaximum());
        if (!first)
            first = hist;
        if (centr==2) legend->AddEntry(hist,"pp 2.76 TeV","l");
        else if (centr==5) legend->AddEntry(hist,"pp 7 TeV","l");
        else legend->AddEntry(hist, label, "p");
        //     break;
    }
//    first->GetYaxis()->SetRangeUser(min * 1.1, max * 1.1);
//     first->GetYaxis()->SetRangeUser(min * 1.1, 0.67/3.6);
    first->GetYaxis()->SetRangeUser(-0.009, 0.2);
    cout << max*1.1 << endl;
    cout << 0.925*first->GetMaximum() << endl;
    legend->Draw();
//    legend2->Draw();
    TLine * li = new TLine(first->GetXaxis()->GetXmin(),0.0,first->GetXaxis()->GetXmax(),0.0);
    li->SetLineStyle(kDashed);
    li->SetLineColor(kBlack);
    li->Draw("same");

    TLatex * tex_Pbp = new TLatex(0.65,0.945*first->GetMaximum(),"p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV");
    tex_Pbp->SetTextFont(62);
    tex_Pbp->SetTextSize(0.05);
    tex_Pbp->Draw();
    
    TLatex * tex_statu = new TLatex(2.3,0.65*first->GetMaximum(),"stat. uncertainties only");
    tex_statu->SetTextFont(62);
    tex_statu->SetTextSize(0.04);
    tex_statu->Draw();
    
    if (eta > 0)
    {
      TLatex * tex_statu = new TLatex(2.3,0.56*first->GetMaximum(),"ridge subtracted");
      tex_statu->SetTextFont(62);
      tex_statu->SetTextSize(0.04);
      tex_statu->Draw();
    }

    TLatex * tex_PtT = new TLatex(2.3,0.85*first->GetMaximum(),sPtTRange->Data());
    tex_PtT->SetTextFont(62);
    tex_PtT->SetTextSize(0.04);
    tex_PtT->Draw();
    
    TLatex * tex_PtA = new TLatex(2.3,0.75*first->GetMaximum(),sPtARange->Data());
    tex_PtA->SetTextFont(62);
    tex_PtA->SetTextSize(0.04);
    tex_PtA->Draw();
    
//     box = new TBox(-1.4, 0.185, -1.3, 0.195);
//     box->SetFillColor(kGray);
//     box->Draw();

    if (i == 2 && j == 2 && eta < 0)
    {
      c->SaveAs("fig2.eps");
      c->SaveAs("fig2.png");
    }
    else if (i == 2 && j == 2 && eta > 0)
    {
      c->SaveAs("fig5.eps");
      c->SaveAs("fig5.png");
    }
    
//     c->SaveAs(Form("%s_%d_%d.png", (etaPhi) ? "eta" : "phi", i, j));
//     c->SaveAs(Form("%s_%d_%d.eps", (etaPhi) ? "eta" : "phi", i, j));
    
    return;
    
    c = new TCanvas;
//    gPad->SetGridx();
//    gPad->SetGridy();
    
    ppHist = GetProjections(i, j, 2, &label);
    
    for (Int_t centr=0; centr<6; centr++)
    {
        if (centr >= 4 || centr == 2)
            continue;
        const char* label = 0;
        TH1* hist = GetProjections(i, j, centr, &label);
        hist->SetStats(0);
        hist->SetLineColor(centr+1);
        hist->Divide(ppHist);
        hist->GetYaxis()->SetRangeUser(1, 3);
        hist->Draw((centr == 0) ? "" : "SAME");
    }
    legend->Draw();
    
    c->SaveAs(Form("phi_%d_%d_ratio.png", i, j));
    c->SaveAs(Form("phi_%d_%d_ratio.eps", i, j));
}


void DrawProjections(const char* fileName, Int_t i, Int_t j, Float_t eta = 1.0, Bool_t etaPhi = kFALSE)
{
  c = new TCanvas;
  gPad->SetGridx();
  gPad->SetGridy();
  
  TFile::Open(fileName);

  legend = new TLegend(0.8, 0.6, 0.99, 0.99);
  legend->SetFillColor(0);
  
  Int_t colors[] = { 1, 2, 1, 4, 6, 2 };
  Int_t markers[] = { 1, 1, 5, 1, 1, 5 };
  
  TH1* first = 0;
  
  Float_t min = 100;
  Float_t max = -100;
  
  for (Int_t centr=0; centr<6; centr++)
  {
/*    if (centr >= 5)
      continue;*/
    const char* label = 0;
    TH1* etaHist = 0;
    TH1* hist = GetProjectionsNew(i, j, centr, &label, eta, &etaHist);
    if (!hist)
      continue;
    if (etaPhi)
      hist = etaHist;
    hist->SetStats(0);
    hist->SetLineColor(colors[centr]);
    hist->SetMarkerColor(colors[centr]);
    hist->SetMarkerStyle(markers[centr]);
    if (etaPhi)
      hist->GetXaxis()->SetRangeUser(-1.79, 1.79);
    c->cd();
    hist->Draw((centr == 0) ? "" : "SAME");
    min = TMath::Min(min, hist->GetMinimum());
    max = TMath::Max(max, hist->GetMaximum());
    if (!first)
      first = hist;
    legend->AddEntry(hist, label);
//     break;
  }
  first->GetYaxis()->SetRangeUser(min * 1.1, max * 1.1);
  legend->Draw();
  
  c->SaveAs(Form("%s_%d_%d.png", (etaPhi) ? "eta" : "phi", i, j));
  
  return;
  
  c = new TCanvas;
  gPad->SetGridx();
  gPad->SetGridy();
  
  ppHist = GetProjections(i, j, 2, &label);
  
  for (Int_t centr=0; centr<6; centr++)
  {
    if (centr >= 4 || centr == 2)
      continue;
    const char* label = 0;
    TH1* hist = GetProjections(i, j, centr, &label);
    hist->SetStats(0);
    hist->SetLineColor(centr+1);
    hist->Divide(ppHist);
    hist->GetYaxis()->SetRangeUser(1, 3);
    hist->Draw((centr == 0) ? "" : "SAME");
  }
  legend->Draw();
  
  c->SaveAs(Form("phi_%d_%d_ratio.png", i, j));
}

void CompareProjections(const char* fileName, Int_t i, Int_t j, Int_t centr)
{
  TFile::Open(fileName);

  const char* label = 0;
  TH1* etaHist = 0;
  TH1* hist = GetProjections(i, j, centr, &label, -1, &etaHist);
  if (!hist)
    continue;
  TH1* hist2 = GetProjections(i, j, centr, &label, 1.2, &etaHist);

  c = new TCanvas;
  gPad->SetGridx();
  gPad->SetGridy();
  
  hist->SetStats(0);
  hist->SetMarkerStyle(20);
  hist->Draw("");

  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(2);
  hist2->SetLineColor(2);
  hist2->Draw("SAME");
}

void DrawProjectionsAll(const char* fileName, Float_t eta = 1, Bool_t etaPhi = kFALSE)
{
  DrawProjections(fileName, 0, 1, eta, etaPhi);
  DrawProjections(fileName, 0, 2, eta, etaPhi);
  DrawProjections(fileName, 1, 1, eta, etaPhi);
  DrawProjections(fileName, 1, 2, eta, etaPhi);
  DrawProjections(fileName, 1, 3, eta, etaPhi);
  DrawProjections(fileName, 2, 3, eta, etaPhi);
}

void DrawProjectionsRidge(const char* fileName, Int_t i, Int_t j, Int_t centr1, Int_t centr2)
{
  Float_t etaLimit = 1.0;
  Float_t outerLimit = 1.6;
  
  TFile::Open(fileName);
  
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, centr1));
  TH2* hist2 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, centr2));
  
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
  hist2->Scale(1.0 / hist2->GetYaxis()->GetBinWidth(1));

  SubtractEtaGapNS(hist1, etaLimit, outerLimit, kFALSE);
  
  proj1y = ((TH2*) hist1)->ProjectionY("proj1y", hist1->GetXaxis()->FindBin(-0.5), hist1->GetXaxis()->FindBin(0.5));
  proj1x = ((TH2*) hist1)->ProjectionX("proj1x", hist1->GetYaxis()->FindBin(-1.79), hist1->GetYaxis()->FindBin(1.79));
  
  Float_t etaPhiScale = 1.0 * (hist1->GetXaxis()->FindBin(0.5) - hist1->GetXaxis()->FindBin(-0.5) + 1) / (hist1->GetYaxis()->FindBin(1.79) - hist1->GetYaxis()->FindBin(-1.79) + 1);
  Printf("%f", etaPhiScale);
  
  proj1y->GetXaxis()->SetTitleOffset(1);
  proj1x->GetXaxis()->SetTitleOffset(1);
  
  clone = (TH1*) proj1x->Clone();
  clone->Fit("pol0", "0", "", 1.2, 1.8);
  Float_t zyam = clone->GetFunction("pol0")->GetParameter(0);

  proj1x->Add(new TF1("func", "-1", -100, 100), zyam);
  proj1y->Add(new TF1("func", "-1", -100, 100), zyam * etaPhiScale);
  
  new TCanvas("c", "c", 800, 800);
  gPad->SetLeftMargin(0.15);
//   hist1->SetTitle("");

  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist1->GetXaxis()->SetTitleOffset(1.5);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->SetStats(kFALSE);
  hist1->DrawCopy("SURF1");
  
  proj2y = ((TH2*) hist2)->ProjectionY("proj2y", hist1->GetXaxis()->FindBin(-0.5), hist1->GetXaxis()->FindBin(0.5));
  proj2x = ((TH2*) hist2)->ProjectionX("proj2x", hist1->GetYaxis()->FindBin(-1.79), hist1->GetYaxis()->FindBin(1.79));

//   proj2y->Scale(1.0 / (hist1->GetXaxis()->FindBin(0.5) - hist1->GetXaxis()->FindBin(-0.5) + 1));
//   proj2x->Scale(1.0 / (hist1->GetYaxis()->FindBin(1.79) - hist1->GetYaxis()->FindBin(-1.79) + 1));
 
  clone = (TH1*) proj2x->Clone();
  clone->Fit("pol0", "0", "", 1.2, 1.8);
  zyam = clone->GetFunction("pol0")->GetParameter(0);

  proj2x->Add(new TF1("func", "-1", -100, 100), zyam);
  proj2y->Add(new TF1("func", "-1", -100, 100), zyam * etaPhiScale);

  proj2y->SetLineColor(2); proj2x->SetLineColor(2);
  
  new TCanvas; proj1y->Draw(); proj2y->Draw("SAME"); gPad->SetGridx(); gPad->SetGridy();
  new TCanvas; proj1x->Draw(); proj2x->Draw("SAME"); gPad->SetGridx(); gPad->SetGridy();
  
  new TCanvas("c2", "c2", 800, 800);
  gPad->SetLeftMargin(0.15);
//   hist2->SetTitle("");
  hist2->GetYaxis()->SetRangeUser(-1.79, 1.79);
  hist2->GetXaxis()->SetTitleOffset(1.5);
  hist2->GetYaxis()->SetTitleOffset(2);
  hist2->SetStats(kFALSE);
  hist2->DrawCopy("SURF1");
}

void CalculateIAA(Int_t i, Int_t j, Int_t centr, Float_t etaBegin, Double_t& nsPeak, Double_t& nsPeakE, Double_t& asPeak, Double_t& asPeakE, Double_t& nsRidge, Double_t& nsRidgeE, char** label, Bool_t subtractRidge)
{
  TH2* hist = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, centr));
  if (!hist)
  {
    nsPeak = 0;
    nsPeakE = 0;
    asPeak = 0;
    asPeakE = 0;
    nsRidge = 0;
    nsRidgeE = 0;
    return;
  }
  
//   new TCanvas; hist->Draw("COLZ");
  
  TH2* hist1 = (TH2*) hist->Clone();
  hist1->Reset();
  
  // copy to one quadrant
  for (Int_t x=1; x<=hist->GetNbinsX(); x++)
    for (Int_t y=1; y<=hist->GetNbinsY(); y++)
    {
      Int_t xTarget = x;
      if (hist->GetXaxis()->GetBinCenter(x) < 0)
	xTarget = hist->GetXaxis()->FindBin(-1.0 * hist->GetXaxis()->GetBinCenter(x));
      else if (hist->GetXaxis()->GetBinCenter(x) > TMath::Pi())
	xTarget = hist->GetXaxis()->FindBin(TMath::TwoPi() - hist->GetXaxis()->GetBinCenter(x));
      
      Int_t yTarget = y;
      if (hist->GetYaxis()->GetBinCenter(y) < 0)
	yTarget = hist->GetYaxis()->FindBin(-1.0 * hist->GetYaxis()->GetBinCenter(y));
      
//       Printf("%d %d --> %d %d", x, y, xTarget, yTarget);
      
      Float_t value = 0;
      Float_t error = 0;
      value += hist->GetBinContent(x, y);
      error += hist->GetBinError(x, y) * hist->GetBinError(x, y);

      value += hist1->GetBinContent(xTarget, yTarget);
      error += hist1->GetBinError(xTarget, yTarget) * hist1->GetBinError(xTarget, yTarget);

      error = TMath::Sqrt(error);
      
      hist1->SetBinContent(xTarget, yTarget, value);
      hist1->SetBinError(xTarget, yTarget, error);
    }
  
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
  new TCanvas; hist1->Draw("COLZ");
//   new TCanvas; hist1->Draw("SURF1");
  
  tokens = TString(hist1->GetTitle()).Tokenize("-");
  centralityStr = new TString;
  if (tokens->GetEntries() > 1)
  {
    *centralityStr = tokens->At(0)->GetName();
    *centralityStr = *centralityStr + "-" + tokens->At(1)->GetName();
  }
  *label = centralityStr->Data();

  Int_t phi1 = hist1->GetXaxis()->FindBin(0.0001);
  Int_t phi3 = hist1->GetXaxis()->FindBin(TMath::Pi()/2 - 0.3);
  Int_t phi2 = phi3 - 1;
  Int_t phi4 = hist1->GetXaxis()->FindBin(TMath::Pi()/2 - 0.1);
  Int_t phi5 = phi4 + 1;
  Int_t phi6 = hist1->GetXaxis()->FindBin(TMath::Pi() - 0.0001);
//   Printf("phi = %d %d %d %d %d %d", phi1, phi2, phi3, phi4, phi5, phi6);
  
  Int_t eta1 = hist1->GetYaxis()->FindBin(0.0001);
  Int_t eta2 = hist1->GetYaxis()->FindBin(etaBegin - 0.0001);
  Int_t eta3 = eta2 + 1;
  Int_t eta4 = hist1->GetYaxis()->FindBin(etaMax - 0.0001);
//   Printf("eta = %d %d %d %d", eta1, eta2, eta3, eta4);
  
  Double_t zyamYieldE, zyamYieldE1, zyamYieldE2;
  nsPeak              = hist1->IntegralAndError(phi1, phi2, eta1, eta2, nsPeakE, "width");
  nsRidge             = hist1->IntegralAndError(phi1, phi2, eta3, eta4, nsRidgeE, "width");
  asPeak              = hist1->IntegralAndError(phi5, phi6, eta1, eta4, asPeakE, "width");
  Double_t zyamYield  = hist1->IntegralAndError(phi3, phi4, eta1, eta4, zyamYieldE, "width");
  Double_t zyamYield1 = hist1->IntegralAndError(phi3, phi4, eta1, eta2, zyamYieldE1, "width");
  Double_t zyamYield2 = hist1->IntegralAndError(phi3, phi4, eta3, eta4, zyamYieldE2, "width");
  
  // factor 4 from folding to one quadrant above
  Double_t zyamYield2Density = zyamYield2 / (hist1->GetXaxis()->GetBinUpEdge(phi4) - hist1->GetXaxis()->GetBinLowEdge(phi3)) / (hist1->GetYaxis()->GetBinUpEdge(eta4) - hist1->GetYaxis()->GetBinLowEdge(eta3)) / 4;
  Double_t zyamYieldE2Density = zyamYieldE2 / (hist1->GetXaxis()->GetBinUpEdge(phi4) - hist1->GetXaxis()->GetBinLowEdge(phi3)) / (hist1->GetYaxis()->GetBinUpEdge(eta4) - hist1->GetYaxis()->GetBinLowEdge(eta3)) / 4;
  
  Double_t nsZyamScaling = 1.0 * (phi2 - phi1 + 1) / (phi4 - phi3 + 1);
  Double_t asZyamScaling = 1.0 * (phi6 - phi5 + 1) / (phi4 - phi3 + 1);
//   Printf("%f %f", nsZyamScaling, asZyamScaling);
  
  nsPeak -= zyamYield1 * nsZyamScaling;
  nsPeakE = TMath::Sqrt(zyamYieldE1 * zyamYieldE1 * nsZyamScaling * nsZyamScaling + nsPeakE * nsPeakE);
  
  nsRidge -= zyamYield2 * nsZyamScaling;
  nsRidgeE = TMath::Sqrt(zyamYieldE2 * zyamYieldE2 * nsZyamScaling * nsZyamScaling + nsRidgeE * nsRidgeE);

  asPeak -= zyamYield * asZyamScaling;
  asPeakE = TMath::Sqrt(zyamYieldE * zyamYieldE * asZyamScaling * asZyamScaling + asPeakE * asPeakE);
  
  if (subtractRidge)
  {
    Double_t nsRidgeScaling = 1.0 * (eta2 - eta1 + 1) / (eta4 - eta3 + 1);
    Double_t asRidgeScaling = 1.0 * (eta4 - eta1 + 1) / (eta4 - eta3 + 1);
    
    nsPeak -= nsRidge * nsRidgeScaling;
    asPeak -= nsRidge * asRidgeScaling;
    nsPeakE = TMath::Sqrt(nsRidgeE * nsRidgeE * nsRidgeScaling * nsRidgeScaling + nsPeakE * nsPeakE);
    asPeakE = TMath::Sqrt(nsRidgeE * nsRidgeE * asRidgeScaling * asRidgeScaling + asPeakE * asPeakE);
  }

  nsRidge /= (etaMax - etaBegin) * 2;
  nsRidgeE /= (etaMax - etaBegin) * 2;
  
  Printf("Peak yields (%d %d %d): %f +- %f; %f +- %f; %f +- %f; %f +- %f", i, j, centr, nsPeak, nsPeakE, asPeak, asPeakE, nsRidge, nsRidgeE, zyamYield2Density, zyamYieldE2Density);

//   Printf("");
}

void PlotIAA(const char* fileName)
{
  Int_t colors[] = { 1, 2, 3, 4, 6, 7 };
  Int_t markers[] = { 20, 21, 22, 23, 24, 25 };
  
  if (0)
  {
    Int_t n = 6;
    Int_t is[] = { 0, 0, 1, 1, 1, 2 };
    Int_t js[] = { 1, 2, 1, 2, 3, 3 };
  
    Int_t centralityBins = 6;
    Float_t centralityX[] = { 10, 30, 110, 50, 80, 120 };
    Float_t centralityEX[] = { 10, 10, 0, 10, 20, 0 };
  }
  else if (0)
  {
    Int_t n = 3;
    Int_t is[] = { 0, 1, 2, 3 };
    Int_t js[] = { 1, 2, 3, 4 };

    Int_t centralityBins = 4;
    Float_t centralityX[] = { 1.5, 6.5, 30, 75 };
    Float_t centralityEX[] = { 1.5, 2.5, 20, 25 };
  }
  else
  {
    Int_t n = 1;
    Int_t is[] = { 0 };
    Int_t js[] = { 1 };
  
    Int_t centralityBins = 6;
    Float_t centralityX[] = { 10, 30, 110, 50, 80, 120 };
    Float_t centralityEX[] = { 10, 10, 0, 10, 20, 0 };
  }
  
  const char* graphTitles[] = { "NS Yield", "AS Yield", "NS Ridge" };

  TCanvas* canvas[3];
  for (Int_t ci=0; ci<3; ci++)
  {
    canvas[ci] = new TCanvas;
    gPad->SetGridx();
    gPad->SetGridy();
  }
  
  TFile::Open(fileName);

  legend = new TLegend(0.4, 0.6, 0.99, 0.99);
  legend->SetFillColor(0);

  for (Int_t i=0; i<n; i++)
  {
    TGraphErrors* graph[5];
    for (Int_t ci=0; ci<5; ci++)
      graph[ci] = new TGraphErrors;

    char* label = 0;
    for (Int_t c=0; c<centralityBins; c++)
    {
      Double_t nsYield, nsError, asYield, asError, nsRidge, nsRidgeE;
      CalculateIAA(is[i], js[i], c, 1.2, nsYield, nsError, asYield, asError, nsRidge, nsRidgeE, &label, kTRUE);

      AddPoint(graph[0], centralityX[c], nsYield, centralityEX[c], nsError);
      AddPoint(graph[1], centralityX[c], asYield, centralityEX[c], asError);
      AddPoint(graph[2], centralityX[c], nsRidge, centralityEX[c], nsRidgeE);

//       if (c != 2 && c != 5)
      {
	CalculateIAA(is[i], js[i], c, 1.2, nsYield, nsError, asYield, asError, nsRidge, nsRidgeE, &label, kFALSE);
	AddPoint(graph[3], centralityX[c], nsYield, 0, nsError);
	AddPoint(graph[4], centralityX[c], asYield, 0, asError);
      }
    }

    for (Int_t ci=0; ci<3; ci++)
    {
      canvas[ci]->cd();
      graph[ci]->SetMarkerStyle(markers[i]);
      graph[ci]->SetMarkerColor(colors[i]);
      graph[ci]->SetLineColor(colors[i]);
      graph[ci]->Draw((i == 0) ? "AP" : "PSAME");
      graph[ci]->GetYaxis()->SetTitle(graphTitles[ci]);
      graph[ci]->GetYaxis()->SetRangeUser(0, 1);
    }
    for (Int_t ci=3; ci<5; ci++)
    {
      canvas[ci-3]->cd();
      graph[ci]->SetLineColor(colors[i]);
      graph[ci]->Sort();
      graph[ci]->Draw("LSAME");
    }
    
    legend->AddEntry(graph[0], label, "P");
  }
  
  canvas[0]->cd();
  legend->Draw();
  canvas[0]->SaveAs("ns_yield.png");

  canvas[1]->cd();
  legend->Draw();
  canvas[1]->SaveAs("as_yield.png");

  canvas[2]->cd();
  legend->Draw();
  canvas[2]->SaveAs("ns_ridge.png");
}

void Draw2D(const char* fileName, Int_t i, Int_t j, Int_t centr)
{
  TFile::Open(fileName);
  
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, centr));
  if (!hist1)
    return 0;
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
  
//   hist1->Rebin2D(2, 2); hist1->Scale(0.25);
  
  hist1->GetYaxis()->SetRangeUser(-1.79, 1.79);
  new TCanvas; 
  hist1->DrawCopy("SURF1");
}

void CorrelationSubtraction(const char* fileName, Int_t i, Int_t j, Int_t centr)
{
  CreateGraphStructure();

  CorrelationSubtraction(fileName, i, j, centr, graphs[0], 0);
}

void CorrelationSubtractionAll(const char* fileName, Int_t centr = 0)
{
  Int_t n = 6;
  Int_t is[] = { 0, 1, 1, 2, 2, 2, 3 };
  Int_t js[] = { 1, 1, 2, 1, 2, 3, 3 };

  CreateGraphStructure();

  for (Int_t i=0; i<n; i++)
    CorrelationSubtraction(fileName, is[i], js[i], centr, graphs[0], 0);
}

Int_t gStudySystematic = 0; // 10 = exclusion zone to 0.5; 11 = exclusion off; 12 = exclusion 0.8 and mirror to AS; 13 = scale peripheral; 14 = exclusion 1.2; 20 = non closure; 30 = baseline; 40 = track cuts; 50/51 = cms comparison; 60 = other peripheral bin

void CorrelationSubtractionHistogram(const char* fileName, Bool_t centralityAxis = kTRUE, const char* outputFileName = "graphs.root")
{
  Int_t n = 6;
  Int_t is[] = { 0, 1, 1, 2, 2, 2, 3 };
  Int_t js[] = { 1, 1, 2, 1, 2, 3, 3 };
  Bool_t symm[] = { 1, 0, 1, 0, 0, 1, 0 };

//   Int_t n = 6;
//   Int_t is[] = { 0, 0, 1, 1, 1, 2 };
//   Int_t js[] = { 1, 2, 1, 2, 3, 3 };

  Int_t centralityBins = 4;

  Int_t colors[] = { 1, 2, 3, 4, 6, 7 };
  Int_t markers[] = { 20, 21, 22, 23, 24, 25, 26 };

  CreateGraphStructure();
  
  Float_t yMax[] = { 0.1, 2, 0.4, 4, 1.5, 50, 50, 4 };

  TCanvas* canvas[8];
  for (Int_t ci=0; ci<8; ci++)
  {
    canvas[ci] = new TCanvas;
    gPad->SetGridx();
    gPad->SetGridy();
    dummy = new TH2F("dummy", Form(";%s;%s", (centralityAxis) ? "Event class" : "<mult>", graphTitles[(ci < 3) ? ci*2 : ci+3]), 100, 0, 60, 100, 0, yMax[ci]);
    if (ci == 0)
      dummy->GetYaxis()->SetTitle("Yield");
    if (ci == 1)
      dummy->GetYaxis()->SetTitle("Width");
    if (ci == 2)
      dummy->GetYaxis()->SetTitle("v2 , v3");
    dummy->SetStats(0);
    dummy->Draw();
  }
  
  legend = new TLegend(0.4, 0.6, 0.99, 0.99);
  legend->SetFillColor(0);

  legend2 = new TLegend(0.4, 0.6, 0.99, 0.99);
  legend2->SetFillColor(0);

  for (Int_t i=0; i<n; i++)
  {
    for (Int_t c=0; c<centralityBins; c++)
    {
      if (c == 2)
	continue;
      
      Double_t nsYield, nsYieldE, asYield, asYieldE, nsWidth, nsWidthE, asWidth, asWidthE, v2, v2E, v3, v3E, centrality;
      CorrelationSubtraction(fileName, is[i], js[i], c, graphs[i], 1.5 * (-n/2 + i), kTRUE, symm[i], centralityAxis);
    }

    for (Int_t ci=0; ci<6; ci++)
    {
      canvas[ci/2]->cd();

      graphs[i][ci]->SetMarkerStyle(markers[i]);
      graphs[i][ci]->SetMarkerColor((ci % 2 == 0) ? 1 : 2);
      graphs[i][ci]->SetLineColor((ci % 2 == 0) ? 1 : 2);
      graphs[i][ci]->GetXaxis()->SetTitle(dummy->GetXaxis()->GetTitle());
      graphs[i][ci]->Draw("PSAME");
    }
    
    for (Int_t ci=6; ci<11; ci++)
    {
      canvas[ci-3]->cd();

      graphs[i][ci]->SetMarkerStyle(markers[i]);
      graphs[i][ci]->SetMarkerColor(1);
      graphs[i][ci]->SetLineColor(1);
      graphs[i][ci]->GetXaxis()->SetTitle(dummy->GetXaxis()->GetTitle());
      graphs[i][ci]->Draw("PSAME");
    }

    legend->AddEntry(graphs[i][0], graphs[i][0]->GetTitle(), "P");
    if (symm[i])
      legend2->AddEntry(graphs[i][0], graphs[i][0]->GetTitle(), "P");
  }
  
  const char* fileNames[] = { "ridge_yield.png", "ridge_width.png", "v2_v3.png", "ridge_yield_ratio.png", "v3_over_v2.png", "remaining_peak.png", "remaining_jet.png", "chi2ndf.png" };
  
  for (Int_t ci=0; ci<8; ci++)
  {
    canvas[ci]->cd();
    if (ci == 2 || ci == 4)
      legend2->Draw();
    else
      legend->Draw();
    canvas[ci]->SaveAs(fileNames[ci]);
  }
  
  if (gStudySystematic != 0)
    Printf(">>>>>>>>>>>> WARNING: gStudySystematic set to %d", gStudySystematic);

  WriteGraphs(outputFileName);
}

void CorrelationSubtraction(const char* fileName, Int_t i, Int_t j, Int_t centr, TGraphErrors** graph, Float_t axisOffset, Bool_t silent = kFALSE, Bool_t symmetricpT = kFALSE, Bool_t centralityAxis = kTRUE)
{
  TFile::Open(fileName);
  
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, centr));
  TH2* hist2 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, ((gStudySystematic == 60) ? 6 : 4)));
  TH1* refMult = (TH1*) gFile->Get("refMult");

  hist1 = (TH2*) hist1->Clone();
  hist2 = (TH2*) hist2->Clone();
  
  if (!hist1 || !hist2)
    return 0;
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
  hist2->Scale(1.0 / hist2->GetYaxis()->GetBinWidth(1));
  
//   hist1->Scale(1.0 / 283550); hist2->Scale(1.0 / 190000); // for 1, 2, 0
//   hist1->Scale(1.0 / 283550); hist2->Scale(1.0 / 75000); // for 2, 3, 0
//   hist1->Scale(1.0 / 283550); hist2->Scale(1.0 / 100000); // for 2, 2, 0
  
  tokens = TString(hist1->GetTitle()).Tokenize("-");
  if (tokens->GetEntries() > 1)
  {
    TString centralityStr;
    centralityStr = tokens->At(0)->GetName();
    centralityStr = centralityStr + "-" + tokens->At(1)->GetName();
//     Printf("%s", centralityStr.Data());
    for (Int_t gID=0; gID<NHists; gID++)
    {
      TString xTitle = graph[gID]->GetXaxis()->GetTitle();
      TString yTitle = graph[gID]->GetYaxis()->GetTitle();
//       Printf("%s", centralityStr.Data());
      graph[gID]->SetTitle(centralityStr);
      graph[gID]->GetXaxis()->SetTitle(xTitle);
      graph[gID]->GetYaxis()->SetTitle(yTitle);
//       Printf("%s %s", yTitle.Data(), graph[gID]->GetYaxis()->GetTitle());
    }
    Double_t centrality = -1;
    if (tokens->GetEntries() == 4)
      centrality = 0.5 * (TString(tokens->At(3)->GetName()).Atoi() + TString(tokens->At(2)->GetName()).Atoi());
  }
  
  hist1->GetYaxis()->SetRangeUser(-1.99, 1.99); hist2->GetYaxis()->SetRangeUser(-1.99, 1.99);
  if (!silent)
  {
    new TCanvas; 
    copy = (TH2*) hist1->DrawCopy("SURF1");
    copy->Rebin2D(2, 2);
    copy->Scale(0.25);
    copy->GetYaxis()->SetRangeUser(-1.99, 1.99);
    
    new TCanvas; 
    copy = (TH2*) hist2->DrawCopy("SURF1");
    copy->Rebin2D(2, 2);
    copy->Scale(0.25);
    copy->GetYaxis()->SetRangeUser(-1.99, 1.99);
  }
  
  if (gStudySystematic == 13)
  {
    Float_t factor[10][10][10];
    
    if (1)
    {
      // VZERO, from dphi_corr_pA_121108_hybrid_corrected.root
      factor[0][1][0] = 0.259;
      factor[0][1][1] = 0.212;
      factor[0][1][3] = 0.150;
      factor[1][1][0] = 0.105;
      factor[1][1][1] = 0.105;
      factor[1][1][3] = 0.081;
      factor[1][2][0] = 0.209;
      factor[1][2][1] = 0.176;
      factor[1][2][3] = 0.121;
      factor[2][1][0] = 0.050;
      factor[2][1][1] = -0.005;
      factor[2][1][3] = 0.031;
      factor[2][2][0] = 0.042;
      factor[2][2][1] = 0.061;
      factor[2][2][3] = 0.022;
      factor[2][3][0] = 0.081;
      factor[2][3][1] = 0.073;
      factor[2][3][3] = 0.058;
    }
    else if (0)
    {
      // ZNA, from dphi_corr_pA_121029_zna.root
      factor[0][1][0] = 0.125;
      factor[0][1][1] = 0.114;
      factor[0][1][3] = 0.052;
      factor[1][1][0] = 0.052;
      factor[1][1][1] = 0.065;
      factor[1][1][3] = 0.024;
      factor[1][2][0] = 0.113;
      factor[1][2][1] = 0.102;
      factor[1][2][3] = 0.060;
      factor[2][1][0] = 0.014;
      factor[2][1][1] = 0.023;
      factor[2][1][3] = -0.017;
      factor[2][2][0] = 0.018;
      factor[2][2][1] = 0.037;
      factor[2][2][3] = -0.005;
      factor[2][3][0] = 0.068;
      factor[2][3][1] = 0.066;
      factor[2][3][3] = 0.035;
    }
    else if (0) //V0, cmsmethod
    {
      factor[0][1][0] = 0.470 +  0.221 + 0.103;
      factor[0][1][1] = 0.432 +  0.186 + 0.080;
      factor[0][1][3] = 0.339 +  0.115 + 0.039;
      factor[1][1][0] = 0.164 +  0.027 + 0.005;
      factor[1][1][1] = 0.160 +  0.025 + 0.004;
      factor[1][1][3] = 0.111 +  0.013 + 0.001;
      factor[1][2][0] = 0.458 +  0.210 + 0.097;
      factor[1][2][1] = 0.385 +  0.148 + 0.057;
      factor[1][2][3] = 0.280 +  0.078 + 0.022;
      factor[2][1][0] = 0.026 +  0.000 + 0.000;
      factor[2][1][1] = -0.018 + 0.001 + -0.000;
      factor[2][1][3] = -0.005 + 0.000 + -0.000;
      factor[2][2][0] = 0.047 +  0.003 + -0.000;
      factor[2][2][1] = 0.044 +  0.002 + 0.000;
      factor[2][2][3] = 0.013 +  0.000 + 0.000;
      factor[2][3][0] = 0.258 +  0.067 + 0.017;
      factor[2][3][1] = 0.176 +  0.031 + 0.005;
      factor[2][3][3] = 0.118 +  0.014 + 0.002;
    }
    else
    {
      factor[1][2][0] = 0.3;
    }
    
    hist2->Scale(1.0 + factor[i][j][centr]);
  }
  
  const Float_t etaFlat = 1.2;
  const Float_t centralEta = 0.5;

  projCentral = hist1->ProjectionY(Form("%s_proj3y", hist1->GetName()), hist1->GetXaxis()->FindBin(-0.7), hist1->GetXaxis()->FindBin(0.7));
  projCentral->Scale(hist1->GetXaxis()->GetBinWidth(1));
  projCentral2 = hist2->ProjectionY(Form("%s_proj4y", hist2->GetName()), hist2->GetXaxis()->FindBin(-0.7), hist2->GetXaxis()->FindBin(0.7));
  projCentral2->Scale(hist2->GetXaxis()->GetBinWidth(1));

  Float_t yEta1 = projCentral->FindBin(-etaMax+0.01);
  Float_t yEta2 = projCentral->FindBin(-etaFlat-0.01);
  Float_t yEta3 = projCentral->FindBin(-centralEta+0.01);
  Float_t yEta4 = projCentral->FindBin(centralEta-0.01);
  Float_t yEta5 = projCentral->FindBin(etaFlat+0.01);
  Float_t yEta6 = projCentral->FindBin(etaMax-0.01);
  Float_t yBaseline = (projCentral->Integral(yEta1, yEta2) + projCentral->Integral(yEta5, yEta6)) / (yEta2 - yEta1 + 1 + yEta6 - yEta5 + 1);
  Float_t yPeak = projCentral->Integral(yEta3, yEta4, "width") - yBaseline * (yEta4 - yEta3 + 1) * projCentral->GetBinWidth(1);
  Printf("Peak (unsubtracted): %f %f", yBaseline, yPeak);

  if (!silent)
  {
    Float_t yBaseline = (projCentral->Integral(yEta1, yEta2) + projCentral->Integral(yEta5, yEta6)) / (yEta2 - yEta1 + 1 + yEta6 - yEta5 + 1);
    Float_t yBaseline2 = (projCentral2->Integral(yEta1, yEta2) + projCentral2->Integral(yEta5, yEta6)) / (yEta2 - yEta1 + 1 + yEta6 - yEta5 + 1);
    Printf("baselines: %f %f", yBaseline, yBaseline2);
    
    projCentral2->Add(new TF1("flat", "1", -5, 5), yBaseline - yBaseline2);
    
    new TCanvas; 
    projCentral->DrawCopy();
    projCentral2->SetLineColor(2);
    projCentral2->DrawCopy("SAME");
  }
  
  Float_t xValue1 = -1;
  Float_t xValue2 = -1;
  Float_t xValue1vn = -1;
  Float_t xValue2vn = -1;
  if (centralityAxis)
  {
    xValue1 = centrality + axisOffset;
    xValue2 = xValue1 + 0.5;
    xValue1vn = centrality + (Int_t) axisOffset;
    xValue2vn = xValue1vn + 0.5;
  }
  else
  {
    if (centrality > 0)
    {
      xValue1 = refMult->GetBinContent(refMult->GetXaxis()->FindBin(centrality));
      xValue2 = xValue1;
      xValue1vn = xValue1;
      xValue2vn = xValue1;
//       Printf("%f %f %f", centrality, xValue1, xValue2);
    }
  }
  
  Float_t exclusion = 0.8;
  if (gStudySystematic == 10)
    exclusion = 0.5;
  else if (gStudySystematic == 11)
    exclusion = 0.0;
  else if (gStudySystematic == 14 || gStudySystematic == 51)
    exclusion = 1.2;
  
  Int_t eta1 = hist1->GetYaxis()->FindBin(-etaMax + 0.001);
  Int_t eta2 = hist1->GetYaxis()->FindBin(-exclusion + 0.001);
  Int_t eta3 = hist1->GetYaxis()->FindBin(exclusion - 0.001);
  Int_t eta6 = hist1->GetYaxis()->FindBin(etaMax - 0.001);
  Int_t phi1Z = hist1->GetXaxis()->FindBin(TMath::Pi()/2 - 0.2);
  Int_t phi2Z = hist1->GetXaxis()->FindBin(TMath::Pi()/2 + 0.2);
  
//   new TCanvas; tmpProj = hist1->ProjectionX("tmp", eta1, eta6); tmpProj->Scale(hist1->GetYaxis()->GetBinWidth(1)); tmpProj->Scale(1.0 / (etaMax * 2)); tmpProj->Draw();
  
  Double_t baseLineE;
  Float_t baseLine = hist1->IntegralAndError(phi1Z, phi2Z, eta1, eta6, baseLineE);
  baseLine /= (phi2Z - phi1Z + 1) * (eta6 - eta1 + 1);
  baseLineE /= (phi2Z - phi1Z + 1) * (eta6 - eta1 + 1);
  
  hist1->Add(hist2, -1);
  
  // phi projection
  // NS
  proj = hist1->ProjectionX(Form("%s_proj1x", hist1->GetName()), eta1, eta2);
  proj2 = hist1->ProjectionX(Form("%s_proj2x", hist1->GetName()), eta3, eta6);
  proj->Add(proj2, 1);
  
  // AS
  projAS = hist1->ProjectionX(Form("%s_proj3x", hist1->GetName()), eta1, eta6);

  // match NS and AS yield
  proj->Scale(1.0 * (eta6 - eta1 + 1) / (eta6 - eta3 + 1 + eta2 - eta1 + 1));
  
  // copy AS
  for (Int_t bin=proj->FindBin(TMath::Pi()/2+0.001); bin<=proj->GetNbinsX(); bin++)
  {
//     Printf("%d %f", bin, projAS->GetBinContent(bin));
    proj->SetBinContent(bin, projAS->GetBinContent(bin));
    proj->SetBinError(bin, projAS->GetBinError(bin));
  }
  
  if (gStudySystematic == 12)
  {
    // get difference between exclusion 0.8 and 0.0, shift to AS, remove from AS
    
    projAll = hist1->ProjectionX(Form("%s_proj4x", hist1->GetName()), eta1, eta6);
    projAll->Add(proj, -1);
    if (!silent)
    {
      new TCanvas; projAll->DrawCopy();
    }

    // From Tim, 13.11.12
    Float_t ratioNSAS[10][10][10];
    ratioNSAS[0][1][4] = 2.060797;
    ratioNSAS[1][1][4] = 1.433904;
    ratioNSAS[1][2][4] = 1.390280;
    ratioNSAS[2][1][4] = 1.013925;
    ratioNSAS[2][2][4] = 1.070422;
    ratioNSAS[2][3][4] = 1.087238;
    ratioNSAS[3][3][4] = 1.023756;
    
    Printf("Using NS/AS ratio %f", ratioNSAS[i][j][4]);

    for (Int_t k=0; k<=projAll->GetNbinsX()/2; k++)
    {
  //     Printf("%d -> %d", i, projAll->GetNbinsX()+1-i);
      projAll->SetBinContent(projAll->GetNbinsX()+1-k, projAll->GetBinContent(k) / ratioNSAS[i][j][4]);
      projAll->SetBinError(projAll->GetNbinsX()+1-k, projAll->GetBinError(k) / ratioNSAS[i][j][4]);
      
      projAll->SetBinContent(k, 0);
      projAll->SetBinError(k, 0);
    }
    if (!silent)
    {
      new TCanvas; projAll->DrawCopy();
    }
    
//     new TCanvas; proj->DrawCopy();
    proj->Add(projAll, -1);
//     new TCanvas; proj->DrawCopy();
  }
  
  proj->Scale(hist1->GetYaxis()->GetBinWidth(1));
  proj->Scale(1.0 / (etaMax * 2));
  proj->Rebin(2); proj->Scale(0.5);

  if (gStudySystematic == 20)
  {
    Printf(">>>>>>>>>>>> Applying non-closure systematics <<<<<<<<<<<<");
    file2 = TFile::Open("non_closure.root");
    non_closure = (TH1*) file2->Get(Form("non_closure_%d_%d_%d", i, j, 0));
    for (Int_t bin=1; bin<=non_closure->GetNbinsX(); bin++)
      non_closure->SetBinError(bin, 0);
    
    proj->Multiply(non_closure);
  }
  
//   proj = hist1->ProjectionX("proj1x", hist1->GetYaxis()->FindBin(0.5), eta6);

  // eta projection
  Int_t binEtaProj1 = 1;
  Int_t binEtaProj2 = hist1->GetXaxis()->FindBin(-TMath::Pi()/3-0.001);
  Int_t binEtaProj3 = binEtaProj2 + 1;
  Int_t binEtaProj4 = hist1->GetXaxis()->FindBin(TMath::Pi()/3-0.001);
  Int_t binEtaProj5 = binEtaProj4+1;
  Int_t binEtaProj6 = hist1->GetXaxis()->FindBin(TMath::Pi()-TMath::Pi()/3-0.001);
  Int_t binEtaProj7 = binEtaProj6+1;
  Int_t binEtaProj8 = hist1->GetXaxis()->FindBin(TMath::Pi()+TMath::Pi()/3-0.001);
  Int_t binEtaProj9 = binEtaProj8+1;
  Int_t binEtaProj10 = hist1->GetNbinsX();
  Printf("%d %d %d %d %d %d %d %d %d %d", binEtaProj1, binEtaProj2, binEtaProj3, binEtaProj4, binEtaProj5, binEtaProj6, binEtaProj7, binEtaProj8, binEtaProj9, binEtaProj10);
  
  etaProj = hist1->ProjectionY(Form("%s_proj1y", hist1->GetName()), binEtaProj3, binEtaProj4);
  etaProj->Scale(1.0 / (binEtaProj4 - binEtaProj3 + 1));
  etaProj2 = hist1->ProjectionY(Form("%s_proj2y", hist1->GetName()), binEtaProj7, binEtaProj8);
  etaProj2->Scale(1.0 / (binEtaProj8 - binEtaProj7 + 1));
  etaProj3 = hist1->ProjectionY(Form("%s_proj5y", hist1->GetName()), binEtaProj1, binEtaProj2);
  etaProj3b = hist1->ProjectionY(Form("%s_proj5by", hist1->GetName()), binEtaProj5, binEtaProj6);
  etaProj3c = hist1->ProjectionY(Form("%s_proj5cy", hist1->GetName()), binEtaProj9, binEtaProj10);
  etaProj3->Add(etaProj3b);
  etaProj3->Add(etaProj3c);
  etaProj3->Scale(1.0 / (binEtaProj2 - binEtaProj1 + 1 + binEtaProj6 - binEtaProj5 + 1 + binEtaProj10 - binEtaProj9 + 1));

  Float_t ySubBaseline = (etaProj->Integral(yEta1, yEta2) + etaProj->Integral(yEta5, yEta6)) / (yEta2 - yEta1 + 1 + yEta6 - yEta5 + 1);
  Float_t ySubPeak = etaProj->Integral(yEta3, yEta4, "width") - ySubBaseline * (yEta4 - yEta3 + 1) * etaProj->GetBinWidth(1);
  Printf("Peak (subtracted): %f %f (%.2f%% of unsubtracted peak)", ySubBaseline, ySubPeak, ySubPeak / yPeak * 100);
  Printf("    factor[%d][%d][%d] = %.3f;", i, j, centr, ySubPeak / yPeak);
  AddPoint(graph[8], xValue1, 100.0 * ySubPeak / yPeak, 0, 0);
  
  hist1->Rebin2D(4, 4); hist1->Scale(1.0 / 16);
  hist1->GetYaxis()->SetRangeUser(-1.99, 1.99);
  
  TString fitOption = "0I";
  if (!silent)
  {
    c = new TCanvas("c", "c", 600, 600);

    PadFor2DCorr();
    
    hist1->GetYaxis()->SetRangeUser(-1.99, 1.99);
    hist1->GetXaxis()->SetTitleOffset(1.7);
    hist1->GetYaxis()->SetTitleOffset(2);
    hist1->GetZaxis()->SetTitle(kCorrFuncTitle);
    hist1->GetZaxis()->SetTitleOffset(2);
    hist1->GetZaxis()->SetNdivisions(504);
    hist1->SetStats(kFALSE);
    hist1->GetZaxis()->CenterTitle(kTRUE);
    hist1->GetYaxis()->CenterTitle(kTRUE);
    hist1->GetXaxis()->CenterTitle(kTRUE);
    hist1->GetXaxis()->SetTitle("#Delta#varphi (rad)");
    hist1->GetYaxis()->SetNdivisions(505);
    
    TString label(hist1->GetTitle());
    hist1->SetTitle("");
    label.ReplaceAll(".00", "");
    label.ReplaceAll(".0", "");
    TObjArray* objArray = label.Tokenize("-");
    TPaveText* paveText = new TPaveText(0.03, 0.86, 0.44, 0.97, "BRNDC");
    paveText->SetTextSize(0.035);
    paveText->SetFillColor(0);
    paveText->SetShadowColor(0);
    paveText->SetBorderSize(0);
    paveText->SetFillStyle(0);
    TString tmpStr(objArray->At(0)->GetName());
    tmpStr.ReplaceAll("p_", "#it{p}_");
    paveText->AddText(tmpStr + "GeV/#it{c}");

    TString tmpStr(objArray->At(1)->GetName());
    tmpStr.ReplaceAll("p_", "#it{p}_");
    paveText->AddText(tmpStr + "GeV/#it{c}");
  
    TString label2(hist2->GetTitle());
    TObjArray* objArray2 = label2.Tokenize("-");

    TPaveText* paveText2 = new TPaveText(0.65, 0.86, 0.98, 0.97, "BRNDC");
    paveText2->SetTextSize(0.035);
    paveText2->SetBorderSize(0);
    paveText2->SetFillColor(0);
    paveText2->SetFillStyle(0);
    paveText2->SetShadowColor(0);
    paveText2->AddText("p-Pb #sqrt{s_{NN}} = 5.02 TeV");
    if (objArray->GetEntries() > 3 && objArray2->GetEntries() > 3)
    {
      TString centrBeginStr;
      centrBeginStr = objArray->At(2)->GetName();
      centrBeginStr.ReplaceAll(" ", "");
      TString centrEndStr;
      centrEndStr = objArray2->At(2)->GetName();
      centrEndStr.ReplaceAll(" ", "");
      paveText2->AddText(Form("(%s-%s) - (%s-%s)", centrBeginStr.Data(), objArray->At(3)->GetName(), centrEndStr.Data(), objArray2->At(3)->GetName()));
    }

    hist1->DrawCopy("SURF1");
    paveText->Draw();
    paveText2->Draw();
    gPad->GetCanvas()->SaveAs(Form("ridge_%d_%d.png", i, j));
    gPad->GetCanvas()->SaveAs(Form("ridge_%d_%d.eps", i, j));
    gPad->GetCanvas()->SaveAs("fig3a.eps");
    
    Float_t fontSize = 0.05;

    c3 = new TCanvas("c3", "c3", 600, 400);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.01);
    etaProj->SetTitle();
//     etaProj->GetYaxis()->SetNdivisions(505);
    etaProj->GetYaxis()->SetLabelSize(fontSize);
    etaProj->GetXaxis()->SetLabelSize(fontSize);
    etaProj->GetXaxis()->SetTitleSize(fontSize);
    etaProj->GetYaxis()->SetTitleSize(fontSize);
    etaProj->GetYaxis()->SetTitle(kProjYieldTitleEta);
    etaProj->GetYaxis()->SetTitleOffset(1.1);
    etaProj->SetStats(0);
    etaProj->GetXaxis()->SetNdivisions(505);
    etaProj->Rebin(2); etaProj->Scale(0.5);
    etaProj2->Rebin(2); etaProj2->Scale(0.5);
    etaProj3->Rebin(2); etaProj3->Scale(0.5);
    etaProj->Rebin(2); etaProj->Scale(0.5);
    etaProj2->Rebin(2); etaProj2->Scale(0.5);
    etaProj3->Rebin(2); etaProj3->Scale(0.5);
    etaProj->GetXaxis()->SetRangeUser(-1.99, 1.99);
//     etaProj->GetYaxis()->SetRangeUser(TMath::Min(etaProj->GetMinimum(), etaProj3->GetMinimum()) * 0.8, TMath::Max(etaProj->GetMaximum(), etaProj3->GetMaximum()) * 1.2);
    etaProj->GetYaxis()->SetRangeUser(proj->GetMinimum() * 0.98, proj->GetMaximum() * 1.065);
    etaProj2->SetLineColor(2); etaProj2->SetMarkerColor(2);
    etaProj3->SetLineColor(4); etaProj3->SetMarkerColor(4);
    etaProj->SetMarkerStyle(24);
    etaProj2->SetMarkerStyle(25);
    etaProj3->SetMarkerStyle(26);
    etaProj->DrawCopy("E0 X0");
    etaProj2->DrawCopy("E0 X0 SAME");
    etaProj3->DrawCopy("E0 X0 SAME");
    
    legend5 = new TLegend(0.48, 0.74, 0.92, 0.97);
    legend5->SetBorderSize(0);
    legend5->SetTextSize(fontSize);
    legend5->SetFillColor(0);
    legend5->AddEntry(etaProj,  "|#Delta#varphi| < #pi/3", "P");
    legend5->AddEntry(etaProj2, "|#Delta#varphi - #pi| < #pi/3", "P");
//     legend5->AddEntry(etaProj3, "#pi/3 < |#Delta#varphi| < 2#pi/3, #Delta#varphi > 4#pi/3", "P");
    legend5->AddEntry(etaProj3, "Remaining #Delta#varphi", "P");
    legend5->Draw();
    
    c = new TCanvas("c2", "c2", 600, 400);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.01);
    
    proj->SetStats(0);
//     proj->GetYaxis()->SetNdivisions(505);
    proj->GetXaxis()->SetTitle("#Delta#varphi (rad)");
    proj->GetYaxis()->SetTitle(kProjYieldTitlePhi);
    proj->GetYaxis()->SetTitleOffset(1.1);
    proj->GetYaxis()->SetLabelSize(fontSize);
    proj->GetXaxis()->SetLabelSize(fontSize);
    proj->GetXaxis()->SetTitleSize(fontSize);
    proj->GetYaxis()->SetTitleSize(fontSize);
    proj->SetTitle();
    proj->SetMarkerStyle(21);
    proj->SetMarkerSize(0.7);
    proj->Draw("E0 X0");
    proj->SetMinimum(proj->GetMinimum() * 0.98); proj->SetMaximum(proj->GetMaximum() * 1.065);
//     proj->SetMinimum(proj->GetMinimum() * 1.2); proj->SetMaximum(proj->GetMaximum() * 0.6);
    fitOption = "I";
    
    if (0)
    {
      Printf("\nPer-trigger yield per unit of delta eta (y) as function of delta phi (x) in the bin %s", centralityStr.Data());
      Printf("Systematic uncertainties are mostly correlated and affect the baseline. Uncorrelated uncertainties are less than 1%% and not indicated explicitly in the table below.");

      for (Int_t k=1; k<=proj->GetNbinsX(); k++)
	Printf("x = %.2f, y = %.4f +- %.4f (stat)", proj->GetXaxis()->GetBinCenter(k), proj->GetBinContent(k), proj->GetBinError(k));
    }
  }
  
  fileProj = TFile::Open("dphi_proj.root", "UPDATE");
  proj->Write(Form("proj_%d_%d_%d_%d", i, j, centr, gStudySystematic));
  fileProj->Close();

  TF1* v2 = new TF1("func", "[0]+2*[1]*cos(2*x)", -5, 5);
  v2->SetLineStyle(2);
//   v2->SetLineWidth(1);
  
  TF1* v2v3 = new TF1("func", "[0]+2*[1]*cos(2*x)+2*[2]*cos(3*x)", -5, 5);
  v2v3->SetLineColor(2);
//   v2v3->FixParameter(2, 0);
  proj->Fit(v2, fitOption);
//   return;
    
  fitOption += "+";
  proj->Fit(v2v3, fitOption, "E0 X0");

  Float_t min = v2v3->GetMinimum();
  Float_t diffMinParam0 = v2v3->GetParameter(0) - min;
  Printf("Chi2: %f ndf: %d chi2/ndf %f:Min: %f %f", v2v3->GetChisquare(), v2v3->GetNDF(), v2v3->GetChisquare() / v2v3->GetNDF(), min, diffMinParam0);
  AddPoint(graph[10], xValue1, v2v3->GetChisquare() / v2v3->GetNDF(), 0, 0);

  AddPoint(graph[11], xValue1, baseLine + diffMinParam0, 0, baseLineE);
//   AddPoint(graph[11], xValue1, baseLine, 0, baseLineE);
  
  Float_t v2value = 0, v2E = 0, v3 = 0, v3E = 0;
  if (v2v3->GetParameter(1) > 0)
  {
    v2value = TMath::Sqrt(v2v3->GetParameter(1) / (baseLine + diffMinParam0));
    v2E = 0.5 * v2value * TMath::Sqrt(v2v3->GetParError(1) * v2v3->GetParError(1) / v2v3->GetParameter(1) / v2v3->GetParameter(1) + baseLineE * baseLineE / baseLine / baseLine);
    if (symmetricpT)
      AddPoint(graph[4], xValue1vn, v2value, 0, v2E);
  }
  if (v2v3->GetParameter(2) > 0)
  {
    v3 = TMath::Sqrt(v2v3->GetParameter(2) / (baseLine + diffMinParam0));
    v3E = 0.5 * v3 * TMath::Sqrt(v2v3->GetParError(2) * v2v3->GetParError(2) / v2v3->GetParameter(2) / v2v3->GetParameter(2) + baseLineE * baseLineE / baseLine / baseLine);
    if (symmetricpT)
      AddPoint(graph[5], xValue2vn, v3, 0, v3E);
  }
  if (v2v3->GetParameter(1) > 0 && v2v3->GetParameter(2) > 0 && symmetricpT)
  {
    Float_t v3v2Ratio = TMath::Sqrt(v2v3->GetParameter(2) / v2v3->GetParameter(1));
    Float_t v3v2RatioE = 0.5 * v3v2Ratio * TMath::Sqrt(v2v3->GetParError(1) * v2v3->GetParError(1) / v2v3->GetParameter(1) / v2v3->GetParameter(1) + v2v3->GetParError(2) * v2v3->GetParError(2) / v2v3->GetParameter(2) / v2v3->GetParameter(2));
    AddPoint(graph[7], xValue1vn, v3v2Ratio, 0, v3v2RatioE);
  }
  
  Printf("Baseline: %f +- %f; v2 = %f +- %f; v3 = %f +- %f", baseLine, baseLineE, v2value, v2E, v3, v3E);
  
  if (gStudySystematic == 30)
  {
    // alternative way for baseline
    
    parabola = new TF1("parabola", "[0] + [1]*(x - [2])**2",  0, TMath::Pi());
    parabola->SetLineColor(3);
    parabola->SetParameters(min, -0.1, TMath::Pi() / 2);
//     parabola->SetParLimits(1, 0, 1);
    proj->Fit(parabola, fitOption, "", TMath::Pi() / 2 - 1, TMath::Pi() / 2 + 1);
//     proj->Fit(parabola, fitOption, "", 0.1, 2);
    
    min = parabola->GetParameter(0);
    Printf("New baseline: %f", min);
  }

  fitOption += "R";
  
  if (0)
  {
    // single gauss fit
    
    TF1* gaus1 = new TF1("gaus1", "[0]+gaus(1)", -0.5 * TMath::Pi(), 0.5 * TMath::Pi());
    gaus1->SetParameters(min, 1, 0, 1);
  //   gaus1->FixParameter(0, min);
    gaus1->FixParameter(2, 0);
    gaus1->SetParLimits(1, 0.001, 10);
    gaus1->SetParLimits(3, 0.1, 2);
    gaus1->SetLineColor(3);
    proj->Fit(gaus1, fitOption);
    
    TF1* gaus2 = new TF1("gaus2", "[0]+gaus(1)", 0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    gaus2->SetParameters(min, 1, TMath::Pi(), 1);
    gaus2->FixParameter(2, TMath::Pi());
    gaus2->SetParLimits(1, 0.001, 10);
    gaus2->SetParLimits(3, 0.1, 2);
    gaus2->SetLineColor(3);
    proj->Fit(gaus2, fitOption);
  
    AddPoint(graph[2], xValue1, gaus1->GetParameter(3), 0, gaus1->GetParError(3));
    AddPoint(graph[3], xValue2, gaus2->GetParameter(3), 0, gaus2->GetParError(3));
  }
  else if (0)
  {
    // combined gauss fit
    
    TF1* gausBoth = new TF1("gaus2", "[0]+[1]*(exp(-0.5*(x/[2])**2)+exp(-0.5*((x-TMath::TwoPi())/[2])**2))+[3]*(exp(-0.5*((x-TMath::Pi())/[4])**2)+exp(-0.5*((x+TMath::Pi())/[4])**2))", -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    gausBoth->SetParameters(min, 1, 1, 1, 1);
  //   gausBoth->FixParameter(0, min);
    gausBoth->SetParLimits(1, 0.0001, 10);
    gausBoth->SetParLimits(3, 0.0001, 10);
    gausBoth->SetParLimits(2, 0.1, 4);
    gausBoth->SetParLimits(4, 0.1, 4);
    gausBoth->SetLineColor(6);
    
    proj->Fit(gausBoth, fitOption);
  
    AddPoint(graph[2], xValue1, gausBoth->GetParameter(2), 0, gausBoth->GetParError(3));
    AddPoint(graph[3], xValue2, gausBoth->GetParameter(4), 0, gausBoth->GetParError(4));
  }
  else
  {
    // RMS without baseline
    projRMS = (TH1*) proj->Clone();
    projRMS->Add(new TF1("f", "1", -10, 10), -min);
    projRMS->GetXaxis()->SetRangeUser(-TMath::Pi()/2 + 0.01, TMath::Pi()/2 - 0.01);
    AddPoint(graph[2], xValue1, projRMS->GetRMS(), 0, projRMS->GetRMSError());
    projRMS->GetXaxis()->SetRangeUser(TMath::Pi()/2 + 0.01, 3 * TMath::Pi()/2 - 0.01);
    AddPoint(graph[3], xValue2, projRMS->GetRMS(), 0, projRMS->GetRMSError());
//     new TCanvas; projRMS->Draw(); return;
  }
  
//   new TCanvas; gausBoth->Draw(); return;

  if (!silent)
  {
    TF1* v2v3_v2 = new TF1("func", "[0]+2*[1]*cos(2*x)", -5, 5);
    v2v3_v2->SetParameters(v2v3->GetParameter(0), v2v3->GetParameter(1));
    v2v3_v2->SetLineStyle(3);
    v2v3_v2->SetLineColor(2);
//     v2v3_v2->Draw("SAME");
    TF1* v2v3_v3 = new TF1("func", "[0]+2*[1]*cos(3*x)", -5, 5);
    v2v3_v3->SetParameters(v2v3->GetParameter(0), v2v3->GetParameter(2));
    v2v3_v3->SetLineStyle(4);
    v2v3_v3->SetLineColor(2);
//     v2v3_v3->Draw("SAME");
    
    line = new TLine(-0.5 * TMath::Pi(), min, 1.5 * TMath::Pi(), min);
    line->SetLineWidth(2);
    line->SetLineColor(4);
    line->Draw();

    legend = new TLegend(0.47, 0.65, 0.88, 0.95);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);

    legend->AddEntry(proj, "Data", "P");
    legend->AddEntry(v2v3, "a_{0} + a_{2} cos(2#Delta#varphi) + a_{3} cos(3#Delta#varphi)", "L");
    legend->AddEntry(v2, "a_{0} + a_{2} cos(2#Delta#varphi)", "L");
    legend->AddEntry(line, "Baseline for yield extraction", "L");

    if (1)
    {
      hijingFile = TFile::Open("dphi_proj_hijing.root");
      if (hijingFile)
      {
	hijing = (TH1*) hijingFile->Get(Form("proj_%d_%d_%d_0", i, j, centr));
	if (hijing)
	{
	  hijing->Add(new TF1("flat", "1", -5, 5), min - ((i == 1) ? 0.192078 : 0.373710));
	  hijing->SetMarkerColor(2);
	  hijing->SetLineColor(2);
	  hijing->SetMarkerStyle(25);
	  hijing->Draw("SAME E0X0");
	  legend->AddEntry(hijing, "HIJING shifted", "P");
	}
      }
    }
    
    legend->SetTextSize(fontSize);
    legend->Draw();
    
//     paveText3 = (TPaveText*) paveText2->Clone();
//     paveText3->SetX1NDC(0.16);
//     paveText3->SetY1NDC(0.77);
//     paveText3->SetX2NDC(0.42);
//     paveText3->SetY2NDC(0.94);
//     paveText3->Draw();
    paveText4 = (TPaveText*) paveText2->Clone();
    paveText4->SetTextSize(fontSize);
    paveText4->SetX1NDC(0.16);
    paveText4->SetY1NDC(0.68);
    paveText4->SetX2NDC(0.42);
    paveText4->SetY2NDC(0.96);

    TString tmpStr(objArray->At(0)->GetName());
    tmpStr.ReplaceAll("p_", "#it{p}_");
    paveText4->AddText(tmpStr + "GeV/#it{c}");

    TString tmpStr(objArray->At(1)->GetName());
    tmpStr.ReplaceAll("p_", "#it{p}_");
    paveText4->AddText(tmpStr + "GeV/#it{c}");
    
    paveText4->Draw();
    
//     DrawLatex(0.27, 0.19, 1, Form("%sGeV/#it{c}       %sGeV/#it{c}", objArray->At(0)->GetName(), objArray->At(1)->GetName()), fontSize);
    
    gPad->GetCanvas()->SaveAs(Form("ridge_fit_%d_%d.png", i, j));
    gPad->GetCanvas()->SaveAs(Form("ridge_fit_%d_%d.eps", i, j));
    gPad->GetCanvas()->SaveAs("fig3b.eps");
    
    c3->cd();
    paveText3 = (TPaveText*) paveText4->Clone();
    paveText3->Draw();
    gPad->GetCanvas()->SaveAs(Form("ridge_eta_%d_%d.png", i, j));
    gPad->GetCanvas()->SaveAs(Form("ridge_eta_%d_%d.eps", i, j));
    gPad->GetCanvas()->SaveAs("fig3c.eps");
  }
  
  if (gStudySystematic != 50 && gStudySystematic != 51)
  {
    Int_t phi1 = 1;
    Int_t phi2 = proj->FindBin(TMath::Pi()/2-0.001);
    Int_t phi3 = phi2 + 1;
    Int_t phi4 = proj->GetNbinsX();
  }
  else
  {
    Printf(">>>> Using |dphi| < 1.2 <<<<");
    Int_t phi1 = proj->FindBin(-1.2);;
    Int_t phi2 = proj->FindBin(1.2);
    Int_t phi3 = proj->FindBin(TMath::Pi() - 1.2);
    Int_t phi4 = proj->FindBin(TMath::Pi() + 1.2);
    Printf("%d %d %d %d", phi1, phi2, phi3, phi4);
  }
  
  Double_t nsYield, asYield, nsYieldE, asYieldE;
  
  nsYield = proj->IntegralAndError(phi1, phi2, nsYieldE, "width");
  nsYield -= min * proj->GetBinWidth(1) * (phi2 - phi1 + 1);
  
  asYield = proj->IntegralAndError(phi3, phi4, asYieldE, "width");
  asYield -= min * proj->GetBinWidth(1) * (phi4 - phi3 + 1);
  
  AddPoint(graph[0], xValue1, nsYield, 0, nsYieldE);
  AddPoint(graph[1], xValue2, asYield, 0, asYieldE);
  
  if (nsYieldE > 0 && asYieldE > 0)
    AddPoint(graph[6], xValue1, nsYield / asYield, 0, nsYield / asYield * TMath::Sqrt(nsYieldE * nsYieldE / nsYield / nsYield + asYieldE * asYieldE / asYield / asYield));
  
  AddPoint(graph[9], xValue1, 100.0 * ySubPeak / (etaMax * 2) / nsYield, 0, 0);
  
  Printf("Yields: %f +- %f ; %f +- %f", nsYield, nsYieldE, asYield, asYieldE);
}

Int_t colors[] = { 1, 2, 3, 4, kGray+1, 6, 7, 8 };

void DrawSeveral(Int_t n, const char** graphFiles, Int_t id)
{
  ReadGraphs(graphFiles[0]);
  TGraphErrors*** base = graphs;

  Float_t yMax[] = { 0.1, 0.1, 2, 2, 0.25, 0.25, 4, 1.5, 50, 50, 4, 4 };
  Int_t markers[] = { 20, 21, 22, 23, 24, 25, 26 };

  TString baseName(graphFiles[0]);
  baseName.ReplaceAll(".", "_");
  TCanvas* canvas = new TCanvas(Form("%s_%d", baseName.Data(), id), Form("%s_%d", baseName.Data(), id), 800, 600);
//   Printf("%p", canvas);
  gPad->SetGridx();
  gPad->SetGridy();
  dummy = new TH2F(Form("hist_%s_%d", graphFiles[0], id), Form(";%s;%s", graphs[0][id]->GetXaxis()->GetTitle(), graphs[0][id]->GetYaxis()->GetTitle()), 100, 0, 60, 100, 0, yMax[id]);
  dummy->SetStats(0);
  dummy->DrawCopy();
  
  TCanvas* canvas2 = new TCanvas(Form("%s_%d_ratio", baseName.Data(), id), Form("%s_%d", baseName.Data(), id), 800, 600);
  gPad->SetGridx();
  gPad->SetGridy();
  dummy = new TH2F(Form("hist_%s_%d_ratio", graphFiles[0], id), Form(";%s;%s ratio", graphs[0][id]->GetXaxis()->GetTitle(), graphs[0][id]->GetYaxis()->GetTitle()), 100, 0, 60, 100, 0, 2);
  dummy->SetStats(0);
  dummy->DrawCopy();

  legend = new TLegend(0.4, 0.6, 0.99, 0.99);
  legend->SetFillColor(0);
  
  for (Int_t fc = 0; fc<n; fc++)
  {
    ReadGraphs(graphFiles[fc]);
    
    for (Int_t i=0; i<NGraphs; i++)
    {
      if (TString(graphFiles[0]).Contains("cms") && i != 2 && i != 5) continue;
    
      canvas->cd();
//       Printf("%p", canvas);
      graphs[i][id]->SetMarkerStyle(markers[i]);
      graphs[i][id]->SetMarkerColor(colors[fc]);
      graphs[i][id]->SetLineColor(colors[fc]);
      GraphShiftX((TGraphErrors*) graphs[i][id]->DrawClone("PSAME"), 0.5 / n * fc);

      if (fc == 0 && graphs[i][id]->GetN() > 0)
	legend->AddEntry(graphs[i][id], graphs[i][id]->GetTitle(), "P");
      
      if (fc > 0)
      {
	canvas2->cd();
	DivideGraphs(graphs[i][id], base[i][id]);
	GraphShiftX((TGraphErrors*) graphs[i][id]->DrawClone("PSAME"), 0.5 / n * fc);
      }
    }
  }

  canvas->cd();
  legend->Draw();
  canvas->SaveAs(Form("%s.eps", canvas->GetName()));

  canvas2->cd();
  legend->Draw();
  canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));
}

void CMSRidge()
{
  CreateGraphStructure();
  
  // scaling done for zyam level difference (note that the values in fig2 have to be multiplied by 2 due to the |dphi| instead of dphi
  
  // Fig3a
  // 0.334516912727  0.0104728088208
  // 0.721360493366  0.0263651653896
  // 1.1900579331  0.0301999626238
  // 1.68965302436  0.024186548724
  // 2.19286114745  0.0123441101352
  // 2.69021366723  0.00994871155963
  // 3.38790257273  0.00434394401877

//   AddPoint(graphs[0][0], 55.3, 0.0263651653896, 0, 0);
  AddPoint(graphs[2][0], 55.383968, 0.0301999626238 * 0.72 * 2 / 1.336, 0, 0);
  AddPoint(graphs[5][0], 55.383968, 0.0123441101352 * 0.154 * 2 / 0.229, 0, 0);
  
  // Fig3b
  // 128.344370861  0.0254622516556
  // 98.1456953642  0.0154940397351
  // 56.6225165563  0.00493774834437
  // 17.6158940397  0.00037880794702
  // Ref multiplicity for centrality 0.000000 to 3.000000: 55.284322
  // Ref multiplicity for centrality 3.000000 to 10.000000: 41.124550
  // Ref multiplicity for centrality 10.000000 to 50.000000: 24.334303
  // Ref multiplicity for centrality 50.000000 to 100.000000: 8.008294

//   AddPoint(graphs[2][0], 55.3, 0.0254622516556, 0, 0);
  AddPoint(graphs[2][0], 40.649288, 0.0154940397351 * 0.5 * 2 / 0.985, 0, 0);
  AddPoint(graphs[2][0], 23.783224, 0.00493774834437 * 0.29 * 2 / 0.506, 0, 0);

  WriteGraphs("graphs_cms.root");
}

void ALICECMSMethod()
{
  CreateGraphStructure();
  
  AddPoint(graphs[2][0], 10, 0.027, 0, 0);
  AddPoint(graphs[4][0], 10, 0.115 / 1.34 * 1, 0, 0);
  AddPoint(graphs[5][0], 10, 0.019 / 0.152 * 0.126, 0, 0);
  
  AddPoint(graphs[2][4], 10, 0.09, 0, 0);
//   AddPoint(graphs[4][4], 10, 0.115, 0, 0);
  AddPoint(graphs[5][4], 10, 0.137, 0, 0);

  AddPoint(graphs[2][5], 10, 0.04, 0, 0);
//   AddPoint(graphs[4][5], 10, 0.046, 0, 0);
  AddPoint(graphs[5][5], 10, 0.06, 0, 0);

  WriteGraphs("graphs_alice_cmsmethod.root");
}

void GetSystematic()
{
  fileProj = TFile::Open("dphi_proj.root", "RECREATE");
  fileProj->Close();

  gROOT->SetBatch(kTRUE);
  
  if (1)
  {
    const char* baseFile = "dphi_corr_pA_121119_3.root";
//     const char* baseFile = "dphi_corr_pA_121108_hybrid_corrected.root";
    TString fileTag = "graphs_121119";
  }
  else if (0)
  {
    const char* baseFile = "dphi_corr_pA_121029_zna.root";
    TString fileTag = "graphs_121113_zna_uncorrected";
  }
  else if (0)
  {
    // needs dphi setting...
    const char* baseFile = "dphi_corr_pA_121119_cms_2.root";
    TString fileTag = "graphs_cms_121119";
  }
  else if (1)
  {
    const char* baseFile = "dphi_corr_pA_121121_cmsmethod.root";
    TString fileTag = "graphs_cmmethod_121121";
  }
  
  CorrelationSubtractionHistogram(baseFile, kTRUE, fileTag + ".root");
  
  gStudySystematic = 10;
  CorrelationSubtractionHistogram(baseFile, kTRUE, fileTag + "_exclusion05.root");
  
  gStudySystematic = 11;
  CorrelationSubtractionHistogram(baseFile, kTRUE, fileTag + "_exclusion00.root");

  gStudySystematic = 12;
  CorrelationSubtractionHistogram(baseFile, kTRUE, fileTag + "_exclusionAS.root");

  gStudySystematic = 13;
  CorrelationSubtractionHistogram(baseFile, kTRUE, fileTag + "_exclusionScale.root");

  gStudySystematic = 14;
  CorrelationSubtractionHistogram(baseFile, kTRUE, fileTag + "_exclusion12.root");
  
//   return;
  
  gStudySystematic = 20;
  CorrelationSubtractionHistogram(baseFile, kTRUE, fileTag + "_nonclosure.root");

  gStudySystematic = 30;
  CorrelationSubtractionHistogram(baseFile, kTRUE, fileTag + "_baseline.root");
  
  gStudySystematic = 40;
  CorrelationSubtractionHistogram("dphi_corr_pA_121116_global.root", kTRUE, fileTag + "_trackcuts.root");

  gStudySystematic = 60;
  CorrelationSubtractionHistogram(baseFile, kTRUE, fileTag + "_otherperipheral.root");
}

void DrawSystematics(TString fileTag = "graphs_121119")
{
  if (1)
  {
    const Int_t n = 6;
    const char* filesBase[] = { "", "_exclusion05", "_exclusion12", "_exclusion00", "_exclusionAS", "_exclusionScale" };
  }
  else
  {
    const Int_t n = 5;
    const char* filesBase[] = { "", "_trackcuts", "_nonclosure", "_baseline", "_otherperipheral" };
  }
  
  const char* files[n];
  for (Int_t i=0; i<n; i++)
  {
    str = new TString;
    str->Form("%s%s.root", fileTag.Data(), filesBase[i]);
    files[i] = str->Data();
  }
  
  const Int_t plots = 6;
  Int_t ids[7] = { 0, 1, 2, 3, 4, 5, 7 };
  
  for (Int_t i=0; i<plots; i++)
  {
    DrawSeveral(n, files, ids[i]);
//     break;
  }
  
  new TCanvas;
  for (Int_t i=0; i<n; i++)
    DrawLatex(0.2, 0.9 - 0.1 * i, colors[i], files[i]);
}

void DrawSeveral(const char* graphFile1, const char* graphFile2, const char* graphFile3, Int_t id)
{
  const char* files[3] = { graphFile1, graphFile2, graphFile3 };
  Int_t n = 1;
  if (graphFile2)
    n++;
  if (graphFile3)
    n++;
  
  DrawSeveral(n, files, id);
}
  
void DrawYield(const char* graphFile)
{
  DrawGraph(graphFile, 0, 1, "Ridge yield per #Delta#eta", "fig4b.eps", kTRUE);
}

void DrawRMS(const char* graphFile)
{
  DrawGraph(graphFile, 2, 3, "#sigma", 0);
}

void Drawv2v3(const char* graphFile)
{
  DrawGraph(graphFile, 4, 5, "#it{v}_{2} , #it{v}_{3}", "fig4a.eps");
}

void AddSystUnc(TGraphErrors* graph, Float_t syst020, Float_t syst2060)
{
  for (Int_t j=0; j<graph->GetN(); j++)
  {
    Float_t syst = syst2060;
    if (graph->GetX()[j] < 20)
      syst = syst020;
  
//     Printf("%d %f", j, syst);
    graph->GetEY()[j] = TMath::Sqrt(graph->GetEY()[j] * graph->GetEY()[j] + syst * syst * graph->GetY()[j] * graph->GetY()[j]);
  }
}

void SetSystUnc(TGraphErrors* graph, Float_t syst020, Float_t syst2060)
{
  for (Int_t j=0; j<graph->GetN(); j++)
  {
    Float_t syst = syst2060;
    if (graph->GetX()[j] < 20)
      syst = syst020;
  
//     Printf("%d %f", j, syst);
    graph->GetEY()[j] = syst * graph->GetY()[j];
  }
}

void SetXError(TGraphErrors* graph, Float_t value)
{
  for (Int_t j=0; j<graph->GetN(); j++)
    graph->GetEX()[j] = value;
}

void DrawGraph(const char* graphFile, Int_t id1, Int_t id2, const char* yLabel = 0, const char* outputFileName, Bool_t corrGraph = kFALSE)
{
  ReadGraphs(graphFile);
  TGraphErrors*** graphsSyst = graphs;
  ReadGraphs(graphFile);
  TGraphErrors*** graphsCombinedError = graphs;
  ReadGraphs(graphFile);
  
  if (yLabel == 0)
    yLabel = graphs[0][id1]->GetYaxis()->GetTitle();
  
  Float_t yMax[] = { 0.12, 0.12, 2, 2, 0.2, 0.25, 4, 1.5, 50, 50, 4 };
  Int_t markers[] = { 20, 21, 22, 23, 29, 33 };
  Int_t markers2[] = { 24, 25, 26, 32, 30, 27 };
  Float_t markerSize[] = { 1.7, 1.7, 1.7, 1.7, 2.0, 2.0 };

  if (1)
  {
    // default systs (for all pT, centr)
    Float_t syst020Array[] = { 0.16, 0.18, 0.15, 0.15, 0.14, 0.23 };
    Float_t syst2060Array[] = { 0.23, 0.25, 0.23, 0.23, 0.18, 0.42 };
    
    // 0.5<1 0.5<1
    Float_t syst020ArrayGr0[] = { 0.23, 0.25, 0.15, 0.15, 0.14, 0.23 };
    Float_t syst2060ArrayGr0[] = { 0.42, 0.43, 0.23, 0.23, 0.18, 0.42 };

    for (Int_t i=0; i<NGraphs; i++)
    {
      Float_t syst020 = syst020Array[id1];
      Float_t syst2060 = syst2060Array[id1];
      
      if (i == 0)
      {
	syst020 = syst020ArrayGr0[id1];
	syst2060 = syst2060ArrayGr0[id1];
      }
      
//       Printf("%d %d %f %f", i, id1, syst020, syst2060);
      SetSystUnc(graphsSyst[i][id1], syst020, syst2060);
      AddSystUnc(graphsCombinedError[i][id1], syst020, syst2060);
      SetXError(graphsSyst[i][id1], 0.3);

      Float_t syst020 = syst020Array[id2];
      Float_t syst2060 = syst2060Array[id2];
      
      if (i == 0)
      {
	syst020 = syst020ArrayGr0[id2];
	syst2060 = syst2060ArrayGr0[id2];
      }

      SetSystUnc(graphsSyst[i][id2], syst020, syst2060);
      AddSystUnc(graphsCombinedError[i][id2], syst020, syst2060);
      SetXError(graphsSyst[i][id2], 0.3);
    }
    
//     graphs[0][id1]->Print();
  }
  
  const char* eventClass[] = { "0-20%", "20-40%", "40-60%" };
  Printf("\n%s:", graphTitles[id1]);
  for (Int_t i=0; i<3; i++)
  {
    Printf("Event class: %s minus 60-100%%", eventClass[i]);
    
    for (Int_t j=0; j<NGraphs; j++)
      if (graphs[j][id1]->GetN() > i)
	Printf("%s: %.4f +- %.4f (stat) +- %.4f (syst)", graphs[j][id1]->GetTitle(), graphs[j][id1]->GetY()[i], graphs[j][id1]->GetEY()[i], graphsSyst[j][id1]->GetEY()[i]);
  }
    
  Printf("\n%s:", graphTitles[id2]);
  for (Int_t i=0; i<3; i++)
  {
    Printf("Event class: %s minus 60-100%%", eventClass[i]);
    
    for (Int_t j=0; j<NGraphs; j++)
      if (graphs[j][id2]->GetN() > i)
	Printf("%s: %.4f +- %.4f (stat) +- %.4f (syst)", graphs[j][id2]->GetTitle(), graphs[j][id2]->GetY()[i], graphs[j][id2]->GetEY()[i], graphsSyst[j][id2]->GetEY()[i]);
  }

  TCanvas* canvas = new TCanvas;
  gPad->SetTopMargin(0.03);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.13);
//   gPad->SetGridx();
//   gPad->SetGridy();
  TH2F* dummy = new TH2F("dummy", Form(";%s;%s", graphs[0][id1]->GetXaxis()->GetTitle(), yLabel), 3, 0, 60, 100, 0, yMax[id1]);
  dummy->GetYaxis()->SetTitleOffset(1.1);
  dummy->SetStats(0);
  dummy->GetYaxis()->SetNdivisions(505);
  dummy->GetXaxis()->SetLabelSize(0.06);
  dummy->GetYaxis()->SetLabelSize(0.06);
  dummy->GetXaxis()->SetTitleSize(0.06);
  dummy->GetYaxis()->SetTitleSize(0.06);
  dummy->Draw();
  
  if (strcmp(graphs[0][id1]->GetXaxis()->GetTitle(), "Event class") == 0)
  {
    dummy->GetXaxis()->SetBinLabel(1, "0-20%");
    dummy->GetXaxis()->SetBinLabel(2, "20-40%");
    dummy->GetXaxis()->SetBinLabel(3, "40-60%");
    dummy->GetXaxis()->SetLabelSize(0.09);
    dummy->GetXaxis()->SetTitleOffset(1);
  }
  
  legend = new TLegend((id1 == 4) ? 0.47 : 0.33, (id1 == 4) ? 0.70 : 0.55, 0.95, 0.92);
  legend->SetNColumns(2);
  if (id1 == 0 || id1 == 2)
    legend->SetHeader("Near side   Away side");
  else if (id1 == 4)
    legend->SetHeader("  #it{v}_{2}       #it{v}_{3}");
    
  legend->SetFillColor(0);
  legend->SetBorderSize(0);

  if (1)
  {
    Int_t fillStyle[11] = { 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008, 3009, 3010, 3011 };
    
    for (Int_t i=0; i<NGraphs; i++)
    {
      graphsSyst[i][id1]->SetMarkerStyle(1);
      graphsSyst[i][id1]->SetMarkerColor(1);
      graphsSyst[i][id1]->SetLineColor(1);
      graphsSyst[i][id1]->SetFillColor(1);
      graphsSyst[i][id1]->SetFillStyle(3001);
      graphsSyst[i][id1]->Draw("2SAME");

      graphsSyst[i][id2]->SetMarkerStyle(1);
      graphsSyst[i][id2]->SetMarkerColor(2);
      graphsSyst[i][id2]->SetLineColor(2);
      graphsSyst[i][id2]->SetFillColor(2);
      graphsSyst[i][id2]->SetFillStyle(3001);
      graphsSyst[i][id2]->Draw("2SAME");
    }
  }
  else
    Printf(">>>>>>>>>>>> SKIPPING SYST");

  
  for (Int_t i=0; i<NGraphs; i++)
  {
    graphs[i][id1]->SetMarkerStyle(markers[i]);
    graphs[i][id1]->SetMarkerSize(markerSize[i]);
    graphs[i][id1]->SetMarkerColor(1);
    graphs[i][id1]->SetLineColor(1);
    graphs[i][id1]->Draw("PSAME");

    graphs[i][id2]->SetMarkerStyle(markers2[i]);
    graphs[i][id2]->SetMarkerSize(markerSize[i]);
    graphs[i][id2]->SetMarkerColor(2);
    graphs[i][id2]->SetLineColor(2);
    graphs[i][id2]->Draw("PSAME");

    TString label(graphs[i][id1]->GetTitle());
    label.ReplaceAll(".00", ".0");
    label.ReplaceAll(".50", ".5");
    
    
/*    label.ReplaceAll(".0", " GeV/#it{c}");
    label.ReplaceAll(".5", ".5 GeV/#it{c}");*/
    TObjArray* objArray = label.Tokenize("-");
    label.Form("%s-%s", objArray->At(0)->GetName(), objArray->At(1)->GetName());
    label.ReplaceAll("-", ";");
    
    if (id1 == 4)
    {
      // reduce label
      label.ReplaceAll(" ", "");
      label.ReplaceAll(";", "<");
      tokens = label.Tokenize("<");
      label.Form("%s < %s < %s < %s ", tokens->At(0)->GetName(), tokens->At(1)->GetName(), tokens->At(4)->GetName(), tokens->At(5)->GetName());
    }
    
    label.ReplaceAll("p_", "#it{p}_");
    label += "GeV/#it{c}";

    if (graphs[i][id1]->GetN() > 0)
    {
      legend->AddEntry(graphs[i][id1], (id1 == 4) ? " " : "    ", "P");
      legend->AddEntry(graphs[i][id2], label, "P");
    }
  }
  
  legend->Draw();
  DrawLatex(0.7, 0.92, 1, "p-Pb #sqrt{s_{NN}} = 5.02 TeV", 0.04);
  
  canvas->SaveAs(outputFileName);
  
  if (!corrGraph)
    return;
  
  corr = new TGraphErrors;
  for (Int_t i=0; i<NGraphs; i++)
  {
    for (Int_t j=0; j<graphs[i][id1]->GetN(); j++)
      AddPoint(corr, graphsCombinedError[i][id1]->GetY()[j], graphsCombinedError[i][id2]->GetY()[j], graphsCombinedError[i][id1]->GetEY()[j], graphsCombinedError[i][id2]->GetEY()[j]);
  }

  new TCanvas;
  gPad->SetTopMargin(0.03);
  gPad->SetRightMargin(0.02);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.13);
  TString titleString;
  titleString.Form(";%s;%s", graphs[0][id1]->GetYaxis()->GetTitle(), graphs[0][id2]->GetYaxis()->GetTitle());
  titleString.ReplaceAll("NS", "Near-side");
  titleString.ReplaceAll("AS", "Away-side");
  titleString.ReplaceAll("Ridge Yield", "ridge yield per #Delta#eta");
  dummy = new TH2F("dummy2", titleString, 100, 0, 0.095, 100, 1e-5, 0.095);
  dummy->GetYaxis()->SetTitleOffset(1.1);
  dummy->SetStats(0);
  dummy->GetXaxis()->SetNdivisions(505);
  dummy->GetYaxis()->SetNdivisions(505);
  dummy->GetXaxis()->SetLabelSize(0.06);
  dummy->GetYaxis()->SetLabelSize(0.06);
  dummy->GetXaxis()->SetTitleSize(0.06);
  dummy->GetYaxis()->SetTitleSize(0.06);
  dummy->Draw();

  corr->Draw("PSAME");
  
  line = new TLine(0, 0, 0.095, 0.095);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw();

  DrawLatex(0.18, 0.92, 1, "p-Pb #sqrt{s_{NN}} = 5.02 TeV", 0.04);
}

void PaperCorrFuncAll(const char* fileName)
{
  Int_t n = 6;
  Int_t is[] = { 0, 1, 1, 2, 2, 2, 3 };
  Int_t js[] = { 1, 1, 2, 1, 2, 3, 3 };
  Int_t centr[] = { 0, 1, 3, 4 };

  for (Int_t i=0; i<n; i++)
    for (Int_t j=0; j<4; j++)
      PaperCorrFunc(fileName, is[i], js[i], centr[j]);
}

void PaperCorrFunc(const char* fileName, Int_t i = 2, Int_t j = 2, Int_t centr = 0)
{
  TFile::Open(fileName);
  
  TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, centr));

  if (!hist1)
    return 0;
  // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
  hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));

  hist1 = (TH2*) hist1->Clone();
  hist1->Rebin2D(2, 2); hist1->Scale(0.25);
//   hist1->Rebin2D(2, 1); hist1->Scale(0.5);
  
  hist1->GetYaxis()->SetRangeUser(-1.99, 1.99);
  hist1->GetXaxis()->SetTitleOffset(1.5);
  hist1->GetYaxis()->SetTitleOffset(2);
  hist1->GetZaxis()->SetTitle(kCorrFuncTitle);
  hist1->SetStats(kFALSE);
  hist1->GetZaxis()->SetTitleOffset(2);
  hist1->GetZaxis()->SetNdivisions(504);
  hist1->SetStats(kFALSE);
  hist1->GetZaxis()->CenterTitle(kTRUE);
  hist1->GetYaxis()->CenterTitle(kTRUE);
  hist1->GetXaxis()->CenterTitle(kTRUE);
  hist1->GetYaxis()->SetNdivisions(505);
  hist1->GetXaxis()->SetTitle("#Delta#varphi (rad)");

  c = new TCanvas("c", "c", 600, 600);
  gPad->SetPad(0, 0, 1, 1);
  gPad->SetLeftMargin(0.2);

  TString label(hist1->GetTitle());
  hist1->SetTitle("");
  label.ReplaceAll(".00", "");
  label.ReplaceAll(".0", "");
//   label.ReplaceAll(".00", " GeV/#it{c}");
//   label.ReplaceAll(".0", " GeV/#it{c}");
  TObjArray* objArray = label.Tokenize("-");
  TPaveText* paveText = new TPaveText(0.03, 0.86, 0.44, 0.97, "BRNDC");
  paveText->SetTextSize(0.035);
  paveText->SetFillColor(0);
  paveText->SetShadowColor(0);
  paveText->SetBorderSize(0);
  paveText->SetFillStyle(0);
  
  TString tmpStr(objArray->At(0)->GetName());
  tmpStr.ReplaceAll("p_", "#it{p}_");
  paveText->AddText(tmpStr + "GeV/#it{c}");

  TString tmpStr(objArray->At(1)->GetName());
  tmpStr.ReplaceAll("p_", "#it{p}_");
  paveText->AddText(tmpStr + "GeV/#it{c}");
//   paveText->AddText(objArray->At(1)->GetName());

  TPaveText* paveText2 = new TPaveText(0.63, 0.86, 0.97, 0.97, "BRNDC");
  paveText2->SetTextSize(0.035);
  paveText2->SetBorderSize(0);
  paveText2->SetFillColor(0);
  paveText2->SetFillStyle(0);
  paveText2->SetShadowColor(0);
  if (objArray->GetEntries() == 4)
  {
    paveText2->AddText(Form("p-Pb #sqrt{s_{NN}} = 5.02 TeV"));
    paveText2->AddText(Form("%s-%s", objArray->At(2)->GetName(), objArray->At(3)->GetName()));
  }
  else
    paveText2->AddText(Form("%s #sqrt{s} = 2.76 TeV", objArray->At(2)->GetName()));

  hist1->Draw("SURF1");
  paveText->Draw();
  paveText2->Draw();
  
  if (i == 2 && j == 2 && centr == 0)
    c->SaveAs("fig1b.eps");
  else if (i == 2 && j == 2 && centr == 4)
    c->SaveAs("fig1a.eps");
  
  c->SaveAs(Form("corr_%d_%d_%d.eps", i, j, centr));
}

void ExtractSystematics(const char* fileName = "dphi_proj.root")
{
  Int_t i = 2; Int_t j = 2;
//   Int_t i = 0; Int_t j = 1;
  Int_t centr = 0;
  
  Int_t nSyst = 8;
  Int_t syst[] = { 0, 10, 11, 12, 13, 20, 30, 40 };
  
  TFile::Open(fileName);
  base = (TH1*) gFile->Get(Form("proj_%d_%d_%d_%d", i, j, centr, syst[0]));
  
  c = new TCanvas;
//   base->Rebin(2); base->Scale(0.5);
  base->Draw();
  
  Float_t baseValue = base->Integral() / base->GetNbinsX();
//   (base->GetBinContent(1) + base->GetBinContent(base->GetNbinsX()/2+1)) / 2;
  
  c2 = new TCanvas;
  c3 = new TCanvas;
  
  Int_t color = 2;

  for (Int_t n = 1; n<nSyst; n++)
  {
    hist = (TH1*) gFile->Get(Form("proj_%d_%d_%d_%d", i, j, centr, syst[n]));
    if (!hist)
      continue;
    
//     hist->Rebin(2); hist->Scale(0.5);
    c->cd();
    hist->SetLineColor(color++);
    Float_t baseValue2 = hist->Integral() / hist->GetNbinsX();
    Printf("%f %f", baseValue, baseValue2);
//     hist->Add(new TF1("func", "1", -5, 5), baseValue - baseValue2);
    hist->Scale(baseValue / baseValue2);
    hist->DrawCopy("SAME");
    
    c2->cd();
    hist->DrawCopy((n == 1) ? "HIST" : "HIST SAME")->Divide(base);

    c3->cd();
    hist->Add(base, -1);
    hist->DrawCopy((n == 1) ? "HIST" : "HIST SAME")->GetYaxis()->SetRangeUser(-0.1, 0.1);
  }
}

void GetProjectionSystematics(Float_t eta)
{
  fileProj = TFile::Open("dphi_proj2.root", "RECREATE");
  fileProj->Close();

  Int_t i = 2; Int_t j = 2;
//   Int_t i = 0; Int_t j = 1;

  Int_t n = 5;
  const char* files[] = { "dphi_corr_pA_121119_3.root", "dphi_corr_pA_121116_global.root", "dphi_corr_pA_121119_3.root", "dphi_corr_pA_121116_hybrid.root" /*pair cuts*/, "dphi_corr_pA_121119_hijing.root" };
  
  for (Int_t centr=0; centr<6; centr++)
  {
    for (Int_t k=0; k<n; k++)
    {
      gStudySystematic = 0;
      if (k == 2)
	gStudySystematic = 20;
      
      TFile::Open(files[k]);
      const char* label = 0;
      TH1* hist = GetProjections(i, j, centr, &label, eta);
      if (!hist)
	continue;
      
      fileProj = TFile::Open("dphi_proj2.root", "UPDATE");
      hist->Write(Form("proj_%d_%d_%d_%d", i, j, centr, k*10));
      fileProj->Close();
    }
  }
}

void CMSPlot()
{
  Int_t centrArr[] = { 0, 1, 3, 4 };
  Int_t nCentr = 4;
  
//   Int_t i = 1; Int_t j = 2;
  Int_t i = 2; Int_t j = 2;

  if (1)
  {
    const char* labels[] = { "Data", "HIJING", "DPMJET" };
  //   const char* files[] = { "dphi_corr_pA_121119_3.root", "dphi_corr_121121_hijing_uncorrected.root", "dphi_corr_121121_dpmjet_uncorrected.root" };
//     const char* files[] = { "dphi_corr_pA_121119_3.root", "dphi_corr_121121_hijing_step0.root", "dphi_corr_121121_dpmjet_step0.root" };
    const char* files[] = { "dphi_corr_pA_121121_cmsmethod.root", "dphi_corr_121122_cmsmethod_hijing_step0.root", "dphi_corr_121122_cmsmethod_dpmjet_step0.root" };
//     const char* files[] = { "dphi_corr_pA_121122_cnd_cmsmethod.root", "dphi_corr_121122_cmsmethod_hijing_cnd_step0.root", 0 };
    Int_t nFiles = 3;
  }
  else
  {
    const char* labels[] = { "ALICE method", "CMS method" };
//     const char* files[] = { "dphi_corr_pA_121119_3.root", "dphi_corr_pA_121121_cmsmethod.root" };
    const char* files[] = { "dphi_corr_121121_hijing_step0.root", "dphi_corr_121122_cmsmethod_hijing_step0.root" };
    Int_t nFiles = 2;
  }
  
  Int_t colors[] = { 1, 2, 4 };
  
  c = new TCanvas("c", "c", 800, 800);
  c->Divide(2, 2);
  
  Float_t maxValue = -1;
  for (Int_t centr = 0; centr < nCentr; centr++)
  {
    c->cd(centr+1);
    
    legend = new TLegend(0.15, 0.55, 0.46, 0.85);
    legend->SetFillColor(0);
    
    for (Int_t fileId = 0; fileId < nFiles; fileId++)
    {
      if (files[fileId] == 0)
	continue;
      
      TFile::Open(files[fileId]);
  
      TH2* hist1 = (TH2*) gFile->Get(Form("dphi_%d_%d_%d", i, j, centrArr[centr]));
      if (!hist1)
	return 0;
      hist1 = (TH2*) hist1->Clone(Form("%s_%.1f", hist1->GetName(), 0.0));

      // NOTE fix normalization. these 2d correlations come out of AliUEHist normalized by dphi bin width, but not deta
      hist1->Scale(1.0 / hist1->GetYaxis()->GetBinWidth(1));
      
      tokens = TString(hist1->GetTitle()).Tokenize("-");
      centralityStr = new TString;
      if (tokens->GetEntries() > 2)
	*centralityStr = tokens->At(2)->GetName();
      if (tokens->GetEntries() > 3)
	*centralityStr = *centralityStr + "-" + tokens->At(3)->GetName();
      
      Float_t etaMin = 1.0;
      
      proj1x = hist1->ProjectionX(Form("proj1x_%d_%d_%d_%d", i, j, centr, fileId), hist1->GetYaxis()->FindBin(-etaMax+0.01), hist1->GetYaxis()->FindBin(-etaMin-0.01));
      proj2x = hist1->ProjectionX(Form("proj2x_%d_%d_%d_%d", i, j, centr, fileId), hist1->GetYaxis()->FindBin(etaMin+0.01), hist1->GetYaxis()->FindBin(etaMax-0.01));
      proj1x->Add(proj2x);
      
      proj1x->Scale(hist1->GetYaxis()->GetBinWidth(1));
      proj1x->Scale(1.0 / (etaMax - etaMin) / 2);

      proj1x->Rebin(2); proj1x->Scale(0.5);
      proj1x->Rebin(2); proj1x->Scale(0.5);
      
      // zyam
      Float_t zyam = proj1x->Integral(proj1x->FindBin(TMath::Pi()/2 - 0.5), proj1x->FindBin(TMath::Pi()/2 + 0.5)) / (proj1x->FindBin(TMath::Pi()/2 + 0.5) - proj1x->FindBin(TMath::Pi()/2 - 0.5) + 1);
      proj1x->Add(new TF1("f", "1", -5, 5), -zyam);

      proj1x->SetStats(kFALSE);
      proj1x->GetXaxis()->SetTitleOffset(1);
      
      proj1x->SetLineColor(colors[fileId]);
      if (maxValue < 0)
	maxValue = proj1x->GetMaximum() * 2;
      proj1x->SetMaximum(maxValue);
//       proj1x->SetMinimum(-0.01);
      proj1x->Draw((fileId == 0) ? "" : "SAME");
      
      legend->AddEntry(proj1x, labels[fileId]);
    }
    legend->Draw();
  }
}

void FourierFactorization(const char* fileName)
{
  ReadGraphs(fileName);
  
  Int_t n = 6;
  Int_t is[] =    { 0, 1, 1, 2, 2, 2, 3 };
  Int_t js[] =    { 1, 1, 2, 1, 2, 3, 3 };
  Bool_t symm[] = { 1, 0, 1, 0, 0, 1, 0 };

  Int_t graphID = 5;
  
  TGraphErrors** symmGraphs[] = { graphs[0], graphs[2], graphs[5] };
  
  for (Int_t i=0; i<symmGraphs[0][graphID]->GetN(); i++)
  {
    Printf("%d", i);
    graphs[1][graphID]->GetY()[i] = TMath::Sqrt(symmGraphs[0][graphID]->GetY()[i] * symmGraphs[1][graphID]->GetY()[i]);
    graphs[3][graphID]->GetY()[i] = TMath::Sqrt(symmGraphs[2][graphID]->GetY()[i] * symmGraphs[1][graphID]->GetY()[i]);
    graphs[4][graphID]->GetY()[i] = TMath::Sqrt(symmGraphs[2][graphID]->GetY()[i] * symmGraphs[2][graphID]->GetY()[i]);
  }
  
//   graphs[1][graphID]->Draw("A*");
  
  WriteGraphs();
}

void CompareATLAS()
{
  atlas = ReadHepdata("/home/jgrosseo/Dropbox/alice-paper-paridge/comparison_atlas/atlas_figaux10.txt", kFALSE, 1, 1);
  atlas->SetMarkerStyle(20);
  atlas->Draw("PA");
  atlas->SetTitle("");
  atlas->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  atlas->GetYaxis()->SetTitle("s_{2} / v_{2}");
  atlas->GetYaxis()->SetRangeUser(0.02, 0.18);

  alice = new TGraphErrors;
  AddPoint(alice, 0.75, 0.0583503, 0.25, 0.00831755);
  AddPoint(alice, 1.5, 0.0953562, 0.5, 0.0134597);
  AddPoint(alice, 3, 0.128309, 1, 0.0189634);
  
  alice->SetMarkerStyle(21);
  alice->SetLineColor(2);
  alice->SetMarkerColor(2);
  alice->Draw("PSAME");
}
