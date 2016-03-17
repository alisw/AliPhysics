//====================================================================
/** 
 * Find bin number correspondig to centrality range 
 * 
 * @param c1 Least centrality 
 * @param c2 Largest centrality 
 * 
 * @return Bin number 
 */
static Int_t PbPbBin(Double_t c1, Double_t c2)
{
  Double_t c = (c1+c2) / 2;
  if      (c <  5) return 0;
  else if (c < 10) return 1;
  else if (c < 20) return 2;
  else if (c < 30) return 3;
  else if (c < 40) return 4;
  else if (c < 50) return 5;
  else if (c < 60) return 6;
  else if (c < 70) return 7;
  else if (c < 80) return 8;
  else if (c < 90) return 9;
  return                  10;
}
  
//____________________________________________________________________
/** 
 * Get the centrality color for PbPb 
 * 
 * @param c1 Lower edge
 * @param c2 Upper edge
 * 
 * @return Color 
 */
static Color_t PbPbColor(Double_t c1, Double_t c2)
{
  const Color_t cc[] = { kMagenta+2,
			 kBlue+2,
			 kAzure-1, // 10,
			 kCyan+2,
			 kGreen+1,
			 kSpring+5,//+10,
			 kYellow+1,
			 kOrange+5,//+10,
			 kRed+1,
			 kPink+5,//+10,
			 kBlack };
  Int_t bin = PbPbBin(c1,c2);
  return cc[bin];
}

//====================================================================
Double_t SysEval(Double_t x, Double_t sMin, Double_t sMax, Double_t xMax)
{
  return sMin + TMath::Power(x/xMax, 2)*(sMax-sMin);
}

Double_t CSysEval(Double_t x, Double_t sMin, Double_t sMax)
{
  return SysEval(x, sMin, sMax, 80);
}
Double_t EtaSysEval(Double_t x, Double_t sMin, Double_t sMax)
{
  return SysEval(x, sMin, sMax, 2);
}

void MakeGSE(std::ostream& o, const TH1* h,
	     Double_t c1, Double_t c2, Bool_t reweigh, Double_t fac=1)
{
  // These are correlated
  Double_t bgMid   = CSysEval(c2, 100*0.02,  100*0.002);
  Double_t cSys    = CSysEval(c2, 100*0.005, 100*0.075);
  Double_t strMid  = 100*0.005;
  Double_t strFwd  = 100*0.05;
  Double_t pidMid  = 100*0;
  Double_t pidFwd  = 100*0.01;
  // Double_t bgSys = (1-c1/100)*(2-0.2) + 0.2;
  // Double_t cSys  = TMath::Power(c1/100,2)*(7.5-0.5) + 0.5;
  
  o << "*dataset:\n"
    << "*dscomment: The pseudo-rapidity density of charged particle\n"
    << "*reackey: PB PB --> CHARGED X\n"
    << "*obskey: DN/DETARAP\n"
    << "*qual: CENTRALITY IN PCT : " << c1 << " TO " << c2 << "\n"
    << "*qual: SQRT(S)/NUCLEON IN GEV : 5023\n";
  if (!reweigh) {
    o << "*dserror: " << strMid << " PCT : Weak decays\n"
      << "*dserror: " << pidMid << " PCT : Particle composition\n"
      << "*dserror: " << bgMid  << " PCT : Background subtraction\n";
  }
  o << "*dserror: 1   PCT : pT extrapolation\n"
    << "*dserror: " << cSys << " PCT : Centrality\n"
    << "*xheader: ETARAP\n"
    << "*yheader: DN/DETARAP\n"
    << "*data: x : y" << std::endl;
  // Define points
  TAxis* xa = h->GetXaxis();
  for (Int_t i = 0; i < h->GetNbinsX(); i++) {
    Int_t    j  = i+1;
    Double_t cc = h->GetBinContent(j)*fac;
    Double_t x  = h->GetXaxis()->GetBinCenter(j);
    Double_t ex = h->GetXaxis()->GetBinWidth(j)/2;
    Double_t xo = TMath::Abs(x)+ex;
    if (cc < 1e-8) continue;
    o << " " << xa->GetBinLowEdge(j) << " TO " << xa->GetBinUpEdge(j)
      << "; " << cc << " +-" << h->GetBinError(j)
      << " (DSYS=" << 0.01+0.01*TMath::Power(xo/2,2) << " PCT:Acceptance";
    if (reweigh) {
      Double_t bg     = EtaSysEval(xo, bgMid, 3*bgMid);
      Double_t pid    = EtaSysEval(xo, pidMid, pidFwd);
      Double_t str    = EtaSysEval(xo, strMid, strFwd);
      o << ",DSYS=" << bg  << " PCT:Background subtraction"
	<< ",DSYS=" << pid << " PCT:Particle composition"
	<< ",DSYS=" << str << " PCT:Weak decay";
      
    }
    o << ");"
      << std::endl;
  }
  o << "*dataend:\n" << std::endl;
}
void FindLeastLargest(TObject* o, Double_t c1, Double_t c2)
{
  Double_t eMin = +10000;
  Double_t eMax = -10000;
  GraphSysErr* g = static_cast<GraphSysErr*>(o);
  for (Int_t i = 0; i < g->GetN(); i++) {
    Double_t eyl, eyh;
    Double_t y  = g->GetYandError(i, true, false, true, false, eyl, eyh);
    Double_t ey = TMath::Max(eyl,eyh)/y;
    eMin        = TMath::Min(ey, eMin);
    eMax        = TMath::Max(ey, eMax);
  }
  Printf("%4.1f - %4.1f%%:  Least: %6.2f%%, Largest %6.2f%%",
	 c1, c2, eMin*100, eMax*100);
}

/** 
 * Extract ALICE PbPb @ 5.02TeV over |eta|<2
 * 
 * @param filename 
 * @param outname 
 */
void
Extract(const char* filename="dndneta.pbpb502.20151124.root",
	const char* outname="TRACKLETS_5023_PbPb.input",
	Bool_t      reweigh=false)
{
  if (filename == 0) return;
  TFile* file = TFile::Open(filename, "READ");
  TObjArray* arr = static_cast<TObjArray*>(file->Get("TObjArray"));
  // Now count number of bins
  Int_t nBins = 0;
  TIter next(arr);
  TObject* obj = 0;
  while ((obj = next())) {
    if (TString(obj->GetName()).Contains("DataCorrSignal")) 
      nBins++;
  }
  Info("ExtractdNdeta", "Defining %d centrality bins", nBins);
  TArrayD c(nBins+1);
  if (nBins == 5) {
    c[0] = 0; c[1] = 10; c[2] = 20; c[3] = 40; c[4] = 60; c[5] = 80;
  }
  else if (nBins >= 9) {
    c[0] =  0; c[1] =  5; c[2] = 10; c[3] = 20; c[4] = 30; c[5] = 40;
    c[6] = 50; c[7] = 60; c[8] = 70; c[9] = 80;
    if (nBins >= 10) c[10] =  90;
    if (nBins >= 11) c[11] = 100;
  }
  
  THStack* all = new THStack("all","all");
  std::ofstream out(outname);
  std::ostream& o = out; // std::cout;
  // std::ostream& o = std::cout;
  
  o << "*author: SHAHOYAN : 2015\n"
    << "*title: Full centrality dependence of the charged "
    << "particle pseudo-rapidity density over the widest "
    << "possible pseudo-rapidity range in Pb-Pb collisions "
    << "at 5.02 TeV\n"
    << "*detector: TRACKLETS\n"
    << "*experiment: CERN-LHC-TRACKLETS\n"
    << "*comment: CERN-LHC: We present the charged particle pseudo-rapidity "
    << "density of charged particles in Pb-Pb collisions at sqrt(s)/nucleon "
    "= 5.02 over the widest possible pseudo-rapidity and centrality range "
    << "possible.\n"  << std::endl;
  
  for (Int_t i = 0; i < nBins; i++) {
    TString hName = Form("bin%d_DataCorrSignal_PbPb",i);
    TH1* h = static_cast<TH1*>(arr->FindObject(hName));
    if (!h) {
      hName.ReplaceAll("PbPb", "PBPB");
      h = static_cast<TH1*>(arr->FindObject(hName));
      if (!h) {
	Warning("", "Histogram (%s) missing for bin %d", hName.Data(), i);
	arr->Print();
	continue;
      }
    }
    
    Color_t      col = PbPbColor(c[i], c[i+1]);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetFillColor(col);
    all->Add(h);
    Info("","Making GSE for %0d%% to %3d%% (%d)",
	 Int_t(c[i]), Int_t(c[i+1]), col);
    
    MakeGSE(o, h, c[i], c[i+1], reweigh);
  }
  // all->Draw("nostack");
  o << "*E" << std::endl;
  out.close();

  TCanvas*        can = new TCanvas("c","C", 1600, 800);
  can->SetRightMargin(0.2);
  can->SetTopMargin(0.01);

  TLegend*        cl = new TLegend(1-can->GetRightMargin(),
				   can->GetBottomMargin(),.99,
				   1-can->GetTopMargin());
  cl->SetFillStyle(0);
  cl->SetBorderSize(0);
  
  gROOT->LoadMacro("$HOME/GraphSysErr/GraphSysErr.C+");
  TList* ll = GraphSysErr::Import(outname);
  // ll->ls();

  TIter next(ll);
  TObject* obj = 0;
  Bool_t first = true;
  TH1* frame = 0;
  Double_t min=100000, max=0;
  Int_t i = 0;
  while ((obj = next())) {
    if (c[i+1] > 80) break;
    GraphSysErr* g = static_cast<GraphSysErr*>(obj);
    Color_t      col = PbPbColor(c[i], c[i+1]);
    TLegendEntry* e =  cl->AddEntry("", Form("%4.1f-%4.1f%%", c[i], c[i+1]),
				    "F");
    e->SetFillColor(col);
    e->SetFillStyle(1001);
    g->SetLineColor(col);
    g->SetMarkerColor(col);
    g->SetFillColor(col);
    // g->Print("qual");
    g->SetDataOption(GraphSysErr::kNoTick);
    g->SetSumOption(GraphSysErr::kBox);
    g->SetSumLineColor(col);
    g->SetSumFillColor(col);
    g->SetCommonSumOption(GraphSysErr::kBox);
    g->SetCommonSumLineColor(col);
    g->SetCommonSumFillColor(col);
    g->SetName(Form("tracklets%03dd%02d_%03dd%02d",
		    Int_t(c[i]),   Int_t(c[i]*100) % 100,
		    Int_t(c[i+1]), Int_t(c[i+1]*100) % 100));
    g->SetTitle(Form("%4.1f - %4.1f%%", c[i], c[i+1]));
    if (first) g->Draw("combine stat quad axis xbase=2.5");
    else       g->Draw("combine stat quad xbase=2.5");
    if (!frame)
      frame = g->GetMulti()->GetHistogram();
    first = false;
    Double_t mn, mx;
    g->GetMinMax("combine stat quad", mn, mx);
    FindLeastLargest(g, c[i], c[i+1]);
    min = TMath::Min(min, mn);
    max = TMath::Max(max, mx);
    i++;
  }
  frame->SetMinimum(min*.9);
  frame->SetMaximum(max*1.1);
  cl->Draw();

  TFile* outFile = TFile::Open(Form("PbPb5023midRapidity%s.root",
				    reweigh ? "Reweighed" : "Normal"),
			       "RECREATE");
  ll->Write("container",TObject::kSingleKey);
  outFile->Write();
  
  can->SaveAs(Form("PbPb5023midRapidity%s.png",
		   reweigh ? "Reweighed" : "Normal"));
}
//
// EOF
// 
