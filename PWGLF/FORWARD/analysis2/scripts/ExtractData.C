// --- Get a single histogram ----------------------------------------
TH1*
GetOne(const TList* list, const char* which, bool mirror)
{
  if (!list) { 
    Error("GetOne", "No list passed");
    return 0;
  }
  TString n(Form("dndeta%s_rebin05", which));
  if (mirror) n.Append("_mirror");
  
  TObject* o = list->FindObject(n);
  if (!o) { 
    Error("GetOne", "Object %s not found in %s", n.Data(), list->GetName());
    return 0;
  }
  TH1* ret = static_cast<TH1*>(o);
  ret->SetLineColor(ret->GetMarkerColor());
  ret->SetDirectory(0);
  return ret;
}

// --- Make point-to-point systematic errors -------------------------
TH1*
MakeSysError(TH1* h, Double_t sysErr)
{
  TString n(h->GetName());
  n.Append("_syserr");
  TH1* ret = static_cast<TH1*>(h->Clone(n));
  ret->SetMarkerStyle(1);
  ret->SetMarkerSize(0);
  ret->SetMarkerColor(kBlue-10);
  ret->SetFillColor(kBlue-10);
  ret->SetFillStyle(1001);
  ret->SetLineColor(kBlue-10);
  ret->SetLineWidth(0);
  ret->SetDirectory(0);
  
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++) { 
    Double_t c = ret->GetBinContent(i);
    if (c < 0.001) { 
      ret->SetBinContent(i,0);
      ret->SetBinError(i,0);
    }
    Double_t e = c * sysErr/100;
    ret->SetBinError(i,e);
  }
  
  return ret;
}

// --- Turn a TGraphAsymmErrors into a histogram ---------------------
TH1* Graph2Hist(const TGraphAsymmErrors* g)
{
  Int_t    nBins = g->GetN();
  TArrayF  bins(nBins+1);
  TArrayF  y(nBins);
  TArrayF  ey(nBins);
  Double_t dx = 0;
  Double_t xmin = 10000;
  Double_t xmax = -10000;
  for (Int_t i = 0; i < nBins; i++) { 
    Double_t x   = g->GetX()[i];
    Double_t exl = g->GetEXlow()[i];
    Double_t exh = g->GetEXhigh()[i];
    xmin             = TMath::Min(x-exl, xmin);
    xmax             = TMath::Max(x+exh, xmax);
    bins.fArray[i]   = x-exl;
    bins.fArray[i+1] = x+exh;
    Double_t dxi = exh+exl;
    if (dxi == 0 && i != 0) dxi = bins.fArray[i]-bins.fArray[i-1];
    if (dx == 0) dx  = dxi;
    else if (dxi != dx) dx = 0;
    
    y.fArray[i]  = g->GetY()[i];
    ey.fArray[i] = TMath::Max(g->GetEYlow()[i],g->GetEYhigh()[i]);

  }
  TString name(g->GetName());
  TString title(g->GetTitle());
  TH1D* h = 0;
  if (dx != 0) {
    h = new TH1D(name.Data(), title.Data(), nBins, 
		 bins[0]-dx/2, bins[nBins]+dx/2);
  }
  else {
    h = new TH1D(name.Data(), title.Data(), nBins, bins.fArray);
  }
  for (Int_t i = 1; i <= nBins; i++) { 
    h->SetBinContent(i, y.fArray[i-1]);
    h->SetBinError(i, ey.fArray[i-1]);
  }
  h->SetMarkerStyle(g->GetMarkerStyle());
  h->SetMarkerColor(g->GetMarkerColor());
  h->SetMarkerSize(g->GetMarkerSize());
  h->SetDirectory(0);
    
  return h;
}

// --- Set Histogram properties --------------------------------------
void SetAttributes(TH1* h, UShort_t sNN, Int_t which, Bool_t mirror)
{
  Int_t marker = (sNN == 900 ? 20     : 21) + (mirror ? 4 : 0);
  Int_t color  = (which == 0 ? kRed+2  :
		  which == 1 ? kMagenta+2 : 
		  kBlue + 2);
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetMarkerStyle(marker);
}

  
// --- Get publlished data via script --------------------------------
TH1* GetPublished(UShort_t sNN, Bool_t isNSD)
{
  Long_t ptr = 0;
  Int_t  typ = (isNSD ? 4 : 1);

  ptr = gROOT->ProcessLine(Form("GetSingle(2,1,%d,%d)", sNN, typ));
  if (!ptr) return 0;

  TGraphAsymmErrors* g = reinterpret_cast<TGraphAsymmErrors*>(ptr);
  TH1* h = Graph2Hist(g);
  h->SetLineColor(h->GetMarkerColor());
  
  return h;
}

// --- Add points to stack -------------------------------------------
Double_t
AddToStack(THStack* stack, UShort_t sNN, Bool_t isNSD, 
	   TList* forward, TList* central, 
	   Double_t fwdSysErr, Double_t cenSysErr,
	   Double_t strangeCorr)
{
  TH1* fd1 = GetOne(forward, "Forward", false);
  TH1* fd2 = GetOne(forward, "Forward", true);
  fd1->Scale(1-strangeCorr/100);
  fd2->Scale(1-strangeCorr/100);
  TH1* fs1 = MakeSysError(fd1, fwdSysErr);
  TH1* fs2 = MakeSysError(fd2, fwdSysErr);

  TH1* cd  = 0;
  TH1* cs  = 0;
  if (central) { 
    cd = GetOne(central, "Central", false);
    cs = MakeSysError(cd, cenSysErr);
  }
  else {
    cd = GetPublished(sNN, isNSD);
  }
  SetAttributes(fd1, sNN, 0, false);
  SetAttributes(fd2, sNN, 0, true);
  if (cd) SetAttributes(cd,  sNN, central ? 1 : 2, false);

  if (cs) stack->Add(cs,  "e2");
  stack->Add(fs1, "e2");
  stack->Add(fs2, "e2");
  if (cd) stack->Add(cd,  "ep");
  stack->Add(fd1, "ep");
  stack->Add(fd2, "ep");

  Double_t mcd = (cd ? cd->GetMaximum() : 0);
  Double_t mfs = fs1->GetMaximum();
  
  return TMath::Max(mcd, mfs);
}

// --- Find a list in directory --------------------------------------
TList*
GetList(const TDirectory* d, const char* what)
{
  if (!d) { 
    Error("GetList", "No diretory passed");
    return 0;
  }
  TList* p = static_cast<TList*>(d->Get(Form("%sResults", what)));
  if (!p) { 
    Error("GetList", "%sResults not found in %s", what, d->GetName());
    return 0;
  }
  TList* r = static_cast<TList*>(p->FindObject("all"));
  if (!r) { 
    Error("GetList", "all not found in %s", p->GetName());
    return 0;
  }
  return r;
}

void
AdjustAxis(TAxis* a)
{
  a->SetTitleFont(132);
  a->SetLabelFont(132);
  a->SetTitleSize(0.08);
  a->SetLabelSize(0.08);
  a->SetTitleOffset(0.5);
  a->SetNdivisions(10);
}
  
// --- Make a stack for a single plot --------------------------------
THStack*
MakeStack(Bool_t isNSD,  Bool_t showClusters, 
	  Double_t strangeCorr, Double_t maxFactor=1.1)
{
  const char* trg = isNSD ? "nsd" : "inel";
  TFile* f0900 =TFile::Open(Form("forward_dndeta_%s%04d.root",trg,900, "READ"));
  TFile* f7000 =TFile::Open(Form("forward_dndeta_%s%04d.root",trg,7000,"READ"));
  if (!f0900 || !f7000) { 
    Error("MakeStack", "Failed to open one or more files (%p, %p)", 
	  f0900, f7000);
    return 0;
  }
  TList* forward0900 = GetList(f0900, "Forward");
  TList* forward7000 = GetList(f7000, "Forward");
  TList* central0900 = 0;
  TList* central7000 = showClusters ? GetList(f7000, "Central") : 0;

  THStack* stack = new THStack("stack", "Stack");
  Double_t sysDen = 5;
  Double_t sysPt  = 2;
  Double_t sysMix = 2;
  Double_t sysNch = 5;
  Double_t sysNor = 2;
  Double_t fwdSys = TMath::Sqrt(sysDen * sysDen +
				sysPt  * sysPt  + 
				sysMix * sysMix + 
				sysNch * sysNch + 
				sysNor * sysNor);
  Info("", "Forward systematic error: %4.1f", fwdSys);
  Double_t m0900 = 
    AddToStack(stack, 900,  isNSD, forward0900, central0900, 
	       fwdSys, 5, strangeCorr);
  Double_t m7000 = 
    AddToStack(stack, 7000, isNSD, forward7000, central7000, 
	       fwdSys, 5, strangeCorr);

  stack->SetMaximum(maxFactor * TMath::Max(m0900, m7000));

  f0900->Close();
  f7000->Close();

  return stack;
}

// --- Draw one stack ------------------------------------------------
Double_t 
DrawStack(Bool_t isNSD, Bool_t showClusters, 
	  Double_t strangeCorr, Double_t maxFactor=1.1)
{
  THStack* stack = MakeStack(isNSD, showClusters, strangeCorr, maxFactor);
  stack->Draw("nostack ep");
  stack->GetHistogram()->SetXTitle("#eta");
  if (!isNSD) 
    stack->GetHistogram()->SetYTitle("#frac{1}{N} #frac{dN_{ch}}{d#eta}");
  if (!isNSD) 
    stack->SetMinimum(.0001);
  AdjustAxis(stack->GetHistogram()->GetXaxis());
  AdjustAxis(stack->GetHistogram()->GetYaxis());
  stack->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  gPad->Clear();
  stack->Draw("nostack ep");

  TLatex* title = new TLatex(1-gPad->GetRightMargin()-.01, 
			     1-gPad->GetTopMargin()-.01, 
			     (isNSD ? "NSD" : "INEL"));
  title->SetNDC();
  title->SetTextSize(0.1);
  title->SetTextFont(132);
  title->SetTextAlign(33);
  title->SetTextColor(kBlue+2);
  title->Draw();

  return stack->GetMaximum();
}

// --- Draw one stack ------------------------------------------------
void
ExtractData(Bool_t showClusters = true)
{
  gROOT->LoadMacro("OtherData.C");
  gStyle->SetOptTitle(0);
  gStyle->SetGridColor(kGray);
  TCanvas* c = new TCanvas("c", "C", 1200, 850);
  c->SetRightMargin(0.01);
  c->SetLeftMargin(0.1);
  c->SetTopMargin(0.02);
  c->SetBottomMargin(0.15);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->cd();

  c->Divide(1, 2, 0, 0);

  TVirtualPad* p = c->cd(1);
  p->SetFillColor(0);
  p->SetRightMargin(0.02);
  p->SetGridx();
  p->SetGridy();
  DrawStack(false, showClusters, 1.5);

  Double_t ms = 1.2;

  TLegend* l = new TLegend(.33, .01, .53, .35);
  l->SetBorderSize(0);
  l->SetFillColor(0);
  l->SetTextFont(132);
  l->SetFillStyle(0);
  TLegendEntry* e = l->AddEntry("", "900GeV", "p");
  e->SetMarkerStyle(20);
  e->SetMarkerSize(ms+.1);
  e = l->AddEntry("", "    7TeV", "p");
  e->SetMarkerStyle(21);
  e->SetMarkerSize(ms);
  e = l->AddEntry("", "Mirrored data", "p");
  e->SetMarkerStyle(24);
  e->SetMarkerSize(ms);
  l->Draw();

  TLegend* l2 = new TLegend(.5, .01, .8, .35);
  l2->SetBorderSize(0);
  l2->SetFillColor(0);
  l2->SetTextFont(132);
  l2->SetFillStyle(0);
  l2->SetColumnSeperation(-.05);
  e = l2->AddEntry("", "Forward", "p");
  e->SetMarkerStyle(20);
  e->SetMarkerColor(kRed+2);
  e->SetLineColor(kRed+2);
  e->SetMarkerSize(ms);
  if (showClusters) {
    e = l2->AddEntry("", "Central", "p");
    e->SetMarkerStyle(20);
    e->SetMarkerColor(kMagenta+2);
    e->SetLineColor(kMagenta+2);
    e->SetMarkerSize(ms);
  }
  e = l2->AddEntry("", "#splitline{Eur.Phys.J.#font[22]{C68}:89-108}"
		   "{Eur.Phys.J.#font[22]{C68}:345--354}", "p");
  e->SetMarkerStyle(20);
  e->SetMarkerColor(kBlue+2);
  e->SetLineColor(kBlue+2);
  e->SetMarkerSize(ms);
  l2->Draw();


  p = c->cd(2);
  p->SetGridx();
  p->SetGridy();
  p->SetFillColor(0);
  p->SetRightMargin(0.02);
  DrawStack(true, showClusters, 1.5);

  // TPad, x1, y1, x2, y2, ts, t1, t2, prel
  gROOT->LoadMacro("AddLogo.C");
  AddLogo((TPad*)p, .4, p->GetBottomMargin()+.05, 
	  .5, .44, 0.08, "", "pp data", true);
  
  c->cd();
  TString out("dndeta_pp_forward");
  if (!showClusters) out.Append("_noclusters");
  c->SaveAs(Form("%s.png", out.Data()));
  c->SaveAs(Form("%s.eps", out.Data()));
}



  
