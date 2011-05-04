//____________________________________________________________________
/** 
 * Get an object with specified name from TCollection @a l 
 * 
 * @param l    Collection
 * @param name Name of object to retrieve 
 * 
 * @return Object, or null 
 */
TObject*
GetObject(const TObject* l, const char* name)
{
  if (!l->IsA()->InheritsFrom(TCollection::Class())) { 
    Error("GetObject", "passed parent %s not a TCollection, but %s",
	  l->GetName(), l->IsA()->GetName());
    return 0;
  }
  TObject* o = static_cast<const TCollection*>(l)->FindObject(name);
  if (!o) { 
    Error("GetObject", "No object '%s' found in '%s'", name, l->GetName());
    return 0;
  }
  return o;
}

//____________________________________________________________________
/** 
 * Get histogram from @a directory/which/sub/hname 
 * 
 * @param dir    Directory 
 * @param which  Name of parent list
 * @param sub    Name of sub-list
 * @param hname  Name of histogram 
 * 
 * @return Pointer to histogram or null
 */
TH1D*
GetHist(TDirectory* dir, 
	const char* which, 
	const char* sub, 
	const char* hname)
{
  if (!dir) return 0;
  TList* parent = static_cast<TList*>(dir->Get(which));
  if (!parent) { 
    Error("GetHist", "List '%s' not found in '%s'", which, dir->GetName());
    return 0;
  }
  TList* child = static_cast<TList*>(parent->FindObject(sub));
  if (!child) { 
    Error("GetHist", "List '%s' not found in '%s'", sub, parent->GetName());
    return 0;
  }
  TObject* o = GetObject(child,hname);
  if (!o) return 0;
  return static_cast<TH1D*>(o);
}


//____________________________________________________________________
/** 
 * Get histogram from 
 * <i>dir/which</i>Results<i>/sub/</i>dndeta<i>which</i>[_rebin<i>rebin</i>]
 * 
 * @param dir    Directory 
 * @param which  Name
 * @param rebin  Optional rebinning 
 * @param sub    Sub-list name
 * 
 * @return Histogram pointer or null
 */
TH1D* 
GetHist(TDirectory* dir, 
	const char* which, 
	UShort_t    rebin,
	const char* sub="all")
{
  TString name(Form("dndeta%s", which));
  if (rebin > 1) name.Append(Form("_rebin%02d", rebin));
  return GetHist(dir, Form("%sResults", which), sub, name);
}

//____________________________________________________________________
/** 
 * Merge two histograms into one 
 * 
 * @param cen    Central part
 * @param fwd    Forward part
 * @param xlow   On return, lower eta bound
 * @param xhigh  On return, upper eta bound
 * 
 * @return Newly allocated histogram or null
 */
TH1* 
Merge(const TH1* cen, const TH1* fwd, Double_t& xlow, Double_t& xhigh)
{
  TH1* tmp = static_cast<TH1*>(fwd->Clone("dndetaMerged"));
  // tmp->SetMarkerStyle(28);
  // tmp->SetMarkerColor(kBlack);
  tmp->SetDirectory(0);
  xlow  = 100;
  xhigh = -100;
  for (Int_t i = 1; i <= tmp->GetNbinsX(); i++) {
    Double_t cc = cen->GetBinContent(i);
    Double_t cf = fwd->GetBinContent(i);
    Double_t ec = cen->GetBinError(i);
    Double_t ef = fwd->GetBinError(i);
    Double_t nc = cf;
    Double_t ne = ef;
    if (cc < 0.001 && cf < 0.01) continue;
    xlow  = TMath::Min(tmp->GetXaxis()->GetBinLowEdge(i),xlow);
    xhigh = TMath::Max(tmp->GetXaxis()->GetBinUpEdge(i),xhigh);
    if (cc > 0.001) {
      nc = cc;
      ne = ec;
      if (cf > 0.001) {
	nc  = (cf + cc) / 2;
	ne  = TMath::Sqrt(ec*ec + ef*ef);
      }
    }
    tmp->SetBinContent(i, nc);
    tmp->SetBinError(i, ne);
  }
  return tmp;
}

//____________________________________________________________________
/** 
 * Function to calculate 
 * @f[
 *  g(x;A_1,A_2,\sigma_1,\sigma_2) = 
 *       A_1\left(\frac{1}{2\pi\sigma_1}e^{(x/\sigma_1)^2} - 
 *           A_2\frac{1}{2\pi\sigma_2}e^{(x/\sigma_2)^2}\right)
 * @f]
 * 
 * @param xp Pointer to x array
 * @param pp Pointer to parameter array 
 * 
 * @return @f$g(x;A_1,A_2,\sigma_1,\sigma_2)@f$
 */
Double_t myFunc(Double_t* xp, Double_t* pp)
{
  Double_t x  = xp[0];
  Double_t a1 = pp[0];
  Double_t a2 = pp[1];
  Double_t s1 = pp[2];
  Double_t s2 = pp[3];
  return a1*(TMath::Gaus(x, 0, s1) - a2 * TMath::Gaus(x, 0, s2));
}
//____________________________________________________________________
/** 
 * Calculate 
 * @f[ 
 *    r(x) = \frac{g(x;A_1,A_2,\sigma_1,\sigma_2)}{
 *                 g(x;A_1',A_2',\sigma'_1,\sigma'_2)}
 * @f] 
 * 
 * @param xp Pointer to X array
 * @param pp Pointer to parameter array (8 entries)
 * 
 * @return @f$r(x)@f$ 
 */
Double_t myRatio(Double_t* xp, Double_t* pp) 
{
  Double_t denom = myFunc(xp, &(pp[4]));
  if (TMath::Abs(denom) < 1.e-6) return 0;
  return myFunc(xp, pp) / denom;
}

//____________________________________________________________________
/** 
 * Fit  @f$g(x;A_1,A_2,\sigma_1,\sigma_2)@f$ to histogram data 
 * 
 * @param tmp    Histogram
 * @param xlow   Lower x bound
 * @param xhigh  Upper x bound 
 *
 * @return Fitted function 
 */
TF1* 
FitMerged(TH1* tmp, Double_t xlow, Double_t xhigh)
{
  TF1* tmpf  = new TF1("tmpf", "gaus", xlow, xhigh);
  tmp->Fit(tmpf, "N", "");
  tmp->SetDirectory(0);

  TF1* fit = new TF1("f", myFunc, xlow, xhigh, 4);
  fit->SetParNames("a_{1}", "a_{2}", "#sigma_{1}", "#sigma_{2}");
  fit->SetParameters(tmpf->GetParameter(0), 
		     .2, 
		     tmpf->GetParameter(2), 
		     tmpf->GetParameter(2)/4);
  fit->SetParLimits(3, 0, 100);
  fit->SetParLimits(4, 0, 100);
  tmp->Fit(fit,"0W","");
  return fit;
}

//____________________________________________________________________
/** 
 * Make band of systematic errors 
 * 
 * @param tmp Histogram
 * @param fit Fit 
 */
void
MakeSysError(TH1* tmp, TF1* fit)
{
  for (Int_t i = 1; i <= tmp->GetNbinsX(); i++) {
    Double_t tc = tmp->GetBinContent(i);
    if (tc < 0.01) continue;
    Double_t x = tmp->GetXaxis()->GetBinCenter(i);
    Double_t y = fit->Eval(x);
    tmp->SetBinContent(i, y);
    tmp->SetBinError(i,sysErr*y);
  }
  return tmp;
}

//____________________________________________________________________
/** 
 * Transform a graph into a histogram 
 * 
 * @param g 
 * 
 * @return 
   */
TH1* 
Graph2Hist(const TGraphAsymmErrors* g)
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

//____________________________________________________________________
/** 
 * Calculate ratio of histogram to function 
 * 
 * @param h     Histogram
 * @param f     Function
 * @param title (Optional) title 
 * 
 * @return Ratio in a histogram 
 */
TH1*
Ratio(TH1* h, TF1* f, const char* title)
{
  TH1* ret = static_cast<TH1*>(h->Clone(Form("%s_%s", 
					       h->GetName(), 
					       f->GetName())));
  ret->SetDirectory(0);
  if (title) ret->SetTitle(title);
  else       ret->SetTitle(Form("%s (data) / %s",
				h->GetTitle(),f->GetTitle()));
  ret->Reset();
  for (Int_t i = 1; i <= ret->GetNbinsX(); i++) { 
    Double_t cc = h->GetBinContent(i);
    if (cc < 0.01) {
      ret->SetBinContent(i,0);
      ret->SetBinError(i,0);
      continue;
    }
    Double_t xx = h->GetBinCenter(i);
    Double_t ee = h->GetBinError(i);
    Double_t ff = f->Eval(xx);
    Double_t yy = cc / ff;
    Double_t ey = ee / ff;
    ret->SetBinContent(i, yy);
    ret->SetBinError(i, ey);
  }
  return ret;
}

//____________________________________________________________________
/** 
 * Get the UA5 data 
 * 
 * @param type   Trigger type (1: INEL, 4: NSD)
 * @param p      On return, positive part or null
 * @param n      On return, negative part or null
 * @param xlow   On return, lower X bound
 * @param xhigh  On return, upper X bound
 * 
 * @return Merged histogram or null 
 */
TH1D* 
GetUA5Data(UShort_t type, TH1*& p, TH1*& n,
	   Double_t& xlow, Double_t& xhigh)
{
  gROOT->SetMacroPath(Form(".:$ALICE_ROOT.trunk/PWG2/FORWARD/analysis2/:%s",
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("OtherData.C");

  p                     = 0;
  n                     = 0;
  TGraphAsymmErrors* gp = GetSingle(UA5,    1, 900, type, 0, 0);
  TGraphAsymmErrors* gn = GetSingle(UA5+10, 1, 900, type, 0, 0);
  if (!gp || !gn) return 0;

  p = Graph2Hist(gp);
  n = Graph2Hist(gn);
  
  Int_t    nn    = p->GetNbinsX();
  xlow           = n->GetXaxis()->GetXmin();
  xhigh          = p->GetXaxis()->GetXmax();
  TH1D*    ret   = new TH1D("ua5", "UA5", 2*nn, xlow, xhigh);
  ret->SetDirectory(0);
  ret->SetMarkerColor(p->GetMarkerColor());
  ret->SetMarkerStyle(p->GetMarkerStyle());

  for (Int_t i = 1; i <= nn; i++) { 
    ret->SetBinContent(nn+i, p->GetBinContent(i));
    ret->SetBinContent(   i, n->GetBinContent(i));
    ret->SetBinError(nn+i, p->GetBinError(i));
    ret->SetBinError(   i, n->GetBinError(i));
  }
  return ret;
}


//____________________________________________________________________
/** 
 * 
 * 
 */
void
DrawUA5Ratios(const char* fname="forward_dndeta.root", UShort_t rebin=5)
{
  TFile* forward_dndeta = TFile::Open(fname, "READ");
  if (!forward_dndeta) { 
    Error("DrawUA5Ratios", "%s not found", fname);
    return;
  }

  UShort_t type = 1;
  TH1D* forward = GetHist(forward_dndeta, "Forward", rebin);
  TH1D* central = GetHist(forward_dndeta, "Central", rebin);

  TObject* sys = GetObject(forward_dndeta->Get("ForwardResults"), "sys");
  TObject* sNN = GetObject(forward_dndeta->Get("ForwardResults"), "sNN");
  TObject* trg = GetObject(forward_dndeta->Get("ForwardResults"), "trigger");
  if (!sys || !sNN || !trg) return;
  if (sys->GetUniqueID() != 1) { 
    Error("DrawUA5Ratios", "Comparison only valid for pp, not %s", 
	  sys->GetTitle());
    return;
  }
  if (sNN->GetUniqueID() != 900) { 
    Error("DrawUA5Ratios", "Comparison only valid for 900GeV, not %dGeV", 
	  sNN->GetUniqueID());
    return;
  }
  if (trg->GetUniqueID() != 1 &&
      trg->GetUniqueID() != 4) { 
    Error("DrawUA5Ratios", 
	  "Comparison only valid for INEL or NSD, not %s (%d)", 
	  trg->GetTitle(), trg->GetUniqueID());
    return;
  }
    

  Double_t alilow, alihigh;
  TH1D* ali     = Merge(central, forward, alilow, alihigh);
  TF1*  alifit  = FitMerged(ali, alilow, alihigh);
  ali->SetTitle("Forward+Central");
  alifit->SetLineColor(kRed+1);
  alifit->SetLineStyle(2);
  alifit->SetName("alifit");
  alifit->SetTitle("Fit to Forward+Central");

  Double_t ua5low, ua5high;
  TH1*  ua5p, *ua5n;
  TH1D* ua5    = GetUA5Data(trg->GetUniqueID(), ua5p, ua5n, ua5low, ua5high);
  TF1*  ua5fit = FitMerged(ua5, ua5low, ua5high);
  ua5fit->SetLineColor(kBlue+1);
  ua5fit->SetLineStyle(3);
  ua5fit->SetName("ua5fit");
  ua5fit->SetTitle("Fit to UA5");

  gStyle->SetOptTitle(0);
  TCanvas* c = new TCanvas("c", "C", 900, 900);
  c->SetFillColor(0);
  c->SetFillStyle(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);

  Double_t yd = .3;
  
  TPad* p1 = new TPad("p1", "p1", 0, yd, 1, 1, 0, 0, 0);
  p1->SetBorderSize(0);
  p1->SetBorderMode(0);
  p1->SetFillColor(0);
  p1->SetFillStyle(0);
  p1->SetRightMargin(0.02);
  p1->SetTopMargin(0.02);
  p1->SetBottomMargin(0.00);
  p1->SetGridx();
  p1->Draw();
  p1->cd();
  
  THStack* all = new THStack("all", "All");
  all->Add(ua5p);
  all->Add(ua5n);
  // all->Add(ua5);
  all->Add(forward);
  all->Add(central);
  // all->Add(merged);
  all->Draw("nostack");
  all->SetMinimum(-.07);
  all->GetXaxis()->SetTitleFont(132);
  all->GetYaxis()->SetTitleFont(132);
  all->GetXaxis()->SetLabelFont(132);
  all->GetYaxis()->SetLabelFont(132);
  all->GetYaxis()->SetDecimals();
  p1->Clear();
  all->Draw("nostack");
  // ua5p->Draw("same p");
  // ua5m->Draw("same p");
  alifit->Draw("same");
  ua5fit->Draw("same");
  
  TLegend* l = new TLegend(.2, .1, .8, .5,
			   Form("pp @ #sqrt{s}=900GeV, %s",trg->GetTitle()));
  l->AddEntry(ua5,     Form("U: %s", ua5->GetTitle()),    "lp");
  l->AddEntry(forward, "F: Forward",                      "lp");
  l->AddEntry(central, "C: Central",                      "lp");
  l->AddEntry(alifit,  
	      Form("f: %s: %4.2f#left[e^{(x/%4.2f)^{2}}-"
		   "%4.2fe^{(x/%4.2f)^{2}}#right]", 
		   alifit->GetTitle(),
		   alifit->GetParameter(0), 
		   alifit->GetParameter(2), 
		   alifit->GetParameter(1), 
		   alifit->GetParameter(3)), "l");
  l->AddEntry(ua5fit,  
	      Form("u: %s: %4.2f#left[e^{(x/%4.2f)^{2}}-"
		   "%4.2fe^{(x/%4.2f)^{2}}#right]", 
		   ua5fit->GetTitle(),
		   ua5fit->GetParameter(0), 
		   ua5fit->GetParameter(2), 
		   ua5fit->GetParameter(1), 
		   ua5fit->GetParameter(3)), "l");
  l->SetTextFont(132);
  l->SetBorderSize(0);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->Draw();

  c->cd();
  TPad* p2 = new TPad("p2", "p2", 0, 0, 1, yd, 0, 0, 0);
  p2->SetBorderSize(0);
  p2->SetBorderMode(0);
  p2->SetFillColor(0);
  p2->SetFillStyle(0);
  p2->SetRightMargin(0.02);
  p2->SetTopMargin(0.00);
  p2->SetBottomMargin(0.15);
  p2->SetGridx();
  p2->Draw();
  p2->cd();

  THStack* ratios = new THStack("Ratios", "Ratios");
  TH1* ua5ali = Ratio(ua5, alifit, 0);
  TH1* aliua5 = Ratio(ali, ua5fit, 0);
  ratios->Add(ua5ali);
  ratios->Add(aliua5);
  ratios->Draw("nostack");
  ratios->SetMinimum(0.4);
  ratios->GetYaxis()->SetTitleSize(2*ratios->GetYaxis()->GetTitleSize());
  ratios->GetYaxis()->SetLabelSize(2*ratios->GetYaxis()->GetLabelSize());
  ratios->GetYaxis()->SetNdivisions(508);
  ratios->GetXaxis()->SetTitleSize(2*ratios->GetXaxis()->GetTitleSize());
  ratios->GetXaxis()->SetLabelSize(2*ratios->GetXaxis()->GetLabelSize());
  ratios->GetXaxis()->SetNdivisions(510);
  ratios->GetXaxis()->SetTitle("#eta");
  ratios->GetXaxis()->SetTitleFont(132);
  ratios->GetYaxis()->SetTitleFont(132);
  ratios->GetXaxis()->SetLabelFont(132);
  ratios->GetYaxis()->SetLabelFont(132);
  ratios->GetYaxis()->SetDecimals();
  p2->Clear();

  TGraphErrors* sysErr = new TGraphErrors(2);
  sysErr->SetPoint(0, all->GetHistogram()->GetXaxis()->GetXmin(),1);
  sysErr->SetPoint(1, all->GetHistogram()->GetXaxis()->GetXmax(),1);
  sysErr->SetPointError(0,0,.07);
  sysErr->SetPointError(1,0,.07);
  sysErr->SetTitle("Systematic error on Forward+Central");
  sysErr->SetFillColor(kYellow+1);
  sysErr->SetFillStyle(3001);  
  sysErr->SetHistogram(ratios->GetHistogram());
  sysErr->DrawClone("ael3");
  ratios->DrawClone("nostack same");

  TF1* fitfit = new TF1("fitfit", myRatio, alilow, alihigh, 8);
  fitfit->SetParameters(ua5fit->GetParameter(0), 
			ua5fit->GetParameter(1), 
			ua5fit->GetParameter(2), 
			ua5fit->GetParameter(3), 
			alifit->GetParameter(0),
			alifit->GetParameter(1),
			alifit->GetParameter(2),
			alifit->GetParameter(3));
  fitfit->Draw("same");
  
  TLegend* ll = new TLegend(.3,.15, .7, .6);
  ll->AddEntry(sysErr, sysErr->GetTitle(), "f");
  ll->AddEntry(ua5ali, ua5ali->GetTitle(), "lp");
  ll->AddEntry(aliua5, aliua5->GetTitle(), "lp");
  ll->AddEntry(fitfit, "UA5 (fit) / Forward+Central (fit)", "lp");
  ll->SetTextFont(132);
  ll->SetBorderSize(0);
  ll->SetFillColor(0);
  ll->SetFillStyle(0);
  ll->Draw();


  c->SaveAs(Form("ua5_ratios_%s_r%02d.png", trg->GetTitle(), rebin));
  gROOT->GetInterpreter()->UnloadFile("OtherData.C");

}

//____________________________________________________________________
//
// EOF
//
