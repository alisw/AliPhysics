//____________________________________________________________________
AliFMDAnaCalibEnergyDistribution*
ViewEDCor(const char* filename)
{
  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gSystem->Load("libPWG0base");
  // gSystem->Load("libPWGLFforward");

  AliFMDAnaParameters* param = AliFMDAnaParameters::Instance();
  param->SetEnergy(900.F);
  param->SetCollisionSystem("pp");
  param->SetMagField(+5.F);
  param->Init(kTRUE,AliFMDAnaParameters::kBackgroundCorrection);

  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("ViewEDCor", "Cannot open file %s", filename);
    return 0;
  }

  AliFMDAnaCalibEnergyDistribution* bg = 
    static_cast<AliFMDAnaCalibEnergyDistribution*>(file->Get("energydistributions"));
  
  return bg;
}


//____________________________________________________________________
struct RingHistos 
{
  TH1D* MakeHisto(const char* name, const char* title, 
		  const TAxis* etaAxis)
  {
    TH1D* ret = new TH1D(Form(name, fD, fR), 
			 Form(title, fD, fR), 
			 etaAxis->GetNbins(), 
			 etaAxis->GetXmin(), 
			 etaAxis->GetXmax());
    ret->SetDirectory(0);
    ret->SetLineColor(fC);
    ret->SetMarkerColor(fC);
    ret->SetMarkerStyle(20+2*fD+(fR == 'I' || fR == 'i' ? 0 : 1));
    ret->SetFillColor(fC);
    ret->SetFillStyle(3001);
    ret->SetXTitle("#eta");
    ret->SetYTitle(title);
    ret->SetStats(0);
    return ret;
  }
  
  
  RingHistos(UShort_t d, Char_t r, const TAxis* etaAxis, Int_t c) 
    : fD(d), fR(r), fC(c) 
  {
    fHists.AddAt(MakeHisto("c%d%c",     "Constant FMD%d%c",      etaAxis),0);
    fHists.AddAt(MakeHisto("mpv%d%c",   "MPV FMD%d%c",           etaAxis),1);
    fHists.AddAt(MakeHisto("w%d%c",     "Width FMD%d%c",         etaAxis),2);
    fHists.AddAt(MakeHisto("alpha%d%c", "#alpha FMD%d%c",        etaAxis),3);
    fHists.AddAt(MakeHisto("beta%d%c",  "#beta FMD%d%c",         etaAxis),4);
    fHists.AddAt(MakeHisto("chi%d%c",   "#chi^{2}/#nu FMD%d%c",  etaAxis),5);
    fHists.AddAt(MakeHisto("acc%d%c",   "Acceptance in FMD%d%c", etaAxis),6);
  }
  void Fill(Int_t iEta, Int_t p, Double_t v, Double_t e=0) 
  {
    TH1D* h = Get(p);
    if (!h) return;
    h->SetBinContent(iEta, v);
    h->SetBinError(iEta, e);
  }
  TH1D* Get(Int_t p) const
  {
    return static_cast<TH1D*>(fHists.At(p));
  }
  UShort_t fD;
  Char_t   fR;
  Int_t    fC;
  TObjArray fHists;
};

//____________________________________________________________________
struct All 
{
  All(const TAxis* etaAxis)
  {
    fFMD1i = new RingHistos(1,'i', etaAxis, kRed+1);
    fFMD2i = new RingHistos(2,'i', etaAxis, kGreen+1);
    fFMD2o = new RingHistos(2,'o', etaAxis, kGreen+7);
    fFMD3i = new RingHistos(3,'i', etaAxis, kBlue+1);
    fFMD3o = new RingHistos(3,'o', etaAxis, kBlue+7);
  }
  RingHistos* Get(UShort_t d, Char_t r) const
  {
    switch (d) { 
    case 1:    return fFMD1i;
    case 2:    return (r == 'I' || r == 'i' ? fFMD2i : fFMD2o);
    case 3:    return (r == 'I' || r == 'i' ? fFMD3i : fFMD3o);
    }
    return 0;
  }
  RingHistos* fFMD1i; 
  RingHistos* fFMD2i; 
  RingHistos* fFMD2o; 
  RingHistos* fFMD3i; 
  RingHistos* fFMD3o; 
};

//____________________________________________________________________
PlotEDCor(const char* filename)
{
  AliFMDAnaCalibEnergyDistribution* ed = ViewEDCor(filename);

  AliFMDAnaParameters* param = AliFMDAnaParameters::Instance();

  TAxis* etaAxis = param->GetBackgroundCorrection(1,'I',4)->GetXaxis();
  All* all = new All(etaAxis);

  for (Int_t iEta = 1; iEta < etaAxis->GetNbins(); iEta++) { 
    Float_t eta = etaAxis->GetBinCenter(iEta);
    Info("PlotEDCor", "eta bin %3d, eta=%f", iEta, eta);

    for (UShort_t d=1; d <= 3; d++) { 
      UShort_t nq = d == 1 ? 1 : 2;
      for (UShort_t q = 0; q < nq; q++) { 
	Char_t r = q == 0 ? 'I' : 'O';
	RingHistos* hists = all->Get(d, r);
	if (!hists) continue;
	
	TH1F* hed = ed->GetEnergyDistribution(d, r, eta);
	if (!hed) continue;

	if (hed->GetEntries() > 0) 
	  hists->Fill(iEta, 6, 1);
	
	TF1* fed = hed->GetFunction("FMDfitFunc");
	if (!fed) continue;

	Int_t npar = fed->GetNpar();
	for (Int_t i = 0; i < npar; i++) { 
	  Double_t par = fed->GetParameter(i);
	  Double_t err = fed->GetParError(i);
	  hists->Fill(iEta, i, par, err);
	}
	Double_t chi2 = fed->GetChisquare();
	Double_t ndf  = fed->GetNDF();
	if (ndf != 0) 
	  hists->Fill(iEta, 5, chi2/ndf);
      }
    }
  }
	  
  TCanvas* c = new TCanvas("c", "c", 800, 800);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  c->Divide(1,7,0,0);
  
  THStack* stack = 0;
  
  const char* yTitles[] = { "Constant", "#Delta_{p}", "w", 
			    "#alpha", "#beta", "#chi^2/NDF", 
			    "Acceptance" };
  for (Int_t i = 0; i < 7; i++) { 
    TVirtualPad* p = c->cd(i+1);
    p->SetFillColor(0);
    p->SetFillStyle(0);
    THStack* stack = new THStack;
    stack->Add(all->Get(1,'i')->Get(i));
    stack->Add(all->Get(2,'i')->Get(i));
    stack->Add(all->Get(2,'o')->Get(i));
    stack->Add(all->Get(3,'i')->Get(i));
    stack->Add(all->Get(3,'o')->Get(i));
    stack->Draw("nostack");
    stack->GetHistogram()->SetYTitle(yTitles[i]);
    stack->GetHistogram()->SetXTitle("#eta");
    TAxis* yaxis = stack->GetHistogram()->GetYaxis();
    yaxis->SetTitleSize(0.15);
    yaxis->SetLabelSize(0.08);
    yaxis->SetTitleOffset(0.35);
    yaxis->SetNdivisions(10);
    TAxis* xaxis = stack->GetHistogram()->GetXaxis();
    xaxis->SetTitleSize(0.15);
    xaxis->SetLabelSize(0.08);
    xaxis->SetTitleOffset(0.35);
    xaxis->SetNdivisions(320);
    p->cd();
  }
  c->cd();
  c->Print("energyDist.png");
}
//____________________________________________________________________
//
// EOF
//

