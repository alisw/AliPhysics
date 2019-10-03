// -*- C++ -*-

const char* momentFileNames[] = {
  "root/2852/lumiRegion_Scan1X%s.root",
  "root/2852/lumiRegion_Scan1Y%s.root",

  "root/2852/lumiRegion_Scan2X%s.root",
  "root/2852/lumiRegion_Scan2Y%s.root"
};

TVectorD MakeStartParameters() {
  TVectorD par(26);

  par[ 0] = 0.010; par[ 1] = 0.010; par[ 2] = 8.0; par[ 3] =  0.0;
  par[ 4] = 1/0.6; par[ 5] = 1/0.8; par[ 6] = 1.1; par[ 7] = -0.0;

  par[ 8] = 0.85;

  par[ 9] = 0.010; par[10] = 0.010; par[11] = 8.0; par[12] = -0.0;
  par[13] = 1/0.6; par[14] = 1/0.6; par[15] = 1.1; par[16] =  0.0;

  par[17] = 0.85;

  par[18] = 0.0; par[19] = 212e-6;

  par[20] = 0.0872; par[21] = 0.3126; par[22] = 0.0;

  par[23] = 1e-7;

  par[24] = 1; par[25] = 1;

  return par;
}

TGraphErrors* ScaleGraph(TGraphErrors *g, Double_t s) {
  Printf("ScaleGraph %f", s);
  if (!g)
    return g;
  Double_t *x = g->GetX();
  Double_t *y = g->GetY();
  Double_t *ex = g->GetEX();
  Double_t *ey = g->GetEY();
  for (Int_t i=0, n=g->GetN(); i<n; ++i) {
    g->SetPoint(i, x[i], s*y[i]);
    g->SetPointError(i, ex[i], s*ey[i]);
  }
  return g;
}

const Double_t wMin[3] = { 0.10, 0.05, 0.50 };
const Double_t wMax[3] = { 0.90, 0.95, 0.90 };

void MakeNonseparationFit_2852(Int_t index_BC=-1, Int_t wRange=0, Bool_t fixCrossingAngles=kTRUE) // -1 -> all BCs, else single BC
{
  gROOT->LoadMacro("AliDoubleGaussianBeamProfile.cxx+");
  gROOT->LoadMacro("AliNonseparationModelFit.cxx+");
  gROOT->LoadMacro("MakePlots.C");
  gROOT->LoadMacro("ExtractFromCanvas.C");

  TTree tBC;
  tBC.ReadFile("txt/2852/BCID.txt");

  Int_t bcid = -1;
  tBC.SetBranchAddress("bcid", &bcid);
  if (index_BC != -1) {
    tBC.GetEntry(index_BC);
  }

  TString txtFileName;
  if (index_BC != -1)
    txtFileName = TString::Format("txt/2852/fit_2852_%d_%d_%d.txt", bcid, wRange, fixCrossingAngles);
  else
    txtFileName = TString::Format("txt/2852/fit_2852_%d_%d.txt",          wRange, fixCrossingAngles);

  std::ofstream ofs(txtFileName.Data());

  TTree tBI;
  tBI.ReadFile("txt/2852/bunchIntensities.txt");
  tBI.Draw("1/(s1b1*s1b2)*(s2b1*s2b2)", "", "GOFF");
  Double_t *relInt = tBI.GetV1();

  TFile::Open("root/2852/T0Rate_vs_Separation.root");
  // TFile::Open("root/2852/Rates_2852.root");

  const char *rnames[4] = {
    "T0_SX1",
    "T0_SY1",
    "T0_SX2",
    "T0_SY2",
  };
  Int_t perm[4] = { 0,1,2,3 };
  TGraphErrors *gRates[4] = { NULL };
  for (Int_t i=0; i<4; ++i) {
    TString nn = TString::Format("gre_T0s%d_BID%d", perm[i], index_BC < 0 ? 0 : index_BC);
    // TString nn = TString::Format("%s;%d", rnames[i], 1+(index_BC < 0 ? 0 : index_BC));
    Printf("%s", nn.Data());
    gRates[i] = (TGraphErrors*)gFile->Get(nn);
    if (i<2)
      ScaleGraph(gRates[i], relInt[0]);
  }
  if (index_BC == -1) {
    for (Int_t i=0; i<4; ++i) {
      for (Int_t j=1; j<tBC.GetEntries(); ++j) {
	tBC.GetEntry(j);
	TString nn = TString::Format("gre_T0s%d_BID%d", perm[i], index_BC < 0 ? 0 : index_BC);
	// TString nn = TString::Format("%s;%d", rnames[i], 1+(index_BC < 0 ? 0 : index_BC));
	TGraphErrors *gr = (TGraphErrors*)gFile->Get(nn);
	if (i<2)
	  ScaleGraph(gr, relInt[j]);
	Double_t *xx0 = gRates[i]->GetX();
	Double_t *ex0 = gRates[i]->GetEX();
	Double_t *yy0 = gRates[i]->GetY();
	Double_t *ey0 = gRates[i]->GetEY();
	Double_t *yy = gr->GetY();
	Double_t *ey = gr->GetEY();
	for (Int_t k=0; k<gRates[i].GetN(); ++k) {
	  gRates[i]->SetPoint(k, xx0[k], yy0[k]+yy[k]);
	  gRates[i]->SetPointError(k, ex0[k], TMath::Sqrt(TMath::Power(ey0[k],2)+
							  TMath::Power(ey[k], 2)));
	}
      }
    }
    bcid=-1;
  }
  Printf("%p %p %p %p",
	 gRates[0],
	 gRates[1],
	 gRates[2],
	 gRates[3]);

  AliNonseparationModelFit f;
  f.SetMom9Limits(1.0, 1.8);
  f.SetWRange(wMin[wRange], wMax[wRange]);
  f.SetFixCrossingAngles(fixCrossingAngles);

  TString s = "";
  if (bcid != -1)
    s = TString::Format("_bcid%d", bcid);

  TF1 *fg = new TF1("gaus", "[0]*TMath::Gaus(x,[1],[2])*(1+[3]*(x-[1])**2+[4]*(x-[1])**4)");
  TF1 *fgg;
  for (Int_t i=0; i<4; ++i) {
    TFile::Open(TString::Format(momentFileNames[i], s.Data()));
    TTree *t = (TTree*)gFile->Get("TBeamSpot");
    fg->SetParameters(8e4,0,0.12,0,0);
    gRates[i]->Fit(fg);
    gRates[i]->Fit(fg);
    gRates[i]->SetLineColor(kBlack);
    gRates[i]->SetLineWidth(1);
    gRates[i]->SetMarkerColor(kBlack);
    gRates[i]->SetMarkerStyle(kFullDotMedium);
    f.Add(i, t, gRates[i]);
  }

  // TVectorD par = MakeStartParameters();
  // f.DoFit(par, TString::Format("root/%d/par_with0TVX%s_%d_%d.root", 2852, s.Data(), wRange, fixCrossingAngles));

  MakePlots(!kTRUE, 2852, kTRUE, bcid, wRange, fixCrossingAngles);
  TString fn = TString::Format("pdf/%d/par_with0TVX%s_%d_%d.pdf_canvas.root", 2852, s.Data(), wRange, fixCrossingAngles);
  TString line = ExtractFromCanvas(fn);
  ofs << line.Data() << std::endl;
}
