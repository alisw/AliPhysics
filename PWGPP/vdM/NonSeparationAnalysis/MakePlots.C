// -*- C++ -*-

Int_t counterCanvas=0;
TFile *fSave = NULL;

TH1* DrawFrame(TVirtualPad* c,
	       Double_t yMin, Double_t yMax,
	       const char* title,
	       const char* yTitle,
	       Bool_t logY=kFALSE,
	       Double_t xRange=0.65,
	       const char* xTitle="Separation [mm]",
	       Double_t bottomMargin=0.15,
	       Bool_t plotLabel=kTRUE) {
  c->SetLeftMargin(0.17);
  c->SetBottomMargin(bottomMargin);
  TH1* hf = c->DrawFrame(-xRange, yMin, +xRange, yMax);
  hf->SetTitle(title);
  hf->GetXaxis()->SetTitle(xTitle);
  hf->GetYaxis()->SetTitle(yTitle);
  hf->SetLabelSize(0.08, "X");
  hf->SetLabelSize(0.06, "Y");
  hf->SetTitleSize(0.08, "XY");
  hf->SetTitleOffset(0.9, "X");
  hf->SetTitleOffset(1.1, "Y");
  if (logY)
    c->SetLogy();

  if (plotLabel) {
    TLatex *tl = new TLatex;
    tl->SetNDC();
    tl->SetTextSize(0.08);
    tl->SetTextAlign(23);
    tl->DrawLatex(0.55, 0.88, "ALICE pp #sqrt{#it{s}}=5 TeV");
  }
  return hf;
}

void SetupErrorPlot(TVirtualPad *pad, Double_t maxSep) {  
  pad->cd()->SetGrid(0,1);
  TH1* hf = DrawFrame(pad, -8, 8, "", "Residual (#sigma)", kFALSE, maxSep, "Separation [mm]", 0.30, kFALSE);
  hf->SetLabelSize(0.17, "X");
  hf->SetLabelSize(0.15, "Y");
  hf->GetYaxis()->SetNdivisions(-4);
  hf->GetYaxis()->CenterTitle(kTRUE);

  hf->SetTitleSize(0.18, "X");
  hf->SetTitleOffset(0.8, "X");

  hf->SetTitleSize(0.13,  "Y");
  hf->SetTitleOffset(0.35, "Y");
}

TGraph *ConfGraph(TGraph *g, Int_t color, Int_t marker, Int_t lineWidth) {
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetLineColor(color);
  g->SetLineWidth(lineWidth);
  return g;
}

TvectorD              par(34);
TMinuit              *gMinuit         =   NULL;
TParameter<Double_t> *pNDF            =   NULL;
TParameter<Double_t> *pChi2           =   NULL;
TParameter<Double_t> *pScaleRateError =   NULL;
TGraphErrors         *gRate[6]        = { NULL };
TGraph               *gRateModel[6]   = { NULL };
TTree                *tBeamSpot[6]    = { NULL };
TVectorD             *muOffsetsX      =   NULL;
TVectorD             *muOffsetsY      =   NULL;
TGraph               *gMomentModel[6][7];
TGraphErrors         *gMoment[6][7];

TGraphErrors* ScaleGraphError(TGraphErrors *g, Double_t s) {
  Printf("ScaleGraphError %p %f", g, s);
  if (!g)
    return g;
  Double_t *ex = g->GetEX();
  Double_t *ey = g->GetEY();
  for (Int_t i=0, n=g->GetN(); i<n; ++i)
    g->SetPointError(i, ex[i], ey[i]/TMath::Sqrt(s));

  return g;
}


Bool_t ReadData(TString rfn) {
  if (!TFile::Open(rfn))
    return kFALSE;

  gMinuit         = (TMinuit*)gFile->Get("m");
  pNDF            = (TParameter<Double_t>*)gFile->Get("ndf");
  pChi2           = (TParameter<Double_t>*)gFile->Get("chi2");
  pScaleRateError = (TParameter<Double_t>*)gFile->Get("scaleRateError");
  muOffsetsX      = (TVectorD*)gFile->Get("muOffsetsX");
  muOffsetsY      = (TVectorD*)gFile->Get("muOffsetsY");
  for (Int_t i=0; i<6; ++i) {
    gRate[i]     = ScaleGraphError((TGraphErrors*)gFile->Get(Form("gRateScan%c_%d", (i%2) ? 'Y' : 'X', i/2)),
				   pScaleRateError->GetVal());
    tBeamSpot[i] = (TTree*)gFile->Get(Form("fMoments%c_%d",  (i%2) ? 'Y' : 'X', i/2));
  }
   
  return gMinuit && pNDF && pChi2 && pScaleRateError && muOffsetsX && muOffsetsY;
}

struct ModelPar {
  ModelPar(const TMinuit *m) {
    Int_t iuint(0);
    TString _name;
    for (Int_t i=0; i<34; ++i) {
      gMinuit->mnpout(i,_name,par[i],err[i],xlolim[i],xuplim[i],iuint);
      name[i] = _name;
    }
  }

  TVectorD GetPar() {
    TVectorD p(34);
    for (Int_t i=0; i<34; ++i)
      p(i) = par[i];
    return p;
  }

  TString name[34];
  Double_t par[34];
  Double_t err[34];
  Double_t xlolim[34];
  Double_t xuplim[34];
} ;


void MakeGraphsModel() {
  ModelPar mp(gMinuit);

  Int_t       scanType(0);
  TVectorD    beamSep(2);
  TVectorD    mom(10);
  TMatrixDSym cov(10);

  TVectorD    profile(8);

  for (Int_t i=0; i<6; ++i) {
    if (!tBeamSpot[i])
      continue;

    for (Int_t k=0; k<7; ++k) {
      gMomentModel[i][k]  = ConfGraph(new TGraph, kRed,   kFullDotMedium, 2);
      gMoment[i][k]       = (TGraphErrors*)ConfGraph(new TGraphErrors, kBlack, kFullDotMedium, 1);
    }
    gRateModel[i] = ConfGraph(new TGraph, kRed, kFullDotMedium, 2);

    tBeamSpot[i]->SetBranchAddress("scanType", &scanType);
    tBeamSpot[i]->SetBranchAddress("beamSep",  beamSep.GetMatrixArray());
    tBeamSpot[i]->SetBranchAddress("modelPar", mom.GetMatrixArray());
    tBeamSpot[i]->SetBranchAddress("modelCov", cov.GetMatrixArray());

    for (Int_t j=0, n=tBeamSpot[i]->GetEntries(); j<n; ++j) {
      tBeamSpot[i]->GetEntry(j);
      if (!AliNonseparationModelFit::IsGoodLumiRegionFit(i, scanType, beamSep, mom, cov))
	  continue;

      for (Int_t k=0; k<7; ++k) {
	gMoment[i][k]->SetPoint(gMoment[i][k]->GetN(), 10*beamSep(scanType), mom(k));
	gMoment[i][k]->SetPointError(gMoment[i][k]->GetN()-1, 0, TMath::Sqrt(cov(k,k)));
      }
    }

    for (Int_t j=0; j<=140; ++j) {
      beamSep(scanType) = 0.001*j-0.07;
      AliDoubleGaussianBeamProfile::Eval(mp.par[24]*beamSep(0) - (*muOffsetsX)[i/2],
					 mp.par[25]*beamSep(1) - (*muOffsetsY)[i/2],
					 mp.GetPar(), profile, 5e-3);
      
      for (Int_t k=0; k<7; ++k)
	gMomentModel[i][k]->SetPoint(gMomentModel[i][k]->GetN(), 10*beamSep(scanType), profile(1+k) + (k<3 ? mp.par[20+k] : 0.0));
      gRateModel[i]->SetPoint(gRateModel[i]->GetN(), 10*beamSep(scanType), mp.par[23]*profile(0));
    }

    tBeamSpot[i]->ResetBranchAddresses();    
  }  
}

Double_t ComputeR(Double_t &eRmin, Double_t &eRmax, Int_t nEval=20) {
  ModelPar mp(gMinuit);

  TF2 fM("fM", AliDoubleGaussianBeamProfile::EvalProfile0, 
	 -.1, +.1, 
	 -.1, +.1, 26);
  
  TVectorD R(nEval);
  Double_t p[26] = { 0 };
  for (Int_t i=0; i<nEval; ++i) {
    for (Int_t j=0; j<20; ++j) {
      p[j] = TMath::Min(mp.xuplim[j],
			TMath::Max(mp.xlolim[j],
				   mp.par[j]+ 3*(i!=0)*(gRandom->Uniform(0.,1.)>0.5 ? 1 : -1)*mp.err[j]));
    }

    fM.SetParameters(p);
    fM.SetNpx(50);
    fM.SetNpy(50);
    const Double_t eps = 1e-10;
    const Double_t x0 = fM.Eval(0,0);
    const Double_t x1 = fM.Integral(-1,1, -eps,+eps)/2/eps;
    const Double_t x2 = fM.Integral(-eps,+eps, -1,1)/2/eps;
    const Double_t x3 = fM.Integral(-1,1, -1,1);
    R[i]  = (x1/x0) * (x2/x0) / (x3/x0);
    Printf("R[%2d]= %f", i, R[i]);
  }
  eRmin = TMath::MinElement(nEval, R.GetMatrixArray()) - R[0];
  eRmax = TMath::MaxElement(nEval, R.GetMatrixArray()) - R[0];

  return R[0];
}

TGraph* MakeGraphDeviation(TGraphErrors *g, TGraph *gm) {
  TGraph *ge = new TGraph();
  Double_t *x  = g->GetX();
  Double_t *y  = g->GetY();
  Double_t *ey = g->GetEY();
  Double_t ndf    = 0.0;
  Double_t chi2   = 0.0;
  Double_t nSigma = 0.0;
  for (Int_t i=0, n=g->GetN(); i<n; ++i) {
    if (TMath::Abs(x[i]) > 0.7)
      continue;
    nSigma = (y[i] - gm->Eval(x[i]))/ey[i];
    ge->SetPoint(ge->GetN(), x[i], nSigma);
    ndf  += 1.0;
    chi2 += nSigma*nSigma;
  }
  ge->SetLineColor(gm->GetLineColor());
  ge->SetMarkerColor(gm->GetLineColor());
  ge->SetLineWidth(2);
  ge->SetMarkerStyle(kFullDotMedium);

  TLatex *tl = new TLatex;
  tl->SetNDC();
  tl->SetTextColor(kMagenta);
  tl->SetTextSize(0.2);
  tl->DrawLatex(0.2, 0.8, Form("<#chi^{2}> = %.2f", chi2/ndf));

  return ge;
}

void SetupPads(TPad **padData, TPad **padError) {
  if (padError) {
    *padData  = new TPad("padData",  "", 0.0, 0.3, 1.0, 1.0);
    (*padData)->SetBottomMargin(0.15);
    (*padData)->Draw();
    *padError = new TPad("padError", "", 0.0, 0.0, 1.0, 0.3);
    (*padError)->SetTopMargin(0.05);
    (*padError)->Draw();
  } else {
    *padData  = new TPad("padDataNoError",  "", 0.0, 0.0, 1.0, 1.0);
    (*padData)->Draw();
  }

}

void MakePlotsMoments(Int_t scanIndex, TString pn, Bool_t drawErrorPlots=kTRUE, Double_t maxSep=0.65) {  
  const char* labelsX[3] = {
    "Horizontal scan 1",
    "Horizontal scan 2",
    "Horizontal scan with offset"
  };
  const char* labelsY[3] = {
    "Vertical scan 1",
    "Vertical scan 2",
    "Vertical scan with offset"
  };
  const char* yLabels[9] = {
    "#LTx#GT [cm]",
    "#LTy#GT [cm]",
    "#LTz#GT [cm]",
    "#sigma_{x} [cm]",
    "#sigma_{y} [cm]",
    "#sigma_{z} [cm]",
    "#rho_{xy}",
    "T0 rate [Hz]",
  };

  const char* frameTitle = ((scanIndex%2) ? labelsY[scanIndex/2] : labelsX[scanIndex/2]);

  TCanvas *c1 = new TCanvas(frameTitle, "", 600, (drawErrorPlots ? 600 : 600*7/10));
  c1->Divide(3,3);
  
  TPad *padData=NULL, *padError=NULL;
  
  Double_t yRange[7][2] = {
    { gMoment[scanIndex][0]->GetMean(2)-0.01, gMoment[scanIndex][0]->GetMean(2)+0.01 },
    { gMoment[scanIndex][1]->GetMean(2)-0.01, gMoment[scanIndex][1]->GetMean(2)+0.01 },
    { -5.0 + drawErrorPlots*0.0001, +5.0   },
    {  0.0 + drawErrorPlots*0.0001,  0.015 },
    {  0.0 + drawErrorPlots*0.0001,  0.015 },
    {  3.0 + drawErrorPlots*0.0001,  7.0,  },
    { -1.0 + drawErrorPlots*0.0001,  1.0,  }
  };

  for (Int_t i=0; i<2; ++i) {
    yRange[i][0] = TMath::Nint(yRange[i][0]*500)/500.0 + drawErrorPlots*0.00001;
    yRange[i][1] = TMath::Nint(yRange[i][1]*500)/500.0;
  }

  // draw all moments
  for (Int_t i=0; i<7; ++i) {
    c1->cd(1+i);

    SetupPads(&padData, drawErrorPlots ? &padError : NULL);
    DrawFrame(padData->cd(), yRange[i][0], yRange[i][1], frameTitle, yLabels[i], kFALSE, maxSep, drawErrorPlots ? "" : "Separation [mm]", drawErrorPlots ? 0.01 : 0.15);

    gMoment[scanIndex][i]->Draw("PE");
    gMomentModel[scanIndex][i]->Draw("C");
    if (drawErrorPlots) {
      SetupErrorPlot(padError->cd(), maxSep);
      MakeGraphDeviation(gMoment[scanIndex][i], gMomentModel[scanIndex][i])->Draw("PC");
    }
  }

  // draw rate
  c1->cd(8);
  SetupPads(&padData, drawErrorPlots ? &padError : NULL);
  DrawFrame(padData->cd(), .11, 1e6, frameTitle, yLabels[7], kTRUE, maxSep, drawErrorPlots ? "" : "Separation [mm]", drawErrorPlots ? 0.01 : 0.15);
  if (gRate[scanIndex]) {
    delete gRate[scanIndex]->FindObject("gaus");
    gRate[scanIndex]->Draw("PE");
    gRateModel[scanIndex]->Draw("C");
    if (drawErrorPlots) {
      SetupErrorPlot(padError->cd(), maxSep);
      MakeGraphDeviation(gRate[scanIndex], gRateModel[scanIndex])->Draw("PC");
    }
  }

  // draw legend
  c1->cd(9);
  TLegend *leg = new TLegend(0.2, 0.35, 0.8, 0.65);
  leg->AddEntry(gMoment[scanIndex][0],      "Data", "PEL");
  leg->AddEntry(gMomentModel[scanIndex][0], "Model", "PL");
  leg->Draw();

  c1->SaveAs(pn);

  fSave->cd();
  c1->Write(Form("canvas_%d", ++counterCanvas));

  delete c1;
}


void MakePlotsPar(TString pn) {
  TCanvas *c1 = new TCanvas;

  c1->Divide(1,2,0,0);
  c1->cd(1);

  TLatex *tl = new TLatex;
  tl->SetNDC();
  tl->SetTextSize(0.06);

  ModelPar mp(gMinuit);

  const Double_t posX[26] = {
    0.1, 0.1, 0.1, 0.1,
    0.3, 0.3, 0.3, 0.3,
    0.25,
    0.6, 0.6, 0.6, 0.6,
    0.8, 0.8, 0.8, 0.8,
    0.75,
    0.4, 0.4,
    0.7, 0.7, 0.7,
    0.0, // not plotted
    0.4
  };
  const Double_t posY[26] = {
    0.9, 0.8, 0.7, 0.6,
    0.9, 0.8, 0.7, 0.6,
    0.5,
    0.9, 0.8, 0.7, 0.6,
    0.9, 0.8, 0.7, 0.6,
    0.5,
    0.35, 0.25,
    0.35, 0.25, 0.15,
    0.0, // not plotted
    0.15
  };

  const Double_t scale[26] = {
    1e4,1e4,1,1,
    1,  1,  1,1,
    1,
    1e4,1e4,1,1,
    1,  1,  1,1,
    1,
    1e6,1e6,

    1,1,1,
    1,
    1
  };

  const char* units[26] = {
    "#mum","#mum","cm","",
    "",    "",    "",  "",
    "",
    "#mum","#mum","cm","",
    "",    "",    "",  "",
    "",
    "#murad","#murad",

    "cm","cm","cm",
    "",
    ""
  };
  const char* fmt[26] = {
    "%s = %.1f#pm%.1f %s","%s = %.1f#pm%.1f %s","%s = %.1f#pm%.1f %s","%s = %.2f#pm%.2f %s",
    "%s = %.2f#pm%.2f %s","%s = %.2f#pm%.2f %s","%s = %.2f#pm%.2f %s","%s = %.2f#pm%.2f %s",
    "%s = %.2f#pm%.2f %s",
    "%s = %.1f#pm%.1f %s","%s = %.1f#pm%.1f %s","%s = %.1f#pm%.1f %s","%s = %.2f#pm%.2f %s",
    "%s = %.2f#pm%.2f %s","%s = %.2f#pm%.2f %s","%s = %.2f#pm%.2f %s","%s = %.2f#pm%.2f %s",
    "%s = %.2f#pm%.2f %s",
    "%s = %.1f#pm%.1f %s","%s = %.1f#pm%.1f %s",

    "%s = %.5f#pm%.5f %s","%s = %.5f#pm%.5f %s","%s = %.3f#pm%.3f %s",
    "",
    "%s = %.3f#pm%.3f %s"
  };
  
  for (Int_t i=0; i<24; ++i) {
    if (i==23)
      continue;

    tl->DrawLatex(posX[i], posY[i], Form(fmt[i],
					 mp.name[i].Data(),
					 scale[i]*mp.par[i],
					 scale[i]*mp.err[i], units[i]));
  }


  Double_t eRmin=0, eRmax=0;
  const Double_t R = ComputeR(eRmin, eRmax, 20);

  tl->DrawLatex(0.1, 0.35, Form("R = %.3f_{%+.3f}^{%+.3f}", R, eRmin,eRmax));
  tl->DrawLatex(0.1, 0.25, Form("#chi^{2}/NDF = %.0f/%.0f=%.1f", pChi2->GetVal(), pNDF->GetVal(), pChi2->GetVal()/pNDF->GetVal()));
  tl->DrawLatex(0.1, 0.15, Form("0TVX rate errors rescaled by %.2f", 1.0/TMath::Sqrt(pScaleRateError->GetVal())));

  c1->SaveAs(pn);
  fSave->cd();
  c1->Write("canvas_par");
  
  delete c1;
}

void MakePlots(Bool_t drawErrorPlots=kTRUE, Int_t fillNumber=4296, Bool_t useRates=kTRUE, Int_t bcid=-1) {

  gROOT->LoadMacro("AliDoubleGaussianBeamProfile.cxx+");
  gROOT->LoadMacro("AliNonseparationModelFit.cxx+");

  TString bcString = (bcid > 0 ? Form("_bcid%d", bcid) : "");
  TString rfn = TString::Format("root/%d/par_%s%s",
				fillNumber,
				useRates ? "with0TVX" : "withoutRates",
				bcString.Data());
  if (!ReadData(rfn+".root"))
    return;
  
  MakeGraphsModel();

  TString pn = rfn; pn.ReplaceAll("root", "pdf");
  pn += ".pdf";

  TString rn = pn;
  rn += "_canvas.root";
  fSave = TFile::Open(rn, "RECREATE");

  MakePlotsPar(pn+"(");
  for (Int_t i=0; i<6; ++i)
    MakePlotsMoments(i, pn + (i==5 ? ")" : ""));

  fSave->Write();
  fSave->Close();
}
