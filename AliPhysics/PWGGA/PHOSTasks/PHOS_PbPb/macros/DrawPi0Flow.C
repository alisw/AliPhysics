DrawPi0Flow(const TString filename = "AnalysisResults.root",
	    const TString eventType = "SemiCentral")
{
  /* $Id$ */

  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  TFile * f = new TFile(filename) ;
  TString listName;
  if (eventType.Contains("Central"))
    listName = Form("AliPHOSPi0Flow/PHOSPi0FlowCoutput1",
		    eventType.Data(),eventType.Data());
  else
    listName = Form("AliPHOSPi0Flow_%s/AliPHOSPi0Flow_%sCoutput1",
		    eventType.Data(),eventType.Data());
  cout << listName << endl;
  TList *histoList = (TList*)f->Get(listName.Data()); // lego train

  TH1F *hev         = (TH1F*)histoList->FindObject("hTotSelEvents") ;
  TH2F* hZvertex    = (TH2F*)histoList->FindObject("hZvertex");
  TH2F* hCenTrack   = (TH2F*)histoList->FindObject("hCenTrack");
  TH2F *hCenPHOS    = (TH2F*)histoList->FindObject("hCenPHOS") ;
  TH2F* phiRP       = (TH2F*)histoList->FindObject("phiRP");
  TH2F* phiRPflat   = (TH2F*)histoList->FindObject("phiRPflat");
  TH2F* phiRPV0A    = (TH2F*)histoList->FindObject("phiRPV0A");
  TH2F* phiRPV0C    = (TH2F*)histoList->FindObject("phiRPV0C");
  TH2F* phiRPV0Aflat= (TH2F*)histoList->FindObject("phiRPV0Aflat");
  TH2F* phiRPV0Cflat= (TH2F*)histoList->FindObject("phiRPV0Cflat");
  TH2F* hCellNXZM1  = (TH2F*)histoList->FindObject("hCellNXZM1");
  TH2F* hCellNXZM2  = (TH2F*)histoList->FindObject("hCellNXZM2");
  TH2F* hCellNXZM3  = (TH2F*)histoList->FindObject("hCellNXZM3");
  TH2F *hPi0All_cen0     = (TH2F*)histoList->FindObject("hPi0All_cen0");
  TH2F *hPi0Allcore_cen0 = (TH2F*)histoList->FindObject("hPi0Allcore_cen0");
  TH2F *hPi0M11          = (TH2F*)histoList->FindObject("hPi0M11");
  TH2F *hPi0M22          = (TH2F*)histoList->FindObject("hPi0M22");
  TH2F *hPi0M33          = (TH2F*)histoList->FindObject("hPi0M33");

  //-----------------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1","Event selection");
  hev->GetXaxis()->SetRangeUser(0,6);
  hev->GetXaxis()->SetBinLabel(1,"Total");
  hev->GetXaxis()->SetBinLabel(2,"Has Vertex");
  hev->GetXaxis()->SetBinLabel(3,"abs(z_vertex) < 10.");
  hev->GetXaxis()->SetBinLabel(4,"Has Centrality");
  hev->GetXaxis()->SetBinLabel(5,"C. upper edge");
  hev->GetXaxis()->SetBinLabel(6,"C. lower edge");
  hev->GetXaxis()->SetBinLabel(7,"Has PClusters");
  hev->SetYTitle("N_{events}");
  hev->GetYaxis()->SetTitleOffset(1.2);
  hev->Draw();
  c1->Print("PHOS_EvSel.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c2 = new TCanvas("c2","Track multiplicity vs centrality");
  gPad->SetLogz();
  hCenTrack->SetXTitle("centrality (%)");
  hCenTrack->SetYTitle("Number of tracks");
  hCenTrack->Draw("colz");
  c2->Print("TrackMultCentrality.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c3 = new TCanvas("c3","PHOS multiplicity vs centrality");
  gPad->SetLogz();
  hCenPHOS->SetXTitle("centrality (%)");
  hCenPHOS->SetYTitle("Number of PHOS clusters");
  hCenPHOS->Draw("colz");
  c3->Print("PHOSMultCentrality.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c3cent = new TCanvas("c3cent","Centrality");
  hCenPHOS->SetXTitle("centrality (%)");
  hCenPHOS->SetYTitle("Number of PHOS clusters");
  TH1D *hCent = hCenPHOS->ProjectionX();
  hCent->SetTitle("Centrality");
  hCent->Draw();
  c3cent->Print("Centrality.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c4V0 = new TCanvas("c4V0","RP phi in VZERO");
  phiRPV0A1 = phiRPV0A->ProjectionX();
  phiRPV0C1 = phiRPV0C->ProjectionX();
  phiRPV0A1->SetLineColor(kRed);
  phiRPV0C1->SetLineColor(kBlue);
  phiRPV0A1->SetXTitle("VZERO RP #phi");
  phiRPV0C1->SetXTitle("VZERO RP #phi");
  phiRPV0A1->SetMinimum(phiRPV0A1->GetMinimum()*0.98);
  phiRPV0A1->SetMaximum(phiRPV0A1->GetMaximum()*1.02);
  phiRPV0A1->Draw();
  phiRPV0C1->Draw("same");
  leg = new TLegend(0.2,0.8,0.39,0.89);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->AddEntry(phiRPV0A1,"V0A","l");
  leg->AddEntry(phiRPV0C1,"V0C","l");
  leg->Draw();
  c4V0->Print("V0RPphi.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c4V0flat = new TCanvas("c4V0flat","RP phi in VZERO flattened");
  phiRPV0Aflat1 = phiRPV0Aflat->ProjectionX();
  phiRPV0Cflat1 = phiRPV0Cflat->ProjectionX();
  phiRPV0Aflat1->SetLineColor(kRed);
  phiRPV0Cflat1->SetLineColor(kBlue);
  phiRPV0Aflat1->SetXTitle("VZERO RP #phi");
  phiRPV0Cflat1->SetXTitle("VZERO RP #phi");
  phiRPV0Aflat1->Draw();
  phiRPV0Cflat1->Draw("same");
  leg = new TLegend(0.2,0.8,0.39,0.89);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->AddEntry(phiRPV0Aflat1,"V0A","l");
  leg->AddEntry(phiRPV0Cflat1,"V0C","l");
  leg->Draw();
  c4V0flat->Print("V0RPphiFlat.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c4TPC = new TCanvas("c4TPC","TPC RP phi");
  phiRP1 = phiRP->ProjectionX();
  phiRP1->SetXTitle("TPC RP #phi");
  phiRP1->Draw();
  c4TPC->Print("TPCRPphi.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c4TPCflat = new TCanvas("c4TPCflat","TPC RP phi flattened");
  phiRPflat1 = phiRPflat->ProjectionX();
  phiRPflat1->SetXTitle("TPC RP #phi");
  phiRPflat1->Draw();
  c4TPCflat->Print("TPCRPphiFlat.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c5 = new TCanvas("c5","XZ in M1");
  c5->SetLogz();
  hCellNXZM1->SetTitle("Cell occupancy in module 4");
  hCellNXZM1->SetXTitle("X (cells)");
  hCellNXZM1->SetYTitle("X (cells)");
  hCellNXZM1->Draw("colz");
  c5->Print("PHOS_XYM4.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c6 = new TCanvas("c6","XZ in M2");
  c6->SetLogz();
  hCellNXZM2->SetTitle("Cell occupancy in module 3");
  hCellNXZM2->SetXTitle("X (cells)");
  hCellNXZM2->SetYTitle("X (cells)");
  hCellNXZM2->Draw("colz");
  c6->Print("PHOS_XYM3.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c7 = new TCanvas("c7","XZ in M3");
  c7->SetLogz();
  hCellNXZM3->SetTitle("Cell occupancy in module 2");
  hCellNXZM3->SetXTitle("X (cells)");
  hCellNXZM3->SetYTitle("X (cells)");
  hCellNXZM3->Draw("colz");
  c7->Print("PHOS_XYM2.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c8 = new TCanvas("c8","gg mass vs pt, PID=All");
  c8->SetLogz();
  hPi0All_cen0->SetTitle("PID: All");
  hPi0All_cen0->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
  hPi0All_cen0->SetYTitle("p_{T} (GeV/c)");
  hPi0All_cen0->Draw("colz");
  c8->Print("PHOS_MggAll.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c9 = new TCanvas("c9","gg mass vs pt, PID=Allcore");
  c9->SetLogz();
  hPi0Allcore_cen0->SetTitle("PID: Allcore");
  hPi0Allcore_cen0->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})"); 
  hPi0Allcore_cen0->SetYTitle("p_{T} (GeV/c)");
  hPi0Allcore_cen0->Draw("colz");
  c9->Print("PHOS_MggAllcore.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c10 = new TCanvas("c10","gg mass in M1");
  hM1 = hPi0M11->ProjectionX("m1",31,400);
  hM1->SetAxisRange(0.,0.29.);
  hM1->SetTitle("#gamma#gamma in module 1");
  Fit1Pi0(hM1,2);
  hM1->Draw();
  TPaveText *txtM1 = new TPaveText(0.6,0.7,0.89,0.89,"NDC");
  txtM1->SetFillColor(kWhite);
  txtM1->SetBorderSize(0);
  txtM1->AddText("p_{T}>3 GeV/c");
  txtM1->AddText(Form("M_{#pi^{0}} = %.2f #pm %.2f MeV/c^{2}",
		      hM1->GetFunction("fitfun")->GetParameter(1)*1000,
		      hM1->GetFunction("fitfun")->GetParError (1)*1000));
  txtM1->AddText(Form("#sigma_{#pi^{0}} = %.2f #pm %.2f MeV/c^{2}",
		      hM1->GetFunction("fitfun")->GetParameter(2)*1000,
		      hM1->GetFunction("fitfun")->GetParError (2)*1000));
  txtM1->Draw();
  c10->Print("PHOS_MggM1.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c11 = new TCanvas("c11","gg mass in M2");
  hM2 = hPi0M22->ProjectionX("m2",31,400);
  hM2->SetAxisRange(0.,0.29.);
  hM2->SetTitle("#gamma#gamma in module 2");
  Fit1Pi0(hM2,2);
  hM2->Draw();
  TPaveText *txtM2 = new TPaveText(0.6,0.7,0.89,0.89,"NDC");
  txtM2->SetFillColor(kWhite);
  txtM2->SetBorderSize(0);
  txtM2->AddText("p_{T}>3 GeV/c");
  txtM2->AddText(Form("M_{#pi^{0}} = %.2f #pm %.2f MeV/c^{2}",
		      hM2->GetFunction("fitfun")->GetParameter(1)*1000,
		      hM2->GetFunction("fitfun")->GetParError (1)*1000));
  txtM2->AddText(Form("#sigma_{#pi^{0}} = %.2f #pm %.2f MeV/c^{2}",
		      hM2->GetFunction("fitfun")->GetParameter(2)*1000,
		      hM2->GetFunction("fitfun")->GetParError (2)*1000));
  txtM2->Draw();
  c11->Print("PHOS_MggM2.eps");

  //-----------------------------------------------------------------------------
  TCanvas *c12 = new TCanvas("c12","gg mass in M3");
  hM3 = hPi0M33->ProjectionX("m3",31,400);
  hM3->SetAxisRange(0.,0.29.);
  hM3->SetTitle("#gamma#gamma in module 3");
  Fit1Pi0(hM3,2);
  hM3->Draw();
  TPaveText *txtM3 = new TPaveText(0.6,0.7,0.89,0.89,"NDC");
  txtM3->SetFillColor(kWhite);
  txtM3->SetBorderSize(0);
  txtM3->AddText("p_{T}>3 GeV/c");
  txtM3->AddText(Form("M_{#pi^{0}} = %.2f #pm %.2f MeV/c^{2}",
		      hM3->GetFunction("fitfun")->GetParameter(1)*1000,
		      hM3->GetFunction("fitfun")->GetParError (1)*1000));
  txtM3->AddText(Form("#sigma_{#pi^{0}} = %.2f #pm %.2f MeV/c^{2}",
		      hM3->GetFunction("fitfun")->GetParameter(2)*1000,
		      hM3->GetFunction("fitfun")->GetParError (2)*1000));
  txtM3->Draw();
  c12->Print("PHOS_MggM3.eps");
}

//=============================================================================
Fit1Pi0(const TH1D *hMass = 0,
       const Int_t polN = 1)
{
  // This script takes a 2D histogram hMassPt with invariant mass vs
  // pt of cluster pairs:
  // - slices them along pt bins,
  // - fits 1D invariant mass specrta by Gaussian+polynomial
  // - calculates the number of pi0's as an integral of Gaussian
  // - calculates background under pi0 as an integral of polynomial

  TF1 *fitfun = 0;
  if      (polN == 1)
    fitfun = new TF1("fitfun",pi0massP1,0.1,0.7,5);
  else if (polN == 2)
    fitfun = new TF1("fitfun",pi0massP2,0.1,0.7,6);
  else if (polN == 101)
    fitfun = new TF1("fitfun",CombiBG,0.0,1.0,6);
  fitfun->SetLineColor(kRed);
  fitfun->SetLineWidth(2);

  if (polN <= 2) {
    fitfun->SetParName(0,"A");
    fitfun->SetParName(1,"m_{0}");
    fitfun->SetParName(2,"#sigma");
    fitfun->SetParName(3,"a_{0}");
    fitfun->SetParName(4,"a_{1}");
    if (polN == 2)
      fitfun->SetParName(5,"a_{2}");
    
    fitfun->SetParLimits(0,  1.000,1.e+5);
    fitfun->SetParLimits(1,  0.12,0.14);
    fitfun->SetParLimits(2,  0.003,0.020);
  }
  else if (polN == 101) {
    fitfun->SetParName(0,"a_{0}");
    fitfun->SetParName(1,"a_{1}");
    fitfun->SetParName(2,"a_{2}");
    fitfun->SetParName(3,"a_{3}");
    fitfun->SetParName(4,"a_{4}");
    fitfun->SetParName(5,"a_{5}");
  }

  Int_t   nM     = hMass->GetNbinsX();
  Float_t mMin   = hMass->GetXaxis()->GetXmin();
  Float_t mMax   = hMass->GetXaxis()->GetXmax();
  Float_t mStep  = (mMax-mMin)/nM;

  hMass->SetXTitle("M_{#gamma#gamma}, GeV/c");
  hMass->SetAxisRange(0.01,1.0);
  Float_t nPi0Max = hMass->Integral(hMass->FindBin(0.30),
				    hMass->FindBin(0.40));
  for (Int_t iM=1; iM<nM; iM++)
    if (hMass->GetBinContent(iM) == 0) hMass->SetBinError(iM,1.0);
  hMass->SetMinimum(0.01);
  if      (polN == 1)
    fitfun->SetParameters(nPi0Max/4,0.135,0.010,1.,0.);
  else if (polN == 2)
    fitfun->SetParameters(nPi0Max/4,0.135,0.010,1.,0.,0.);
//   else if (polN == 101)
//     fitfun->SetParameters(1.,1.,100.,1.,1.,10.);
  Double_t mFitMin, mFitMax;
  mFitMin=0.10;
  mFitMax=0.18;
  hMass->Fit("fitfun","Q","",mFitMin,mFitMax);
  hMass->SetAxisRange(mFitMin,mFitMax);
  fitfun->SetNpx(1./mStep);
//   TH1F *bgFun = fitfun->GetHistogram();
//   hMass->Add(bgFun,-1.);
}

//-----------------------------------------------------------------------------
Double_t pi0massP2(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                       (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0];
  return gaus+back;
}

//-----------------------------------------------------------------------------
Double_t pi0massP1(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                       (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0];
  return gaus+back;
}

//-----------------------------------------------------------------------------
Double_t CombiBG(Double_t *x, Double_t *par)
{
  if (x[0] > 0.12 && x[0] < 0.15) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t back = par[0] + par[1]*x[0] 
    + par[2]*TMath::Power(x[0],2)
    + par[3]*TMath::Power(x[0],3) 
    + par[4]*TMath::Power(x[0],4)
    + par[5]*TMath::Power(x[0],5);
  return back;
}

//-----------------------------------------------------------------------------
Double_t CombiBGExp(Double_t *x, Double_t *par)
{
  if (x[0] > 0.12 && x[0] < 0.15) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t back = par[0] * TMath::Power(x[0],1.5) * TMath::Exp(-par[1]*x[0]);
  return back;
}

