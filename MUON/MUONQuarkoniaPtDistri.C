//
// Macro for analysis of 2d InvMass v Pt distribution for the generation of:
// - Pt distribution of quarkonia, 
// - Mass and sigma of the quarkonia resonance as a function of pT
// - and invariant mass distributions ( all, background and subtracted)
//   as a function of pT
// Usual utilisation:
// QuarkoniaPtDistri("ProdMuonCocktailPer1/MUONmassPlot.root",4,2.9,3.3)
//
//  MUONmassPlotFile  root file form MUONmassPlot_ESD.C
//  NRes=1 bining of the pT distribution 
//    ==1, equal to MUONmassPlotFile 2D histo, 
//    ==2, bin twice larger
//  InvMass1, InvMass2 Range of the quarkonia resonance
//
// Gines MARTINEZ, Subatech, May 04
//

void MUONQuarkoniaPtDistri(char * MUONmassPlotFile="MUONmassPlot.root", Int_t NRes=4, Float_t InvMass1=2.9, Float_t InvMass2=3.3)
{
  TFile MUONmassPlot(MUONmassPlotFile);

  TH2F * hInvMassAll_vs_Pt = (TH2F*) MUONmassPlot.Get("hInvMassAll_vs_Pt");
  TH2F * hInvMassBgk_vs_Pt = (TH2F*) MUONmassPlot.Get("hInvMassBgk_vs_Pt");
  hInvMassAll_vs_Pt->Sumw2();
  hInvMassBgk_vs_Pt->Sumw2();

  gROOT->cd();

  Int_t   Nbins = hInvMassAll_vs_Pt->GetYaxis()->GetNbins();
  Float_t ptMin = hInvMassAll_vs_Pt->GetYaxis()->GetXmin();
  Float_t ptMax = hInvMassAll_vs_Pt->GetYaxis()->GetXmax();

  const Int_t NumberOfPoints = Nbins/NRes;
  Float_t *pT = new Float_t[NumberOfPoints];
  Float_t *e_pT = new Float_t[NumberOfPoints];
  Float_t *dNdpT = new Float_t[NumberOfPoints];
  Float_t *e_dNdpT = new Float_t[NumberOfPoints];

  Float_t *pT1 = new Float_t[NumberOfPoints];
  Float_t *e_pT1 = new Float_t[NumberOfPoints];
  Float_t *Cent = new Float_t[NumberOfPoints];
  Float_t *e_Cent = new Float_t[NumberOfPoints];

  Float_t *pT2 = new Float_t[NumberOfPoints];
  Float_t *e_pT2 = new Float_t[NumberOfPoints];
  Float_t *Sigma = new Float_t[NumberOfPoints];
  Float_t *e_Sigma = new Float_t[NumberOfPoints];

  TH1D * hInvMassAll[NumberOfPoints+1];
  TH1D * hInvMassBgk[NumberOfPoints+1];
  TH1D * hInvMassSub[NumberOfPoints+1];

  char histoname[60];

  Int_t ibin, ibin1,ibin2,  ipoint;
  Int_t integral1, integral2;
  for(ipoint=0; ipoint<=NumberOfPoints; ipoint++) {

    if (ipoint!=NumberOfPoints) {
      ibin1 = ipoint * NRes;
      ibin2 = ibin1+  NRes;
      pT[ipoint]=0.;
      e_pT[ipoint]=0.;
      dNdpT[ipoint]=0.;
      e_dNdpT[ipoint]=0.;
      for(ibin=ibin1; ibin<ibin2; ibin++) {
	pT[ipoint]+=hInvMassAll_vs_Pt->GetYaxis()->GetBinCenter(ibin);
      }
      pT[ipoint]/=NRes;
      pT1[ipoint]=pT[ipoint];
      pT2[ipoint]=pT[ipoint];
      
      e_pT[ipoint] = (hInvMassAll_vs_Pt->GetYaxis()->GetBinCenter(ibin2)-hInvMassAll_vs_Pt->GetYaxis()->GetBinCenter(ibin1))/2.;
      e_pT1[ipoint]=e_pT[ipoint];
      e_pT2[ipoint]=e_pT[ipoint];
    } 
    else {
      ibin1 = 0;
      ibin2 = Nbins;
    }


    if (ipoint!=NumberOfPoints) sprintf(histoname,"hInvMassAll_%2.0f",100*pT[ipoint]);
    else sprintf(histoname,"hInvMassAll_full");
    TH1D * nume = hInvMassAll_vs_Pt->ProjectionX("nume",ibin1,ibin2-1);
    hInvMassAll[ipoint] =  new TH1D(*nume);
    hInvMassAll[ipoint]->SetName(histoname); hInvMassAll[ipoint]->SetTitle(histoname);
    hInvMassAll[ipoint]->Sumw2();

    if (ipoint!=NumberOfPoints)  sprintf(histoname,"hInvMassSub_%2.0f",100*pT[ipoint]);
    else sprintf(histoname,"hInvMassSub_full");
    hInvMassSub[ipoint] =  new TH1D(*nume);
    hInvMassSub[ipoint]->SetName(histoname); hInvMassSub[ipoint]->SetTitle(histoname);
    hInvMassSub[ipoint]->Sumw2();

    if (ipoint!=NumberOfPoints)  sprintf(histoname,"hInvMassBgk_%2.0f",100*pT[ipoint]);
    else sprintf(histoname,"hInvMassBgk_full");
    TH1D * deno = hInvMassBgk_vs_Pt->ProjectionX("deno",ibin1,ibin2-1);
    hInvMassBgk[ipoint]=  new TH1D(*deno);
    hInvMassBgk[ipoint]->SetName(histoname);hInvMassBgk[ipoint]->SetTitle(histoname);
    hInvMassBgk[ipoint]->Sumw2();

    hInvMassSub[ipoint]->Add(hInvMassBgk[ipoint],-1.);
    hInvMassSub[ipoint]->Fit("gaus","","",InvMass1,InvMass2);

    integral1 = hInvMassSub[ipoint]->FindBin( (hInvMassSub[ipoint]->GetFunction("gaus")->GetParameter(1)-2.*hInvMassSub[ipoint]->GetFunction("gaus")->GetParameter(2) ) );
    integral2 = hInvMassSub[ipoint]->FindBin( (hInvMassSub[ipoint]->GetFunction("gaus")->GetParameter(1)+2.*hInvMassSub[ipoint]->GetFunction("gaus")->GetParameter(2) ) );

    // Renormalisation to take into account the resonance statistics
    if ( hInvMassBgk[ipoint]->Integral()>0.) {
      Float_t renormalisation = (hInvMassBgk[ipoint]->Integral() - hInvMassSub[ipoint]->Integral(integral1,integral2))/hInvMassBgk[ipoint]->Integral();
      hInvMassBgk[ipoint]->Scale(renormalisation);
      hInvMassSub[ipoint]->Reset("ICE");
      hInvMassSub[ipoint]->Add(hInvMassAll[ipoint],+1.);
      hInvMassSub[ipoint]->Add(hInvMassBgk[ipoint],-1.);
    }
    if (ipoint!=NumberOfPoints) { 
      Cent[ipoint] = hInvMassSub[ipoint]->GetFunction("gaus")->GetParameter(1);
      e_Cent[ipoint] = hInvMassSub[ipoint]->GetFunction("gaus")->GetParError(1);
      Sigma[ipoint] = hInvMassSub[ipoint]->GetFunction("gaus")->GetParameter(2);
      e_Sigma[ipoint] = hInvMassSub[ipoint]->GetFunction("gaus")->GetParError(2);
      dNdpT[ipoint] = hInvMassSub[ipoint]->Integral(integral1,integral2);
      for(ibin=integral1; ibin<=integral2; ibin++) {
	e_dNdpT[ipoint]+=hInvMassSub[ipoint]->GetBinContent(ibin)*hInvMassSub[ipoint]->GetBinContent(ibin);
      }
      e_dNdpT[ipoint] = TMath::Sqrt( e_dNdpT[ipoint]);
    }
  }

  MUONmassPlot.Close();
    
  TFile out("QuarkoniaPtDistri.root","RECREATE");
  TGraphErrors * gePtDistri = new TGraphErrors(NumberOfPoints,pT,dNdpT,e_pT, e_dNdpT);
  TGraphErrors * geCenter   = new TGraphErrors(NumberOfPoints,pT1,Cent,e_pT1, e_Cent);
  TGraphErrors * geSigma   = new TGraphErrors(NumberOfPoints,pT2,Sigma,e_pT2, e_Sigma);
  gePtDistri->Write();
  geCenter->Write();
  geSigma->Write();
  for(ipoint=0;ipoint<NumberOfPoints; ipoint++) {
    hInvMassAll[ipoint]->Write();
    hInvMassBgk[ipoint]->Write();
    hInvMassSub[ipoint]->Write();
  }
  out.Write();
  out.Close();

  TCanvas *c1 =new TCanvas("c1","c1");
  TH2D * hframe1 = new TH2D("hframe1","hframe1",10,0,12,10,InvMass1,InvMass2);
  hframe1->SetXTitle("Transverse Momentum (GeV/c)");
  hframe1->SetYTitle("Resonance Centroide (GeV/C^{2})");
  hframe1->Draw();
  geCenter->SetMarkerStyle(27);
  geCenter->Draw("p");
  
  TCanvas *c2 =new TCanvas("c2","c2");
  TH2D * hframe2 = new TH2D("hframe2","hframe2",10,0,12,10,0.03,0.200); 
  hframe2->SetXTitle("Transverse Momentum (GeV/c)");
  hframe2->SetYTitle("Resonance Sigma (GeV/C^2)");
  hframe2->Draw();
  geSigma->SetMarkerStyle(28);
  geSigma->Draw("p");
  
  TCanvas *c3 =new TCanvas("c3","c3");
  TH2D * hframe3 = new TH2D("hframe3","hframe3",10,0,12,10,0.1,10000);
  hframe3->SetXTitle("Transverse Momentum (GeV/c)");
  hframe3->SetYTitle("dN/dpT (GeV^{-1})");
  hframe3->Draw();
  c3->SetLogy();
  hframe3->Draw();
  gePtDistri->SetMarkerStyle(29);
  gePtDistri->Draw("p");

  TCanvas *c4 =new TCanvas("c4","c4");
  c4->SetLogy();
  c4->Divide(3,3);
  Int_t i; 
  Float_t izoom1;
  Float_t izoom2; 
  for (i=1;i<10;i++) {
    izoom1 = hInvMassSub[i]->FindBin( (hInvMassSub[i]->GetFunction("gaus")->GetParameter(1)-15.*hInvMassSub[i]->GetFunction("gaus")->GetParameter(2) ) );
    izoom2 = hInvMassSub[i]->FindBin( (hInvMassSub[i]->GetFunction("gaus")->GetParameter(1)+15.*hInvMassSub[i]->GetFunction("gaus")->GetParameter(2) ) );
    c4->cd(i);
    hInvMassAll[i]->SetXTitle("Invariant Mass (GeV/C^{2})");
    hInvMassAll[i]->SetYTitle("Counts (1/25 MeV)");
    hInvMassAll[i]->GetXaxis()->SetRange(izoom1,izoom2);
    hInvMassAll[i]->Draw();
    hInvMassBgk[i]->SetLineColor(4);
    hInvMassBgk[i]->Draw("histo,same");
  }
 
  TCanvas *c5 =new TCanvas("c5","c5"); 
  c5->Divide(3,3);
  Int_t i;
  for (i=1;i<10;i++) {
    c5->cd(i);
    izoom1 = hInvMassSub[i]->FindBin( (hInvMassSub[i]->GetFunction("gaus")->GetParameter(1)-15.*hInvMassSub[i]->GetFunction("gaus")->GetParameter(2) ) );
    izoom2 = hInvMassSub[i]->FindBin( (hInvMassSub[i]->GetFunction("gaus")->GetParameter(1)+15.*hInvMassSub[i]->GetFunction("gaus")->GetParameter(2) ) );
    hInvMassSub[i]->GetXaxis()->SetRange(izoom1,izoom2);
    hInvMassSub[i]->SetXTitle("Invariant Mass (GeV/C^{2})");
    hInvMassSub[i]->SetYTitle("Counts (1/25 MeV)");
    hInvMassSub[i]->Draw();
  }

  TCanvas *c6 =new TCanvas("c6","c6"); 
  c6->Divide(2,1);
  c6->cd(1);
  izoom1 = hInvMassSub[NumberOfPoints]->FindBin( (hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(1)-15.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(2) ) );
  izoom2 = hInvMassSub[NumberOfPoints]->FindBin( (hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(1)+15.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(2) ) );
  c4->cd(i);
  hInvMassAll[NumberOfPoints]->SetXTitle("Invariant Mass (GeV/C^{2})");
  hInvMassAll[NumberOfPoints]->SetYTitle("Counts (1/25 MeV)");
  hInvMassAll[NumberOfPoints]->GetXaxis()->SetRange(izoom1,izoom2);
  hInvMassAll[NumberOfPoints]->Draw();
  hInvMassBgk[NumberOfPoints]->SetLineColor(4);
  hInvMassBgk[NumberOfPoints]->Draw("histo,same");
  c6->cd(2);
  izoom1 = hInvMassSub[NumberOfPoints]->FindBin( (hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(1)-15.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(2) ) );
  izoom2 = hInvMassSub[NumberOfPoints]->FindBin( (hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(1)+15.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(2) ) );
  hInvMassSub[NumberOfPoints]->GetXaxis()->SetRange(izoom1,izoom2);
  hInvMassSub[NumberOfPoints]->SetXTitle("Invariant Mass (GeV/C^{2})");
  hInvMassSub[NumberOfPoints]->SetYTitle("Counts (1/25 MeV)");
  hInvMassSub[NumberOfPoints]->Draw(); 

  //Number of resonances
  
    integral1 = hInvMassSub[NumberOfPoints]->FindBin( (hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(1)-2.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(2) ) );
    integral2 = hInvMassSub[NumberOfPoints]->FindBin( (hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(1)+2.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(2) ) );
    
    Int_t   NumberOfResonances = hInvMassSub[NumberOfPoints]->Integral(integral1,integral2);
    Int_t e_NumberOfResonances = 0;
    for(ibin=integral1; ibin<=integral2; ibin++) {
      e_NumberOfResonances+=hInvMassSub[NumberOfPoints]->GetBinContent(ibin)*hInvMassSub[NumberOfPoints]->GetBinContent(ibin);
    }
    e_NumberOfResonances= TMath::Sqrt( e_NumberOfResonances);
    printf(">>> Number of resonances is %d+-%d in the inv mass range (2sigma) %f, %f \n",NumberOfResonances,e_NumberOfResonances,
	   hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(1)-2.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(2) ,
	   hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(1)+2.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(2) );;

    printf(">>> Centroide is %f +-%f GeV\n", hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(1), hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParError(1));
    printf(">>> Sigma is %f +-%f MeV \n", 1000.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParameter(2), 1000.*hInvMassSub[NumberOfPoints]->GetFunction("gaus")->GetParError(2));

}
