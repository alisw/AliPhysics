void VertexResolutionsFromCmpHF(Int_t pdgSel=421,Int_t nprongsSel=2) 
{
  //
  // Computes secondary vertex resolutions from the ntuple
  // written to CmpHF.root by AliAnalysisTaskSECompareHF
  // A.Dainese
  //

  gStyle->SetOptStat(0);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);

  // input
  TFile *in = new TFile("CmpHF.root");
  TNtuple *nt = (TNtuple*)in->Get("fNtupleCmp");
  Float_t pdg,nprongs,ptRec,vxRec,vxTrue,errVx,vyRec,vyTrue,errVy,vzRec,vzTrue,errVz,mRec;
  nt->SetBranchAddress("pdg",&pdg);
  nt->SetBranchAddress("nprongs",&nprongs);
  nt->SetBranchAddress("PtRec",&ptRec);
  nt->SetBranchAddress("VxRec",&vxRec);
  nt->SetBranchAddress("VxTrue",&vxTrue);
  nt->SetBranchAddress("ErrVx",&errVx);
  nt->SetBranchAddress("VyRec",&vyRec);
  nt->SetBranchAddress("VyTrue",&vyTrue);
  nt->SetBranchAddress("ErrVy",&errVy);
  nt->SetBranchAddress("VzRec",&vzRec);
  nt->SetBranchAddress("VzTrue",&vzTrue);
  nt->SetBranchAddress("ErrVz",&errVz);
  nt->SetBranchAddress("Mrec",&mRec);

  Int_t entries = (Int_t)nt->GetEntries();

  // histograms for 15 pt bins
  TH1F *hV2x_Reco_ = new TH1F[15];
  TH1F *hV2y_Reco_ = new TH1F[15];
  TH1F *hV2z_Reco_ = new TH1F[15];

  TH1F *hV2pullx_Reco_ = new TH1F[15];
  TH1F *hV2pully_Reco_ = new TH1F[15];
  TH1F *hV2pullz_Reco_ = new TH1F[15];

  TH1F *hChi2_ = new TH1F[15];

  TH1F *hDum = new TH1F("hDum","",100,-1000,1000);
  TH1F *hDum2 = new TH1F("hDum2","",100,-5,5);
  for(Int_t j=0; j<15; j++) {
    hV2x_Reco_[j] = *hDum;
    hV2y_Reco_[j] = *hDum;
    hV2z_Reco_[j] = *hDum;
    hV2pullx_Reco_[j] = *hDum2;
    hV2pully_Reco_[j] = *hDum2;
    hV2pullz_Reco_[j] = *hDum2;
  }
  delete hDum; hDum=0;
  delete hDum2; hDum2=0;



  Int_t bin;
  // loop on candidates
  for(Int_t i=0; i<entries; i++) {
    nt->GetEvent(i);
    if(i%10000==0) cout<<" Processing entry "<<i<<" of "<<entries<<endl;
    if(TMath::Abs(pdg)!=pdgSel || nprongs!=nprongsSel) continue;
    bin = GetPtBin(ptRec);
    if(bin==-1) continue;
    hV2x_Reco_[bin].Fill(10000.*(vxRec-vxTrue));
    hV2y_Reco_[bin].Fill(10000.*(vyRec-vyTrue));
    hV2z_Reco_[bin].Fill(10000.*(vzRec-vzTrue));
    hV2pullx_Reco_[bin].Fill((vxRec-vxTrue)/errVx);
    hV2pully_Reco_[bin].Fill((vyRec-vyTrue)/errVy);
    hV2pullz_Reco_[bin].Fill((vzRec-vzTrue)/errVz);
  }


  // output histograms
  Int_t nbins = 15;
  Float_t xbins[16]={0,0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8,10,12,14};
  TH1F *hResV2x_Reco = new TH1F("hResV2x_Reco","",nbins,xbins);
  hResV2x_Reco->SetMaximum(300);
  hResV2x_Reco->SetMinimum(0);
  hResV2x_Reco->SetXTitle("p_{t} D^{0} [GeV/c]");
  hResV2x_Reco->SetYTitle("resolution D^{0}#rightarrow K#pi vertex [#mu m]");
  hResV2x_Reco->SetLineWidth(1);
  hResV2x_Reco->SetMarkerStyle(20);

  TH1F *hResV2y_Reco = (TH1F*)hResV2x_Reco->Clone("hResV2y_Reco");
  hResV2y_Reco->SetYTitle("resolution y  D^{0} vertex [#mu m]");

  TH1F *hResV2z_Reco = (TH1F*)hResV2x_Reco->Clone("hResV2z_Reco");
  hResV2z_Reco->SetYTitle("resolution z  D^{0} vertex [#mu m]");
  hResV2z_Reco->SetMaximum(200);


  TH1F *hPullV2x_Reco = new TH1F("hPullV2x_Reco","",nbins,xbins);
  hPullV2x_Reco->SetMaximum(4);
  hPullV2x_Reco->SetMinimum(0);
  hPullV2x_Reco->SetXTitle("p_{t} D^{0} [GeV/c]");
  hPullV2x_Reco->SetYTitle("pull D^{0}#rightarrow K#pi vertex");
  hPullV2x_Reco->SetLineWidth(1);
  hPullV2x_Reco->SetMarkerStyle(20);

  TH1F *hPullV2y_Reco = (TH1F*)hPullV2x_Reco->Clone("hPullV2y_Reco");
  hPullV2y_Reco->SetYTitle("pull y  D^{0} vertex");
  hPullV2y_Reco->SetMaximum(4);

  TH1F *hPullV2z_Reco = (TH1F*)hPullV2x_Reco->Clone("hPullV2z_Reco");
  hPullV2z_Reco->SetYTitle("pull z  D^{0} vertex");
  hPullV2z_Reco->SetMaximum(4);


  // fit gaussian
  TF1 *g = new TF1("g","gaus",-10000.,10000.);
  TCanvas *cccx = new TCanvas("cccx","cccx",0,0,800,800);
  cccx->Divide(5,3);
  TCanvas *cccz = new TCanvas("cccz","cccz",0,0,800,800);
  cccz->Divide(5,3);
  // fit
  for(Int_t j=0;j<15;j++) {
    cccz->cd(j+1);
    //  Secondary vertex   RECONSTRUCTED

    g->SetRange(-3.*hV2y_Reco_[j].GetRMS(),+3.*hV2y_Reco_[j].GetRMS());
    if(j>-9) hV2y_Reco_[j].Rebin(4);
    hV2y_Reco_[j].Fit("g","R,Q");
    hResV2y_Reco->SetBinContent(j+1,hV2y_Reco_[j].GetRMS());
    hResV2y_Reco->SetBinContent(j+1,g->GetParameter(2));
    hResV2y_Reco->SetBinError(j+1,g->GetParError(2));

    g->SetRange(-3.*hV2z_Reco_[j].GetRMS(),+3.*hV2z_Reco_[j].GetRMS());
    if(j>-9) hV2z_Reco_[j].Rebin(4);
    hV2z_Reco_[j].Fit("g","R,Q");
    hResV2z_Reco->SetBinContent(j+1,hV2z_Reco_[j].GetRMS());
    hResV2z_Reco->SetBinContent(j+1,g->GetParameter(2));
    hResV2z_Reco->SetBinError(j+1,g->GetParError(2));

    cccx->cd(j+1);
    g->SetRange(-3.*hV2x_Reco_[j].GetRMS(),+3.*hV2x_Reco_[j].GetRMS());
    if(j>-9) hV2x_Reco_[j].Rebin(4);
    hV2x_Reco_[j].Draw();
    hV2x_Reco_[j].Fit("g","R,Q");
    hResV2x_Reco->SetBinContent(j+1,hV2x_Reco_[j].GetRMS());
    hResV2x_Reco->SetBinContent(j+1,g->GetParameter(2));
    hResV2x_Reco->SetBinError(j+1,g->GetParError(2));

  }



  TCanvas *ccccx = new TCanvas("ccccx","ccccx",0,0,800,800);
  ccccx->Divide(5,3);
  TCanvas *ccccz = new TCanvas("ccccz","ccccz",0,0,800,800);
  ccccz->Divide(5,3);
  // fit
  for(Int_t j=0;j<15;j++) {
    ccccz->cd(j+1);
    g->SetRange(-3.*hV2pully_Reco_[j].GetRMS(),+3.*hV2pully_Reco_[j].GetRMS());
    if(j>9) hV2pully_Reco_[j].Rebin(4);
    hV2pully_Reco_[j].Fit("g","R,Q");
    hPullV2y_Reco->SetBinContent(j+1,g->GetParameter(2));
    hPullV2y_Reco->SetBinError(j+1,g->GetParError(2));

    g->SetRange(-3.*hV2pullz_Reco_[j].GetRMS(),+3.*hV2pullz_Reco_[j].GetRMS());
    if(j>9) hV2pullz_Reco_[j].Rebin(4);
    hV2pullz_Reco_[j].Fit("g","R,Q");
    hPullV2z_Reco->SetBinContent(j+1,g->GetParameter(2));
    hPullV2z_Reco->SetBinError(j+1,g->GetParError(2));

    ccccx->cd(j+1);
    g->SetRange(-3.*hV2pullx_Reco_[j].GetRMS(),+3.*hV2pullx_Reco_[j].GetRMS());
    if(j>9) hV2pullx_Reco_[j].Rebin(4);
    hV2pullx_Reco_[j].Draw();
    hV2pullx_Reco_[j].Fit("g","R,Q");
    hPullV2x_Reco->SetBinContent(j+1,g->GetParameter(2));
    hPullV2x_Reco->SetBinError(j+1,g->GetParError(2));
  }

  
  // Draw Secondary vertex resolution
  TCanvas *cV2 = new TCanvas("cV2","cV2",0,0,700,900);
  cV2->Divide(1,3);
  cV2->cd(1);
  hResV2x_Reco->Draw("p,e");
  cV2->cd(2);
  hResV2y_Reco->Draw("p,e");
  cV2->cd(3);
  hResV2z_Reco->Draw("p,e");

  // Draw Secondary vertex pulls
  TCanvas *cV2pull = new TCanvas("cV2pull","cV2pull",0,0,700,900);
  cV2pull->Divide(1,3);
  cV2pull->cd(1);
  hPullV2x_Reco->Draw("p,e");
  cV2pull->cd(2);
  hPullV2y_Reco->Draw("p,e");
  cV2pull->cd(3);
  hPullV2z_Reco->Draw("p,e");
  

  // Draw Secondary vertex resolution
  TCanvas *cSummary = new TCanvas("cSummary","cSummary",0,0,900,500);
  cSummary->Divide(2,1);
  cSummary->cd(1);
  hResV2x_Reco->SetMarkerStyle(20);
  hResV2x_Reco->SetMarkerColor(1);
  hResV2x_Reco->Draw("p,e");
  hResV2y_Reco->SetMarkerStyle(21);
  hResV2y_Reco->SetMarkerColor(2);
  hResV2y_Reco->Draw("p,e,same");
  hResV2z_Reco->SetMarkerStyle(22);
  hResV2z_Reco->SetMarkerColor(4);
  hResV2z_Reco->Draw("p,e,same");
  TLegend *legg = new TLegend(0.5,0.5,0.9,0.9);
  legg->SetFillStyle(0);
  legg->SetBorderSize(0); 
  legg->AddEntry(hResV2x_Reco,"x coordinate","p");
  legg->AddEntry(hResV2y_Reco,"y coordinate","p");
  legg->AddEntry(hResV2z_Reco,"z coordinate","p");
  legg->Draw();
  cSummary->cd(2);
  hPullV2x_Reco->SetMarkerStyle(20);
  hPullV2x_Reco->SetMarkerColor(1);
  hPullV2x_Reco->Draw("p,e");
  hPullV2y_Reco->SetMarkerStyle(21);
  hPullV2y_Reco->SetMarkerColor(2);
  hPullV2y_Reco->Draw("p,e,same");
  hPullV2z_Reco->SetMarkerStyle(22);
  hPullV2z_Reco->SetMarkerColor(4);
  hPullV2z_Reco->Draw("p,e,same");
  legg->Draw();

  return;
}
//---------------------------------------------------------------------------
Int_t GetPtBin(Double_t pt) {

  if(pt<0.5) return 0; 
  if(pt<1) return 1; 
  if(pt<1.5) return 2; 
  if(pt<2) return 3; 
  if(pt<2.5) return 4; 
  if(pt<3) return 5; 
  if(pt<3.5) return 6; 
  if(pt<4) return 7; 
  if(pt<5) return 8; 
  if(pt<6) return 9; 
  if(pt<7) return 10; 
  if(pt<8) return 11; 
  if(pt<10) return 12; 
  if(pt<12) return 13; 
  if(pt<14) return 14; 
  return -1; 

}












