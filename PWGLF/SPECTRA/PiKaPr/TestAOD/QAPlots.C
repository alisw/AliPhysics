gROOT->LoadMacro("~/works/AOD/AliSpectraAODTrackCuts.cxx+g");
gROOT->LoadMacro("~/works/AOD/AliSpectraAODEventCuts.cxx+g");
gROOT->LoadMacro("~/works/AOD/AliSpectraAODHistoManager.cxx+g");
gROOT->LoadMacro("~/works/AOD/AliSpectraAODPID.cxx+g");
gROOT->LoadMacro("~/works/AOD/AliAnalysisTaskSpectraAOD.cxx+g");

void QAPlots( AliSpectraAODHistoManager* hman_data, AliSpectraAODHistoManager* hman_mc,
	      AliSpectraAODEventCuts* ecuts_data, AliSpectraAODEventCuts* ecuts_mc,
	      AliSpectraAODTrackCuts* tcuts_data, AliSpectraAODTrackCuts* tcuts_mc){
   
  //vtx distr in data and MC before and after event selection
  TCanvas *cVtx=new TCanvas("Vtxdistr","Vtxdistr",700,500);
  TH1F *hVtxBef_data=ecuts_data->GetHistoVtxBefSel();
  hVtxBef_data->Scale(1./ecuts_data->NumberOfProcessedEvents());
  hVtxBef_data->SetTitle(Form("%s - data",hVtxBef_data->GetTitle()));
  TH1F *hVtxBef_mc=ecuts_mc->GetHistoVtxBefSel();
  hVtxBef_mc->Scale(1./ecuts_mc->NumberOfProcessedEvents());
  hVtxBef_mc->SetTitle(Form("%s - mc",hVtxBef_mc->GetTitle()));
  hVtxBef_mc->SetLineColor(2);
  cVtx->Divide(1,2);
  cVtx->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hVtxBef_data->DrawClone("lhist");
  hVtxBef_mc->DrawClone("lhistsame");
  gPad->BuildLegend();
  TH1F *hVtxAft_data=ecuts_data->GetHistoVtxAftSel();
  hVtxAft_data->Scale(1./ecuts_data->NumberOfEvents());
  hVtxAft_data->SetTitle(Form("%s - data",hVtxAft_data->GetTitle()));
  TH1F *hVtxAft_mc=ecuts_mc->GetHistoVtxAftSel();
  hVtxAft_mc->Scale(1./ecuts_mc->NumberOfEvents());
  hVtxAft_mc->SetTitle(Form("%s - mc",hVtxAft_mc->GetTitle()));
  hVtxAft_mc->SetLineColor(2);
  cVtx->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hVtxAft_data->DrawClone("lhist");
  hVtxAft_mc->DrawClone("lhistsame");
  gPad->BuildLegend();
  
  //eta distr in data and MC before and after event selection
  TCanvas *cEta=new TCanvas("Etadistr","Etadistr",700,500);
  TH1F *hEtaBef_data=ecuts_data->GetHistoEtaBefSel();
  hEtaBef_data->Scale(1./ecuts_data->NumberOfProcessedEvents());
  hEtaBef_data->SetTitle(Form("%s - data",hEtaBef_data->GetTitle()));
  TH1F *hEtaBef_mc=ecuts_mc->GetHistoEtaBefSel();
  hEtaBef_mc->Scale(1./ecuts_mc->NumberOfProcessedEvents());
  hEtaBef_mc->SetTitle(Form("%s - mc",hEtaBef_mc->GetTitle()));
  hEtaBef_mc->SetLineColor(2);
  cEta->Divide(1,2);
  cEta->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hEtaBef_data->DrawClone("lhist");
  hEtaBef_mc->DrawClone("lhistsame");
  gPad->BuildLegend();
  TH1F *hEtaAft_data=ecuts_data->GetHistoEtaAftSel();
  hEtaAft_data->Scale(1./ecuts_data->NumberOfEvents());
  hEtaAft_data->SetTitle(Form("%s - data",hEtaAft_data->GetTitle()));
  TH1F *hEtaAft_mc=ecuts_mc->GetHistoEtaAftSel();
  hEtaAft_mc->Scale(1./ecuts_mc->NumberOfEvents());
  hEtaAft_mc->SetTitle(Form("%s - mc",hEtaAft_mc->GetTitle()));
  hEtaAft_mc->SetLineColor(2);
  cEta->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hEtaAft_data->DrawClone("lhist");
  hEtaAft_mc->DrawClone("lhistsame");
  gPad->BuildLegend();

  //Nch distr in data and MC before and after event selection
  TCanvas *cNCh=new TCanvas("NChdistr","NChdistr",700,500);
  gPad->SetGridy();
  gPad->SetGridx();
  TH1F *hNChAft_data=ecuts_data->GetHistoNChAftSel();
  hNChAft_data->Scale(1./hNChAft_data->GetEntries());
  hNChAft_data->SetTitle(Form("%s - data",hNChAft_data->GetTitle()));
  TH1F *hNChAft_mc=ecuts_mc->GetHistoNChAftSel();
  hNChAft_mc->Scale(1./hNChAft_mc->GetEntries());
  hNChAft_mc->SetTitle(Form("%s - mc",hNChAft_mc->GetTitle()));
  hNChAft_mc->SetLineColor(2);
  hNChAft_data->DrawClone("lhist");
  hNChAft_mc->DrawClone("lhistsame");
  gPad->BuildLegend();

  //Eta Phi at high Pt in data and Monte Carlo
  TCanvas *cEtaPhi=new TCanvas("EtaPhi","EtaPhi",700,500);
  cEtaPhi->Divide(2,1);
  TH2F *hEtaPhi_data=tcuts_data->GetHistoEtaPhiHighPt()->Clone("hEtaPhi_data");
  hEtaPhi_data->Rebin2D(4);
  hEtaPhi_data->Scale(1./hEtaPhi_data->GetEntries());
  cEtaPhi->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hEtaPhi_data->DrawClone("colz");
  TH2F *hEtaPhi_mc=tcuts_mc->GetHistoEtaPhiHighPt()->Clone("hEtaPhi_mc");
  hEtaPhi_mc->Rebin2D(4);
  hEtaPhi_mc->Scale(1./hEtaPhi_mc->GetEntries());
  cEtaPhi->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hEtaPhi_mc->DrawClone("colz");
  TCanvas *cRatioEtaPhi=new TCanvas("RatioEtaPhi","RatioEtaPhi",700,500);
  gPad->SetGridy();
  gPad->SetGridx();
  hEtaPhi_data->Divide(hEtaPhi_mc);
  hEtaPhi_data->SetMaximum(2);
  hEtaPhi_data->SetMinimum(-0.00000001);
  hEtaPhi_data->DrawClone("colz");

  //Track selection in data and Monte Carlo
  TCanvas *cTrackCuts=new TCanvas("TrackCuts","TrackCuts",700,500);
  cTrackCuts->Divide(1,2);
  TH1F *hTrCuts_data=new TH1F("hTrCuts_data","hTrCuts_data",11,0,11);
  for(Int_t ibin=1;ibin<=tcuts_data->GetHistoCuts()->GetNbinsX();ibin++){
    hTrCuts_data->SetBinContent(ibin,tcuts_data->GetHistoCuts()->GetBinContent(ibin)/tcuts_data->GetHistoCuts()->GetBinContent(6));
    hTrCuts_data->GetXaxis()->SetBinLabel(ibin,tcuts_data->GetHistoCuts()->GetXaxis()->GetBinLabel(ibin));
  }
  TH1F *hTrCuts_mc=new TH1F("hTrCuts_mc","hTrCuts_mc",11,0,11);
  for(Int_t ibin=1;ibin<=tcuts_mc->GetHistoCuts()->GetNbinsX();ibin++){
    hTrCuts_mc->SetBinContent(ibin,tcuts_mc->GetHistoCuts()->GetBinContent(ibin)/tcuts_mc->GetHistoCuts()->GetBinContent(6));
    hTrCuts_mc->GetXaxis()->SetBinLabel(ibin,tcuts_mc->GetHistoCuts()->GetXaxis()->GetBinLabel(ibin));
  }
  hTrCuts_mc->SetLineColor(2);
  cTrackCuts->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hTrCuts_data->SetMinimum(0);
  hTrCuts_data->DrawClone();
  hTrCuts_mc->DrawClone("same");
  gPad->BuildLegend();
  cTrackCuts->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hTrCuts_data->Divide(hTrCuts_mc);
  hTrCuts_data->SetTitle("DATA/MC");
  hTrCuts_data->DrawClone();
  gPad->BuildLegend();


 
  //dedx in data and MC
  TCanvas *cPIDSig=new TCanvas("cPIDSig","cPIDSig",700,500);
  cPIDSig->Divide(2,2);
  cPIDSig->cd(1);
  TH2F *PIDSig_data = (TH2F*)((TH2F*)hman_data->GetPIDHistogram("hHistPIDTPC"))->Clone();
  PIDSig_data->SetYTitle("TPC signal");
  gPad->SetLogz();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_data->DrawClone("colz");
  cPIDSig->cd(2);
  TH2F *PIDSig_mc = (TH2F*)((TH2F*)hman_mc->GetPIDHistogram("hHistPIDTPC"))->Clone();
  PIDSig_mc->SetYTitle("TPC signal");
  gPad->SetLogz();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_mc->DrawClone("colz");
  cPIDSig->cd(3);
  TH2F *PIDSig_data = (TH2F*)((TH2F*)hman_data->GetPIDHistogram("hHistPIDTOF"))->Clone();
  PIDSig_data->SetYTitle("TOF signal");
  gPad->SetLogz();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_data->DrawClone("colz");
  cPIDSig->cd(4);
  TH2F *PIDSig_mc = (TH2F*)((TH2F*)hman_mc->GetPIDHistogram("hHistPIDTOF"))->Clone();
  PIDSig_mc->SetYTitle("TOF signal/100");
  gPad->SetLogz();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_mc->DrawClone("colz");

  //dedx projection in data and MC
  Double_t Proj1[2]={0.6,0.7};
  Double_t Proj2[2]={1.1,1.2};
  TCanvas *cPIDSigProjection=new TCanvas("cPIDSigProjection","cPIDSigProjection",700,500);
  cPIDSigProjection->Divide(2,2);
  //TPC
  TH2F *PIDSig_data = (TH2F*)((TH2F*)hman_data->GetPIDHistogram("hHistPIDTPC"))->Clone();
  TH1F *PIDSig_data_Proj1=(TH1F*)PIDSig_data->ProjectionY(Form("TPC, data [%.1f,%.1f]",Proj1[0],Proj1[1]),
  							  PIDSig_data->GetXaxis()->FindBin(Proj1[0]),PIDSig_data->GetXaxis()->FindBin(Proj1[1]));
  TH1F *PIDSig_data_Proj2=(TH1F*)PIDSig_data->ProjectionY(Form("TPC, data [%.1f,%.1f]",Proj2[0],Proj2[1]),
  							  PIDSig_data->GetXaxis()->FindBin(Proj2[0]),PIDSig_data->GetXaxis()->FindBin(Proj2[1]));
  PIDSig_data_Proj1->SetTitle(Form("TPC, data [%.1f,%.1f]",Proj1[0],Proj1[1]));
  PIDSig_data_Proj2->SetTitle(Form("TPC, data [%.1f,%.1f]",Proj2[0],Proj2[1]));
  TH2F *PIDSig_mc = (TH2F*)((TH2F*)hman_mc->GetPIDHistogram("hHistPIDTPC"))->Clone();
  TH1F *PIDSig_mc_Proj1=(TH1F*)PIDSig_mc->ProjectionY(Form("TPC, mc [%.1f,%.1f]",Proj1[0],Proj1[1]),
						      PIDSig_mc->GetXaxis()->FindBin(Proj1[0]),PIDSig_mc->GetXaxis()->FindBin(Proj1[1]));
  TH1F *PIDSig_mc_Proj2=(TH1F*)PIDSig_mc->ProjectionY(Form("TPC, mc [%.1f,%.1f]",Proj2[0],Proj2[1]),
						      PIDSig_mc->GetXaxis()->FindBin(Proj2[0]),PIDSig_mc->GetXaxis()->FindBin(Proj2[1]));
  PIDSig_mc_Proj1->SetTitle(Form("TPC, mc [%.1f,%.1f]",Proj1[0],Proj1[1]));
  PIDSig_mc_Proj2->SetTitle(Form("TPC, mc [%.1f,%.1f]",Proj2[0],Proj2[1]));
  PIDSig_mc_Proj1->SetLineColor(2);
  PIDSig_mc_Proj2->SetLineColor(2);
  cPIDSigProjection->cd(1);
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetGridx();
  PIDSig_data_Proj1->DrawNormalized("lhist");
  PIDSig_mc_Proj1->DrawNormalized("lhistsame");
  gPad->BuildLegend();
  cPIDSigProjection->cd(2);
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_data_Proj2->DrawNormalized("lhist");
  PIDSig_mc_Proj2->DrawNormalized("lhistsame");
  gPad->BuildLegend();
  //TOF
  TH2F *PIDSig_data = (TH2F*)((TH2F*)hman_data->GetPIDHistogram("hHistPIDTOF"))->Clone();
  TH1F *PIDSig_data_Proj1=(TH1F*)PIDSig_data->ProjectionY(Form("TOF, data [%.1f,%.1f]",Proj1[0],Proj1[1]),
  							  PIDSig_data->GetXaxis()->FindBin(Proj1[0]),PIDSig_data->GetXaxis()->FindBin(Proj1[1]));
  TH1F *PIDSig_data_Proj2=(TH1F*)PIDSig_data->ProjectionY(Form("TOF, data [%.1f,%.1f]",Proj2[0],Proj2[1]),
  							  PIDSig_data->GetXaxis()->FindBin(Proj2[0]),PIDSig_data->GetXaxis()->FindBin(Proj2[1]));
  PIDSig_data_Proj1->SetTitle(Form("TOF, data [%.1f,%.1f]",Proj1[0],Proj1[1]));
  PIDSig_data_Proj2->SetTitle(Form("TOF, data [%.1f,%.1f]",Proj2[0],Proj2[1]));
  TH2F *PIDSig_mc = (TH2F*)((TH2F*)hman_mc->GetPIDHistogram("hHistPIDTOF"))->Clone();
  TH1F *PIDSig_mc_Proj1=(TH1F*)PIDSig_mc->ProjectionY(Form("TOF, mc [%.1f,%.1f]",Proj1[0],Proj1[1]),
						      PIDSig_mc->GetXaxis()->FindBin(Proj1[0]),PIDSig_mc->GetXaxis()->FindBin(Proj1[1]));
  TH1F *PIDSig_mc_Proj2=(TH1F*)PIDSig_mc->ProjectionY(Form("TOF, mc [%.1f,%.1f]",Proj2[0],Proj2[1]),
						      PIDSig_mc->GetXaxis()->FindBin(Proj2[0]),PIDSig_mc->GetXaxis()->FindBin(Proj2[1]));
  PIDSig_mc_Proj1->SetTitle(Form("TOF, mc [%.1f,%.1f]",Proj1[0],Proj1[1]));
  PIDSig_mc_Proj2->SetTitle(Form("TOF, mc [%.1f,%.1f]",Proj2[0],Proj2[1]));
  PIDSig_mc_Proj1->SetLineColor(2);
  PIDSig_mc_Proj2->SetLineColor(2);
  cPIDSigProjection->cd(3);
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_data_Proj1->DrawNormalized("lhist");
  PIDSig_mc_Proj1->DrawNormalized("lhistsame");
  gPad->BuildLegend();
  cPIDSigProjection->cd(4);
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetGridx();
  PIDSig_data_Proj2->DrawNormalized("lhist");
  PIDSig_mc_Proj2->DrawNormalized("lhistsame");
  gPad->BuildLegend();
  
  //nsig in data and MC
  for(Int_t ipart=0;ipart<3;ipart++){
    TCanvas *cnsig=new TCanvas(Form("cnsig%s",Particle[ipart].Data()),Form("cnsig%s",Particle[ipart].Data()),700,500);
    cnsig->Divide(2,3);
    cnsig->cd(1);
    TH2F *nsig_data = (TH2F*)((TH2F*)hman_data->GetNSigHistogram(Form("hHistNSig%sPtTPC",Particle[ipart].Data())))->Clone();
    nsig_data->GetYaxis()->SetRangeUser(-5,5);
    nsig_data->Scale(1/nsig_data->GetMaximum());
    nsig_data->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_data->DrawClone("colz");
    cnsig->cd(2);
    TH2F *nsig_mc = (TH2F*)((TH2F*)hman_mc->GetNSigHistogram(Form("hHistNSig%sPtTPC",Particle[ipart].Data())))->Clone();
    nsig_mc->GetYaxis()->SetRangeUser(-5,5);
    nsig_mc->Scale(1/nsig_mc->GetMaximum());
    nsig_mc->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_mc->DrawClone("colz");
    cnsig->cd(3);
    TH2F *nsig_data = (TH2F*)((TH2F*)hman_data->GetNSigHistogram(Form("hHistNSig%sPtTOF",Particle[ipart].Data())))->Clone();
    nsig_data->GetYaxis()->SetRangeUser(-5,5);
    nsig_data->Scale(1/nsig_data->GetMaximum());
    nsig_data->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_data->DrawClone("colz");
    cnsig->cd(4);
    TH2F *nsig_mc = (TH2F*)((TH2F*)hman_mc->GetNSigHistogram(Form("hHistNSig%sPtTOF",Particle[ipart].Data())))->Clone();
    nsig_mc->GetYaxis()->SetRangeUser(-5,5);
    nsig_mc->Scale(1/nsig_mc->GetMaximum());
    nsig_mc->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_mc->DrawClone("colz");
    cnsig->cd(5);
    TH2F *nsig_data = (TH2F*)((TH2F*)hman_data->GetNSigHistogram(Form("hHistNSig%sPtTPCTOF",Particle[ipart].Data())))->Clone();
    nsig_data->GetYaxis()->SetRangeUser(0,5);
    nsig_data->Scale(1/nsig_data->GetMaximum());
    nsig_data->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_data->DrawClone("colz");
    cnsig->cd(6);
    TH2F *nsig_mc = (TH2F*)((TH2F*)hman_mc->GetNSigHistogram(Form("hHistNSig%sPtTPCTOF",Particle[ipart].Data())))->Clone();
    nsig_mc->GetYaxis()->SetRangeUser(0,5);
    nsig_mc->Scale(1/nsig_mc->GetMaximum());
    nsig_mc->SetXTitle("Pt (GeV/c)");
    gPad->SetLogz();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_mc->DrawClone("colz");
  }
  
  //NSigma projection in data and MC
  Double_t Proj1[2]={0.4,0.5};
  Double_t Proj2[2]={1.1,1.2};
  TCanvas *cnsigProjection=new TCanvas("cnsigProjection","cnsigProjection",700,500);
  cnsigProjection->Divide(2,3);
  //TPC
  for(Int_t ipart=0;ipart<3;ipart++){
    TH2F *nsig_data = (TH2F*)((TH2F*)hman_data->GetNSigHistogram(Form("hHistNSig%sPtTPC",Particle[ipart].Data())))->Clone();
    TH1F *nsig_data_Proj1=(TH1F*)nsig_data->ProjectionY(Form("TPC NsigProjectionv %s, data [%.1f,%.1f]",Particle[ipart].Data(),Proj1[0],Proj1[1]),
							nsig_data->GetXaxis()->FindBin(Proj1[0]),nsig_data->GetXaxis()->FindBin(Proj1[1]));
    TH1F *nsig_data_Proj2=(TH1F*)nsig_data->ProjectionY(Form("TPC NsigProjection %s, data [%.1f,%.1f]",Particle[ipart].Data(),Proj2[0],Proj2[1]),
							nsig_data->GetXaxis()->FindBin(Proj2[0]),nsig_data->GetXaxis()->FindBin(Proj2[1]));
    nsig_data_Proj1->SetTitle(Form("TPC NsigProjection %s, data [%.1f,%.1f]",Particle[ipart].Data(),Proj1[0],Proj1[1]));
    nsig_data_Proj2->SetTitle(Form("TPC NsigProjection %s, data [%.1f,%.1f]",Particle[ipart].Data(),Proj2[0],Proj2[1]));
    TH2F *nsig_mc = (TH2F*)((TH2F*)hman_mc->GetNSigHistogram(Form("hHistNSig%sPtTPC",Particle[ipart].Data())))->Clone();
    TH1F *nsig_mc_Proj1=(TH1F*)nsig_mc->ProjectionY(Form("TPC NsigProjection %s, mc [%.1f,%.1f]",Particle[ipart].Data(),Proj1[0],Proj1[1]),
						    nsig_mc->GetXaxis()->FindBin(Proj1[0]),nsig_mc->GetXaxis()->FindBin(Proj1[1]));
    TH1F *nsig_mc_Proj2=(TH1F*)nsig_mc->ProjectionY(Form("TPC NsigProjection %s, mc [%.1f,%.1f]",Particle[ipart].Data(),Proj2[0],Proj2[1]),
						    nsig_mc->GetXaxis()->FindBin(Proj2[0]),nsig_mc->GetXaxis()->FindBin(Proj2[1]));
    nsig_mc_Proj1->SetTitle(Form("TPC NsigProjection %s, mc [%.1f,%.1f]",Particle[ipart].Data(),Proj1[0],Proj1[1]));
    nsig_mc_Proj2->SetTitle(Form("TPC NsigProjection %s, mc [%.1f,%.1f]",Particle[ipart].Data(),Proj2[0],Proj2[1]));
    nsig_mc_Proj1->SetLineColor(2);
    nsig_mc_Proj2->SetLineColor(2);
    cnsigProjection->cd(2*ipart+1);
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_data_Proj1->DrawNormalized("lhist");
    nsig_mc_Proj1->DrawNormalized("lhistsame");
    gPad->BuildLegend();
    cnsigProjection->cd(2*ipart+2);
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->SetGridx();
    nsig_data_Proj2->DrawNormalized("lhist");
    nsig_mc_Proj2->DrawNormalized("lhistsame");
    gPad->BuildLegend();
  }
  return;

  //Muon over Pion Ratio
  Printf("\n\n-> Muon Over Pion");
  TCanvas *cMu=new TCanvas("cMu","cMu");
  TH1F *hMuOverPi[2];
  TH1F *hMuOverPi_bis[2];
  for(Int_t icharge=0;icharge<2;icharge++){
    TString hname=Form("hHistPtRecTruePrimaryMuon%s",Sign[icharge].Data());
    hMuOverPi[icharge]=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
    hname=Form("hHistPtRecTruePrimaryPion%s",Sign[icharge].Data());
    hMuOverPi[icharge]->Divide((TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone());
    hname=Form("hHistPtRecTrueMuon%s",Sign[icharge].Data());
    hMuOverPi_bis[icharge]=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
    hname=Form("hHistPtRecTruePrimaryPion%s",Sign[icharge].Data());
    hMuOverPi_bis[icharge]->Divide((TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone());
    hMuOverPi_bis[icharge]->SetMarkerColor(2);
    hMuOverPi_bis[icharge]->SetLineColor(2);
    if(icharge==0)hMuOverPi_bis[icharge]->DrawClone();
    else hMuOverPi_bis[icharge]->DrawClone("same");
    hMuOverPi[icharge]->DrawClone("same");
  }
  
  //Contamination
  Printf("\n\n-> Contamination from MC");
  TH1F *Cont[6];
  TCanvas *ccont=new TCanvas("ccont","ccont",700,500);
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      TString hname=Form("hHistPtRecTrue%s%s",Particle[ipart].Data(),Sign[icharge].Data());
      Printf("Getting %s",hname.Data());
      Cont[index] =(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
      Cont[index]->SetName(Form("Cont_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      Cont[index]->SetTitle(Form("RecTrue/RecSigma %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
      Cont[index]->SetMarkerStyle(Marker[index]);
      Cont[index]->SetMarkerColor(Color[ipart]);
      Cont[index]->SetLineColor(Color[ipart]);
      hname=Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
      Printf("... and divide it by %s",hname.Data());
      Cont[index]->Divide((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1));
      if(index==0)Cont[index]->DrawClone();
      else Cont[index]->DrawClone("same");
      //Spectra[index]->Multiply(Cont[index]);
    }
  } 
  gPad->BuildLegend();


  //Raw yield
  TH1F *hRaw_data_allCh=(TH1F*)((TH1F*) hman_data->GetPtHistogram1D("hHistPtRec",-1,-1))->Clone();
  hRaw_data_allCh->Scale(1./ecuts_data->NumberOfEvents(),"width");
  TH1F *hRaw_mc_allCh=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D("hHistPtRec",-1,-1))->Clone();
  hRaw_mc_allCh->Scale(1./ecuts_mc->NumberOfEvents(),"width");
      
  TCanvas *cRaw=new TCanvas("cRaw","cRaw",700,500);
  cRaw->Divide(2,2);
  cRaw->cd(1);
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetGridx();
  hRaw_data_allCh->DrawClone();
  hRaw_mc_allCh->DrawClone("lhistsame");
  cRaw->cd(2);
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetGridx();
  hRaw_data_allCh->DrawClone();
  hRaw_mc_allCh->DrawClone("lhistsame"); 
  hRaw_data_allCh->Divide(hRaw_mc_allCh);
  hRaw_data_allCh->SetLineStyle(2);
  cRaw->cd(3);
  gPad->SetGridy();
  gPad->SetGridx();
  hRaw_data_allCh->DrawClone("lhist");
  cRaw->cd(4);
  gPad->SetGridy();
  gPad->SetGridx();
  hRaw_data_allCh->DrawClone("lhist");
  TH1F *hRaw_data[6];
  TH1F *hRaw_mc[6];
  for(Int_t icharge=0;icharge<2;icharge++){
    for(Int_t ipart=0;ipart<3;ipart++){
      Int_t index=ipart+3*icharge;
      TString hname=Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()); //MC correction for Prim+Cont+Eff
      hRaw_data[index]=(TH1F*)((TH1F*) hman_data->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
      hRaw_data[index]->SetMarkerStyle(Marker[index]);
      hRaw_data[index]->SetMarkerColor(Color[ipart]);
      hRaw_data[index]->SetLineColor(Color[ipart]);
      hRaw_data[index]->Scale(1./ecuts_data->NumberOfEvents(),"width");
      hRaw_mc[index]=(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
      hRaw_mc[index]->SetMarkerStyle(Marker[index]);
      hRaw_mc[index]->SetMarkerColor(Color[ipart]);
      hRaw_mc[index]->SetLineColor(Color[ipart]);
      hRaw_mc[index]->Scale(1./ecuts_mc->NumberOfEvents(),"width");
      for(Int_t ibin=0;ibin<hRaw_data[index]->GetNbinsX();ibin++){
	if(hRaw_data[index]->GetBinCenter(ibin)<Range[ipart]){
	  hRaw_data[index]->SetBinContent(ibin,0);
	  hRaw_data[index]->SetBinError(ibin,0);
	  hRaw_mc[index]->SetBinContent(ibin,0);
	  hRaw_mc[index]->SetBinError(ibin,0);
	}
      }
      cRaw->cd(icharge+1);
      hRaw_data[index]->DrawClone("same");
      hRaw_mc[index]->DrawClone("lhistsame");
      hRaw_data[index]->Divide(hRaw_mc[index]);
      cRaw->cd(icharge+3);
      hRaw_data[index]->DrawClone("lhistsame");
    }
  } 
  

  
  //PID Efficiency //TO BE IMPLEMENTED, need to add TrueSigma in the Task
  // Printf("\n\n-> PID efficiency from MC");
  // TH1F *PIDEff[6];
  // TCanvas *ccont=new TCanvas("ccont","ccont",700,500);
  // for(Int_t icharge=0;icharge<2;icharge++){
  //   for(Int_t ipart=0;ipart<3;ipart++){
  //     Int_t index=ipart+3*icharge;
  //     TString hname=Form("hHistPtRecTrue%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  //     Printf("Getting %s",hname.Data());
  //     PIDEff[index] =(TH1F*)((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
  //     PIDEff[index]->SetName(Form("PIDEff_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  //     PIDEff[index]->SetTitle(Form("RecTrue/RecSigma %s%s",Particle[ipart].Data(),Sign[icharge].Data()));
  //     PIDEff[index]->SetMarkerStyle(Marker[index]);
  //     PIDEff[index]->SetMarkerColor(Color[ipart]);
  //     PIDEff[index]->SetLineColor(Color[ipart]);
  //     hname=Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
  //     Printf("... and divide it by %s",hname.Data());
  //     PIDEff[index]->Divide((TH1F*) hman_mc->GetPtHistogram1D(hname.Data(),-1,-1));
  //     if(index==0)PIDEff[index]->DrawClone();
  //     else PIDEff[index]->DrawClone("same");
  //     //Spectra[index]->Multiply(PIDEff[index]);
  //   }
  // } 
  // gPad->BuildLegend();

  //qvec in data and MC
  // for (Int_t icharge=0;icharge<2;icharge++){
  //   TCanvas *cqvec=new TCanvas(Form("cqvec%s",Charge[icharge].Data()),Form("cqvec%s",Charge[icharge].Data()),700,500);
  //   cqvec->Divide(2,1);
  //   cqvec->cd(1);
  //   TString hname=Form("hHistqVec%s",Charge[icharge].Data());
  //   TH2F *qvec_data = (TH2F*)((TH2F*)hman_data->GetqVecHistogram(hname.Data()))->Clone();
  //   gPad->SetLogz();
  //   gPad->SetGridy();
  //   gPad->SetGridx();
  //   qvec_data->DrawClone("colz");
  //   cqvec->cd(2);
  //   TH2F *qvec_mc = (TH2F*)((TH2F*)hman_mc->GetqVecHistogram(hname.Data()))->Clone();
  //   gPad->SetLogz();
  //   gPad->SetGridy();
  //   gPad->SetGridx();
  //   qvec_mc->DrawClone("colz");
    
  //   TCanvas *cProjqvec=new TCanvas(Form("cProjqvec%s",Charge[icharge].Data()),Form("cProjqvec%s",Charge[icharge].Data()),700,500);
  //   cProjqvec->Divide(3,2);
  //   for(Int_t iproj=0;iproj<3;iproj++){
  //     TH1F *proj_data=(TH1F*)qvec_data->ProjectionX("data",qvec_data->GetYaxis()->FindBin(projl[iproj]),qvec_data->GetYaxis()->FindBin(proju[iproj]));
  //     if(proj_data->GetEntries()==0)continue;
  //     proj_data->Scale(1/proj_data->GetEntries());
  //     proj_data->SetTitle(Form("data q%s [%.0f,%.0f]",Charge[icharge].Data(),projl[iproj],proju[iproj]));
  //     TH1F *proj_mc=(TH1F*)qvec_mc->ProjectionX("mc",qvec_mc->GetYaxis()->FindBin(projl[iproj]),qvec_mc->GetYaxis()->FindBin(proju[iproj]));
  //     proj_mc->Scale(1/proj_mc->GetEntries());
  //     proj_mc->SetTitle(Form("mc q%s [%.0f,%.0f]",Charge[icharge].Data(),projl[iproj],proju[iproj]));
  //     proj_mc->SetLineColor(2);
  //     cProjqvec->cd(iproj+1);
  //     gPad->SetGridy();
  //     gPad->SetGridx();
  //     proj_data->DrawClone();
  //     proj_mc->DrawClone("same");
  //     gPad->BuildLegend();
  //     proj_data->Divide(proj_mc);
  //     cProjqvec->cd(iproj+4);
  //     proj_data->DrawClone();
  //   }
  // }
  
  
  
} 
