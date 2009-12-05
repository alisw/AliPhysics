void PlotITSTrackingHists(TString fname="ITS.Performance.root") {
  //
  // Macro to plot the histos from the task AliAnalysisTaskITSTrackingCheck
  // A. Dainese 28.11.09
  // 

  gStyle->SetOptStat(0);

  TFile *f= new TFile(fname.Data());

  TList *list=(TList*)f->Get("cOutput");
  TH1F *fHistNclsITSSA = (TH1F*)list->FindObject("fHistNclsITSSA");
  TH1F *fHistClusterMapITSSA = (TH1F*)list->FindObject("fHistClusterMapITSSA");
  TH1F *fHistClusterMapITSSAok = (TH1F*)list->FindObject("fHistClusterMapITSSAok");
  TH1F *fHistClusterMapITSSAbad = (TH1F*)list->FindObject("fHistClusterMapITSSAbad");
  TH1F *fHistClusterMapITSSAskipped = (TH1F*)list->FindObject("fHistClusterMapITSSAskipped");
  TH1F *fHistClusterMapITSSAoutinz = (TH1F*)list->FindObject("fHistClusterMapITSSAoutinz");
  TH1F *fHistClusterMapITSSAokoutinzbad = (TH1F*)list->FindObject("fHistClusterMapITSSAokoutinzbad");
  TH1F *fHistClusterMapITSSAnorefit = (TH1F*)list->FindObject("fHistClusterMapITSSAnorefit");
  TH1F *fHistClusterMapITSSAnocls = (TH1F*)list->FindObject("fHistClusterMapITSSAnocls");
  TH1F *fHistNclsITSSAInAcc = (TH1F*)list->FindObject("fHistNclsITSSAInAcc");
  TH1F *fHistClusterMapITSSAInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAInAcc");
  TH1F *fHistClusterMapITSSAokInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAokInAcc");
  TH1F *fHistClusterMapITSSAbadInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAbadInAcc");
  TH1F *fHistClusterMapModuleITSSAokInAcc = (TH1F*)list->FindObject("fHistClusterMapModuleITSSAokInAcc");
  TH1F *fHistClusterMapModuleITSSAbadInAcc = (TH1F*)list->FindObject("fHistClusterMapModuleITSSAbadInAcc");
  TH1F *fHistClusterMapITSSAskippedInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAskippedInAcc");
  TH1F *fHistClusterMapITSSAoutinzInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAoutinzInAcc");
  TH1F *fHistClusterMapITSSAokoutinzbadInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAokoutinzbadInAcc");
  TH1F *fHistClusterMapITSSAnorefitInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAnorefitInAcc");
  TH1F *fHistClusterMapITSSAnoclsInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAnoclsInAcc");
  TH1F *fHistClusterMapModuleITSSAnoclsInAcc = (TH1F*)list->FindObject("fHistClusterMapModuleITSSAnoclsInAcc");
  TH1F *fHistNclsITSMI = (TH1F*)list->FindObject("fHistNclsITSMI");
  TH1F *fHistClusterMapITSMI = (TH1F*)list->FindObject("fHistClusterMapITSMI");
  TH1F *fHistClusterMapITSMIok = (TH1F*)list->FindObject("fHistClusterMapITSMIok");
  TH1F *fHistClusterMapITSMIbad = (TH1F*)list->FindObject("fHistClusterMapITSMIbad");
  TH1F *fHistClusterMapITSMIskipped = (TH1F*)list->FindObject("fHistClusterMapITSMIskipped");
  TH1F *fHistClusterMapITSMIoutinz = (TH1F*)list->FindObject("fHistClusterMapITSMIoutinz");
  TH1F *fHistClusterMapITSMIokoutinzbad = (TH1F*)list->FindObject("fHistClusterMapITSMIokoutinzbad");
  TH1F *fHistClusterMapITSMInorefit = (TH1F*)list->FindObject("fHistClusterMapITSMInorefit");
  TH1F *fHistClusterMapITSMInocls = (TH1F*)list->FindObject("fHistClusterMapITSMInocls");

  TH1F *fHistPhiTPCInAcc = (TH1F*)list->FindObject("fHistPhiTPCInAcc");
  TH1F *fHistPhiITSMIokbadoutinz6InAcc = (TH1F*)list->FindObject("fHistPhiITSMIokbadoutinz6InAcc");

  TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIokbadoutinz6InAcc = (TH1F*)list->FindObject("fHistPtITSMIokbadoutinz6InAcc");
  TH1F *fHistPtITSMIokbadoutinz5InAcc = (TH1F*)list->FindObject("fHistPtITSMIokbadoutinz5InAcc");
  TH1F *fHistPtITSMIokbadoutinz4InAcc = (TH1F*)list->FindObject("fHistPtITSMIokbadoutinz4InAcc");

  //---------------------------------------------------------------

  TLegend *l1=new TLegend(0.5,0.5,0.9,0.9);
  TLegend *l2=new TLegend(0.5,0.5,0.9,0.9);

  TCanvas *c1= new TCanvas("c1","c1",10,10,600,500);
  fHistNclsITSSA->SetMinimum(0);
  fHistNclsITSSA->SetLineColor(1);
  l1->AddEntry(fHistNclsITSSA,"ITS-SA","l");
  fHistNclsITSSA->Draw();
  fHistNclsITSSAInAcc->SetLineColor(4);
  l1->AddEntry(fHistNclsITSSAInAcc,"ITS-SA in acc.","l");
  fHistNclsITSSAInAcc->Draw("same");
  fHistNclsITSMI->SetLineColor(2);
  l1->AddEntry(fHistNclsITSMI,"ITS from TPC","l");
  fHistNclsITSMI->Draw("same");
  l1->Draw();


  TCanvas *c2 =new TCanvas("c2","c2",10,10,1200,800);
  c2->Divide(3,2);
  c2->cd(1);
  //
  fHistClusterMapITSSAokoutinzbad->SetLineColor(1);
  fHistClusterMapITSSAokoutinzbad->SetMarkerColor(1);
  fHistClusterMapITSSAokoutinzbad->SetMarkerStyle(20);
  fHistClusterMapITSSAokoutinzbad->Draw();
  fHistClusterMapITSMIokoutinzbad->SetLineColor(2);
  fHistClusterMapITSMIokoutinzbad->SetMarkerColor(2);
  fHistClusterMapITSMIokoutinzbad->SetMarkerStyle(20);
  fHistClusterMapITSMIokoutinzbad->Draw("same");
  l1->Draw();
  //
  c2->cd(2);
  fHistClusterMapITSSAok->SetLineColor(1);
  fHistClusterMapITSSAok->SetMarkerColor(1);
  fHistClusterMapITSSAok->SetMarkerStyle(21);
  fHistClusterMapITSSAok->Draw();
  fHistClusterMapITSMIok->SetLineColor(2);
  fHistClusterMapITSMIok->SetMarkerColor(2);
  fHistClusterMapITSMIok->SetMarkerStyle(21);
  fHistClusterMapITSMIok->Draw("same");
  //
  c2->cd(3);
  fHistClusterMapITSSAoutinz->SetLineColor(1);
  fHistClusterMapITSSAoutinz->SetMarkerColor(1);
  fHistClusterMapITSSAoutinz->SetMarkerStyle(22);
  fHistClusterMapITSSAoutinz->Draw();
  fHistClusterMapITSMIoutinz->SetLineColor(2);
  fHistClusterMapITSMIoutinz->SetMarkerColor(2);
  fHistClusterMapITSMIoutinz->SetMarkerStyle(22);
  fHistClusterMapITSMIoutinz->Draw("same");
  //
  c2->cd(4);
  fHistClusterMapITSSAbad->SetLineColor(1);
  fHistClusterMapITSSAbad->SetMarkerColor(1);
  fHistClusterMapITSSAbad->SetMarkerStyle(23);
  fHistClusterMapITSSAbad->Draw();
  fHistClusterMapITSMIbad->SetLineColor(2);
  fHistClusterMapITSMIbad->SetMarkerColor(2);
  fHistClusterMapITSMIbad->SetMarkerStyle(23);
  fHistClusterMapITSMIbad->Draw("same");
  //
  c2->cd(5);
  fHistClusterMapITSSAnocls->SetLineColor(1);
  fHistClusterMapITSSAnocls->SetMarkerColor(1);
  fHistClusterMapITSSAnocls->SetMarkerStyle(24);
  fHistClusterMapITSSAnocls->Draw();
  fHistClusterMapITSMInocls->SetLineColor(2);
  fHistClusterMapITSMInocls->SetMarkerColor(2);
  fHistClusterMapITSMInocls->SetMarkerStyle(24);
  fHistClusterMapITSMInocls->Draw("same");

  TCanvas *c3 =new TCanvas("c3","c3",10,10,1200,400);
  fHistClusterMapModuleITSSAokInAcc->SetLineColor(1);
  fHistClusterMapModuleITSSAokInAcc->Draw();
  l2->AddEntry(fHistClusterMapModuleITSSAokInAcc,"ok","l");
  fHistClusterMapModuleITSSAbadInAcc->SetLineColor(2);
  fHistClusterMapModuleITSSAbadInAcc->Draw("same");
  l2->AddEntry(fHistClusterMapModuleITSSAbadInAcc,"bad","l");
  fHistClusterMapModuleITSSAnoclsInAcc->SetLineColor(3);
  fHistClusterMapModuleITSSAnoclsInAcc->Draw("same");
  l2->AddEntry(fHistClusterMapModuleITSSAnoclsInAcc,"no cls","l");
  l2->Draw();


  TCanvas *c4 =new TCanvas("c4","c4",10,10,500,500);
  c4->SetGridy();
  fHistPhiITSMIokbadoutinz6InAcc->Divide(fHistPhiTPCInAcc);
  fHistPhiITSMIokbadoutinz6InAcc->SetMinimum(0);
  fHistPhiITSMIokbadoutinz6InAcc->SetMaximum(1.5);
  fHistPhiITSMIokbadoutinz6InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPhiITSMIokbadoutinz6InAcc->SetTitle("Fraction of tracks with 6 layers ok");
  fHistPhiITSMIokbadoutinz6InAcc->Draw();

  TLegend *l3=new TLegend(0.5,0.5,0.9,0.9);
  TLegend *l4=new TLegend(0.5,0.5,0.9,0.9);
  TCanvas *c5 =new TCanvas("c5","c5",10,10,1200,600);
  c5->Divide(2,1);
  c5_1->SetGridy();
  c5_2->SetGridy();
  c5->cd(1);
  TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);
  fHistPtITSMIge2InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMIge2InAcc->SetMaximum(1.5);
  fHistPtITSMIge2InAcc->SetMinimum(0);
  fHistPtITSMIge2InAcc->SetTitle("Fraction of prolonged tracks with N ITS points");
  fHistPtITSMIge2InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIge2InAcc->Draw();
  l3->AddEntry(fHistPtITSMIge2InAcc,">=2 cls","l");
  fHistPtITSMI6InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMI6InAcc->SetLineColor(2);
  l3->AddEntry(fHistPtITSMI6InAcc,"6 cls","l");
  fHistPtITSMI6InAcc->Draw("same");
  fHistPtITSMI5InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMI5InAcc->SetLineColor(3);
  l3->AddEntry(fHistPtITSMI5InAcc,"5 cls","l");
  fHistPtITSMI5InAcc->Draw("same");
  fHistPtITSMI4InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMI4InAcc->SetLineColor(4);
  l3->AddEntry(fHistPtITSMI4InAcc,"4 cls","l");
  fHistPtITSMI4InAcc->Draw("same");
  fHistPtITSMI3InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMI3InAcc->SetLineColor(6);
  l3->AddEntry(fHistPtITSMI3InAcc,"3 cls","l");
  fHistPtITSMI3InAcc->Draw("same");
  fHistPtITSMI2InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMI2InAcc->SetLineColor(7);
  l3->AddEntry(fHistPtITSMI2InAcc,"2 cls","l");
  fHistPtITSMI2InAcc->Draw("same");
  fHistPtITSMISPDInAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMISPDInAcc->SetLineColor(9);
  l3->AddEntry(fHistPtITSMISPDInAcc,"2SPD + any","l");
  fHistPtITSMISPDInAcc->Draw("same");
  fHistPtITSMIge2InAcc->Draw("same");
  l3->Draw();
  c5->cd(2);
  TH1F *fHistPtITSMIokbadoutinzge4InAcc = (TH1F*)fHistPtITSMIokbadoutinz6InAcc->Clone("fHistPtITSMIokbadoutinzge4InAcc");
  fHistPtITSMIokbadoutinzge4InAcc->Add(fHistPtITSMIokbadoutinz5InAcc);
  fHistPtITSMIokbadoutinzge4InAcc->Add(fHistPtITSMIokbadoutinz4InAcc);
  fHistPtITSMIokbadoutinzge4InAcc->SetMaximum(1.5);
  fHistPtITSMIokbadoutinzge4InAcc->SetMinimum(0);
  fHistPtITSMIokbadoutinzge4InAcc->SetTitle("Fraction of prolonged tracks with N ITS layers \"ok\"");
  fHistPtITSMIokbadoutinzge4InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIokbadoutinzge4InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMIokbadoutinzge4InAcc->SetLineColor(1);
  fHistPtITSMIokbadoutinzge4InAcc->Draw();
  fHistPtITSMIokbadoutinz6InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMIokbadoutinz6InAcc->SetLineColor(2);
  fHistPtITSMIokbadoutinz6InAcc->Draw("same");
  fHistPtITSMIokbadoutinz5InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMIokbadoutinz5InAcc->SetLineColor(3);
  fHistPtITSMIokbadoutinz5InAcc->Draw("same");
  fHistPtITSMIokbadoutinz4InAcc->Divide(fHistPtTPCInAcc);
  fHistPtITSMIokbadoutinz4InAcc->SetLineColor(4);
  fHistPtITSMIokbadoutinz4InAcc->Draw("same");
  fHistPtITSMIokbadoutinzge4InAcc->Draw("same");
  l4->AddEntry(fHistPtITSMIokbadoutinzge4InAcc,">=4 layers","l");
  l4->AddEntry(fHistPtITSMIokbadoutinz6InAcc,"6 layers","l");
  l4->AddEntry(fHistPtITSMIokbadoutinz5InAcc,"5 layers","l");
  l4->AddEntry(fHistPtITSMIokbadoutinz4InAcc,"4 layers","l");
  l4->Draw();


  // PLOT ALIGNMENT CHECKS
  //
  TH1F *hSPDTrackletsTBdxy = new TH1F("hSPDTrackletsTBdxy","SPD tracklets; SPD tracklet to SPD vertex distance in (x,y) [cm]; entries",100,-0.1,0.1);
  TH1F *hSPDTrackletsLRdxy = new TH1F("hSPDTrackletsLRdxy","SPD tracklets; SPD tracklet to SPD vertex distance in (x,y) [cm]; entries",100,-0.1,0.1);
  TH1F *hSPDTrackletsTBdz = new TH1F("hSPDTrackletsTBdz","SPD tracklets; SPD tracklet to SPD vertex distance in z [cm]; entries",100,-0.1,0.1);
  TH1F *hSPDTrackletsLRdz = new TH1F("hSPDTrackletsLRdz","SPD tracklets; SPD tracklet to SPD vertex distance in z [cm]; entries",100,-0.1,0.1);

  TNtuple *fNtupleITSAlignSPDTracklets = (TNtuple*)list->FindObject("fNtupleITSAlignSPDTracklets");
  Float_t dxy,dz,phi,pt;
  fNtupleITSAlignSPDTracklets->SetBranchAddress("pt",&pt);
  fNtupleITSAlignSPDTracklets->SetBranchAddress("phi",&phi);
  fNtupleITSAlignSPDTracklets->SetBranchAddress("dxy",&dxy);
  fNtupleITSAlignSPDTracklets->SetBranchAddress("dz",&dz);

  for(Int_t i=0;i<fNtupleITSAlignSPDTracklets->GetEntries();i++) {
    fNtupleITSAlignSPDTracklets->GetEvent(i);
    if(pt<1.3) continue;
    if(TMath::Abs(TMath::Abs(phi)-0.5*TMath::Pi())<0.25*TMath::Pi()) {
      // top-bottom
      hSPDTrackletsTBdxy->Fill(dxy);
      hSPDTrackletsTBdz->Fill(dz);
    } else {
      // left-right
      hSPDTrackletsLRdxy->Fill(dxy);
      hSPDTrackletsLRdz->Fill(dz);
    }
  }

  TLegend *l6=new TLegend(0.5,0.5,0.9,0.9);

  TCanvas *c6 = new TCanvas("c6","c6",0,0,1000,500);
  c6->Divide(2,1);
  c6->cd(1);
  hSPDTrackletsTBdxy->SetLineColor(4);
  hSPDTrackletsTBdxy->Draw();
  hSPDTrackletsLRdxy->SetLineColor(2);
  hSPDTrackletsLRdxy->Draw("same");
  l6->AddEntry(hSPDTrackletsTBdxy,"top-bottom","l");
  l6->AddEntry(hSPDTrackletsLRdxy,"left-right","l");
  l6->Draw();
  c6->cd(2);
  hSPDTrackletsTBdz->SetLineColor(4);
  hSPDTrackletsTBdz->Draw();
  hSPDTrackletsLRdz->SetLineColor(2);
  hSPDTrackletsLRdz->Draw("same");
  l6->Draw();


  TH1F *hSPDExtraClsTBdxy = new TH1F("hSPDExtraClsTBdxy","SPD extra clusters; track-to-point distance in (x,y) [cm]; entries",100,-0.1,0.1);
  TH1F *hSPDExtraClsLRdxy = new TH1F("hSPDExtraClsLRdxy","SPD extra clusters; track-to-point distance in (x,y) [cm]; entries",100,-0.1,0.1);
  TH1F *hSPDExtraClsTBdz = new TH1F("hSPDExtraClsTBdz","SPD extra clusters; track-to-point distance in z [cm]; entries",100,-0.1,0.1);
  TH1F *hSPDExtraClsLRdz = new TH1F("hSPDExtraClsLRdz","SPD extra clusters; track-to-point distance in z [cm]; entries",100,-0.1,0.1);

  TNtuple *fNtupleITSAlignExtra = (TNtuple*)list->FindObject("fNtupleITSAlignExtra");
  Float_t layer,npoints,x,y,z;
  fNtupleITSAlignExtra->SetBranchAddress("layer",&layer);
  fNtupleITSAlignExtra->SetBranchAddress("npoints",&npoints);
  fNtupleITSAlignExtra->SetBranchAddress("x",&x);
  fNtupleITSAlignExtra->SetBranchAddress("y",&y);
  fNtupleITSAlignExtra->SetBranchAddress("z",&z);
  fNtupleITSAlignExtra->SetBranchAddress("dxy",&dxy);
  fNtupleITSAlignExtra->SetBranchAddress("dz",&dz);
  fNtupleITSAlignExtra->SetBranchAddress("pt",&pt);

  for(Int_t i=0;i<fNtupleITSAlignExtra->GetEntries();i++) {
    fNtupleITSAlignExtra->GetEvent(i);
    if(pt<0.5) continue;
    if(layer!=0 && layer!=1) continue;
    if(npoints<4) continue;
    phi = TMath::ATan2(y,x);
    if(TMath::Abs(TMath::Abs(phi)-0.5*TMath::Pi())<0.25*TMath::Pi()) {
      // top-bottom
      hSPDExtraClsTBdxy->Fill(dxy);
      hSPDExtraClsTBdz->Fill(dz);
    } else {
      // left-right
      hSPDExtraClsLRdxy->Fill(dxy);
      hSPDExtraClsLRdz->Fill(dz);
    }
  }

  TCanvas *c7 = new TCanvas("c7","c7",0,0,1000,500);
  c7->Divide(2,1);
  c7->cd(1);
  hSPDExtraClsTBdxy->SetLineColor(4);
  hSPDExtraClsTBdxy->Draw();
  hSPDExtraClsLRdxy->SetLineColor(2);
  hSPDExtraClsLRdxy->Draw("same");
  l6->Draw();
  c7->cd(2);
  hSPDExtraClsTBdz->SetLineColor(4);
  hSPDExtraClsTBdz->Draw();
  hSPDExtraClsLRdz->SetLineColor(2);
  hSPDExtraClsLRdz->Draw("same");
  l6->Draw();


  return;
}
