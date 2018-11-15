Int_t PlotTimestamp(TCanvas* c1)
{
  double rightlegend = 1 - 10./c1->GetWw();
  //the function plots a timestamp, the used Aliroot version, and the number of runs
  TString sTimestamp  = gSystem->GetFromPipe("date");
  TString sAlirootVer;
  if (gSystem->GetFromPipe("echo $ALICE_VER") == "master"){
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("wdir=`pwd`; cd $ALICE_ROOT/../src; git describe; cd $wdir;");
  }
  else {
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("echo $ALIROOT_VERSION");
  }

  TString sAliphysicsVer;
  if (gSystem->GetFromPipe("echo $ALIPHYSICS_VER") == "master"){
    sAliphysicsVer = "AliPhysics: " + gSystem->GetFromPipe("wdir=`pwd`; cd $ALICE_PHYSICS/../src; git describe; cd $wdir;");
  }
  else {
    sAliphysicsVer = "AliPhysics: " + gSystem->GetFromPipe("echo $ALIPHYSICS_VERSION");
  }
  TLatex* latTime = new TLatex(rightlegend,0.99,sTimestamp.Data());
  latTime->SetTextSize(0.03);
  latTime->SetTextAlign(31);
  latTime->SetNDC();
  latTime->Draw("same");
  TLatex* latAliroot = new TLatex(rightlegend,0.95,sAlirootVer.Data());
  latAliroot->SetTextSize(0.03);
  latAliroot->SetTextAlign(31);
  latAliroot->SetNDC();
  latAliroot->Draw("same");
  TLatex* latAliphysics = new TLatex(rightlegend,0.91,sAliphysicsVer.Data());
  latAliphysics->SetTextSize(0.03);
  latAliphysics->SetTextAlign(31);
  latAliphysics->SetNDC();
  latAliphysics->Draw("same");
return 1;
}


int drawPerformanceTPCQAMatch(const char* inFile = "perf.root") {
  //
  // Draw control histograms
  // and generate output pictures
  //

  gSystem->Load("libSTAT");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTPCcalib");
  gSystem->Load("libTRDcalib");
  gSystem->Load("libT0calib");
  gSystem->Load("libTOFcalib");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.025);
  TH1::AddDirectory(kFALSE);

  //
  // set criteria
  //
  Float_t mineta = -0.8;
  Float_t maxeta = 0.799;
  Float_t minNclust = 70.0;
  Double_t ptMax = 10;
  //
  // open input file
  //
  TFile *file = TFile::Open(inFile);
  if(!file)
    return -9;
  cout<<"QA file opened"<<endl;
  file->cd();
  // get the TPC list
  if(gROOT->FindObject("TPC_PerformanceQA")) TPC_PerformanceQA->cd();
  cout<<"TPC_PerformanceQA opened"<<endl;
  TList *TPC = (TList*)gROOT->FindObject("TPCQA");
  if(TPC==NULL) TPC = (TList*)gROOT->FindObject("TPCQA_v0_c0");
  if(TPC==NULL) TPC = (TList*)gROOT->FindObject("TPCQA_v0_c30");
  if(TPC==NULL) TPC = (TList*)gROOT->FindObject("TPCQA_v0_c70");
  if(TPC==NULL) TPC = (TList*)gROOT->FindObject("TPC");
  if(TPC==NULL) return(0);
  cout<<"TPCQA list found"<<endl;
  // get TPC performance object
  AliPerformanceTPC *obj = (AliPerformanceTPC*)TPC->FindObject("AliPerformanceTPC");
  if(obj==NULL) return(0);
  cout<<"what about here after obj  "<<endl;
  // get folder with histograms
  TFolder *fold = obj->GetAnalysisFolder();
  if(!fold) return(0);
  cout<<"what about here after folder  "<<endl;
  //
  // get the HLT list
  //	file->cd();
  //	if(gROOT->FindObject("HLT_PerformanceQA")) HLT_PerformanceQA->cd();
  //	TList *HLT = (TList*)gROOT->FindObject("HLTQA");


  //
  // Draw histograms
  //

  //
  // event level
  //

  TH1 *h1D = 0;
  TH2 *h2D = 0;
  TH3 *h3D = 0;

  // default gaus fit. Used for setting proper fit ranges
  TF1* defaultFit = 0x0;

  h1D = (TH1*)fold->FindObject("h_tpc_event_1");
  Double_t NEvents = h1D->GetEntries();

  cout<<"number of events    "<<NEvents<<endl;

  //
  // ===| Event info |==========================================================
  //
  {
  TCanvas *can1 = new TCanvas("can1","TPC event information",1200,800);
  can1->Divide(3,2);

  can1->cd(1);
  fold->FindObject("h_tpc_event_6")->Draw("histe");

  can1->cd(2);
  gPad->SetLogy();
  h1D = (TH1*)fold->FindObject("h_tpc_event_recvertex_0");
  h1D->GetXaxis()->SetRangeUser(-1.,1.);
  h1D->Draw("histe");

  can1->cd(3);
  gPad->SetLogy();
  h1D = (TH1*)fold->FindObject("h_tpc_event_recvertex_1");
  h1D->GetXaxis()->SetRangeUser(-1.,1.);
  h1D->Draw("histe");
  PlotTimestamp(can1);

  can1->cd(4);
  gPad->SetLogy();
  fold->FindObject("h_tpc_event_recvertex_2")->Draw("histe");

  can1->cd(5);
  gPad->SetLogy();
  TH1 *hp = fold->FindObject("h_tpc_event_recvertex_3");
  hp->SetTitle("Track.Multi., ncl>70, |dcar|<3 cm, |dcaz|<3 cm");
  hp->Draw("histe");

  can1->cd(6);
  gPad->SetLogy();
  hp = (TH1*)fold->FindObject("h_tpc_event_recvertex_4");
  hp->SetTitle("Pos/neg(red) Track.Multi. ncl>70, |dcar|<3 cm, |dcaz|<3 cm");
  hp->Draw("histe");
  TH1* he = fold->FindObject("h_tpc_event_recvertex_5");
  he->SetLineColor(kRed);
  he->Draw("histesame");

  can1->SaveAs("TPC_event_info.png");
  }

  //
  // ===| Track distribution eta phi pt |=======================================
  //
  {
  TCanvas *can2 = new TCanvas("can2","#eta , #phi and p_{t}",1200,800);

  can2->Divide(3,2);

  can2->cd(1);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_0_5_6");
  h3D->GetXaxis()->SetRangeUser(minNclust,160);
  h3D->Project3D("yz")->Draw("colz");
  h3D->Project3D("yz")->SetTitle("#eta vs #phi, positive tracks");
  can2->Update();

  can2->cd(4);
  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_0_5_6");
  h3D->GetXaxis()->SetRangeUser(minNclust,160);
  h3D->Project3D("yz")->Draw("colz");
  h3D->Project3D("yz")->SetTitle("#eta vs #phi, negative tracks");
  can2->Update();

  can2->cd(2);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_0_5_6");
  h3D->GetXaxis()->SetRangeUser(minNclust,160);
  h1D = h3D->Project3D("y");
  h1D->SetTitle("#eta of pos/neg(red) charged tracks, ncl>70");
  h1D->Draw("histe");

  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_0_5_6");
  h3D->GetXaxis()->SetRangeUser(minNclust,160);
  h1D = h3D->Project3D("y");
  h1D->SetLineColor(kRed);
  h1D->Draw("histesame");

  can2->cd(5);
  gPad->SetLogy();
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_0_5_7");
  h3D->GetXaxis()->SetRangeUser(minNclust,160);
  h3D->GetYaxis()->SetRangeUser(mineta,maxeta);
  h1D = h3D->Project3D("z");
  h1D->Scale(1,"width");
  h1D->SetTitle("p_{T} of pos/neg(red) charged tracks, ncl>70, |#eta|<0.8");
  h1D->Draw("histe");

  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_0_5_7");
  h3D->GetXaxis()->SetRangeUser(minNclust,160);
  h3D->GetYaxis()->SetRangeUser(mineta,maxeta);
  h1D = h3D->Project3D("z");
  h1D->Scale(1,"width");
  h1D->SetLineColor(kRed);
  h1D->Draw("histesame");

  can2->cd(3);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_5_6_7");
  h3D->GetZaxis()->SetRangeUser(2,20);
  h3D->GetXaxis()->SetRangeUser(mineta,-0.00001);
  TH1 *h1D = h3D->Project3D("y");
  h1D->SetTitle("#phi of pos/neg(red) charged tracks, pt>1GeV/c, -0.8<#eta<0.0");
  h1D->Draw("histe");

  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_5_6_7");
  h3D->GetZaxis()->SetRangeUser(2,20);
  h3D->GetXaxis()->SetRangeUser(mineta,-0.00001);
  h1D = h3D->Project3D("y");
  h1D->SetLineColor(kRed);
  h1D->Draw("histesame");
  PlotTimestamp(can2);

  can2->cd(6);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_5_6_7");
  TH3 *h3D1 = h3D->Clone("h3D1");
  h3D1->GetXaxis()->SetRangeUser(0.0,maxeta);
  h3D1->GetZaxis()->SetRangeUser(2,20);
  h1D = h3D1->Project3D("y");
  h1D->SetTitle("#phi of pos/neg(red) charged tracks, pt>1GeV/c, 0.0<#eta<0.8");
  h1D->Draw("histe");

  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_5_6_7");
  TH3 *h3D2 = h3D->Clone("h3D2");
  h3D2->GetXaxis()->SetRangeUser(0.0,maxeta);
  h3D2->GetZaxis()->SetRangeUser(2,20);
  h1D = h3D2->Project3D("y");
  h1D->SetLineColor(kRed);
  h1D->Draw("histesame");

  can2->SaveAs("eta_phi_pt.png");
  }

  //
  // ===| Cluster Occupancy |===================================================
  //
  {
  TCanvas *can3 = new TCanvas("can3","Cluster Occupancy",700,700);
  can3->Divide(1,2);

  can3->cd(1);
  TH3 *h3D_1 = (TH3*)fold->FindObject("h_tpc_clust_0_1_2");
  TH3 *h3D_2 = (TH3*) h3D_1->Clone("h3D_2");
  h3D_1->GetZaxis()->SetRangeUser(0,0.99);
  h3D_1->Project3D("xy")->Draw("colz");
  h3D_1->Project3D("xy")->SetTitle("Cluster Occupancy A Side");
  if(NEvents > 0)
    h3D_1->Project3D("xy")->Scale(1.0/NEvents);
  can3->Update();
  PlotTimestamp(can3);

  can3->cd(2);
  h3D_2->GetZaxis()->SetRangeUser(1,2) ;
  h3D_2->Project3D("xy")->Draw("colz");
  h3D_2->Project3D("xy")->SetTitle("Cluster Occupancy C Side");
  if(NEvents>0)
    h3D_2->Project3D("xy")->Scale(1.0/NEvents);

  can3->SaveAs("cluster_occupancy.png");
  }

  //
  // ===| Cluster in detail |===================================================
  //
  {
  TCanvas *can4 = new TCanvas("can4","Clusters in Detail",1200,800);
  can4->Divide(3,2);

  // set proper fit range for # cluster plots
  //   the fits go either up to 160 clusters or as fraction up to 1.1
  //   a range of 200 should be more than fine
  defaultFit = (TF1*)gROOT->GetFunction("gaus");
  if (!defaultFit) defaultFit = new TF1("gaus","gaus");
  defaultFit->SetParameters(1,1,0.01);
  defaultFit->SetParLimits(1, 0.01, 200);
  defaultFit->SetParLimits(2, 0.001, 100);

  TObjArray *arr1 = new TObjArray();
  TObjArray *arr2 = new TObjArray();
  TObjArray *arr3 = new TObjArray();
  TObjArray *arr4 = new TObjArray();
  TObjArray *arr5 = new TObjArray();
  TObjArray *arr6 = new TObjArray();

  can4->cd(1);
  h3D = (TH3*)fold->FindObject("h_tpc_track_all_recvertex_0_5_7");
  h2D = (TH2*)h3D->Project3D("xy");
  h2D->SetTitle("nCluster vs #eta, |dcar|<3 cm, |dcaz|<3 cm");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr1);
  h2D->Draw("colz");
  arr1->At(1)->Draw("same");

  can4->cd(4);
  h3D = (TH3*)fold->FindObject("h_tpc_track_all_recvertex_2_5_7");
  h2D = (TH2*)h3D->Project3D("xy");
  h2D->SetTitle("Findable clusters, |dcar|<3 cm, |dcaz|<3 cm");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr2);
  h2D->Draw("colz");
  arr2->At(1)->Draw("same");

  can4->cd(2);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_0_5_6");
  TH3 *h3D11 = h3D->Clone("h3D11");
  h3D11->Add(((TH3*)fold->FindObject("h_tpc_track_neg_recvertex_0_5_6")),1);
  h3D11->GetYaxis()->SetRangeUser(mineta,-0.00001);
  h2D = (TH2*)h3D11->Project3D("xz");
  h2D->SetTitle("nCluster vs #phi, -0.8<#eta<0.0, |dcar|<3 cm, |dcaz|<3 cm");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr3);
  h2D->Draw("colz");
  arr3->At(1)->Draw("same");

  can4->cd(5);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_0_5_6");
  TH3 *h3D22 = h3D->Clone("h3D22");
  h3D22->Add(((TH3*)fold->FindObject("h_tpc_track_neg_recvertex_0_5_6")),1);
  h3D22->GetYaxis()->SetRangeUser(0.0,maxeta);
  h2D->SetTitle("nCluster vs #phi, 0.0<#eta<0.8, |dcar|<3 cm, |dcaz|<3 cm");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr4);
  h2D->Draw("colz");
  arr4->At(1)->Draw("same");

  can4->cd(3);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_2_5_6");
  TH3 *h3D33 = h3D->Clone("h3D33");
  h3D33->Add(((TH3*)fold->FindObject("h_tpc_track_neg_recvertex_2_5_6")),1);
  h3D33->GetYaxis()->SetRangeUser(mineta,-0.00001);
  h2D = (TH2*)h3D33->Project3D("xz");
  h2D->SetTitle("Findable clusters vs #phi, -0.8<#eta<0.0, |dcar|<3 cm, |dcaz|<3 cm");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr5);
  h2D->Draw("colz");
  arr5->At(1)->Draw("same");
  PlotTimestamp(can4);


  can4->cd(6);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_2_5_6");
  TH3 *h3D44 = h3D->Clone("h3D44");
  h3D44->Add(((TH3*)fold->FindObject("h_tpc_track_neg_recvertex_2_5_6")),1);
  h3D44->GetYaxis()->SetRangeUser(0.0,maxeta);
  h2D = (TH2*)h3D44->Project3D("xz");
  h2D->SetTitle("Findalbe clusters vs #phi, 0.0<#eta<0.8, |dcar|<3 cm, |dcaz|<3 cm");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr6);
  h2D->Draw("colz");
  arr6->At(1)->Draw("same");

  can4->Update();
  can4->SaveAs("cluster_in_detail.png");
  }

  //
  // ===| DCA in detail |=======================================================
  //
  {
  TCanvas *can5 = new TCanvas("can5","DCA In Detail",1200,800);
  can5->Divide(3,2);

  // set proper fit range for # cluster plots
  //   the range of the plot is -3 to 3
  defaultFit = (TF1*)gROOT->GetFunction("gaus");
  if (!defaultFit) defaultFit = new TF1("gaus","gaus");
  defaultFit->SetParameters(1,0,0.1);
  defaultFit->SetParLimits(1, -3, 3);
  defaultFit->SetParLimits(2, 0.001, 1);

  TObjArray *arr7 = new TObjArray();
  TObjArray *arr8 = new TObjArray();

  can5->cd(1);
  h3D = (TH3*)fold->FindObject("h_tpc_track_all_recvertex_3_5_7");
  h3D->GetYaxis()->SetRangeUser(-1,1);
  h2D = (TH2*)h3D->Project3D("xy");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr7);
  h2D->SetTitle("DCAR vs #eta");
  h2D->Draw("colz");
  arr7->At(1)->Draw("same");

  can5->cd(2);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_3_5_7");
  h3D->GetYaxis()->SetRangeUser(-1,1);
  h3D->Project3D("xy")->Draw("colz");
  h3D->Project3D("xy")->SetTitle("DCAR vs #eta of pos. charged tracks");

  can5->cd(3);
  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_3_5_7");
  h3D->GetYaxis()->SetRangeUser(-1,1);
  h3D->Project3D("xy")->Draw("colz");
  h3D->Project3D("xy")->SetTitle("DCAR vs #eta of neg. charged tracks");
  PlotTimestamp(can5);

  can5->cd(4);
  h3D = (TH3*)fold->FindObject("h_tpc_track_all_recvertex_4_5_7");
  h3D->GetYaxis()->SetRangeUser(-1,1);
  h2D = (TH2*)h3D->Project3D("xy");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr8);
  h2D->SetTitle("DCAZ vs #eta");
  h2D->Draw("colz");
  arr8->At(1)->Draw("same");

  can5->cd(5);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_4_5_7");
  h3D->GetYaxis()->SetRangeUser(-1,1);
  h3D->Project3D("xy")->Draw("colz");
  h3D->Project3D("xy")->SetTitle("DCAZ vs #eta of pos. charged tracks");

  can5->cd(6);
  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_4_5_7");
  h3D->GetYaxis()->SetRangeUser(-1,1);
  h3D->Project3D("xy")->Draw("colz");
  h3D->Project3D("xy")->SetTitle("DCAZ vs #eta of neg. charged tracks");

  can5->SaveAs("dca_in_detail.png");
  }

  //
  // ===| DCAr vs pT |==========================================================
  //
  {
  TCanvas *can51 = new TCanvas("can51","DCAr versus pT",700,800);
  can51->Divide(2,2);

  // set proper fit range for # cluster plots
  //   the range of the plot is -3 to 3
  defaultFit = (TF1*)gROOT->GetFunction("gaus");
  if (!defaultFit) defaultFit = new TF1("gaus","gaus");
  defaultFit->SetParameters(1,0,0.1);
  defaultFit->SetParLimits(1, -3, 3);
  defaultFit->SetParLimits(2, 0.001, 3);

  TObjArray *arr9 = new TObjArray();
  TObjArray *arr10 = new TObjArray();
  TObjArray *arr11 = new TObjArray();
  TObjArray *arr12 = new TObjArray();

  can51->cd(1);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_3_5_7");
  TH3 *h3Dp = h3D->Clone("h3Dp");
  h3D->SetAxisRange(0.25,ptMax,"Z");
  h3D->GetYaxis()->SetRangeUser(0.0,maxeta);
  h2D  = (TH2*)h3D->Project3D("xz");
  h2D->Draw("colz");
  h2D->SetTitle("DCAR vs pT of pos. charged tracks(A Side)");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr9);
  TH1 *width1 = (TH1*)arr9->At(2);
  width1->Draw("same");
  width1->SetMarkerColor(2);
  width1->SetLineColor(2);
  arr9->At(1)->Draw("same");

  can51->cd(2);
  //h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_3_5_7");
  h3Dp->SetAxisRange(0.25,ptMax,"Z");
  h3Dp->GetYaxis()->SetRangeUser(mineta,-0.00001);
  h2D  = (TH2*)h3Dp->Project3D("xz");
  h2D->Draw("colz");
  h2D->SetTitle("DCAR vs pT of pos. charged tracks(C Side)");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr10);
  TH1 *width2 = (TH1*)arr10->At(2);
  width2->Draw("same");
  width2->SetMarkerColor(2);
  width2->SetLineColor(2);
  arr10->At(1)->Draw("same");
  PlotTimestamp(can51);

  can51->cd(3);
  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_3_5_7");
  TH3 *h3Dn = h3D->Clone("h3Dn");
  h3D->SetAxisRange(0.25,ptMax,"Z");
  h3D->GetYaxis()->SetRangeUser(0.0,maxeta);
  h2D  = (TH2*)h3D->Project3D("xz");
  h2D->Draw("colz");
  h2D->SetTitle("DCAR vs pT of neg. charged tracks(A Side)");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr11);
  TH1 *width3 = (TH1*)arr11->At(2);
  width3->SetMarkerColor(2);
  width3->SetLineColor(2);
  width3->Draw("same");
  arr11->At(1)->Draw("same");

  can51->cd(4);
  //h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_3_5_7");
  h3Dn->SetAxisRange(0.25,ptMax,"Z");
  h3Dn->GetYaxis()->SetRangeUser(mineta,-0.00001);
  h2D  = (TH2*)h3Dn->Project3D("xz");
  h2D->Draw("colz");
  h2D->SetTitle("DCAR vs pT of neg. charged tracks(C Side)");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr12);
  TH1 *width4 = (TH1*)arr12->At(2);
  width4->SetMarkerColor(2);
  width4->SetLineColor(2);
  width4->Draw("same");
  arr12->At(1)->Draw("same");

  can51->SaveAs("dcar_pT.png");
  }

  // ===========================================================================
  // ===| get TPC dEdx performance object |=====================================
  // ===========================================================================
  AliPerformanceDEdx *obj1 = TPC->FindObject("AliPerformanceDEdxTPCInner");
  if(obj1==NULL) return(0);

  // get folder with histograms
  TFolder *fold1 = obj1->GetAnalysisFolder();
  if(!fold1) return(0);

  //
  // ===| dE/dx plot |==========================================================
  //
  {
  TCanvas *can6 = new TCanvas("can6","TPC dEdX",1200,800);
  can6->Divide(3,2);

  can6->cd(1);
  gPad->SetLogz();
  fold1->FindObject("h_tpc_dedx_0_1")->Draw("colz");

  can6->cd(2);
  gPad->SetLogz();
  fold1->FindObject("h_tpc_dedx_0_5")->Draw("colz");

  can6->cd(3);
  gPad->SetLogz();
  fold1->FindObject("h_tpc_dedx_0_6")->Draw("colz");
  PlotTimestamp(can6);

  can6->cd(4);
  gPad->SetLogx();
  gPad->SetLogz();
  TH2 *h2 = fold1->FindObject("h_tpc_dedx_0_7");
  h2->GetXaxis()->SetRangeUser(0.1,10);
  h2->Draw("colz");
  ////////////////////////////////////////////////////////////////////
  can6->cd(5);
  gPad->SetLogz();
  //fold1->FindObject("h_tpc_dedx_mips_a_0_1")->Draw("colz");
  TH2 *htest = fold1->FindObject("h_tpc_dedx_mips_a_0_1");
  htest->GetYaxis()->SetRangeUser(30,60);
  htest->Draw("colz");
  can6->cd(6);
  gPad->SetLogz();
  //fold1->FindObject("h_tpc_dedx_mips_c_0_1")->Draw("colz");
  TH2 *htest1 = fold1->FindObject("h_tpc_dedx_mips_c_0_1");
  htest1->GetYaxis()->SetRangeUser(30,60);
  htest1->Draw("colz");

  /////////////////////////////////////////////////////////////////////
  can6->SaveAs("TPC_dEdx_track_info.png");
  }

  //
  // ===| DCA vs phi |==========================================================
  //
  {
  TCanvas *can7 = new TCanvas("can7","DCA vs #phi",1200,800);
  can7->Divide(4,2);

  can7->cd(1);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_3_5_6");
  TH3 *h3D71 = h3D->Clone("h3D71");
  h3D->GetYaxis()->SetRangeUser(0.0,maxeta);
  h3D->Project3D("xz")->Draw("colz");
  h3D->Project3D("xz")->SetTitle("DCAR vs #phi of pos. charged tracks(A)");

  can7->cd(2);
  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_3_5_6");
  TH3 *h3D72 = h3D->Clone("h3D72");
  h3D->GetYaxis()->SetRangeUser(0.0,maxeta);
  h3D->Project3D("xz")->Draw("colz");
  h3D->Project3D("xz")->SetTitle("DCAR vs #phi of neg. charged tracks(A)");

  can7->cd(3);
  h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_4_5_6");
  TH3 *h3D73 = h3D->Clone("h3D73");
  h3D->GetYaxis()->SetRangeUser(0.0,maxeta);
  h3D->Project3D("xz")->Draw("colz");
  h3D->Project3D("xz")->SetTitle("DCAZ vs #phi of pos. charged tracks(A)");

  can7->cd(4);
  h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_4_5_6");
  TH3 *h3D74 = h3D->Clone("h3D74");
  h3D->GetYaxis()->SetRangeUser(0.0,maxeta);
  h3D->Project3D("xz")->Draw("colz");
  h3D->Project3D("xz")->SetTitle("DCAZ vs #phi of neg. charged tracks(A)");
  PlotTimestamp(can7);

  can7->cd(5);
  //h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_3_5_6");
  h3D71->GetYaxis()->SetRangeUser(mineta,-0.00001);
  h3D71->Project3D("xz")->Draw("colz");
  h3D71->Project3D("xz")->SetTitle("DCAR vs #phi of pos. charged tracks(C)");

  can7->cd(6);
  //h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_3_5_6");
  h3D72->GetYaxis()->SetRangeUser(mineta,-0.00001);
  h3D72->Project3D("xz")->Draw("colz");
  h3D72->Project3D("xz")->SetTitle("DCAR vs #phi of neg. charged tracks(C)");

  can7->cd(7);
  //h3D = (TH3*)fold->FindObject("h_tpc_track_pos_recvertex_4_5_6");
  h3D73->GetYaxis()->SetRangeUser(mineta,-0.00001);
  h3D73->Project3D("xz")->Draw("colz");
  h3D73->Project3D("xz")->SetTitle("DCAZ vs #phi of pos. charged tracks(C)");

  can7->cd(8);
  //h3D = (TH3*)fold->FindObject("h_tpc_track_neg_recvertex_4_5_6");
  h3D74->GetYaxis()->SetRangeUser(mineta,-0.00001);
  h3D74->Project3D("xz")->Draw("colz");
  h3D74->Project3D("xz")->SetTitle("DCAZ vs #phi of neg. charged tracks(C)");


  can7->SaveAs("dca_and_phi.png");
  }


  // ===========================================================================
  // ===| Get performance objects ==============================================
  // ===========================================================================
  //
  AliPerformanceMatch *obj2 = (AliPerformanceMatch*)TPC->FindObject("AliPerformanceMatchTPCITS");
  TFolder *pMatch = obj2->GetAnalysisFolder();

  AliPerformanceMatch *obj3 = (AliPerformanceMatch*)TPC->FindObject("AliPerformanceMatchITSTPC");
  TFolder *pPull = obj3->GetAnalysisFolder();

  //
  // event level
  //

  TH1 *h1D = 0;
  TH1 *h1D1 = 0;
  TH1 *h1D2 = 0;
  TH1 *h1D3 = 0;
  TH2 *h2D = 0;
  TH2 *h2D1 = 0;

  //
  // ===| TPC-ITS matching efficiency |=========================================
  //
  {
  TCanvas *can8 = new TCanvas("can8","TPC-ITS Matching Efficiency",800,800);
  can8->Divide(2,2);

  can8->cd(1);
  h2D = (TH2*)(pMatch->FindObject("h_tpc_match_trackingeff_all_2_3"));
  h2D1 = (TH2*)(pMatch->FindObject("h_tpc_match_trackingeff_tpc_2_3"));
  TH2 *h2D2 = h2D->Clone("h2D2");
  TH2 *h2D3 = h2D1->Clone("h2D3");

  h2D->GetXaxis()->SetRangeUser(0,1.5);
  h2D1->GetXaxis()->SetRangeUser(0,1.5);
  h1D = h2D->ProjectionY();
  h1D1 = h2D1->ProjectionY();
  h1D1->Divide(h1D);
  h1D1->GetYaxis()->SetRangeUser(0,1.05);
  h1D1->SetTitle("TPC-ITS Matching Efficiency (A)");
  h1D1->Draw("e0");

  can8->cd(2);
  h2D2->GetXaxis()->SetRangeUser(-1.5,0);
  h2D3->GetXaxis()->SetRangeUser(-1.5,0);
  h1D2 = h2D2->ProjectionY();
  h1D3 = h2D3->ProjectionY();
  h1D3->Divide(h1D2);
  h1D3->SetLineColor(2);
  h1D3->GetYaxis()->SetRangeUser(0,1.05);
  h1D3->SetTitle("TPC-ITS Matching Efficiency (C)");
  h1D3->Draw("e0");
  PlotTimestamp(can8);

  can8->cd(3);
  h2D = (TH2*)(pMatch->FindObject("h_tpc_match_trackingeff_all_1_3"));
  h2D1 = (TH2*)(pMatch->FindObject("h_tpc_match_trackingeff_tpc_1_3"));
  TH2 *h2D4 = h2D->Clone("h2D4");
  TH2 *h2D5 = h2D1->Clone("h2D5");

  h2D->GetXaxis()->SetRangeUser(0,1.5);
  h2D1->GetXaxis()->SetRangeUser(0,1.5);
  h1D = h2D->ProjectionY();
  h1D1 = h2D1->ProjectionY();
  h1D1->Divide(h1D);
  h1D1->SetTitle("TPC-ITS Matching Efficiency (A)");
  h1D1->Draw("e0");

  can8->cd(4);
  h2D4->GetXaxis()->SetRangeUser(-1.5,0);
  h2D5->GetXaxis()->SetRangeUser(-1.5,0);
  h1D2 = h2D4->ProjectionY();
  h1D3 = h2D5->ProjectionY();
  h1D3->Divide(h1D2);
  h1D3->SetLineColor(2);
  h1D3->SetTitle("TPC-ITS Matching Efficiency (C)");
  h1D3->Draw("e0");

  can8->SaveAs("TPC-ITS.png");
  //  TH2 *h2D = 0;
  }

  //
  // ===| Pulls tracking parameters vs 1/pT |===================================
  //
  {
  TCanvas *can9 = new TCanvas("can9","Pulls of TPC Tracks vs 1/pT",1200,800);
  can9->Divide(3,2);

  // set proper fit range for # cluster plots
  //   the range of the plot is -5 to 5
  defaultFit = (TF1*)gROOT->GetFunction("gaus");
  if (!defaultFit) defaultFit = new TF1("gaus","gaus");
  defaultFit->SetParameters(1,0,0.1);
  defaultFit->SetParLimits(1, -5, 5);
  defaultFit->SetParLimits(2, 0.001, 3);

  TObjArray *arr1 = new TObjArray();
  TObjArray *arr2 = new TObjArray();
  TObjArray *arr3 = new TObjArray();
  TObjArray *arr4 = new TObjArray();
  TObjArray *arr5 = new TObjArray();

  can9->cd(1);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_0_7"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr1);
  h2D->Draw("colz");
  arr1->At(1)->Draw("same");

  can9->cd(2);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_1_7"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr2);
  h2D->Draw("colz");
  arr2->At(1)->Draw("same");

  can9->cd(3);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_2_7"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr3);
  h2D->Draw("colz");
  arr3->At(1)->Draw("same");

  can9->cd(4);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_3_7"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr4);
  h2D->Draw("colz");
  arr4->At(1)->Draw("same");

  can9->cd(5);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_4_7"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr5);
  h2D->Draw("colz");
  arr5->At(1)->Draw("same");

  can9->cd(6);
  PlotTimestamp(can9);

  can9->SaveAs("pull-pt.png");
  }

  //
  // ===| Pulls tracking parameters vs eta |====================================
  //
  {
  TCanvas *can10 = new TCanvas("can10","Pulls of TPC Tracks vs Eta",1200,800);
  can10->Divide(3,2);

  // set proper fit range for # cluster plots
  //   the range of the plot is -5 to 5
  defaultFit = (TF1*)gROOT->GetFunction("gaus");
  if (!defaultFit) defaultFit = new TF1("gaus","gaus");
  defaultFit->SetParameters(1,0,0.1);
  defaultFit->SetParLimits(1, -5, 5);
  defaultFit->SetParLimits(2, 0.001, 3);

  TObjArray *arr6 = new TObjArray();
  TObjArray *arr7 = new TObjArray();
  TObjArray *arr8 = new TObjArray();
  TObjArray *arr9 = new TObjArray();
  TObjArray *arr10 = new TObjArray();

  can10->cd(1);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_0_6"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr6);
  h2D->Draw("colz");
  arr6->At(1)->Draw("same");

  can10->cd(2);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_1_6"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr7);
  h2D->Draw("colz");
  arr7->At(1)->Draw("same");

  can10->cd(3);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_2_6"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr8);
  h2D->Draw("colz");
  arr8->At(1)->Draw("same");

  can10->cd(4);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_3_6"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr9);
  h2D->Draw("colz");
  arr9->At(1)->Draw("same");

  can10->cd(5);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_4_6"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr10);
  h2D->Draw("colz");
  arr10->At(1)->Draw("same");

  can10->cd(6);
  PlotTimestamp(can10);

  can10->SaveAs("pull-eta.png");
  }

  //
  // ===| Pulls tracking parameters vs phi |====================================
  //
  {
  TCanvas *can11 = new TCanvas("can11","Pulls of TPC Tracks vs Phi",1200,800);
  can11->Divide(3,2);

  // set proper fit range for # cluster plots
  //   the range of the plot is -5 to 5
  defaultFit = (TF1*)gROOT->GetFunction("gaus");
  if (!defaultFit) defaultFit = new TF1("gaus","gaus");
  defaultFit->SetParameters(1,0,0.1);
  defaultFit->SetParLimits(1, -5, 5);
  defaultFit->SetParLimits(2, 0.001, 3);

  TObjArray *arr11 = new TObjArray();
  TObjArray *arr12 = new TObjArray();
  TObjArray *arr13 = new TObjArray();
  TObjArray *arr14 = new TObjArray();
  TObjArray *arr15 = new TObjArray();

  can11->cd(1);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_0_5"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr11);
  h2D->Draw("colz");
  arr11->At(1)->Draw("same");

  can11->cd(2);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_1_5"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr12);
  h2D->Draw("colz");
  arr12->At(1)->Draw("same");

  can11->cd(3);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_2_5"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr13);
  h2D->Draw("colz");
  arr13->At(1)->Draw("same");

  can11->cd(4);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_3_5"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr14);
  h2D->Draw("colz");
  arr14->At(1)->Draw("same");

  can11->cd(5);
  h2D = (TH2*)(pPull->FindObject("h_tpc_match_pull_4_5"));
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr15);
  h2D->Draw("colz");
  arr15->At(1)->Draw("same");

  can11->cd(6);
  PlotTimestamp(can11);

  can11->SaveAs("pull-phi.png");
  }

  AliPerformanceMatch *obj4 = (AliPerformanceMatch*)TPC->FindObject("AliPerformanceMatchTPCConstrain");
  TFolder *pConstrain = obj4->GetAnalysisFolder();

  //
  // ===| phi contrained pulls |================================================
  //
  {
  TCanvas *can12 = new TCanvas("can12","#delta_{sin#phi}/#sigma_{sin#phi}",800,800);
  can12->Divide(2,2);

  h3D = (TH3*)pConstrain->FindObject("h_tpc_constrain_tpc_0_2_3");
  TH3 *h31 = h3D->Clone("h31");

  can12->cd(1);
  h3D->GetZaxis()->SetRangeUser(0,maxeta);
  //  h3D->GetYaxis()->SetRangeUser(0.25,10);
  h2D = (TH2*)h3D->Project3D("xy");
  h2D->Draw("colz");
  h2D->SetTitle("A Side");
  h2D->SetYTitle("(sin#phi_{TPC} - sin#phi_{Global})/#sigma");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr11);
  arr11->At(1)->Draw("same");
  TH1 *width1 = (TH1*)arr11->At(2);
  width1->SetMarkerColor(2);
  width1->SetLineColor(2);
  width1->Draw("same");

  /*  h3D->Project3D("xy")->Draw("colz");
      h3D->Project3D("xy")->SetTitle("A Side");
      h3D->Project3D("xy")->SetYTitle("#delta_{sin#phi}/#sigma_{sin#phi}");
      h3D->Project3D("xy")->FitSlicesY(0,0,-1,0,"QNR",arr11);
      arr11->At(1)->Draw("same");  */

  can12->cd(2);
  h31->GetZaxis()->SetRangeUser(mineta,-0.001);
  //  h31->GetYaxis()->SetRangeUser(0.25,10);
  h2D = (TH2*)h31->Project3D("xy");
  h2D->Draw("colz");
  h2D->SetTitle("C Side");
  h2D->SetYTitle("(sin#phi_{TPC} - sin#phi_{Global})/#sigma");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr12);
  arr12->At(1)->Draw("same");
  TH1 *width2 = (TH1*)arr12->At(2);
  width2->Draw("same");
  width2->SetMarkerColor(2);
  width2->SetLineColor(2);
  PlotTimestamp(can1);

  /*  h31->Project3D("xy")->Draw("colz");
      h31->Project3D("xy")->SetTitle("C Side");
      h31->Project3D("xy")->SetYTitle("#delta_{sin#phi}/#sigma_{sin#phi}");
      h31->Project3D("xy")->FitSlicesY(0,0,-1,0,"QNR",arr12);
      arr12->At(1)->Draw("same");  */

  can12->cd(3);
  h3D = (TH3*)pConstrain->FindObject("h_tpc_constrain_tpc_0_1_3");
  h3D->GetZaxis()->SetRangeUser(0,maxeta);
  TH3 *h32 = h3D->Clone("h32");
  h2D = (TH2*)h3D->Project3D("xy");
  h2D->Draw("colz");
  h2D->SetTitle("A Side");
  h2D->SetYTitle("(sin#phi_{TPC} - sin#phi_{Global})/#sigma");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr13);
  arr13->At(1)->Draw("same");
  TH1 *width3 = (TH1*)arr13->At(2);
  width3->SetMarkerColor(2);
  width3->SetLineColor(2);
  width3->Draw("same");

  /*  h3D->Project3D("xy")->Draw("colz");
      h3D->Project3D("xy")->SetTitle("A Side");
      h3D->Project3D("xy")->SetYTitle("#delta_{sin#phi}/#sigma_{sin#phi}");
      h3D->Project3D("xy")->FitSlicesY(0,0,-1,0,"QNR",arr13);
      arr13->At(1)->Draw("same");  */

  can12->cd(4);
  h32->GetZaxis()->SetRangeUser(mineta,-0.001);
  h2D = (TH2*)h32->Project3D("xy");
  h2D->Draw("colz");
  h2D->SetTitle("C Side");
  h2D->SetYTitle("(sin#phi_{TPC} - sin#phi_{Global})/#sigma");
  h2D->FitSlicesY(0,0,-1,0,"QNR",arr14);
  arr14->At(1)->Draw("same");
  TH1 *width4 = (TH1*)arr14->At(2);
  width4->SetMarkerColor(2);
  width4->SetLineColor(2);
  width4->Draw("same");

  /*  h32->Project3D("xy")->Draw("colz");
      h32->Project3D("xy")->SetTitle("C Side");
      h32->Project3D("xy")->SetYTitle("#delta_{sin#phi}/#sigma_{sin#phi}");
      h32->Project3D("xy")->FitSlicesY(0,0,-1,0,"QNR",arr14);
      arr14->At(1)->Draw("same");  */

  can12->SaveAs("pullPhiConstrain.png");
  }

  //
  // resolution and efficiency plots from David - added by Patrick
  //
  /*
     AliPerformanceRes *objPerfRes = (AliPerformanceRes*) TPC->FindObject("AliPerformanceRes");
     if (objPerfRes == NULL) {printf("Error getting AliPerformanceRes\n");}
     TFolder *folderRes = objPerfRes->GetAnalysisFolder();
     TH1F* h_resPt_vs_Pt       = (TH1F*)folderRes->FindObject("h_res_4_vs_9");
     TH1F* h_mean_resPt_vs_Pt  = (TH1F*)folderRes->FindObject("h_mean_res_4_vs_9");
     TH1F* h_pullPt_vs_Pt      = (TH1F*)folderRes->FindObject("h_pull_4_vs_9");
     TH1F* h_mean_pullPt_vs_Pt = (TH1F*)folderRes->FindObject("h_mean_pull_4_vs_9");

     TCanvas *can13 = new TCanvas("can13","Resolution p_{T}",800,800);
     can13->Divide(2,2);
     can13->cd(1);
     h_resPt_vs_Pt      ->Draw();
     can13->cd(2);
     h_mean_resPt_vs_Pt ->Draw();
     can13->cd(3);
     h_pullPt_vs_Pt     ->Draw();
     can13->cd(4);
     h_mean_pullPt_vs_Pt->Draw();

     can13->SaveAs("res_pT_1overpT.png");

     AliPerformanceEff *objPerfEff = (AliPerformanceEff*) TPC->FindObject("AliPerformanceEff");
     if (objPerfEff == NULL) {printf("Error getting AliPerformanceEff\n");}
     TFolder *folderEff = objPerfEff->GetAnalysisFolder();
     TH1F* eta_all              = (TH1F*)folderEff->FindObject("etaRecEff");
     TH1F* eta_all_neg          = (TH1F*)folderEff->FindObject("etaRecEffNeg");
     TH1F* eta_all_pos          = (TH1F*)folderEff->FindObject("etaRecEffPos");
     TH1F* phi_all              = (TH1F*)folderEff->FindObject("phiRecEff");
     TH1F* phi_all_neg          = (TH1F*)folderEff->FindObject("phiRecEffNeg");
     TH1F* phi_all_pos          = (TH1F*)folderEff->FindObject("phiRecEffPos");
     TH1F* pt_all               = (TH1F*)folderEff->FindObject("ptRecEff");
     TH1F* pt_all_neg           = (TH1F*)folderEff->FindObject("ptRecEffNeg");
     TH1F* pt_all_pos           = (TH1F*)folderEff->FindObject("ptRecEffPos");
     TH1F* eta_all_findable     = (TH1F*)folderEff->FindObject("etaRecEffF");
     TH1F* eta_all_findable_neg = (TH1F*)folderEff->FindObject("etaRecEffFNeg");
     TH1F* eta_all_findable_pos = (TH1F*)folderEff->FindObject("etaRecEffFPos");
     TH1F* phi_all_findable     = (TH1F*)folderEff->FindObject("phiRecEffF");
     TH1F* phi_all_findable_neg = (TH1F*)folderEff->FindObject("phiRecEffFNeg");
     TH1F* phi_all_findable_pos = (TH1F*)folderEff->FindObject("phiRecEffFPos");
     TH1F* pt_all_findable      = (TH1F*)folderEff->FindObject("ptRecEffF");
     TH1F* pt_all_findable_neg  = (TH1F*)folderEff->FindObject("ptRecEffFNeg");
     TH1F* pt_all_findable_pos  = (TH1F*)folderEff->FindObject("ptRecEffFPos");
     TH1F* eta_Pi               = (TH1F*)folderEff->FindObject("etaRecEffPi");
     TH1F* eta_Pi_neg           = (TH1F*)folderEff->FindObject("etaRecEffPiNeg");
     TH1F* eta_Pi_pos           = (TH1F*)folderEff->FindObject("etaRecEffPiPos");
     TH1F* phi_Pi               = (TH1F*)folderEff->FindObject("phiRecEffPi");
     TH1F* phi_Pi_neg           = (TH1F*)folderEff->FindObject("phiRecEffPiNeg");
     TH1F* phi_Pi_pos           = (TH1F*)folderEff->FindObject("phiRecEffPiPos");
     TH1F* pt_Pi                = (TH1F*)folderEff->FindObject("ptRecEffPi");
     TH1F* pt_Pi_neg            = (TH1F*)folderEff->FindObject("ptRecEffPiNeg");
     TH1F* pt_Pi_pos            = (TH1F*)folderEff->FindObject("ptRecEffPiPos");
     TH1F* eta_K                = (TH1F*)folderEff->FindObject("etaRecEffK");
     TH1F* eta_K_neg            = (TH1F*)folderEff->FindObject("etaRecEffKNeg");
     TH1F* eta_K_pos            = (TH1F*)folderEff->FindObject("etaRecEffKPos");
     TH1F* phi_K                = (TH1F*)folderEff->FindObject("phiRecEffK");
     TH1F* phi_K_neg            = (TH1F*)folderEff->FindObject("phiRecEffKNeg");
     TH1F* phi_K_pos            = (TH1F*)folderEff->FindObject("phiRecEffKPos");
     TH1F* pt_K                 = (TH1F*)folderEff->FindObject("ptRecEffK");
     TH1F* pt_K_neg             = (TH1F*)folderEff->FindObject("ptRecEffKNeg");
     TH1F* pt_K_pos             = (TH1F*)folderEff->FindObject("ptRecEffKPos");
     TH1F* eta_P                = (TH1F*)folderEff->FindObject("etaRecEffP");
     TH1F* eta_P_neg            = (TH1F*)folderEff->FindObject("etaRecEffPNeg");
     TH1F* eta_P_pos            = (TH1F*)folderEff->FindObject("etaRecEffPPos");
     TH1F* phi_P                = (TH1F*)folderEff->FindObject("phiRecEffP");
     TH1F* phi_P_neg            = (TH1F*)folderEff->FindObject("phiRecEffPNeg");
     TH1F* phi_P_pos            = (TH1F*)folderEff->FindObject("phiRecEffPPos");
     TH1F* pt_P                 = (TH1F*)folderEff->FindObject("ptRecEffP");
     TH1F* pt_P_neg             = (TH1F*)folderEff->FindObject("ptRecEffPNeg");
     TH1F* pt_P_pos             = (TH1F*)folderEff->FindObject("ptRecEffPPos");
     eta_all             ->SetLineWidth(2);
  eta_all_neg         ->SetLineColor(kRed);
  eta_all_pos         ->SetLineColor(kBlue);
  phi_all             ->SetLineWidth(2);
  phi_all_neg         ->SetLineColor(kRed);
  phi_all_pos         ->SetLineColor(kBlue);
  pt_all              ->SetLineWidth(2);
  pt_all_neg          ->SetLineColor(kRed);
  pt_all_pos          ->SetLineColor(kBlue);
  eta_all_findable    ->SetLineWidth(2);
  eta_all_findable_neg->SetLineColor(kRed);
  eta_all_findable_pos->SetLineColor(kBlue);
  phi_all_findable    ->SetLineWidth(2);
  phi_all_findable_neg->SetLineColor(kRed);
  phi_all_findable_pos->SetLineColor(kBlue);
  pt_all_findable     ->SetLineWidth(2);
  pt_all_findable_neg ->SetLineColor(kRed);
  pt_all_findable_pos ->SetLineColor(kBlue);
  eta_Pi              ->SetLineWidth(2);
  eta_Pi_neg          ->SetLineColor(kRed);
  eta_Pi_pos          ->SetLineColor(kBlue);
  phi_Pi              ->SetLineWidth(2);
  phi_Pi_neg          ->SetLineColor(kRed);
  phi_Pi_pos          ->SetLineColor(kBlue);
  pt_Pi               ->SetLineWidth(2);
  pt_Pi_neg           ->SetLineColor(kRed);
  pt_Pi_pos           ->SetLineColor(kBlue);
  eta_K               ->SetLineWidth(2);
  eta_K_neg           ->SetLineColor(kRed);
  eta_K_pos           ->SetLineColor(kBlue);
  phi_K               ->SetLineWidth(2);
  phi_K_neg           ->SetLineColor(kRed);
  phi_K_pos           ->SetLineColor(kBlue);
  pt_K                ->SetLineWidth(2);
  pt_K_neg            ->SetLineColor(kRed);
  pt_K_pos            ->SetLineColor(kBlue);
  eta_P               ->SetLineWidth(2);
  eta_P_neg           ->SetLineColor(kRed);
  eta_P_pos           ->SetLineColor(kBlue);
  phi_P               ->SetLineWidth(2);
  phi_P_neg           ->SetLineColor(kRed);
  phi_P_pos           ->SetLineColor(kBlue);
  pt_P                ->SetLineWidth(2);
  pt_P_neg            ->SetLineColor(kRed);
  pt_P_pos            ->SetLineColor(kBlue);

  TCanvas *can14 = new TCanvas("can14","Efficiency All",1000,800);
  can14->Divide(3, 2);
  can14->cd(1);
  eta_all             ->Draw();
  eta_all_neg         ->Draw("same");
  eta_all_pos         ->Draw("same");
  can14->cd(2);
  phi_all             ->Draw();
  phi_all_neg         ->Draw("same");
  phi_all_pos         ->Draw("same");
  can14->cd(3);
  pt_all              ->Draw();
  pt_all_neg          ->Draw("same");
  pt_all_pos          ->Draw("same");
  can14->cd(4);
  eta_all_findable    ->Draw();
  eta_all_findable_neg->Draw("same");
  eta_all_findable_pos->Draw("same");
  can14->cd(5);
  phi_all_findable    ->Draw();
  phi_all_findable_neg->Draw("same");
  phi_all_findable_pos->Draw("same");
  can14->cd(6);
  pt_all_findable     ->Draw();
  pt_all_findable_neg ->Draw("same");
  pt_all_findable_pos ->Draw("same");

  can14->SaveAs("eff_all+all_findable.png");

  TCanvas *can15 = new TCanvas("can15","Efficiency Pi K P",1000,1000);
  can15->Divide(3, 3);
  can15->cd(1);
  eta_Pi              ->Draw();
  eta_Pi_neg          ->Draw("same");
  eta_Pi_pos          ->Draw("same");
  can15->cd(2);
  phi_Pi              ->Draw();
  phi_Pi_neg          ->Draw("same");
  phi_Pi_pos          ->Draw("same");
  can15->cd(3);
  pt_Pi               ->Draw();
  pt_Pi_neg           ->Draw("same");
  pt_Pi_pos           ->Draw("same");
  can15->cd(4);
  eta_K               ->Draw();
  eta_K_neg           ->Draw("same");
  eta_K_pos           ->Draw("same");
  can15->cd(5);
  phi_K               ->Draw();
  phi_K_neg           ->Draw("same");
  phi_K_pos           ->Draw("same");
  can15->cd(6);
  pt_K                ->Draw();
  pt_K_neg            ->Draw("same");
  pt_K_pos            ->Draw("same");
  can15->cd(7);
  eta_P               ->Draw();
  eta_P_neg           ->Draw("same");
  eta_P_pos           ->Draw("same");
  can15->cd(8);
  phi_P               ->Draw();
  phi_P_neg           ->Draw("same");
  phi_P_pos           ->Draw("same");
  can15->cd(9);
  pt_P                ->Draw();
  pt_P_neg            ->Draw("same");
  pt_P_pos            ->Draw("same");

  can15->SaveAs("eff_Pi_K_P.png");

  //get more histos from THnSparse...
  THnSparseF* EffHisto = (THnSparseF*) objPerfEff->GetEffHisto();
  //mceta:mcphi:mcpt:pid:recStatus:findable:charge:nclones:nfakes
  //pid:e-=0,mu-=1,pi+=2,K+=3,p+=4
  TH1* h_pid    = (TH1*) EffHisto->Projection(3);
  TH1* h_charge = (TH1*) EffHisto->Projection(6);
  //TH1* h_find  = (TH1*) EffHisto->Projection(5);
  //TH1* eta_All = (TH1*) EffHisto->Projection(0);

  cout<<"before setrange"<<endl;
  EffHisto->GetAxis(5)->SetRangeUser(1, 1.999); //set findable
  EffHisto->GetAxis(3)->SetRangeUser(2, 2.999); //set pion
  cout<<"after setrange"<<endl;
  TH1* h_charge_sel = (TH1*) EffHisto->Projection(6);
  cout<<"after projection"<<endl;
  EffHisto->GetAxis(6)->UnZoom(); //set all charges
  TH1* eta_Pi_findable = (TH1*) EffHisto->Projection(0);
  TH1* phi_Pi_findable = (TH1*) EffHisto->Projection(1);
  TH1* pt_Pi_findable  = (TH1*) EffHisto->Projection(2);
  EffHisto->GetAxis(6)->SetRangeUser(-1.5, -0.499); //set negative
  TH1* eta_Pi_findable_neg = (TH1*) EffHisto->Projection(0);
  TH1* phi_Pi_findable_neg = (TH1*) EffHisto->Projection(1);
  TH1* pt_Pi_findable_neg  = (TH1*) EffHisto->Projection(2);
  EffHisto->GetAxis(6)->SetRangeUser(0.5, 1.499); //set positive
  TH1* eta_Pi_findable_pos = (TH1*) EffHisto->Projection(0);
  TH1* phi_Pi_findable_pos = (TH1*) EffHisto->Projection(1);
  TH1* pt_Pi_findable_pos  = (TH1*) EffHisto->Projection(2);
  //EffHisto->GetAxis(3)->SetRangeUser(3, 3.999); //set kaon
  //EffHisto->GetAxis(6)->UnZoom(); //set all charges

  //EffHisto->GetAxis(3)->SetRangeUser(4, 4.999); //set proton

  eta_Pi_findable     ->SetLineWidth(2);
  eta_Pi_findable_neg ->SetLineColor(kRed);
  eta_Pi_findable_pos ->SetLineColor(kBlue);
  phi_Pi_findable     ->SetLineWidth(2);
  phi_Pi_findable_neg ->SetLineColor(kRed);
  phi_Pi_findable_pos ->SetLineColor(kBlue);
  pt_Pi_findable      ->SetLineWidth(2);
  pt_Pi_findable_neg  ->SetLineColor(kRed);
  pt_Pi_findable_pos  ->SetLineColor(kBlue);
  //eta_K_findable      ->SetLineWidth(2);
  //eta_K_findable_neg  ->SetLineColor(kRed);
  //eta_K_findable_pos  ->SetLineColor(kBlue);
  //phi_K_findable      ->SetLineWidth(2);
  //phi_K_findable_neg  ->SetLineColor(kRed);
  //phi_K_findable_pos  ->SetLineColor(kBlue);
  //pt_K_findable       ->SetLineWidth(2);
  //pt_K_findable_neg   ->SetLineColor(kRed);
  //pt_K_findable_pos   ->SetLineColor(kBlue);
  //eta_P_findable      ->SetLineWidth(2);
  //eta_P_findable_neg  ->SetLineColor(kRed);
  //eta_P_findable_pos  ->SetLineColor(kBlue);
  //phi_P_findable      ->SetLineWidth(2);
  //phi_P_findable_neg  ->SetLineColor(kRed);
  //phi_P_findable_pos  ->SetLineColor(kBlue);
  //pt_P_findable       ->SetLineWidth(2);
  //pt_P_findable_neg   ->SetLineColor(kRed);
  //pt_P_findable_pos   ->SetLineColor(kBlue);

  TCanvas *can16 = new TCanvas("can16","Efficiency Pi K P findable",1000,1000);
  can16->Divide(3,3);
  can16->cd(7);
  h_pid->Draw();
  can16->cd(8);
  h_charge->Draw();
  can16->cd(9);
  h_charge_sel->Draw();
  can16->cd(1);
  eta_Pi_findable->Draw();
  eta_Pi_findable_neg->Draw("same");
  eta_Pi_findable_pos->Draw("same");
  can16->cd(2);
  phi_Pi_findable->Draw();
  phi_Pi_findable_neg->Draw("same");
  phi_Pi_findable_pos->Draw("same");
  can16->cd(3);
  pt_Pi_findable->Draw();
  pt_Pi_findable_neg->Draw("same");
  pt_Pi_findable_pos->Draw("same");

  can16->SaveAs("eff_Pi_K_P_findable.png");
  */

    ofstream fout("runqa_exist");
  //if(NEvents>0)
  fout.precision(10);
  fout<<NEvents<<endl;
  fout.close();

}
