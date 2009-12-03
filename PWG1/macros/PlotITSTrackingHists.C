void PlotITSTrackingHists() {
  //
  // Macro to plot the histos from the task AliAnalysisTaskITSTrackingCheck
  // A. Dainese 28.11.09
  // 

  gStyle->SetOptStat(0);

  TFile *f= new TFile("ITS.Performance.root");

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
  fHistClusterMapITSMIokoutinzbad->SetLineColor(1);
  fHistClusterMapITSMIokoutinzbad->SetMarkerColor(1);
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


  return;
}
