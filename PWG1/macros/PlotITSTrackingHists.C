void Test(TString fname="QAresults.root") {

  TCanvas *clist;
  Int_t cnum=1;
  TFile *f=new TFile(fname.Data());
  PlotITSTPCMatchingEff(f,clist,cnum);


  return;
}
//--------------------------------------------------------------------------
Bool_t PlotITSTPCMatchingEff(TFile *f, TCanvas*& clist,Int_t& cnum) {

  Int_t cnum=1;
  //clist = new TCanvas[2];


  TCanvas* cITSTPCmatch = new TCanvas("cITSTPCmatch","cITSTPCmatch",10,10,1200,600);
  cITSTPCmatch->Divide(2,1);
  cITSTPCmatch_1->SetGridy();
  cITSTPCmatch_2->SetGridy();
  cITSTPCmatch_1->SetLogx();
  cITSTPCmatch_2->SetLogx();

  clist = cITSTPCmatch;

  if(!f) return kFALSE;

  TList *list=0;
  TList *listSPD=0;
  TDirectoryFile *dir=0;

  // count active SPD HSs
  dir=(TDirectoryFile*)f->GetDirectory("SPD_Performance");
  if(dir) listSPD = (TList*)dir->Get("coutput1");
  if(!dir) return kFALSE;

  Float_t spdFrac[2];
  TH1F *hnHSsSPD=new TH1F("hnHSsSPD","Active HSs in SPD layers 1 and 2; layer; HSs",2,0.5,2.5);
  if(listSPD) {
    //listSPD->Print();
    TH1F *hFiredChip = (TH1F*)listSPD->FindObject("hFiredChip");
    Int_t nHSsInner=0,nHSsOuter=0;
    for(Int_t i=0;i<400;i++) if(hFiredChip->GetBinContent(i)>0) nHSsInner++;
    for(Int_t i=400;i<1200;i++) if(hFiredChip->GetBinContent(i)>0) nHSsOuter++;
    nHSsInner = (Int_t)(nHSsInner/10);
    nHSsOuter = (Int_t)(nHSsOuter/10);
    hnHSsSPD->SetBinContent(1,nHSsInner);
    hnHSsSPD->SetBinContent(2,nHSsOuter);
    spdFrac[0]=(Float_t)nHSsInner/40.;
    spdFrac[1]=(Float_t)nHSsOuter/80.;
  }
  TGraph *spdFrac0=new TGraph(1);
  spdFrac0->SetPoint(0,0.08,spdFrac[0]);
  spdFrac0->SetMarkerColor(1); spdFrac0->SetMarkerStyle(20);
  TGraph *spdFrac1=new TGraph(1);
  spdFrac1->SetPoint(0,0.08,spdFrac[1]);
  spdFrac1->SetMarkerColor(1); spdFrac1->SetMarkerStyle(24);
  TLegend *l2=new TLegend(0.1,0.62,0.5,0.93);
  l2->SetBorderSize(1);
  l2->AddEntry(spdFrac0,"Frac. active SPD0","p");
  l2->AddEntry(spdFrac1,"Frac. active SPD1","p");

  //
  // Efficiencies for CENTRAL
  //
  dir=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
  if(dir) list = (TList*)dir->Get("cOutputITS_3500_10000");
  if(!list) return kFALSE;

  TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
  TH1F *fHistPtITSTPCsel = (TH1F*)list->FindObject("fHistPtITSTPCsel");
  TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);


  TLegend *l3=new TLegend(0.5,0.62,0.95,0.93);
  l3->SetBorderSize(1);
  cITSTPCmatch->cd(1);
  fHistPtITSMIge2InAcc->SetTitle("Fraction of prolonged tracks with N ITS points: central");
  fHistPtITSMIge2InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIge2InAcc->SetMaximum(1.6);
  fHistPtITSMIge2InAcc->SetMinimum(0);
  fHistPtITSMIge2InAcc->GetXaxis()->SetRangeUser(0.1,30);
  fHistPtITSMIge2InAcc->Draw();
  l3->AddEntry(fHistPtITSMIge2InAcc,">=2 cls","l");
  fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI6InAcc->SetLineColor(2);
  l3->AddEntry(fHistPtITSMI6InAcc,"6 cls","l");
  fHistPtITSMI6InAcc->Draw("same");
  fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI5InAcc->SetLineColor(3);
  l3->AddEntry(fHistPtITSMI5InAcc,"5 cls","l");
  fHistPtITSMI5InAcc->Draw("same");
  fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI4InAcc->SetLineColor(4);
  l3->AddEntry(fHistPtITSMI4InAcc,"4 cls","l");
  fHistPtITSMI4InAcc->Draw("same");
  fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI3InAcc->SetLineColor(6);
  l3->AddEntry(fHistPtITSMI3InAcc,"3 cls","l");
  fHistPtITSMI3InAcc->Draw("same");
  fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI2InAcc->SetLineColor(7);
  l3->AddEntry(fHistPtITSMI2InAcc,"2 cls","l");
  fHistPtITSMI2InAcc->Draw("same");
  fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMISPDInAcc->SetLineColor(9);
  l3->AddEntry(fHistPtITSMISPDInAcc,"2SPD + any","l");
  fHistPtITSMISPDInAcc->Draw("same");
  fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIoneSPDInAcc->SetLineColor(15);
  l3->AddEntry(fHistPtITSMIoneSPDInAcc,">=1SPD + any","l");
  fHistPtITSMIoneSPDInAcc->Draw("same");
  fHistPtITSTPCsel->Divide(fHistPtITSTPCsel,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSTPCsel->SetLineColor(kAzure+1);
  l3->AddEntry(fHistPtITSTPCsel,">=1SPD + any + d_{0} cut","l");
  fHistPtITSTPCsel->Draw("same");
  fHistPtITSMIge2InAcc->Draw("same");
  l3->Draw();
  l2->Draw();
  spdFrac0->Draw("p");
  spdFrac1->Draw("p");
  /*
  if(ioValues) {
    Int_t index,ptbin;
    ptbin=fHistPtITSMIge2InAcc->FindBin(0.201);
    index=2;
    ioValues[index]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(1.001);
    index=3;
    ioValues[index]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(10.001);
    index=4;
    ioValues[index]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(0.201);
    index=5;
    ioValues[index]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(1.001);
    index=6;
    ioValues[index]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(10.001);
    index=7;
    ioValues[index]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(0.201);
    index=8;
    ioValues[index]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(1.001);
    index=9;
    ioValues[index]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(10.001);
    index=10;
    ioValues[index]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(0.201);
    index=11;
    ioValues[index]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(1.001);
    index=12;
    ioValues[index]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(10.001);
    index=13;
    ioValues[index]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
  }
  */



  //
  // Efficiencies for PERIPHERAL
  //
  dir=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
  if(dir) list = (TList*)dir->Get("cOutputITS_70_310");
  if(!list) return kFALSE;

  TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
  TH1F *fHistPtITSTPCsel = (TH1F*)list->FindObject("fHistPtITSTPCsel");
  TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);


  cITSTPCmatch->cd(2);
  fHistPtITSMIge2InAcc->SetTitle("Fraction of prolonged tracks with N ITS points: peripheral");
  fHistPtITSMIge2InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIge2InAcc->SetMaximum(1.6);
  fHistPtITSMIge2InAcc->SetMinimum(0);
  fHistPtITSMIge2InAcc->GetXaxis()->SetRangeUser(0.1,30);
  fHistPtITSMIge2InAcc->Draw();
  fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI6InAcc->SetLineColor(2);
  fHistPtITSMI6InAcc->Draw("same");
  fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI5InAcc->SetLineColor(3);
  fHistPtITSMI5InAcc->Draw("same");
  fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI4InAcc->SetLineColor(4);
  fHistPtITSMI4InAcc->Draw("same");
  fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI3InAcc->SetLineColor(6);
  fHistPtITSMI3InAcc->Draw("same");
  fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI2InAcc->SetLineColor(7);
  fHistPtITSMI2InAcc->Draw("same");
  fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMISPDInAcc->SetLineColor(9);
  fHistPtITSMISPDInAcc->Draw("same");
  fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIoneSPDInAcc->SetLineColor(15);
  fHistPtITSMIoneSPDInAcc->Draw("same");
  fHistPtITSTPCsel->Divide(fHistPtITSTPCsel,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSTPCsel->SetLineColor(kAzure+1);
  fHistPtITSTPCsel->Draw("same");
  fHistPtITSMIge2InAcc->Draw("same");
  l3->Draw();
  l2->Draw();
  spdFrac0->Draw("p");
  spdFrac1->Draw("p");
  /*
  if(ioValues) {
    Int_t index,ptbin;
    ptbin=fHistPtITSMIge2InAcc->FindBin(0.201);
    index=2;
    ioValues[index]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(1.001);
    index=3;
    ioValues[index]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(10.001);
    index=4;
    ioValues[index]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(0.201);
    index=5;
    ioValues[index]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(1.001);
    index=6;
    ioValues[index]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(10.001);
    index=7;
    ioValues[index]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(0.201);
    index=8;
    ioValues[index]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(1.001);
    index=9;
    ioValues[index]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(10.001);
    index=10;
    ioValues[index]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(0.201);
    index=11;
    ioValues[index]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(1.001);
    index=12;
    ioValues[index]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(10.001);
    index=13;
    ioValues[index]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
  }
  */



  return kTRUE;
}

Bool_t PlotITSTrackingHists(TString fname="ITS.Performance.root",
			    Float_t *ioValues=0,Float_t *ioErrors=0) 
{
  //
  // Macro to plot the histos from the task AliAnalysisTaskITSTrackingCheck
  // A. Dainese 28.11.09
  // 

  Bool_t plotAlignmentChecks=kFALSE;
  //gStyle->SetOptStat(0);

  //if(fname.Contains("alien")) TGrid::Connect("alien://");

  TFile *f= TFile::Open(fname.Data());
  if(!f) return kFALSE;

  TList *list=(TList*)f->Get("cOutputITS");
  TList *listSPD=0;
  TDirectoryFile *dir=0;
  if(!list) {
    dir=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
    if(dir) list = (TList*)dir->Get("cOutputITS_70_310");
    // count active SPD HSs
    dir=(TDirectoryFile*)f->GetDirectory("SPD_Performance");
    if(dir) listSPD = (TList*)dir->Get("coutput1");
  }

  if(!list) return kFALSE;

  TH1F *hnHSsSPD=new TH1F("hnHSsSPD","Active HSs in SPD layers 1 and 2; layer; HSs",2,0.5,2.5);
  if(listSPD) {
    //listSPD->Print();
    TH1F *hFiredChip = (TH1F*)listSPD->FindObject("hFiredChip");
    Int_t nHSsInner=0,nHSsOuter=0;
    for(Int_t i=0;i<400;i++) if(hFiredChip->GetBinContent(i)>0) nHSsInner++;
    for(Int_t i=400;i<1200;i++) if(hFiredChip->GetBinContent(i)>0) nHSsOuter++;
    nHSsInner = (Int_t)(nHSsInner/10);
    nHSsOuter = (Int_t)(nHSsOuter/10);
    hnHSsSPD->SetBinContent(1,nHSsInner);
    hnHSsSPD->SetBinContent(2,nHSsOuter);
    if(ioValues) ioValues[0]=(Float_t)nHSsInner/40.;
    if(ioValues) ioValues[1]=(Float_t)nHSsOuter/80.;
  }

  TCanvas *cSPD=new TCanvas("cSPD","cSPD");
  hnHSsSPD->SetMaximum(90);
  hnHSsSPD->SetMinimum(0);
  hnHSsSPD->Draw();

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
  TH1F *fHistClusterMapModuleITSMIokInAcc = (TH1F*)list->FindObject("fHistClusterMapModuleITSMIokInAcc");
  TH1F *fHistClusterMapModuleITSMIbadInAcc = (TH1F*)list->FindObject("fHistClusterMapModuleITSMIbadInAcc");
  TH1F *fHistClusterMapITSSAskippedInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAskippedInAcc");
  TH1F *fHistClusterMapITSSAoutinzInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAoutinzInAcc");
  TH1F *fHistClusterMapITSSAokoutinzbadInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAokoutinzbadInAcc");
  TH1F *fHistClusterMapITSSAnorefitInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAnorefitInAcc");
  TH1F *fHistClusterMapITSSAnoclsInAcc = (TH1F*)list->FindObject("fHistClusterMapITSSAnoclsInAcc");
  TH1F *fHistClusterMapModuleITSMInoclsInAcc = (TH1F*)list->FindObject("fHistClusterMapModuleITSMInoclsInAcc");
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
  TH1F *fHistPtTPCInAccSfromStrange = (TH1F*)list->FindObject("fHistPtTPCInAccSfromStrange");
  TH1F *fHistPtTPCInAccPfromStrange = (TH1F*)list->FindObject("fHistPtTPCInAccPfromStrange");
  TH1F *fHistPtTPCInAccMCtwoSPD = (TH1F*)list->FindObject("fHistPtTPCInAccMCtwoSPD");
  TH1F *fHistPtTPCInAccMConeSPD = (TH1F*)list->FindObject("fHistPtTPCInAccMConeSPD");
  TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
  TH1F *fHistPtITSMI6InAccFake = (TH1F*)list->FindObject("fHistPtITSMI6InAccFake");
  TH1F *fHistPtITSMI5InAccFake = (TH1F*)list->FindObject("fHistPtITSMI5InAccFake");
  TH1F *fHistPtITSMI4InAccFake = (TH1F*)list->FindObject("fHistPtITSMI4InAccFake");
  TH1F *fHistPtITSMI3InAccFake = (TH1F*)list->FindObject("fHistPtITSMI3InAccFake");
  TH1F *fHistPtITSMI2InAccFake = (TH1F*)list->FindObject("fHistPtITSMI2InAccFake");
  TH1F *fHistPtITSMISPDInAccFake = (TH1F*)list->FindObject("fHistPtITSMISPDInAccFake");
  TH1F *fHistPtITSMIoneSPDInAccFake = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAccFake");
  TH1F *fHistPtITSMIoneSPDthreeSDDSSDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDthreeSDDSSDInAcc");
  TH1F *fHistPtITSMIokbadoutinz6InAcc = (TH1F*)list->FindObject("fHistPtITSMIokbadoutinz6InAcc");
  TH1F *fHistPtITSMIokbadoutinz5InAcc = (TH1F*)list->FindObject("fHistPtITSMIokbadoutinz5InAcc");
  TH1F *fHistPtITSMIokbadoutinz4InAcc = (TH1F*)list->FindObject("fHistPtITSMIokbadoutinz4InAcc");
  TH1F *fHistPtTPCInAccP = (TH1F*)list->FindObject("fHistPtTPCInAccP");
  TH1F *fHistPtITSMI6InAccP = (TH1F*)list->FindObject("fHistPtITSMI6InAccP");
  TH1F *fHistPtITSMI5InAccP = (TH1F*)list->FindObject("fHistPtITSMI5InAccP");
  TH1F *fHistPtITSMI4InAccP = (TH1F*)list->FindObject("fHistPtITSMI4InAccP");
  TH1F *fHistPtITSMI3InAccP = (TH1F*)list->FindObject("fHistPtITSMI3InAccP");
  TH1F *fHistPtITSMI2InAccP = (TH1F*)list->FindObject("fHistPtITSMI2InAccP");
  TH1F *fHistPtITSMISPDInAccP = (TH1F*)list->FindObject("fHistPtITSMISPDInAccP");
  TH1F *fHistPtITSMIoneSPDInAccP = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAccP");
  TH1F *fHistPtTPCInAccS = (TH1F*)list->FindObject("fHistPtTPCInAccS");
  TH1F *fHistPtITSMI6InAccS = (TH1F*)list->FindObject("fHistPtITSMI6InAccS");
  TH1F *fHistPtITSMI5InAccS = (TH1F*)list->FindObject("fHistPtITSMI5InAccS");
  TH1F *fHistPtITSMI4InAccS = (TH1F*)list->FindObject("fHistPtITSMI4InAccS");
  TH1F *fHistPtITSMI3InAccS = (TH1F*)list->FindObject("fHistPtITSMI3InAccS");
  TH1F *fHistPtITSMI2InAccS = (TH1F*)list->FindObject("fHistPtITSMI2InAccS");
  TH1F *fHistPtITSMISPDInAccS = (TH1F*)list->FindObject("fHistPtITSMISPDInAccS");
  TH1F *fHistPtITSMIoneSPDInAccS = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAccS");

  TH1F *fHistRProdVtxInAccP = (TH1F*)list->FindObject("fHistRProdVtxInAccP");
  TH1F *fHistRProdVtxInAccS = (TH1F*)list->FindObject("fHistRProdVtxInAccS");

  TH1F *fHistxlocSDDok = (TH1F*)list->FindObject("fHistxlocSDDok");
  TH1F *fHistxlocSDDall = (TH1F*)list->FindObject("fHistxlocSDDall");
  TH1F *fHistzlocSDDok = (TH1F*)list->FindObject("fHistzlocSDDok");
  TH1F *fHistzlocSDDall = (TH1F*)list->FindObject("fHistzlocSDDall");

  TH1F *fHistPtITSTPCsel = (TH1F*)list->FindObject("fHistPtITSTPCsel");
  TH1F *fHistPtITSTPCselP = (TH1F*)list->FindObject("fHistPtITSTPCselP");
  TH1F *fHistPtITSTPCselS = (TH1F*)list->FindObject("fHistPtITSTPCselS");
  TH1F *fHistPtITSTPCselFake = (TH1F*)list->FindObject("fHistPtITSTPCselFake");
  TH1F *fHistPtITSTPCselPfromStrange = (TH1F*)list->FindObject("fHistPtITSTPCselSfromStrange");
  TH1F *fHistPtITSTPCselSfromStrange = (TH1F*)list->FindObject("fHistPtITSTPCselSfromStrange");
  TH1F *fHistPtITSTPCselSfromMat = (TH1F*)list->FindObject("fHistPtITSTPCselSfromMat");

  TH1F *fHistPtTPCInAccSfromMat = (TH1F*)list->FindObject("fHistPtTPCInAccSfromMat");


  //ReweightStrange(fHistPtTPCInAcc,fHistPtTPCInAccPfromStrange,fHistPtTPCInAccSfromStrange);
  //ReweightStrange(fHistPtITSTPCsel,fHistPtITSTPCselPfromStrange,fHistPtITSTPCselSfromStrange);

  //---------------------------------------------------------------
  TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
  fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);
  TH1F *fHistPtITSMIge2InAccFake = 0;
  if(fHistPtITSMI6InAccFake) {
    fHistPtITSMIge2InAccFake=(TH1F*)fHistPtITSMI6InAccFake->Clone("fHistPtITSMIge2InAccFake");
    fHistPtITSMIge2InAccFake->Add(fHistPtITSMI5InAccFake);
    fHistPtITSMIge2InAccFake->Add(fHistPtITSMI4InAccFake);
    fHistPtITSMIge2InAccFake->Add(fHistPtITSMI3InAccFake);
    fHistPtITSMIge2InAccFake->Add(fHistPtITSMI2InAccFake);
  }
  TH1F *fHistPtITSMIge2InAccP = (TH1F*)fHistPtITSMI6InAccP->Clone("fHistPtITSMIge2InAccP");
  fHistPtITSMIge2InAccP->Add(fHistPtITSMI5InAccP);
  fHistPtITSMIge2InAccP->Add(fHistPtITSMI4InAccP);
  fHistPtITSMIge2InAccP->Add(fHistPtITSMI3InAccP);
  fHistPtITSMIge2InAccP->Add(fHistPtITSMI2InAccP);
  TH1F *fHistPtITSMIge2InAccS = (TH1F*)fHistPtITSMI6InAccS->Clone("fHistPtITSMIge2InAccS");
  fHistPtITSMIge2InAccS->Add(fHistPtITSMI5InAccS);
  fHistPtITSMIge2InAccS->Add(fHistPtITSMI4InAccS);
  fHistPtITSMIge2InAccS->Add(fHistPtITSMI3InAccS);
  fHistPtITSMIge2InAccS->Add(fHistPtITSMI2InAccS);

  // fake fraction
  if(fHistPtITSMIge2InAccFake) {
    fHistPtITSMIge2InAccFake->Divide(fHistPtITSMIge2InAccFake,fHistPtITSMIge2InAcc,1,1,"B");
    fHistPtITSMI2InAccFake->Divide(fHistPtITSMI2InAccFake,fHistPtITSMI2InAcc,1,1,"B");
    fHistPtITSMI3InAccFake->Divide(fHistPtITSMI3InAccFake,fHistPtITSMI3InAcc,1,1,"B");
    fHistPtITSMI4InAccFake->Divide(fHistPtITSMI4InAccFake,fHistPtITSMI4InAcc,1,1,"B");
    fHistPtITSMI5InAccFake->Divide(fHistPtITSMI5InAccFake,fHistPtITSMI5InAcc,1,1,"B");
    fHistPtITSMI6InAccFake->Divide(fHistPtITSMI6InAccFake,fHistPtITSMI6InAcc,1,1,"B");
    fHistPtITSMISPDInAccFake->Divide(fHistPtITSMISPDInAccFake,fHistPtITSMISPDInAcc,1,1,"B");
    fHistPtITSMIoneSPDInAccFake->Divide(fHistPtITSMIoneSPDInAccFake,fHistPtITSMIoneSPDInAcc,1,1,"B");
    if (fHistPtITSTPCselFake)fHistPtITSTPCselFake->Divide(fHistPtITSTPCselFake,fHistPtITSTPCsel,1,1,"B");
  }

  TH1F* fHistPtITSMISPDInAccMC=(TH1F*)fHistPtITSMISPDInAcc->Clone("fHistPtITSMISPDInAccMC");
  TH1F* fHistPtITSMIoneSPDInAccMC=(TH1F*)fHistPtITSMIoneSPDInAcc->Clone("fHistPtITSMIoneSPDInAccMC");


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

  TCanvas *c3 =new TCanvas("c3","c3",10,10,1200,800);
  c3->Divide(1,2);
  c3->cd(1);
  fHistClusterMapModuleITSMIokInAcc->SetLineColor(1);
  fHistClusterMapModuleITSMIokInAcc->Draw();
  l2->AddEntry(fHistClusterMapModuleITSMIokInAcc,"ok","l");
  fHistClusterMapModuleITSMIbadInAcc->SetLineColor(3);
  fHistClusterMapModuleITSMIbadInAcc->Draw("same");
  l2->AddEntry(fHistClusterMapModuleITSMIbadInAcc,"bad","l");
  fHistClusterMapModuleITSMInoclsInAcc->SetLineColor(2);
  fHistClusterMapModuleITSMInoclsInAcc->Draw("same");
  l2->AddEntry(fHistClusterMapModuleITSMInoclsInAcc,"no cls","l");
  l2->Draw();
  c3->cd(2);
  TH1F *fHistClusterMapModuleITSMIallInAcc=(TH1F*)fHistClusterMapModuleITSMIokInAcc->Clone("fHistClusterMapModuleITSMIallInAcc");
  fHistClusterMapModuleITSMIallInAcc->Add(fHistClusterMapModuleITSMIbadInAcc);
  fHistClusterMapModuleITSMIallInAcc->Add(fHistClusterMapModuleITSMInoclsInAcc);
  TH1F *fHistClusterMapModuleITSMIokRatioInAcc=(TH1F*)fHistClusterMapModuleITSMIokInAcc->Clone("fHistClusterMapModuleITSMIokRatioInAcc");
  fHistClusterMapModuleITSMIokRatioInAcc->Divide(fHistClusterMapModuleITSMIallInAcc);
  for(Int_t ib=1;ib<=fHistClusterMapModuleITSMIokRatioInAcc->GetNbinsX();ib++) {
    fHistClusterMapModuleITSMIokRatioInAcc->SetBinError(ib,0);
    if(fHistClusterMapModuleITSMIallInAcc->GetBinContent(ib)) fHistClusterMapModuleITSMIokRatioInAcc->SetBinError(ib,TMath::Sqrt(fHistClusterMapModuleITSMIokRatioInAcc->GetBinContent(ib)*(1.-fHistClusterMapModuleITSMIokRatioInAcc->GetBinContent(ib))/fHistClusterMapModuleITSMIallInAcc->GetBinContent(ib))); 
  }
  fHistClusterMapModuleITSMIokRatioInAcc->Draw();
  TH1F *fHistClusterMapModuleITSMIbadRatioInAcc=(TH1F*)fHistClusterMapModuleITSMIbadInAcc->Clone("fHistClusterMapModuleITSMIbadRatioInAcc");
  fHistClusterMapModuleITSMIbadRatioInAcc->Divide(fHistClusterMapModuleITSMIallInAcc);
  for(Int_t ib=1;ib<=fHistClusterMapModuleITSMIokRatioInAcc->GetNbinsX();ib++) {
    fHistClusterMapModuleITSMIbadRatioInAcc->SetBinError(ib,0);
    if(fHistClusterMapModuleITSMIallInAcc->GetBinContent(ib)) fHistClusterMapModuleITSMIbadRatioInAcc->SetBinError(ib,TMath::Sqrt(fHistClusterMapModuleITSMIbadRatioInAcc->GetBinContent(ib)*(1.-fHistClusterMapModuleITSMIbadRatioInAcc->GetBinContent(ib))/fHistClusterMapModuleITSMIallInAcc->GetBinContent(ib))); 
  }
  fHistClusterMapModuleITSMIbadRatioInAcc->Draw("same");
  TH1F *fHistClusterMapModuleITSMInoclsRatioInAcc=(TH1F*)fHistClusterMapModuleITSMInoclsInAcc->Clone("fHistClusterMapModuleITSMInoclsRatioInAcc");
  fHistClusterMapModuleITSMInoclsRatioInAcc->Divide(fHistClusterMapModuleITSMIallInAcc);
  for(Int_t ib=1;ib<=fHistClusterMapModuleITSMIokRatioInAcc->GetNbinsX();ib++) {
    fHistClusterMapModuleITSMInoclsRatioInAcc->SetBinError(ib,0);
    if(fHistClusterMapModuleITSMIallInAcc->GetBinContent(ib)) fHistClusterMapModuleITSMInoclsRatioInAcc->SetBinError(ib,TMath::Sqrt(fHistClusterMapModuleITSMInoclsRatioInAcc->GetBinContent(ib)*(1.-fHistClusterMapModuleITSMInoclsRatioInAcc->GetBinContent(ib))/fHistClusterMapModuleITSMIallInAcc->GetBinContent(ib))); 
  }
  fHistClusterMapModuleITSMInoclsRatioInAcc->Draw("same");

  TCanvas *c3b = new TCanvas("c3b","c3b");
  c3b->Divide(2,1);
  c3b->cd(1);
  fHistxlocSDDok->Divide(fHistxlocSDDall);
  fHistxlocSDDok->Draw();
  c3b->cd(2);
  fHistzlocSDDok->Divide(fHistzlocSDDall);
  fHistzlocSDDok->Draw();

  TCanvas *c4 =new TCanvas("c4","c4",10,10,500,500);
  c4->SetGridy();
  fHistPhiITSMIokbadoutinz6InAcc->Divide(fHistPhiTPCInAcc);
  fHistPhiITSMIokbadoutinz6InAcc->SetMinimum(0);
  fHistPhiITSMIokbadoutinz6InAcc->SetMaximum(1.5);
  fHistPhiITSMIokbadoutinz6InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPhiITSMIokbadoutinz6InAcc->SetTitle("Fraction of tracks with 6 layers ok");
  fHistPhiITSMIokbadoutinz6InAcc->Draw();

  TCanvas *c5c =new TCanvas("c5c","c5c",10,10,600,600);
  c5c->SetGridy();
  TH1F *fHistPtTPCInAccSecFrac = (TH1F*)fHistPtTPCInAccS->Clone("fHistPtTPCInAccSecFrac");  
  fHistPtTPCInAccSecFrac->Divide(fHistPtTPCInAcc);
  fHistPtTPCInAccSecFrac->SetLineColor(kOrange+7);
  fHistPtTPCInAccSecFrac->SetTitle("Fraction of secondaries");
  fHistPtTPCInAccSecFrac->SetYTitle("sec/all");
  fHistPtTPCInAccSecFrac->Draw();
  TH1F *fHistPtITSMI2InAccSecFrac = (TH1F*)fHistPtITSMI2InAccS->Clone("fHistPtITSMI2InAccSecFrac");  
  fHistPtITSMI2InAccSecFrac->Divide(fHistPtITSMI2InAcc);
  fHistPtITSMI2InAccSecFrac->SetLineColor(7);
  fHistPtITSMI2InAccSecFrac->Draw("same"); 
  TH1F *fHistPtITSMIge2InAccSecFrac = (TH1F*)fHistPtITSMIge2InAccS->Clone("fHistPtITSMIge2InAccSecFrac");  
  fHistPtITSMIge2InAccSecFrac->Divide(fHistPtITSMIge2InAcc);
  fHistPtITSMIge2InAccSecFrac->Draw("same");
  TH1F *fHistPtITSMI3InAccSecFrac = (TH1F*)fHistPtITSMI3InAccS->Clone("fHistPtITSMI3InAccSecFrac");  
  fHistPtITSMI3InAccSecFrac->Divide(fHistPtITSMI3InAcc);
  fHistPtITSMI3InAccSecFrac->SetLineColor(6);
  fHistPtITSMI3InAccSecFrac->Draw("same");
  TH1F *fHistPtITSMI4InAccSecFrac = (TH1F*)fHistPtITSMI4InAccS->Clone("fHistPtITSMI4InAccSecFrac");  
  fHistPtITSMI4InAccSecFrac->Divide(fHistPtITSMI4InAcc);
  fHistPtITSMI4InAccSecFrac->SetLineColor(4);
  fHistPtITSMI4InAccSecFrac->Draw("same");
  TH1F *fHistPtITSMI5InAccSecFrac = (TH1F*)fHistPtITSMI5InAccS->Clone("fHistPtITSMI5InAccSecFrac");  
  fHistPtITSMI5InAccSecFrac->Divide(fHistPtITSMI5InAcc);
  fHistPtITSMI5InAccSecFrac->SetLineColor(3);
  fHistPtITSMI5InAccSecFrac->Draw("same");
  TH1F *fHistPtITSMI6InAccSecFrac = (TH1F*)fHistPtITSMI6InAccS->Clone("fHistPtITSMI6InAccSecFrac");  
  fHistPtITSMI6InAccSecFrac->Divide(fHistPtITSMI6InAcc);
  fHistPtITSMI6InAccSecFrac->SetLineColor(2);
  fHistPtITSMI6InAccSecFrac->SetLineColor(2);
  fHistPtITSMI6InAccSecFrac->Draw("same");
  TH1F *fHistPtITSMISPDInAccSecFrac = (TH1F*)fHistPtITSMISPDInAccS->Clone("fHistPtITSMISPDInAccSecFrac");  
  fHistPtITSMISPDInAccSecFrac->Divide(fHistPtITSMISPDInAcc);
  fHistPtITSMISPDInAccSecFrac->SetLineColor(9);
  fHistPtITSMISPDInAccSecFrac->Draw("same");
  TH1F *fHistPtITSMIoneSPDInAccSecFrac = (TH1F*)fHistPtITSMIoneSPDInAccS->Clone("fHistPtITSMIoneSPDInAccSecFrac");  
  fHistPtITSMIoneSPDInAccSecFrac->Divide(fHistPtITSMIoneSPDInAcc);
  fHistPtITSMIoneSPDInAccSecFrac->SetLineColor(15);
  fHistPtITSMIoneSPDInAccSecFrac->Draw("same");
  TH1F *fHistPtITSTPCselSecFrac = (TH1F*)fHistPtITSTPCselS->Clone("fHistPtITSTPCselSecFrac");  
  fHistPtITSTPCselSecFrac->Divide(fHistPtITSTPCsel);
  fHistPtITSTPCselSecFrac->SetLineColor(15);
  fHistPtITSTPCselSecFrac->Draw("same");


  TLegend *l3=new TLegend(0.5,0.5,0.9,0.9);
  TLegend *l4=new TLegend(0.5,0.5,0.9,0.9);
  TCanvas *c5 =new TCanvas("c5","c5",10,10,1200,600);
  c5->Divide(2,1);
  c5_1->SetGridy();
  c5_2->SetGridy();
  c5_1->SetLogx();
  c5_2->SetLogx();
  c5->cd(1);
  fHistPtITSMIge2InAcc->SetMaximum(1.5);
  fHistPtITSMIge2InAcc->SetMinimum(0);
  fHistPtITSMIge2InAcc->SetTitle("Fraction of prolonged tracks with N ITS points");
  fHistPtITSMIge2InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIge2InAcc->Draw();
  l3->AddEntry(fHistPtITSMIge2InAcc,">=2 cls","l");
  fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI6InAcc->SetLineColor(2);
  l3->AddEntry(fHistPtITSMI6InAcc,"6 cls","l");
  fHistPtITSMI6InAcc->Draw("same");
  fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI5InAcc->SetLineColor(3);
  l3->AddEntry(fHistPtITSMI5InAcc,"5 cls","l");
  fHistPtITSMI5InAcc->Draw("same");
  fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI4InAcc->SetLineColor(4);
  l3->AddEntry(fHistPtITSMI4InAcc,"4 cls","l");
  fHistPtITSMI4InAcc->Draw("same");
  fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI3InAcc->SetLineColor(6);
  l3->AddEntry(fHistPtITSMI3InAcc,"3 cls","l");
  fHistPtITSMI3InAcc->Draw("same");
  fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMI2InAcc->SetLineColor(7);
  l3->AddEntry(fHistPtITSMI2InAcc,"2 cls","l");
  fHistPtITSMI2InAcc->Draw("same");
  fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMISPDInAcc->SetLineColor(9);
  l3->AddEntry(fHistPtITSMISPDInAcc,"2SPD + any","l");
  fHistPtITSMISPDInAcc->Draw("same");
  fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIoneSPDInAcc->SetLineColor(15);
  l3->AddEntry(fHistPtITSMIoneSPDInAcc,">=1SPD + any","l");
  fHistPtITSMIoneSPDInAcc->Draw("same");
  //fHistPtITSMIoneSPDthreeSDDSSDInAcc->Divide(fHistPtITSMIoneSPDthreeSDDSSDInAcc,fHistPtTPCInAcc,1,1,"B");
  //fHistPtITSMIoneSPDthreeSDDSSDInAcc->SetLineColor(kOrange);
  //l3->AddEntry(fHistPtITSMIoneSPDthreeSDDSSDInAcc,">=1SPD + >=3SDDSSD","l");
  //fHistPtITSMIoneSPDthreeSDDSSDInAcc->Draw("same");
  fHistPtITSTPCsel->Divide(fHistPtITSTPCsel,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSTPCsel->SetLineColor(kAzure+1);
  l3->AddEntry(fHistPtITSTPCsel,">=1SPD + any + d_{0} cut","l");
  fHistPtITSTPCsel->Draw("same");
  fHistPtITSMIge2InAcc->Draw("same");
  l3->Draw();
  if(ioValues) {
    Int_t index,ptbin;
    ptbin=fHistPtITSMIge2InAcc->FindBin(0.201);
    index=2;
    ioValues[index]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(1.001);
    index=3;
    ioValues[index]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(10.001);
    index=4;
    ioValues[index]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(0.201);
    index=5;
    ioValues[index]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(1.001);
    index=6;
    ioValues[index]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(10.001);
    index=7;
    ioValues[index]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(0.201);
    index=8;
    ioValues[index]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(1.001);
    index=9;
    ioValues[index]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(10.001);
    index=10;
    ioValues[index]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(0.201);
    index=11;
    ioValues[index]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(1.001);
    index=12;
    ioValues[index]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(10.001);
    index=13;
    ioValues[index]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[index]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
  }
  c5->cd(2);
  TH1F *fHistPtITSMIokbadoutinzge4InAcc = (TH1F*)fHistPtITSMIokbadoutinz6InAcc->Clone("fHistPtITSMIokbadoutinzge4InAcc");
  fHistPtITSMIokbadoutinzge4InAcc->Add(fHistPtITSMIokbadoutinz5InAcc);
  fHistPtITSMIokbadoutinzge4InAcc->Add(fHistPtITSMIokbadoutinz4InAcc);
  fHistPtITSMIokbadoutinzge4InAcc->SetMaximum(1.5);
  fHistPtITSMIokbadoutinzge4InAcc->SetMinimum(0);
  fHistPtITSMIokbadoutinzge4InAcc->SetTitle("Fraction of prolonged tracks with N ITS layers \"ok\"");
  fHistPtITSMIokbadoutinzge4InAcc->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIokbadoutinzge4InAcc->Divide(fHistPtITSMIokbadoutinzge4InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIokbadoutinzge4InAcc->SetLineColor(1);
  fHistPtITSMIokbadoutinzge4InAcc->Draw();
  fHistPtITSMIokbadoutinz6InAcc->Divide(fHistPtITSMIokbadoutinz6InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIokbadoutinz6InAcc->SetLineColor(2);
  fHistPtITSMIokbadoutinz6InAcc->Draw("same");
  fHistPtITSMIokbadoutinz5InAcc->Divide(fHistPtITSMIokbadoutinz5InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIokbadoutinz5InAcc->SetLineColor(3);
  fHistPtITSMIokbadoutinz5InAcc->Draw("same");
  fHistPtITSMIokbadoutinz4InAcc->Divide(fHistPtITSMIokbadoutinz4InAcc,fHistPtTPCInAcc,1,1,"B");
  fHistPtITSMIokbadoutinz4InAcc->SetLineColor(4);
  fHistPtITSMIokbadoutinz4InAcc->Draw("same");
  fHistPtITSMIokbadoutinzge4InAcc->Draw("same");
  l4->AddEntry(fHistPtITSMIokbadoutinzge4InAcc,">=4 layers","l");
  l4->AddEntry(fHistPtITSMIokbadoutinz6InAcc,"6 layers","l");
  l4->AddEntry(fHistPtITSMIokbadoutinz5InAcc,"5 layers","l");
  l4->AddEntry(fHistPtITSMIokbadoutinz4InAcc,"4 layers","l");
  l4->Draw();

  if(fHistPtITSMIge2InAccFake) {
    TCanvas *c5f =new TCanvas("c5f","c5f",10,10,600,600);
    c5f->SetGridy();
    c5f->SetLogx();
    fHistPtITSMIge2InAccFake->SetMaximum(1.5);
    fHistPtITSMIge2InAccFake->SetMinimum(0);
    fHistPtITSMIge2InAccFake->SetTitle("Fraction of fake tracks with N ITS points");
    fHistPtITSMIge2InAccFake->SetYTitle("Fraction of fakes");
    fHistPtITSMIge2InAccFake->Draw();
    fHistPtITSMI6InAccFake->SetLineColor(2);
    fHistPtITSMI6InAccFake->Draw("same");
    fHistPtITSMI5InAccFake->SetLineColor(3);
    //fHistPtITSMI5InAccFake->Draw("same");
    fHistPtITSMI4InAccFake->SetLineColor(4);
    //fHistPtITSMI4InAccFake->Draw("same");
    fHistPtITSMI3InAccFake->SetLineColor(6);
    //fHistPtITSMI3InAccFake->Draw("same");
    fHistPtITSMI2InAccFake->SetLineColor(7);
    //fHistPtITSMI2InAccFake->Draw("same");
    fHistPtITSMISPDInAccFake->SetLineColor(9);
    fHistPtITSMISPDInAccFake->Draw("same");
    fHistPtITSMIoneSPDInAccFake->SetLineColor(15);
    fHistPtITSMIoneSPDInAccFake->Draw("same");
    if(fHistPtITSTPCselFake) {
      fHistPtITSTPCselFake->SetLineColor(kAzure+1);
      fHistPtITSTPCselFake->Draw("same");
    }
    fHistPtITSMIge2InAccFake->Draw("same");
    l3->Draw();
  }


  TLegend *l4d=new TLegend(0.5,0.5,0.9,0.9);
  TCanvas *c5d =new TCanvas("c5d","c5d",10,10,600,600);
  c5d->SetGridy();
  c5d->SetGridy();
  fHistPtITSMISPDInAccMC->SetMaximum(1.5);
  fHistPtITSMISPDInAccMC->SetMinimum(0);
  fHistPtITSMISPDInAccMC->SetTitle("Fraction of prolonged tracks with N ITS points");
  fHistPtITSMISPDInAccMC->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMISPDInAccMC->Divide(fHistPtITSMISPDInAccMC,fHistPtTPCInAccMCtwoSPD,1,1,"B");
  fHistPtITSMISPDInAccMC->SetLineColor(9);
  l4d->AddEntry(fHistPtITSMISPDInAccMC,"2SPD + any","l");
  fHistPtITSMISPDInAccMC->Draw();
  fHistPtITSMIoneSPDInAccMC->Divide(fHistPtITSMIoneSPDInAccMC,fHistPtTPCInAccMConeSPD,1,1,"B");
  fHistPtITSMIoneSPDInAccMC->SetLineColor(15);
  l4d->AddEntry(fHistPtITSMIoneSPDInAccMC,">=1SPD + any","l");
  fHistPtITSMIoneSPDInAccMC->Draw("same");
  l4d->Draw();

  TFile *eff=new TFile("eff.root","recreate");
  fHistPtITSMIge2InAcc->Write();
  fHistPtITSMI6InAcc->Write();
  fHistPtITSMI5InAcc->Write();
  fHistPtITSMI4InAcc->Write();
  fHistPtITSMI3InAcc->Write();
  fHistPtITSMI2InAcc->Write();
  fHistPtITSMISPDInAcc->Write();
  fHistPtITSMIoneSPDInAcc->Write();
  //fHistPtITSMIoneSPDthreeSDDSSDInAcc->Write();
  fHistPtITSTPCsel->Write();
  fHistNclsITSMI->Write();
  fHistClusterMapModuleITSMIokRatioInAcc->Write();
  fHistClusterMapModuleITSMIbadRatioInAcc->Write();
  fHistClusterMapModuleITSMInoclsRatioInAcc->Write();
  eff->Close();

  c5c->cd();
  l3->Draw();


  TCanvas *c5b =new TCanvas("c5b","c5b",10,10,600,600);
  c5b->SetLogy();
  fHistRProdVtxInAccP->SetLineColor(2);
  fHistRProdVtxInAccP->SetMinimum(1);
  fHistRProdVtxInAccS->SetMinimum(1);
  fHistRProdVtxInAccP->Draw();
  fHistRProdVtxInAccS->Draw("sames");

  TCanvas *c5a =new TCanvas("c5a","c5a",10,10,1200,600);
  c5a->Divide(2,1);
  c5a_1->SetLogx();
  c5a_2->SetLogx();
  c5a_1->SetGridy();
  c5a_2->SetGridy();
  c5a->cd(1);
  fHistPtITSMIge2InAccP->SetMaximum(1.5);
  fHistPtITSMIge2InAccP->SetMinimum(0);
  fHistPtITSMIge2InAccP->SetTitle("Fraction of prolonged tracks with N ITS points");
  fHistPtITSMIge2InAccP->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIge2InAccP->Draw();
  fHistPtITSMIge2InAccP->Divide(fHistPtITSMIge2InAccP,fHistPtTPCInAccP,1,1,"B");
  fHistPtITSMI6InAccP->Divide(fHistPtITSMI6InAccP,fHistPtTPCInAccP,1,1,"B");
  fHistPtITSMI6InAccP->SetLineColor(2);
  fHistPtITSMI6InAccP->Draw("same");
  fHistPtITSMI5InAccP->Divide(fHistPtITSMI5InAccP,fHistPtTPCInAccP,1,1,"B");
  fHistPtITSMI5InAccP->SetLineColor(3);
  fHistPtITSMI5InAccP->Draw("same");
  fHistPtITSMI4InAccP->Divide(fHistPtITSMI4InAccP,fHistPtTPCInAccP,1,1,"B");
  fHistPtITSMI4InAccP->SetLineColor(4);
  fHistPtITSMI4InAccP->Draw("same");
  fHistPtITSMI3InAccP->Divide(fHistPtITSMI3InAccP,fHistPtTPCInAccP,1,1,"B");
  fHistPtITSMI3InAccP->SetLineColor(6);
  fHistPtITSMI3InAccP->Draw("same");
  fHistPtITSMI2InAccP->Divide(fHistPtITSMI2InAccP,fHistPtTPCInAccP,1,1,"B");
  fHistPtITSMI2InAccP->SetLineColor(7);
  fHistPtITSMI2InAccP->Draw("same");
  fHistPtITSMISPDInAccP->Divide(fHistPtITSMISPDInAccP,fHistPtTPCInAccP,1,1,"B");
  fHistPtITSMISPDInAccP->SetLineColor(9);
  fHistPtITSMISPDInAccP->Draw("same");
  fHistPtITSMIoneSPDInAccP->Divide(fHistPtITSMIoneSPDInAccP,fHistPtTPCInAccP,1,1,"B");
  fHistPtITSMIoneSPDInAccP->SetLineColor(15);
  fHistPtITSMIoneSPDInAccP->Draw("same");
  fHistPtITSTPCselP->Divide(fHistPtITSTPCselP,fHistPtTPCInAccP,1,1,"B");
  fHistPtITSTPCselP->SetLineColor(kAzure+1);
  fHistPtITSTPCselP->Draw("same");
  fHistPtITSMIge2InAccP->Draw("same");
  l3->Draw();
  c5a->cd(2);
  fHistPtITSMIge2InAccS->SetMaximum(1.5);
  fHistPtITSMIge2InAccS->SetMinimum(0);
  fHistPtITSMIge2InAccS->SetTitle("Fraction of prolonged tracks with N ITS points");
  fHistPtITSMIge2InAccS->SetYTitle("ITS+TPC / TPC");
  fHistPtITSMIge2InAccS->Divide(fHistPtITSMIge2InAccS,fHistPtTPCInAccS,1,1,"B");
  fHistPtITSMIge2InAccS->Draw();
  fHistPtITSMI6InAccS->Divide(fHistPtITSMI6InAccS,fHistPtTPCInAccS,1,1,"B");
  fHistPtITSMI6InAccS->SetLineColor(2);
  fHistPtITSMI6InAccS->Draw("same");
  fHistPtITSMI5InAccS->Divide(fHistPtITSMI5InAccS,fHistPtTPCInAccS,1,1,"B");
  fHistPtITSMI5InAccS->SetLineColor(3);
  fHistPtITSMI5InAccS->Draw("same");
  fHistPtITSMI4InAccS->Divide(fHistPtITSMI4InAccS,fHistPtTPCInAccS,1,1,"B");
  fHistPtITSMI4InAccS->SetLineColor(4);
  fHistPtITSMI4InAccS->Draw("same");
  fHistPtITSMI3InAccS->Divide(fHistPtITSMI3InAccS,fHistPtTPCInAccS,1,1,"B");
  fHistPtITSMI3InAccS->SetLineColor(6);
  fHistPtITSMI3InAccS->Draw("same");
  fHistPtITSMI2InAccS->Divide(fHistPtITSMI2InAccS,fHistPtTPCInAccS,1,1,"B");
  fHistPtITSMI2InAccS->SetLineColor(7);
  fHistPtITSMI2InAccS->Draw("same");
  fHistPtITSMISPDInAccS->Divide(fHistPtITSMISPDInAccS,fHistPtTPCInAccS,1,1,"B");
  fHistPtITSMISPDInAccS->SetLineColor(9);
  fHistPtITSMISPDInAccS->Draw("same");
  fHistPtITSMIoneSPDInAccS->Divide(fHistPtITSMIoneSPDInAccS,fHistPtTPCInAccS,1,1,"B");
  fHistPtITSMIoneSPDInAccS->SetLineColor(15);
  fHistPtITSMIoneSPDInAccS->Draw("same");
  fHistPtITSTPCselS->Divide(fHistPtITSTPCselS,fHistPtTPCInAccS,1,1,"B");
  fHistPtITSTPCselS->SetLineColor(kAzure+1);
  fHistPtITSTPCselS->Draw("same");
  fHistPtITSMIge2InAccS->Draw("same");
  l3->Draw();


  if(!plotAlignmentChecks) return kTRUE;

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
    if(pt<.01) continue;
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
  hSPDTrackletsLRdxy->Draw("sames");
  l6->AddEntry(hSPDTrackletsTBdxy,"top-bottom","l");
  l6->AddEntry(hSPDTrackletsLRdxy,"left-right","l");
  l6->Draw();
  c6->cd(2);
  hSPDTrackletsTBdz->SetLineColor(4);
  hSPDTrackletsTBdz->Draw();
  hSPDTrackletsLRdz->SetLineColor(2);
  hSPDTrackletsLRdz->Draw("sames");
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
    if(pt<0.1) continue;
    if(layer!=0 && layer!=1) continue;
    if(npoints<3) continue;
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


  return kTRUE;
}
//-----------------------------------------------------------------------------
void PlotEffRatio() {

  TFile *f= new TFile("eff.root");

  TH1F *fHistPtITSMI6InAcc = (TH1F*)f->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAcc = (TH1F*)f->Get("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAcc = (TH1F*)f->Get("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAcc = (TH1F*)f->Get("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAcc = (TH1F*)f->Get("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMIge2InAcc = (TH1F*)f->Get("fHistPtITSMIge2InAcc");
  TH1F *fHistPtITSMISPDInAcc = (TH1F*)f->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)f->Get("fHistPtITSMIoneSPDInAcc");
  TH1F *fHistPtITSMIoneSPDthreeSDDSSDInAcc = (TH1F*)f->Get("fHistPtITSMIoneSPDthreeSDDSSDInAcc");
  TH1F *fHistPtITSTPCsel = (TH1F*)f->Get("fHistPtITSTPCsel");
  TH1F *fHistClusterMapModuleITSMIokRatioInAcc = (TH1F*)f->Get("fHistClusterMapModuleITSMIokRatioInAcc");
  TH1F *fHistNclsITSMI = (TH1F*)f->Get("fHistNclsITSMI");


  TFile *fMC= new TFile("../MC/MCLHC11b10a/eff.root");

  TH1F *fHistPtITSMI6InAccMC = (TH1F*)fMC->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAccMC = (TH1F*)fMC->Get("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAccMC = (TH1F*)fMC->Get("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAccMC = (TH1F*)fMC->Get("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAccMC = (TH1F*)fMC->Get("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMIge2InAccMC = (TH1F*)fMC->Get("fHistPtITSMIge2InAcc");
  TH1F *fHistPtITSMISPDInAccMC = (TH1F*)fMC->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAccMC = (TH1F*)fMC->Get("fHistPtITSMIoneSPDInAcc");
  TH1F *fHistPtITSMIoneSPDthreeSDDSSDInAccMC = (TH1F*)fMC->Get("fHistPtITSMIoneSPDthreeSDDSSDInAcc");
  TH1F *fHistPtITSTPCselMC = (TH1F*)fMC->Get("fHistPtITSTPCsel");
  TH1F *fHistClusterMapModuleITSMIokRatioInAccMC = (TH1F*)fMC->Get("fHistClusterMapModuleITSMIokRatioInAcc");
  TH1F *fHistNclsITSMIMC = (TH1F*)fMC->Get("fHistNclsITSMIMC");

  /*
  TCanvas *c0 = new TCanvas("c0","c0",0,0,600,500);
  c0->SetGridy();
  fHistNclsITSMI->SetLineColor(4);
  fHistNclsITSMI->SetMarkerColor(4);
  fHistNclsITSMI->SetLineStyle(1);
  fHistNclsITSMI->SetMarkerStyle(20);
  fHistNclsITSMI->Draw();
  fHistNclsITSMIMC->SetLineColor(4);
  fHistNclsITSMIMC->SetMarkerColor(4);
  fHistNclsITSMIMC->SetLineStyle(2);
  fHistNclsITSMIMC->SetMarkerStyle(24);
  fHistNclsITSMIMC->Draw("same");
  */
  TCanvas *c = new TCanvas("c","c",0,0,600,500);
  c->SetGridy();
  c->SetLogx();
  
  fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAccMC);
  fHistPtITSMI6InAcc->Draw();
  fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAccMC);
  fHistPtITSMI5InAcc->Draw("same");
  fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAccMC);
  fHistPtITSMI4InAcc->Draw("same");
  fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAccMC);
  fHistPtITSMI3InAcc->Draw("same");
  fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAccMC);
  fHistPtITSMI2InAcc->Draw("same");
  fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAccMC);
  fHistPtITSMIge2InAcc->Draw("same");
  fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAccMC);
  fHistPtITSMISPDInAcc->Draw("same");
  fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAccMC);
  fHistPtITSMIoneSPDInAcc->Draw("same");
  //fHistPtITSMIoneSPDthreeSDDSSDInAcc->Divide(fHistPtITSMIoneSPDthreeSDDSSDInAccMC);
  //fHistPtITSMIoneSPDthreeSDDSSDInAcc->Draw("same");
  
  fHistPtITSTPCsel->Divide(fHistPtITSTPCselMC);
  fHistPtITSTPCsel->Draw("same");
  
  TLegend *l3=new TLegend(0.5,0.5,0.9,0.9);
  l3->AddEntry(fHistPtITSMIge2InAcc,">=2 cls","l");
  l3->AddEntry(fHistPtITSMI6InAcc,"6 cls","l");
  l3->AddEntry(fHistPtITSMI5InAcc,"5 cls","l");
  l3->AddEntry(fHistPtITSMI4InAcc,"4 cls","l");
  l3->AddEntry(fHistPtITSMI3InAcc,"3 cls","l");
  l3->AddEntry(fHistPtITSMI2InAcc,"2 cls","l");
  l3->AddEntry(fHistPtITSMISPDInAcc,"2SPD + any","l");
  l3->AddEntry(fHistPtITSMIoneSPDInAcc,">=1SPD + any","l");
  //l3->AddEntry(fHistPtITSMIoneSPDthreeSDDSSDInAcc,">=1SPD + >=3SDDSSD","l");
  l3->AddEntry(fHistPtITSTPCsel,">=1SPD + any + d_{0} cut","l");
  l3->Draw();
  
  TCanvas *cc = new TCanvas("cc","cc",0,0,600,500);
  cc->Divide(1,2);
  cc->cd(1);
  fHistClusterMapModuleITSMIokRatioInAccMC->SetLineColor(4);
  fHistClusterMapModuleITSMIokRatioInAccMC->Draw();
  fHistClusterMapModuleITSMIokRatioInAcc->Draw("same");
  cc->cd(2);
  for(Int_t i=1;i<=fHistClusterMapModuleITSMIokRatioInAccMC->GetNbinsX();i++) {
    if(fHistClusterMapModuleITSMIokRatioInAccMC->GetBinContent(i)<0.02) {
      fHistClusterMapModuleITSMIokRatioInAccMC->SetBinContent(i,1); 
      if(fHistClusterMapModuleITSMIokRatioInAcc->GetBinContent(i)<0.02) {
	fHistClusterMapModuleITSMIokRatioInAcc->SetBinContent(i,1); 
      }
    }
  }
  TH1F* fHistClusterMapDataOverMC=(TH1F*)fHistClusterMapModuleITSMIokRatioInAcc->Clone("fHistClusterMapDataOverMC");
  fHistClusterMapDataOverMC->Divide(fHistClusterMapModuleITSMIokRatioInAccMC);
  fHistClusterMapDataOverMC->Draw();

  return;
}

//-----------------------------------------------------------------------------
void PlotEffOfficial(Bool_t drawRatio=kTRUE) {

  gStyle->SetOptStat(0);

  gROOT->LoadMacro("/Users/dainesea/ALICEWorkInProgress.C");

  TFile *f= new TFile("eff_115328_115393_minipass.root");

  TH1F *fHistPtITSMI6InAcc = (TH1F*)f->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAcc = (TH1F*)f->Get("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAcc = (TH1F*)f->Get("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAcc = (TH1F*)f->Get("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAcc = (TH1F*)f->Get("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMIge2InAcc = (TH1F*)f->Get("fHistPtITSMIge2InAcc");
  fHistPtITSMIge2InAcc->SetLineColor(1);
  fHistPtITSMIge2InAcc->SetLineStyle(1);
  fHistPtITSMIge2InAcc->SetMarkerColor(1);
  fHistPtITSMIge2InAcc->SetMarkerStyle(21);
  TH1F *fHistPtITSMISPDInAcc = (TH1F*)f->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)f->Get("fHistPtITSMIoneSPDInAcc");
  fHistPtITSMIoneSPDInAcc->SetLineColor(2);
  fHistPtITSMIoneSPDInAcc->SetLineStyle(1);
  fHistPtITSMIoneSPDInAcc->SetMarkerColor(2);
  fHistPtITSMIoneSPDInAcc->SetMarkerStyle(20);
  TH1F *fHistPtITSTPCsel = (TH1F*)f->Get("fHistPtITSTPCsel");
  TH1F *fHistClusterMapModuleITSMIokRatioInAcc = (TH1F*)f->Get("fHistClusterMapModuleITSMIokRatioInAcc");
  TH1F *fHistNclsITSMI = (TH1F*)f->Get("fHistNclsITSMI");


  TFile *fMC= new TFile("eff_lhc10b2.root");

  TH1F *fHistPtITSMI6InAccMC = (TH1F*)fMC->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMI5InAccMC = (TH1F*)fMC->Get("fHistPtITSMI5InAcc");
  TH1F *fHistPtITSMI4InAccMC = (TH1F*)fMC->Get("fHistPtITSMI4InAcc");
  TH1F *fHistPtITSMI3InAccMC = (TH1F*)fMC->Get("fHistPtITSMI3InAcc");
  TH1F *fHistPtITSMI2InAccMC = (TH1F*)fMC->Get("fHistPtITSMI2InAcc");
  TH1F *fHistPtITSMIge2InAccMC = (TH1F*)fMC->Get("fHistPtITSMIge2InAcc");
  fHistPtITSMIge2InAccMC->SetLineColor(1);
  fHistPtITSMIge2InAccMC->SetLineStyle(2);
  fHistPtITSMIge2InAccMC->SetMarkerColor(1);
  fHistPtITSMIge2InAccMC->SetMarkerStyle(25);
  TH1F *fHistPtITSMISPDInAccMC = (TH1F*)fMC->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAccMC = (TH1F*)fMC->Get("fHistPtITSMIoneSPDInAcc");
  fHistPtITSMIoneSPDInAccMC->SetLineColor(2);
  fHistPtITSMIoneSPDInAccMC->SetLineStyle(2);
  fHistPtITSMIoneSPDInAccMC->SetMarkerColor(2);
  fHistPtITSMIoneSPDInAccMC->SetMarkerStyle(24);
  TH1F *fHistPtITSTPCselMC = (TH1F*)fMC->Get("fHistPtITSTPCsel");
  TH1F *fHistClusterMapModuleITSMIokRatioInAccMC = (TH1F*)fMC->Get("fHistClusterMapModuleITSMIokRatioInAcc");
  TH1F *fHistNclsITSMIMC = (TH1F*)fMC->Get("fHistNclsITSMI");


  TCanvas *c0 = new TCanvas("c0","c0",0,0,600,500);
  c0->SetGridy();
  fHistNclsITSMI->SetLineColor(4);
  fHistNclsITSMI->SetMarkerColor(4);
  fHistNclsITSMI->SetLineStyle(1);
  fHistNclsITSMI->SetMarkerStyle(20);
  fHistNclsITSMI->Scale(1./fHistNclsITSMI->Integral());
  fHistNclsITSMI->Draw();
  fHistNclsITSMIMC->SetLineColor(4);
  fHistNclsITSMIMC->SetMarkerColor(4);
  fHistNclsITSMIMC->SetLineStyle(2);
  fHistNclsITSMIMC->SetMarkerStyle(24);
  fHistNclsITSMIMC->Scale(1./fHistNclsITSMIMC->Integral());
  fHistNclsITSMIMC->Draw("same");

  TLegend *l0=new TLegend(0.5,0.5,0.9,0.9);
  l0->SetFillStyle(0);
  l0->SetBorderSize(0);
  l0->AddEntry(fHistNclsITSMI,"Data","lp");
  l0->AddEntry(fHistNclsITSMIMC,"MC","lp");
  l0->Draw();
 
  ALICEWorkInProgress(c0,"20/05/2010");

  TCanvas *c = new TCanvas("c","c",0,0,600,500);
  c->SetGridy();
  c->SetLogx();

  TLegend *l3=new TLegend(0.5,0.5,0.9,0.9);
  l3->SetFillStyle(0);
  l3->SetBorderSize(0);
  
  if(drawRatio) {
    //fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAccMC);
    //fHistPtITSMI6InAcc->Draw();
    //fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAccMC);
    //fHistPtITSMI5InAcc->Draw("same");
    //fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAccMC);
    //fHistPtITSMI4InAcc->Draw("same");
    //fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAccMC);
    //fHistPtITSMI3InAcc->Draw("same");
    //fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAccMC);
    //fHistPtITSMI2InAcc->Draw("same");
    fHistPtITSMIge2InAcc->SetYTitle("Data efficiency / MC efficiency");
    fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAccMC);
    fHistPtITSMIge2InAcc->Draw();
    //fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAccMC);
    //fHistPtITSMISPDInAcc->Draw("same");
    fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAccMC);
    fHistPtITSMIoneSPDInAcc->Draw("same");
    
    //fHistPtITSTPCsel->Divide(fHistPtITSTPCselMC);
    //fHistPtITSTPCsel->Draw("same");
    l3->AddEntry(fHistPtITSMIge2InAcc,"at least 2 ITS hits","lp");
    //l3->AddEntry(fHistPtITSMI6InAcc,"6 cls","l");
    //l3->AddEntry(fHistPtITSMI5InAcc,"5 cls","l");
    //l3->AddEntry(fHistPtITSMI4InAcc,"4 cls","l");
    //l3->AddEntry(fHistPtITSMI3InAcc,"3 cls","l");
    //l3->AddEntry(fHistPtITSMI2InAcc,"2 cls","l");
    //l3->AddEntry(fHistPtITSMISPDInAcc,"2 SPD hits","lp");
    l3->AddEntry(fHistPtITSMIoneSPDInAcc,"at least 1 SPD hit","lp");
    //l3->AddEntry(fHistPtITSTPCsel,">=1SPD + any + d_{0} cut","l");

  } else {
    fHistPtITSMIge2InAcc->SetYTitle("ITS prolongation efficiency");
    fHistPtITSMIge2InAcc->Draw();
    fHistPtITSMIge2InAccMC->Draw("same");
    fHistPtITSMIoneSPDInAcc->Draw("same");
    fHistPtITSMIoneSPDInAccMC->Draw("same");
    l3->AddEntry(fHistPtITSMIge2InAcc,"at least 2 ITS hits (Data)","lp");
    l3->AddEntry(fHistPtITSMIoneSPDInAcc,"at least 1 SPD hit (Data)","lp");
    l3->AddEntry(fHistPtITSMIge2InAccMC,"at least 2 ITS hits (MC)","lp");
    l3->AddEntry(fHistPtITSMIoneSPDInAccMC,"at least 1 SPD hit (MC)","lp");
  } 

  l3->Draw();
  
  ALICEWorkInProgress(c,"20/05/2010");

  return;
}

void PlotImpPar_rphi(Int_t rebin=1) {


  TFile *fMC= new TFile("AnalysisResults_onlynonfakes.root");

  TList *list=(TList*)fMC->Get("cOutputITS");
  TDirectoryFile *dir=0;
  if(!list) {
    dir=(TDirectoryFile*)fMC->GetDirectory("ITS_Performance");
    if(dir) list = (TList*)dir->Get("cOutputITS");
  }

  if(!list) return kFALSE;

  TH1F *fHistd0rphiITSMIoneSPDInAccP150200MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP150200");
  TH1F *fHistd0rphiITSMIoneSPDInAccS150200MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS150200");
  TH1F *fHistd0rphiITSMIoneSPDInAccS150200fromStrangeMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS150200fromStrange");
  TH1F *fHistd0rphiITSMIoneSPDInAccS150200fromMatMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS150200fromMat");
  TH1F *fHistd0rphiITSMIoneSPDInAccP350450MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP350450");
  TH1F *fHistd0rphiITSMIoneSPDInAccS350450MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS350450");
  TH1F *fHistd0rphiITSMIoneSPDInAccS350450fromStrangeMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS350450fromStrange");
  TH1F *fHistd0rphiITSMIoneSPDInAccS350450fromMatMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS350450fromMat");
  TH1F *fHistd0rphiITSMIoneSPDInAccP500700MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP500700");
  TH1F *fHistd0rphiITSMIoneSPDInAccS500700MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS500700");
  TH1F *fHistd0rphiITSMIoneSPDInAccS500700fromStrangeMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS500700fromStrange");
  TH1F *fHistd0rphiITSMIoneSPDInAccS500700from211MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS500700from211");
  TH1F *fHistd0rphiITSMIoneSPDInAccS500700from22MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS500700from22");
  TH1F *fHistd0rphiITSMIoneSPDInAccS500700from310MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS500700from310");
  TH1F *fHistd0rphiITSMIoneSPDInAccS500700from321MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS500700from321");
  TH1F *fHistd0rphiITSMIoneSPDInAccS500700from3122MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS500700from3122");
  TH1F *fHistd0rphiITSMIoneSPDInAccS500700fromMatMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS500700fromMat");
  TH1F *fHistd0rphiITSMIoneSPDInAccP10001500MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP10001500");
  TH1F *fHistd0rphiITSMIoneSPDInAccS10001500MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS10001500");
  TH1F *fHistd0rphiITSMIoneSPDInAccS10001500fromStrangeMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS10001500fromStrange");
  TH1F *fHistd0rphiITSMIoneSPDInAccS10001500fromMatMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS10001500fromMat");
  TH1F *fHistd0rphiITSMIoneSPDInAccP25004000MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP25004000");
  TH1F *fHistd0rphiITSMIoneSPDInAccS25004000MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS25004000");
  TH1F *fHistd0rphiITSMIoneSPDInAccS25004000fromStrangeMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS25004000fromStrange");
  TH1F *fHistd0rphiITSMIoneSPDInAccS25004000fromMatMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS25004000fromMat");
  TH1F *fHistd0rphiITSMIoneSPDInAccP40008000MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP40008000");
  TH1F *fHistd0rphiITSMIoneSPDInAccS40008000MC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS40008000");
  TH1F *fHistd0rphiITSMIoneSPDInAccS40008000fromStrangeMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS40008000fromStrange");
  TH1F *fHistd0rphiITSMIoneSPDInAccS40008000fromMatMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS40008000fromMat"); 

  TH1F *fHistd0rphiITSMIoneSPDInAcc150200MC=(TH1F*)fHistd0rphiITSMIoneSPDInAccP150200MC->Clone("fHistd0rphiITSMIoneSPDInAcc150200MC");
  fHistd0rphiITSMIoneSPDInAcc150200MC->Add(fHistd0rphiITSMIoneSPDInAccS150200MC);
  fHistd0rphiITSMIoneSPDInAcc150200MC->Scale(1./fHistd0rphiITSMIoneSPDInAcc150200MC->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc350450MC=(TH1F*)fHistd0rphiITSMIoneSPDInAccP350450MC->Clone("fHistd0rphiITSMIoneSPDInAcc350450MC");
  fHistd0rphiITSMIoneSPDInAcc350450MC->Add(fHistd0rphiITSMIoneSPDInAccS350450MC); 
  fHistd0rphiITSMIoneSPDInAcc350450MC->Scale(1./fHistd0rphiITSMIoneSPDInAcc350450MC->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc500700MC=(TH1F*)fHistd0rphiITSMIoneSPDInAccP500700MC->Clone("fHistd0rphiITSMIoneSPDInAcc500700MC");
  fHistd0rphiITSMIoneSPDInAcc500700MC->Add(fHistd0rphiITSMIoneSPDInAccS500700MC);
  fHistd0rphiITSMIoneSPDInAcc500700MC->Scale(1./fHistd0rphiITSMIoneSPDInAcc500700MC->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc10001500MC=(TH1F*)fHistd0rphiITSMIoneSPDInAccP10001500MC->Clone("fHistd0rphiITSMIoneSPDInAcc10001500MC");
  fHistd0rphiITSMIoneSPDInAcc10001500MC->Add(fHistd0rphiITSMIoneSPDInAccS10001500MC);
  fHistd0rphiITSMIoneSPDInAcc10001500MC->Scale(1./fHistd0rphiITSMIoneSPDInAcc10001500MC->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc25004000MC=(TH1F*)fHistd0rphiITSMIoneSPDInAccP25004000MC->Clone("fHistd0rphiITSMIoneSPDInAcc25004000MC");
  fHistd0rphiITSMIoneSPDInAcc25004000MC->Add(fHistd0rphiITSMIoneSPDInAccS25004000MC);
  fHistd0rphiITSMIoneSPDInAcc25004000MC->Scale(1./fHistd0rphiITSMIoneSPDInAcc25004000MC->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc40008000MC=(TH1F*)fHistd0rphiITSMIoneSPDInAccP40008000MC->Clone("fHistd0rphiITSMIoneSPDInAcc40008000MC");
  fHistd0rphiITSMIoneSPDInAcc40008000MC->Add(fHistd0rphiITSMIoneSPDInAccS40008000MC);
  fHistd0rphiITSMIoneSPDInAcc40008000MC->Scale(1./fHistd0rphiITSMIoneSPDInAcc40008000MC->GetEntries());
  

  TFile *f= new TFile("AnalysisResults_onlynonfakes.root");
  list=(TList*)f->Get("cOutputITS");
  if(!list) {
    dir=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
    if(dir) list = (TList*)dir->Get("cOutputITS");
  }

  if(!list) return kFALSE;

  TH1F *fHistd0rphiITSMIoneSPDInAccP150200 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP150200");
  TH1F *fHistd0rphiITSMIoneSPDInAccS150200 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS150200");
  TH1F *fHistd0rphiITSMIoneSPDInAccP350450 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP350450");
  TH1F *fHistd0rphiITSMIoneSPDInAccS350450 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS350450");
  TH1F *fHistd0rphiITSMIoneSPDInAccP500700 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP500700");
  TH1F *fHistd0rphiITSMIoneSPDInAccS500700 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS500700");
  TH1F *fHistd0rphiITSMIoneSPDInAccP10001500 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP10001500");
  TH1F *fHistd0rphiITSMIoneSPDInAccS10001500 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS10001500");
  TH1F *fHistd0rphiITSMIoneSPDInAccP25004000 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP25004000");
  TH1F *fHistd0rphiITSMIoneSPDInAccS25004000 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS25004000");
  TH1F *fHistd0rphiITSMIoneSPDInAccP40008000 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP40008000");
  TH1F *fHistd0rphiITSMIoneSPDInAccS40008000 = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS40008000");
  TH1F *fHistd0rphiITSMIoneSPDInAcc150200=(TH1F*)fHistd0rphiITSMIoneSPDInAccP150200->Clone("fHistd0rphiITSMIoneSPDInAcc150200");
  fHistd0rphiITSMIoneSPDInAcc150200->Add(fHistd0rphiITSMIoneSPDInAccS150200);
  fHistd0rphiITSMIoneSPDInAcc150200->Scale(1./fHistd0rphiITSMIoneSPDInAcc150200->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc350450=(TH1F*)fHistd0rphiITSMIoneSPDInAccP350450->Clone("fHistd0rphiITSMIoneSPDInAcc350450");
  fHistd0rphiITSMIoneSPDInAcc350450->Add(fHistd0rphiITSMIoneSPDInAccS350450);
  fHistd0rphiITSMIoneSPDInAcc350450->Scale(1./fHistd0rphiITSMIoneSPDInAcc350450->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc500700=(TH1F*)fHistd0rphiITSMIoneSPDInAccP500700->Clone("fHistd0rphiITSMIoneSPDInAcc500700");
  fHistd0rphiITSMIoneSPDInAcc500700->Add(fHistd0rphiITSMIoneSPDInAccS500700);
  fHistd0rphiITSMIoneSPDInAcc500700->Scale(1./fHistd0rphiITSMIoneSPDInAcc500700->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc10001500=(TH1F*)fHistd0rphiITSMIoneSPDInAccP10001500->Clone("fHistd0rphiITSMIoneSPDInAcc10001500");
  fHistd0rphiITSMIoneSPDInAcc10001500->Add(fHistd0rphiITSMIoneSPDInAccS10001500);
  fHistd0rphiITSMIoneSPDInAcc10001500->Scale(1./fHistd0rphiITSMIoneSPDInAcc10001500->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc25004000=(TH1F*)fHistd0rphiITSMIoneSPDInAccP25004000->Clone("fHistd0rphiITSMIoneSPDInAcc25004000");
  fHistd0rphiITSMIoneSPDInAcc25004000->Add(fHistd0rphiITSMIoneSPDInAccS25004000);
  fHistd0rphiITSMIoneSPDInAcc25004000->Scale(1./fHistd0rphiITSMIoneSPDInAcc25004000->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAcc40008000=(TH1F*)fHistd0rphiITSMIoneSPDInAccP40008000->Clone("fHistd0rphiITSMIoneSPDInAcc40008000");
  fHistd0rphiITSMIoneSPDInAcc40008000->Add(fHistd0rphiITSMIoneSPDInAccS40008000);
  fHistd0rphiITSMIoneSPDInAcc40008000->Scale(1./fHistd0rphiITSMIoneSPDInAcc40008000->GetEntries()); 


  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(3,2);
  c1->cd(1);
  fHistd0rphiITSMIoneSPDInAcc150200MC->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc150200->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc150200MC->SetLineColor(2);
  fHistd0rphiITSMIoneSPDInAcc150200MC->Draw();
  fHistd0rphiITSMIoneSPDInAcc150200->SetLineColor(4);
  fHistd0rphiITSMIoneSPDInAcc150200->Draw("same");
  c1->cd(2);
  fHistd0rphiITSMIoneSPDInAcc350450MC->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc350450->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc350450MC->SetLineColor(2);
  fHistd0rphiITSMIoneSPDInAcc350450MC->Draw();
  fHistd0rphiITSMIoneSPDInAcc350450->SetLineColor(4);
  fHistd0rphiITSMIoneSPDInAcc350450->Draw("same"); 
  c1->cd(3);
  fHistd0rphiITSMIoneSPDInAcc500700MC->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc500700->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc500700MC->SetLineColor(2);
  fHistd0rphiITSMIoneSPDInAcc500700MC->Draw();
  fHistd0rphiITSMIoneSPDInAcc500700->SetLineColor(4);
  fHistd0rphiITSMIoneSPDInAcc500700->Draw("same");
  c1->cd(4);
  fHistd0rphiITSMIoneSPDInAcc10001500MC->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc10001500->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc10001500MC->SetLineColor(2);
  fHistd0rphiITSMIoneSPDInAcc10001500MC->Draw();
  fHistd0rphiITSMIoneSPDInAcc10001500->SetLineColor(4);
  fHistd0rphiITSMIoneSPDInAcc10001500->Draw("same");
   c1->cd(5);
  fHistd0rphiITSMIoneSPDInAcc25004000MC->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc25004000->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc25004000MC->SetLineColor(2);
  fHistd0rphiITSMIoneSPDInAcc25004000MC->Draw();
  fHistd0rphiITSMIoneSPDInAcc25004000->SetLineColor(4);
  fHistd0rphiITSMIoneSPDInAcc25004000->Draw("same");
  c1->cd(6);
  fHistd0rphiITSMIoneSPDInAcc40008000MC->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc40008000->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAcc40008000MC->SetLineColor(2);
  fHistd0rphiITSMIoneSPDInAcc40008000MC->Draw();
  fHistd0rphiITSMIoneSPDInAcc40008000->SetLineColor(4);
  fHistd0rphiITSMIoneSPDInAcc40008000->Draw("same");
  
  
  TCanvas *c1a = new TCanvas("c1a","c1a");
  c1a->Divide(3,1);
  c1a->cd(1);
  TH1F* fHistd0rphiITSMIoneSPDInAcc150200Ratio=(TH1F*)fHistd0rphiITSMIoneSPDInAcc150200->Clone("fHistd0rphiITSMIoneSPDInAcc150200Ratio");
  fHistd0rphiITSMIoneSPDInAcc150200Ratio->Divide(fHistd0rphiITSMIoneSPDInAcc150200MC);
  fHistd0rphiITSMIoneSPDInAcc150200Ratio->Draw();
  c1a->cd(2);
  TH1F* fHistd0rphiITSMIoneSPDInAcc500700Ratio=(TH1F*)fHistd0rphiITSMIoneSPDInAcc500700->Clone("fHistd0rphiITSMIoneSPDInAcc500700Ratio");
  fHistd0rphiITSMIoneSPDInAcc500700Ratio->Divide(fHistd0rphiITSMIoneSPDInAcc500700MC);
  fHistd0rphiITSMIoneSPDInAcc500700Ratio->Draw();
  c1a->cd(3);
  TH1F* fHistd0rphiITSMIoneSPDInAcc10001500Ratio=(TH1F*)fHistd0rphiITSMIoneSPDInAcc10001500->Clone("fHistd0rphiITSMIoneSPDInAcc10001500Ratio");
  fHistd0rphiITSMIoneSPDInAcc10001500Ratio->Divide(fHistd0rphiITSMIoneSPDInAcc10001500MC);
  fHistd0rphiITSMIoneSPDInAcc10001500Ratio->Draw();



  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(3,2);
  c2_1->SetLogy();
  c2_2->SetLogy();
  c2_3->SetLogy();
  c2_4->SetLogy();
  c2_5->SetLogy();
  c2_6->SetLogy();
  c2->cd(1);
  fHistd0rphiITSMIoneSPDInAccP150200MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS150200MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS150200fromStrangeMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS150200fromMatMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccP150200MC->SetLineColor(1);
  fHistd0rphiITSMIoneSPDInAccP150200MC->Draw();
  fHistd0rphiITSMIoneSPDInAccS150200MC->SetLineColor(6);
  fHistd0rphiITSMIoneSPDInAccS150200MC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS150200fromStrangeMC->SetLineColor(8);
  fHistd0rphiITSMIoneSPDInAccS150200fromStrangeMC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS150200fromMatMC->SetLineColor(9);
  fHistd0rphiITSMIoneSPDInAccS150200fromMatMC->Draw("same");
  c2->cd(2);
  fHistd0rphiITSMIoneSPDInAccP350450MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS350450MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS350450fromStrangeMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS350450fromMatMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccP350450MC->SetLineColor(1);
  fHistd0rphiITSMIoneSPDInAccP350450MC->Draw();
  fHistd0rphiITSMIoneSPDInAccS350450MC->SetLineColor(6);
  fHistd0rphiITSMIoneSPDInAccS350450MC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS350450fromStrangeMC->SetLineColor(8);
  fHistd0rphiITSMIoneSPDInAccS350450fromStrangeMC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS350450fromMatMC->SetLineColor(9);
  fHistd0rphiITSMIoneSPDInAccS350450fromMatMC->Draw("same"); 
 c2->cd(3);
  fHistd0rphiITSMIoneSPDInAccP500700MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700fromStrangeMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700fromMatMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccP500700MC->SetLineColor(1);
  fHistd0rphiITSMIoneSPDInAccP500700MC->Draw();
  fHistd0rphiITSMIoneSPDInAccS500700MC->SetLineColor(6);
  fHistd0rphiITSMIoneSPDInAccS500700MC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS500700fromStrangeMC->SetLineColor(8);
  fHistd0rphiITSMIoneSPDInAccS500700fromStrangeMC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS500700fromMatMC->SetLineColor(9);
  fHistd0rphiITSMIoneSPDInAccS500700fromMatMC->Draw("same");
  c2->cd(4);
  fHistd0rphiITSMIoneSPDInAccP10001500MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS10001500MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS10001500fromStrangeMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS10001500fromMatMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccP10001500MC->SetLineColor(1);
  fHistd0rphiITSMIoneSPDInAccP10001500MC->Draw();
  fHistd0rphiITSMIoneSPDInAccS10001500MC->SetLineColor(6);
  fHistd0rphiITSMIoneSPDInAccS10001500MC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS10001500fromStrangeMC->SetLineColor(8);
  fHistd0rphiITSMIoneSPDInAccS10001500fromStrangeMC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS10001500fromMatMC->SetLineColor(9);
  fHistd0rphiITSMIoneSPDInAccS10001500fromMatMC->Draw("same");
  c2->cd(5);
  fHistd0rphiITSMIoneSPDInAccP25004000MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS25004000MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS25004000fromStrangeMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS25004000fromMatMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccP25004000MC->SetLineColor(1);
  fHistd0rphiITSMIoneSPDInAccP25004000MC->Draw();
  fHistd0rphiITSMIoneSPDInAccS25004000MC->SetLineColor(6);
  fHistd0rphiITSMIoneSPDInAccS25004000MC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS25004000fromStrangeMC->SetLineColor(8);
  fHistd0rphiITSMIoneSPDInAccS25004000fromStrangeMC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS25004000fromMatMC->SetLineColor(9);
  fHistd0rphiITSMIoneSPDInAccS25004000fromMatMC->Draw("same");
  c2->cd(6);
  fHistd0rphiITSMIoneSPDInAccP40008000MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS40008000MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS40008000fromStrangeMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS40008000fromMatMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccP40008000MC->SetLineColor(1);
  fHistd0rphiITSMIoneSPDInAccP40008000MC->Draw();
  fHistd0rphiITSMIoneSPDInAccS40008000MC->SetLineColor(6);
  fHistd0rphiITSMIoneSPDInAccS40008000MC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS40008000fromStrangeMC->SetLineColor(8);
  fHistd0rphiITSMIoneSPDInAccS40008000fromStrangeMC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS40008000fromMatMC->SetLineColor(9);
  fHistd0rphiITSMIoneSPDInAccS40008000fromMatMC->Draw("same");

  TCanvas *c2b = new TCanvas("c2b","c2b");
  fHistd0rphiITSMIoneSPDInAccS500700fromStrangeMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700fromStrangeMC->SetLineColor(8);
  fHistd0rphiITSMIoneSPDInAccS500700fromStrangeMC->Draw();
  fHistd0rphiITSMIoneSPDInAccS500700from310MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700from310MC->SetLineColor(2);
  fHistd0rphiITSMIoneSPDInAccS500700from310MC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS500700from321MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700from321MC->SetLineColor(5);
  fHistd0rphiITSMIoneSPDInAccS500700from321MC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS500700from3122MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700from3122MC->SetLineColor(1);
  fHistd0rphiITSMIoneSPDInAccS500700from3122MC->Draw("same");

  TCanvas *c2c = new TCanvas("c2c","c2c");
  fHistd0rphiITSMIoneSPDInAccS500700fromMatMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700fromMatMC->Draw();
  fHistd0rphiITSMIoneSPDInAccS500700from211MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700from211MC->SetLineColor(2);
  fHistd0rphiITSMIoneSPDInAccS500700from211MC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccS500700from22MC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccS500700from22MC->SetLineColor(1);
  fHistd0rphiITSMIoneSPDInAccS500700from22MC->Draw("same");



  TCanvas *c3 = new TCanvas("c3","c3");
  c3->Divide(3,2);
  c3_1->SetLogx();
  c3_2->SetLogx();
  c3_3->SetLogx();
  c3_4->SetLogx();
  c3_5->SetLogx();
  c3_6->SetLogx();
  Float_t d0cut[15]={0.0301,0.0401,0.0501,0.0601,0.0801,0.101,0.151,0.201,0.251,0.301,0.401,0.601,0.801,1.001,1.401};
  
  c3->cd(1);
  Float_t fracP150200[15],fracS150200[15],fracSfromStrange150200[15],fracS150200Strp30[15],fracS150200Strm30[15],fracS150200Matm10[15],fracS150200Matp10[15];
  Float_t intPtot150200=fHistd0rphiITSMIoneSPDInAccP150200MC->Integral(1,fHistd0rphiITSMIoneSPDInAccP150200MC->GetNbinsX());
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAccP150200MC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAccP150200MC->FindBin(+d0cut[i]);
    Float_t intPcut150200=fHistd0rphiITSMIoneSPDInAccP150200MC->Integral(bin1,bin2);
    Float_t intScut150200=fHistd0rphiITSMIoneSPDInAccS150200MC->Integral(bin1,bin2);
    Float_t intSfromStrangecut150200=fHistd0rphiITSMIoneSPDInAccS150200fromStrangeMC->Integral(bin1,bin2);
    Float_t intScut150200Strp30 = intScut150200 + 0.3*intSfromStrangecut150200; 
    Float_t intScut150200Strm30 = intScut150200 - 0.3*intSfromStrangecut150200; 
    Float_t intScut150200Matp10 = intScut150200 + 0.1*(intScut150200-intSfromStrangecut150200); 
    Float_t intScut150200Matm10 = intScut150200 - 0.1*(intScut150200-intSfromStrangecut150200); 
    fracP150200[i]=intPcut150200/intPtot150200;
    fracS150200[i]=1.-intScut150200/(intPcut150200+intScut150200);
    fracS150200Strp30[i]=1.-intScut150200Strp30/(intPcut150200+intScut150200Strp30);
    fracS150200Strm30[i]=1.-intScut150200Strm30/(intPcut150200+intScut150200Strm30);
    fracS150200Matp10[i]=1.-intScut150200Matp10/(intPcut150200+intScut150200Matp10);
    fracS150200Matm10[i]=1.-intScut150200Matm10/(intPcut150200+intScut150200Matm10);
    fracSfromStrange150200[i]=1.-intSfromStrangecut150200/(intPcut150200+intScut150200);
  }
  TGraph *gfracP150200=new TGraph(15,d0cut,fracP150200);
  gfracP150200->SetMarkerColor(2);
  gfracP150200->SetMarkerStyle(20);
  gfracP150200->Draw("ap");
  TGraph *gfracS150200=new TGraph(15,d0cut,fracS150200);
  gfracS150200->SetMarkerColor(4);
  gfracS150200->SetMarkerStyle(21);
  gfracS150200->Draw("p");
  TGraph *gfracSfromStrange150200=new TGraph(15,d0cut,fracSfromStrange150200);
  gfracSfromStrange150200->SetMarkerColor(8);
  gfracSfromStrange150200->SetMarkerStyle(22);
  gfracSfromStrange150200->Draw("p");  
   
  c3->cd(2);
  Float_t fracP350450[15],fracS350450[15],fracSfromStrange350450[15],fracS350450Strp30[15],fracS350450Strm30[15],fracS350450Matm10[15],fracS350450Matp10[15];
  Float_t intPtot350450=fHistd0rphiITSMIoneSPDInAccP350450MC->Integral(1,fHistd0rphiITSMIoneSPDInAccP350450MC->GetNbinsX());
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAccP350450MC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAccP350450MC->FindBin(+d0cut[i]);
    Float_t intPcut350450=fHistd0rphiITSMIoneSPDInAccP350450MC->Integral(bin1,bin2);
    Float_t intScut350450=fHistd0rphiITSMIoneSPDInAccS350450MC->Integral(bin1,bin2);
    Float_t intSfromStrangecut350450=fHistd0rphiITSMIoneSPDInAccS350450fromStrangeMC->Integral(bin1,bin2);
    Float_t intScut350450Strp30 = intScut350450 + 0.3*intSfromStrangecut350450; 
    Float_t intScut350450Strm30 = intScut350450 - 0.3*intSfromStrangecut350450; 
    Float_t intScut350450Matp10 = intScut350450 + 0.1*(intScut350450-intSfromStrangecut350450); 
    Float_t intScut350450Matm10 = intScut350450 - 0.1*(intScut350450-intSfromStrangecut350450); 
    fracP350450[i]=intPcut350450/intPtot350450;
    fracS350450[i]=1.-intScut350450/(intPcut350450+intScut350450);
    fracS350450Strp30[i]=1.-intScut350450Strp30/(intPcut350450+intScut350450Strp30);
    fracS350450Strm30[i]=1.-intScut350450Strm30/(intPcut350450+intScut350450Strm30);
    fracS350450Matp10[i]=1.-intScut350450Matp10/(intPcut350450+intScut350450Matp10);
    fracS350450Matm10[i]=1.-intScut350450Matm10/(intPcut350450+intScut350450Matm10);
    fracSfromStrange350450[i]=1.-intSfromStrangecut350450/(intPcut350450+intScut350450);
  }
  TGraph *gfracP350450=new TGraph(15,d0cut,fracP350450);
  gfracP350450->SetMarkerColor(2);
  gfracP350450->SetMarkerStyle(20);
  gfracP350450->Draw("ap");
  TGraph *gfracS350450=new TGraph(15,d0cut,fracS350450);
  gfracS350450->SetMarkerColor(4);
  gfracS350450->SetMarkerStyle(21);
  gfracS350450->Draw("p");
  TGraph *gfracSfromStrange350450=new TGraph(15,d0cut,fracSfromStrange350450);
  gfracSfromStrange350450->SetMarkerColor(8);
  gfracSfromStrange350450->SetMarkerStyle(22);
  gfracSfromStrange350450->Draw("p");

  c3->cd(3);
  Float_t fracP500700[15],fracS500700[15],fracSfromStrange500700[15],fracS500700Strp30[15],fracS500700Strm30[15],fracS500700Matm10[15],fracS500700Matp10[15];
  Float_t intPtot500700=fHistd0rphiITSMIoneSPDInAccP500700MC->Integral(1,fHistd0rphiITSMIoneSPDInAccP500700MC->GetNbinsX());
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAccP500700MC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAccP500700MC->FindBin(+d0cut[i]);
    Float_t intPcut500700=fHistd0rphiITSMIoneSPDInAccP500700MC->Integral(bin1,bin2);
    Float_t intScut500700=fHistd0rphiITSMIoneSPDInAccS500700MC->Integral(bin1,bin2);
    Float_t intSfromStrangecut500700=fHistd0rphiITSMIoneSPDInAccS500700fromStrangeMC->Integral(bin1,bin2);
    Float_t intScut500700Strp30 = intScut500700 + 0.3*intSfromStrangecut500700; 
    Float_t intScut500700Strm30 = intScut500700 - 0.3*intSfromStrangecut500700; 
    Float_t intScut500700Matp10 = intScut500700 + 0.1*(intScut500700-intSfromStrangecut500700); 
    Float_t intScut500700Matm10 = intScut500700 - 0.1*(intScut500700-intSfromStrangecut500700); 
    fracP500700[i]=intPcut500700/intPtot500700;
    fracS500700[i]=1.-intScut500700/(intPcut500700+intScut500700);
    fracS500700Strp30[i]=1.-intScut500700Strp30/(intPcut500700+intScut500700Strp30);
    fracS500700Strm30[i]=1.-intScut500700Strm30/(intPcut500700+intScut500700Strm30);
    fracS500700Matp10[i]=1.-intScut500700Matp10/(intPcut500700+intScut500700Matp10);
    fracS500700Matm10[i]=1.-intScut500700Matm10/(intPcut500700+intScut500700Matm10);
    fracSfromStrange500700[i]=1.-intSfromStrangecut500700/(intPcut500700+intScut500700);
  }
  TGraph *gfracP500700=new TGraph(15,d0cut,fracP500700);
  gfracP500700->SetMarkerColor(2);
  gfracP500700->SetMarkerStyle(20);
  gfracP500700->Draw("ap");
  TGraph *gfracS500700=new TGraph(15,d0cut,fracS500700);
  gfracS500700->SetMarkerColor(4);
  gfracS500700->SetMarkerStyle(21);
  gfracS500700->Draw("p");
  TGraph *gfracSfromStrange500700=new TGraph(15,d0cut,fracSfromStrange500700);
  gfracSfromStrange500700->SetMarkerColor(8);
  gfracSfromStrange500700->SetMarkerStyle(22);
  gfracSfromStrange500700->Draw("p");

  c3->cd(4);
  Float_t fracP10001500[15],fracS10001500[15],fracSfromStrange10001500[15],fracS10001500Strp30[15],fracS10001500Strm30[15],fracS10001500Matm10[15],fracS10001500Matp10[15];
  Float_t intPtot10001500=fHistd0rphiITSMIoneSPDInAccP10001500MC->Integral(1,fHistd0rphiITSMIoneSPDInAccP10001500MC->GetNbinsX());
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAccP10001500MC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAccP10001500MC->FindBin(+d0cut[i]);
    Float_t intPcut10001500=fHistd0rphiITSMIoneSPDInAccP10001500MC->Integral(bin1,bin2);
    Float_t intScut10001500=fHistd0rphiITSMIoneSPDInAccS10001500MC->Integral(bin1,bin2);
    Float_t intSfromStrangecut10001500=fHistd0rphiITSMIoneSPDInAccS10001500fromStrangeMC->Integral(bin1,bin2);
    Float_t intScut10001500Strp30 = intScut10001500 + 0.3*intSfromStrangecut10001500; 
    Float_t intScut10001500Strm30 = intScut10001500 - 0.3*intSfromStrangecut10001500; 
    Float_t intScut10001500Matp10 = intScut10001500 + 0.1*(intScut10001500-intSfromStrangecut10001500); 
    Float_t intScut10001500Matm10 = intScut10001500 - 0.1*(intScut10001500-intSfromStrangecut10001500); 
    fracP10001500[i]=intPcut10001500/intPtot10001500;
    fracS10001500[i]=1.-intScut10001500/(intPcut10001500+intScut10001500);
    fracS10001500Strp30[i]=1.-intScut10001500Strp30/(intPcut10001500+intScut10001500Strp30);
    fracS10001500Strm30[i]=1.-intScut10001500Strm30/(intPcut10001500+intScut10001500Strm30);
    fracS10001500Matp10[i]=1.-intScut10001500Matp10/(intPcut10001500+intScut10001500Matp10);
    fracS10001500Matm10[i]=1.-intScut10001500Matm10/(intPcut10001500+intScut10001500Matm10);
    fracSfromStrange10001500[i]=1.-intSfromStrangecut10001500/(intPcut10001500+intScut10001500);
  }
  TGraph *gfracP10001500=new TGraph(15,d0cut,fracP10001500);
  gfracP10001500->SetMarkerColor(2);
  gfracP10001500->SetMarkerStyle(20);
  gfracP10001500->Draw("ap");
  TGraph *gfracS10001500=new TGraph(15,d0cut,fracS10001500);
  gfracS10001500->SetMarkerColor(4);
  gfracS10001500->SetMarkerStyle(21);
  gfracS10001500->Draw("p");
  TGraph *gfracSfromStrange10001500=new TGraph(15,d0cut,fracSfromStrange10001500);
  gfracSfromStrange10001500->SetMarkerColor(8);
  gfracSfromStrange10001500->SetMarkerStyle(22);
  gfracSfromStrange10001500->Draw("p");
  
  c3->cd(5);
  Float_t fracP25004000[15],fracS25004000[15],fracSfromStrange25004000[15],fracS25004000Strp30[15],fracS25004000Strm30[15],fracS25004000Matm10[15],fracS25004000Matp10[15];
  Float_t intPtot25004000=fHistd0rphiITSMIoneSPDInAccP25004000MC->Integral(1,fHistd0rphiITSMIoneSPDInAccP25004000MC->GetNbinsX());
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAccP25004000MC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAccP25004000MC->FindBin(+d0cut[i]);
    Float_t intPcut25004000=fHistd0rphiITSMIoneSPDInAccP25004000MC->Integral(bin1,bin2);
    Float_t intScut25004000=fHistd0rphiITSMIoneSPDInAccS25004000MC->Integral(bin1,bin2);
    Float_t intSfromStrangecut25004000=fHistd0rphiITSMIoneSPDInAccS25004000fromStrangeMC->Integral(bin1,bin2);
    Float_t intScut25004000Strp30 = intScut25004000 + 0.3*intSfromStrangecut25004000; 
    Float_t intScut25004000Strm30 = intScut25004000 - 0.3*intSfromStrangecut25004000; 
    Float_t intScut25004000Matp10 = intScut25004000 + 0.1*(intScut25004000-intSfromStrangecut25004000); 
    Float_t intScut25004000Matm10 = intScut25004000 - 0.1*(intScut25004000-intSfromStrangecut25004000); 
    fracP25004000[i]=intPcut25004000/intPtot25004000;
    fracS25004000[i]=1.-intScut25004000/(intPcut25004000+intScut25004000);
    fracS25004000Strp30[i]=1.-intScut25004000Strp30/(intPcut25004000+intScut25004000Strp30);
    fracS25004000Strm30[i]=1.-intScut25004000Strm30/(intPcut25004000+intScut25004000Strm30);
    fracS25004000Matp10[i]=1.-intScut25004000Matp10/(intPcut25004000+intScut25004000Matp10);
    fracS25004000Matm10[i]=1.-intScut25004000Matm10/(intPcut25004000+intScut25004000Matm10);
    fracSfromStrange25004000[i]=1.-intSfromStrangecut25004000/(intPcut25004000+intScut25004000);
  }
  TGraph *gfracP25004000=new TGraph(15,d0cut,fracP25004000);
  gfracP25004000->SetMarkerColor(2);
  gfracP25004000->SetMarkerStyle(20);
  gfracP25004000->Draw("ap");
  TGraph *gfracS25004000=new TGraph(15,d0cut,fracS25004000);
  gfracS25004000->SetMarkerColor(4);
  gfracS25004000->SetMarkerStyle(21);
  gfracS25004000->Draw("p");
  TGraph *gfracSfromStrange25004000=new TGraph(15,d0cut,fracSfromStrange25004000);
  gfracSfromStrange25004000->SetMarkerColor(8);
  gfracSfromStrange25004000->SetMarkerStyle(22);
  gfracSfromStrange25004000->Draw("p");

  c3->cd(6);
  Float_t fracP40008000[15],fracS40008000[15],fracSfromStrange40008000[15],fracS40008000Strp30[15],fracS40008000Strm30[15],fracS40008000Matm10[15],fracS40008000Matp10[15];
  Float_t intPtot40008000=fHistd0rphiITSMIoneSPDInAccP40008000MC->Integral(1,fHistd0rphiITSMIoneSPDInAccP40008000MC->GetNbinsX());
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAccP40008000MC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAccP40008000MC->FindBin(+d0cut[i]);
    Float_t intPcut40008000=fHistd0rphiITSMIoneSPDInAccP40008000MC->Integral(bin1,bin2);
    Float_t intScut40008000=fHistd0rphiITSMIoneSPDInAccS40008000MC->Integral(bin1,bin2);
    Float_t intSfromStrangecut40008000=fHistd0rphiITSMIoneSPDInAccS40008000fromStrangeMC->Integral(bin1,bin2);
    Float_t intScut40008000Strp30 = intScut40008000 + 0.3*intSfromStrangecut40008000; 
    Float_t intScut40008000Strm30 = intScut40008000 - 0.3*intSfromStrangecut40008000; 
    Float_t intScut40008000Matp10 = intScut40008000 + 0.1*(intScut40008000-intSfromStrangecut40008000); 
    Float_t intScut40008000Matm10 = intScut40008000 - 0.1*(intScut40008000-intSfromStrangecut40008000); 
    fracP40008000[i]=intPcut40008000/intPtot40008000;
    fracS40008000[i]=1.-intScut40008000/(intPcut40008000+intScut40008000);
    fracS40008000Strp30[i]=1.-intScut40008000Strp30/(intPcut40008000+intScut40008000Strp30);
    fracS40008000Strm30[i]=1.-intScut40008000Strm30/(intPcut40008000+intScut40008000Strm30);
    fracS40008000Matp10[i]=1.-intScut40008000Matp10/(intPcut40008000+intScut40008000Matp10);
    fracS40008000Matm10[i]=1.-intScut40008000Matm10/(intPcut40008000+intScut40008000Matm10);
    fracSfromStrange40008000[i]=1.-intSfromStrangecut40008000/(intPcut40008000+intScut40008000);
  }
  TGraph *gfracP40008000=new TGraph(15,d0cut,fracP40008000);
  gfracP40008000->SetMarkerColor(2);
  gfracP40008000->SetMarkerStyle(20);
  gfracP40008000->Draw("ap");
  TGraph *gfracS40008000=new TGraph(15,d0cut,fracS40008000);
  gfracS40008000->SetMarkerColor(4);
  gfracS40008000->SetMarkerStyle(21);
  gfracS40008000->Draw("p");
  TGraph *gfracSfromStrange40008000=new TGraph(15,d0cut,fracSfromStrange40008000);
  gfracSfromStrange40008000->SetMarkerColor(8);
  gfracSfromStrange40008000->SetMarkerStyle(22);
  gfracSfromStrange40008000->Draw("p");

  
  TCanvas *c4 = new TCanvas("c4","c4");
  c4->Divide(3,1);
  c4_1->SetLogx();
  c4_2->SetLogx();
  c4_3->SetLogx();

  c4->cd(1);
  Float_t intDatacut150200[15],intDataPcut150200[15],intDataPall150200[15];
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAcc150200->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAcc150200->FindBin(+d0cut[i]);
    intDatacut150200[i]=fHistd0rphiITSMIoneSPDInAcc150200->Integral(bin1,bin2);
    intDataPcut150200[i]=intDatacut150200[i]*fracS150200[i];
    intDataPall150200[i]=intDataPcut150200[i]/fracP150200[i];
  }
  TGraph *gintDatacut150200=new TGraph(15,d0cut,intDatacut150200);
  gintDatacut150200->SetMarkerColor(1);
  gintDatacut150200->SetMarkerStyle(20);
  gintDatacut150200->Draw("ap");
  TGraph *gintDataPcut150200=new TGraph(15,d0cut,intDataPcut150200);
  gintDataPcut150200->SetMarkerColor(1);
  gintDataPcut150200->SetMarkerStyle(24);
  gintDataPcut150200->Draw("p");
  TGraph *gintDataPall150200=new TGraph(15,d0cut,intDataPall150200);
  gintDataPall150200->SetMarkerColor(2);
  gintDataPall150200->SetMarkerStyle(22);
  gintDataPall150200->Draw("p");

  c4->cd(2);
  Float_t intDatacut500700[15],intDataPcut500700[15],intDataPall500700[15];
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAcc500700->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAcc500700->FindBin(+d0cut[i]);
    intDatacut500700[i]=fHistd0rphiITSMIoneSPDInAcc500700->Integral(bin1,bin2);
    intDataPcut500700[i]=intDatacut500700[i]*fracS500700[i];
    intDataPall500700[i]=intDataPcut500700[i]/fracP500700[i];
  }
  TGraph *gintDatacut500700=new TGraph(15,d0cut,intDatacut500700);
  gintDatacut500700->SetMarkerColor(1);
  gintDatacut500700->SetMarkerStyle(20);
  gintDatacut500700->Draw("ap");
  TGraph *gintDataPcut500700=new TGraph(15,d0cut,intDataPcut500700);
  gintDataPcut500700->SetMarkerColor(1);
  gintDataPcut500700->SetMarkerStyle(24);
  gintDataPcut500700->Draw("p");
  TGraph *gintDataPall500700=new TGraph(15,d0cut,intDataPall500700);
  gintDataPall500700->SetMarkerColor(2);
  gintDataPall500700->SetMarkerStyle(22);
  gintDataPall500700->Draw("p");

  c4->cd(3);
  Float_t intDatacut10001500[15],intDataPcut10001500[15],intDataPall10001500[15];
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAcc10001500->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAcc10001500->FindBin(+d0cut[i]);
    intDatacut10001500[i]=fHistd0rphiITSMIoneSPDInAcc10001500->Integral(bin1,bin2);
    intDataPcut10001500[i]=intDatacut10001500[i]*fracS10001500[i];
    intDataPall10001500[i]=intDataPcut10001500[i]/fracP10001500[i];
  }
  TGraph *gintDatacut10001500=new TGraph(15,d0cut,intDatacut10001500);
  gintDatacut10001500->SetMarkerColor(1);
  gintDatacut10001500->SetMarkerStyle(20);
  gintDatacut10001500->Draw("ap");
  TGraph *gintDataPcut10001500=new TGraph(15,d0cut,intDataPcut10001500);
  gintDataPcut10001500->SetMarkerColor(1);
  gintDataPcut10001500->SetMarkerStyle(24);
  gintDataPcut10001500->Draw("p");
  TGraph *gintDataPall10001500=new TGraph(15,d0cut,intDataPall10001500);
  gintDataPall10001500->SetMarkerColor(2);
  gintDataPall10001500->SetMarkerStyle(22);
  gintDataPall10001500->Draw("p");

  
  TCanvas *c5 = new TCanvas("c5","c5");
  c5->Divide(3,1);
  c5_1->SetLogx();
  c5_2->SetLogx();
  c5_3->SetLogx();

  c5->cd(1);
  Float_t intDataPall150200Strp30[15],intDataPall150200Strm30[15],intDataPall150200Matp10[15],intDataPall150200Matm10[15];
  for(Int_t i=0;i<15;i++) {
    intDataPall150200Strp30[i]=intDatacut150200[i]*fracS150200Strp30[i]/fracP150200[i];
    intDataPall150200Strm30[i]=intDatacut150200[i]*fracS150200Strm30[i]/fracP150200[i];
    intDataPall150200Matp10[i]=intDatacut150200[i]*fracS150200Matp10[i]/fracP150200[i];
    intDataPall150200Matm10[i]=intDatacut150200[i]*fracS150200Matm10[i]/fracP150200[i];
  }
  gintDataPall150200->Draw("ap");

  TGraph *gintDataPall150200Strp30=new TGraph(15,d0cut,intDataPall150200Strp30);
  gintDataPall150200Strp30->SetMarkerColor(1);
  gintDataPall150200Strp30->SetMarkerStyle(22);
  gintDataPall150200Strp30->Draw("p");
  TGraph *gintDataPall150200Strm30=new TGraph(15,d0cut,intDataPall150200Strm30);
  gintDataPall150200Strm30->SetMarkerColor(1);
  gintDataPall150200Strm30->SetMarkerStyle(22);
  gintDataPall150200Strm30->Draw("p");
  TGraph *gintDataPall150200Matp10=new TGraph(15,d0cut,intDataPall150200Matp10);
  gintDataPall150200Matp10->SetMarkerColor(3);
  gintDataPall150200Matp10->SetMarkerStyle(22);
  gintDataPall150200Matp10->Draw("p");
  TGraph *gintDataPall150200Matm10=new TGraph(15,d0cut,intDataPall150200Matm10);
  gintDataPall150200Matm10->SetMarkerColor(3);
  gintDataPall150200Matm10->SetMarkerStyle(22);
  gintDataPall150200Matm10->Draw("p");

  c5->cd(2);
  Float_t intDataPall500700Strp30[15],intDataPall500700Strm30[15],intDataPall500700Matp10[15],intDataPall500700Matm10[15];
  for(Int_t i=0;i<15;i++) {
    intDataPall500700Strp30[i]=intDatacut500700[i]*fracS500700Strp30[i]/fracP500700[i];
    intDataPall500700Strm30[i]=intDatacut500700[i]*fracS500700Strm30[i]/fracP500700[i];
    intDataPall500700Matp10[i]=intDatacut500700[i]*fracS500700Matp10[i]/fracP500700[i];
    intDataPall500700Matm10[i]=intDatacut500700[i]*fracS500700Matm10[i]/fracP500700[i];
  }
  gintDataPall500700->Draw("ap");

  TGraph *gintDataPall500700Strp30=new TGraph(15,d0cut,intDataPall500700Strp30);
  gintDataPall500700Strp30->SetMarkerColor(1);
  gintDataPall500700Strp30->SetMarkerStyle(22);
  gintDataPall500700Strp30->Draw("p");
  TGraph *gintDataPall500700Strm30=new TGraph(15,d0cut,intDataPall500700Strm30);
  gintDataPall500700Strm30->SetMarkerColor(1);
  gintDataPall500700Strm30->SetMarkerStyle(22);
  gintDataPall500700Strm30->Draw("p");
  TGraph *gintDataPall500700Matp10=new TGraph(15,d0cut,intDataPall500700Matp10);
  gintDataPall500700Matp10->SetMarkerColor(3);
  gintDataPall500700Matp10->SetMarkerStyle(22);
  gintDataPall500700Matp10->Draw("p");
  TGraph *gintDataPall500700Matm10=new TGraph(15,d0cut,intDataPall500700Matm10);
  gintDataPall500700Matm10->SetMarkerColor(3);
  gintDataPall500700Matm10->SetMarkerStyle(22);
  gintDataPall500700Matm10->Draw("p");

  c5->cd(3);
  Float_t intDataPall10001500Strp30[15],intDataPall10001500Strm30[15],intDataPall10001500Matp10[15],intDataPall10001500Matm10[15];
  for(Int_t i=0;i<15;i++) {
    intDataPall10001500Strp30[i]=intDatacut10001500[i]*fracS10001500Strp30[i]/fracP10001500[i];
    intDataPall10001500Strm30[i]=intDatacut10001500[i]*fracS10001500Strm30[i]/fracP10001500[i];
    intDataPall10001500Matp10[i]=intDatacut10001500[i]*fracS10001500Matp10[i]/fracP10001500[i];
    intDataPall10001500Matm10[i]=intDatacut10001500[i]*fracS10001500Matm10[i]/fracP10001500[i];
  }
  gintDataPall10001500->Draw("ap");

  TGraph *gintDataPall10001500Strp30=new TGraph(15,d0cut,intDataPall10001500Strp30);
  gintDataPall10001500Strp30->SetMarkerColor(1);
  gintDataPall10001500Strp30->SetMarkerStyle(22);
  gintDataPall10001500Strp30->Draw("p");
  TGraph *gintDataPall10001500Strm30=new TGraph(15,d0cut,intDataPall10001500Strm30);
  gintDataPall10001500Strm30->SetMarkerColor(1);
  gintDataPall10001500Strm30->SetMarkerStyle(22);
  gintDataPall10001500Strm30->Draw("p");
  TGraph *gintDataPall10001500Matp10=new TGraph(15,d0cut,intDataPall10001500Matp10);
  gintDataPall10001500Matp10->SetMarkerColor(3);
  gintDataPall10001500Matp10->SetMarkerStyle(22);
  gintDataPall10001500Matp10->Draw("p");
  TGraph *gintDataPall10001500Matm10=new TGraph(15,d0cut,intDataPall10001500Matm10);
  gintDataPall10001500Matm10->SetMarkerColor(3);
  gintDataPall10001500Matm10->SetMarkerStyle(22);
  gintDataPall10001500Matm10->Draw("p");
  
  return;
}
//---------------------------------------------------------------------------
void PlotImpPar_z() {

  TFile *fMC= new TFile("ITS.Performance_lhc10a8_100k.root");

  TList *list=(TList*)fMC->Get("cOutputITS");
  TH1F *fHistd0zITSMIoneSPDInAccP150200MC = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccP150200");
  TH1F *fHistd0zITSMIoneSPDInAccS150200MC = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccS150200");
  TH1F *fHistd0zITSMIoneSPDInAccP500700MC = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccP500700");
  TH1F *fHistd0zITSMIoneSPDInAccS500700MC = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccS500700");
  TH1F *fHistd0zITSMIoneSPDInAccP10001500MC = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccP10001500");
  TH1F *fHistd0zITSMIoneSPDInAccS10001500MC = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccS10001500");
  TH1F *fHistd0zITSMIoneSPDInAcc150200MC=(TH1F*)fHistd0zITSMIoneSPDInAccP150200MC->Clone("fHistd0zITSMIoneSPDInAcc150200MC");
  fHistd0zITSMIoneSPDInAcc150200MC->Add(fHistd0zITSMIoneSPDInAccS150200MC);
  fHistd0zITSMIoneSPDInAcc150200MC->Scale(1./fHistd0zITSMIoneSPDInAcc150200MC->GetEntries());
  TH1F *fHistd0zITSMIoneSPDInAcc500700MC=(TH1F*)fHistd0zITSMIoneSPDInAccP500700MC->Clone("fHistd0zITSMIoneSPDInAcc500700MC");
  fHistd0zITSMIoneSPDInAcc500700MC->Add(fHistd0zITSMIoneSPDInAccS500700MC);
  fHistd0zITSMIoneSPDInAcc500700MC->Scale(1./fHistd0zITSMIoneSPDInAcc500700MC->GetEntries());
  TH1F *fHistd0zITSMIoneSPDInAcc10001500MC=(TH1F*)fHistd0zITSMIoneSPDInAccP10001500MC->Clone("fHistd0zITSMIoneSPDInAcc10001500MC");
  fHistd0zITSMIoneSPDInAcc10001500MC->Add(fHistd0zITSMIoneSPDInAccS10001500MC);
  fHistd0zITSMIoneSPDInAcc10001500MC->Scale(1./fHistd0zITSMIoneSPDInAcc10001500MC->GetEntries());


  TFile *f= new TFile("ITS.Performance_104892.root");

  TList *list=(TList*)f->Get("cOutputITS");
  TH1F *fHistd0zITSMIoneSPDInAccP150200 = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccP150200");
  TH1F *fHistd0zITSMIoneSPDInAccS150200 = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccS150200");
  TH1F *fHistd0zITSMIoneSPDInAccP500700 = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccP500700");
  TH1F *fHistd0zITSMIoneSPDInAccS500700 = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccS500700");
  TH1F *fHistd0zITSMIoneSPDInAccP10001500 = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccP10001500");
  TH1F *fHistd0zITSMIoneSPDInAccS10001500 = (TH1F*)list->FindObject("fHistd0zITSMIoneSPDInAccS10001500");
  TH1F *fHistd0zITSMIoneSPDInAcc150200=(TH1F*)fHistd0zITSMIoneSPDInAccP150200->Clone("fHistd0zITSMIoneSPDInAcc150200");
  fHistd0zITSMIoneSPDInAcc150200->Add(fHistd0zITSMIoneSPDInAccS150200);
  fHistd0zITSMIoneSPDInAcc150200->Scale(1./fHistd0zITSMIoneSPDInAcc150200->GetEntries());
  TH1F *fHistd0zITSMIoneSPDInAcc500700=(TH1F*)fHistd0zITSMIoneSPDInAccP500700->Clone("fHistd0zITSMIoneSPDInAcc500700");
  fHistd0zITSMIoneSPDInAcc500700->Add(fHistd0zITSMIoneSPDInAccS500700);
  fHistd0zITSMIoneSPDInAcc500700->Scale(1./fHistd0zITSMIoneSPDInAcc500700->GetEntries());
  TH1F *fHistd0zITSMIoneSPDInAcc10001500=(TH1F*)fHistd0zITSMIoneSPDInAccP10001500->Clone("fHistd0zITSMIoneSPDInAcc10001500");
  fHistd0zITSMIoneSPDInAcc10001500->Add(fHistd0zITSMIoneSPDInAccS10001500);
  fHistd0zITSMIoneSPDInAcc10001500->Scale(1./fHistd0zITSMIoneSPDInAcc10001500->GetEntries());


  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(3,1);
  c1->cd(1);
  fHistd0zITSMIoneSPDInAcc150200MC->SetLineColor(2);
  fHistd0zITSMIoneSPDInAcc150200MC->Draw();
  fHistd0zITSMIoneSPDInAcc150200->SetLineColor(4);
  fHistd0zITSMIoneSPDInAcc150200->Draw("same");
  c1->cd(2);
  fHistd0zITSMIoneSPDInAcc500700MC->SetLineColor(2);
  fHistd0zITSMIoneSPDInAcc500700MC->Draw();
  fHistd0zITSMIoneSPDInAcc500700->SetLineColor(4);
  fHistd0zITSMIoneSPDInAcc500700->Draw("same");
  c1->cd(3);
  fHistd0zITSMIoneSPDInAcc10001500MC->SetLineColor(2);
  fHistd0zITSMIoneSPDInAcc10001500MC->Draw();
  fHistd0zITSMIoneSPDInAcc10001500->SetLineColor(4);
  fHistd0zITSMIoneSPDInAcc10001500->Draw("same");

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(3,1);
  c2_1->SetLogy();
  c2_2->SetLogy();
  c2_3->SetLogy();
  c2->cd(1);
  fHistd0zITSMIoneSPDInAccP150200MC->SetMinimum(1);
  fHistd0zITSMIoneSPDInAccS150200MC->SetMinimum(1);
  fHistd0zITSMIoneSPDInAccP150200MC->SetLineColor(1);
  fHistd0zITSMIoneSPDInAccP150200MC->Draw();
  fHistd0zITSMIoneSPDInAccS150200MC->SetLineColor(6);
  fHistd0zITSMIoneSPDInAccS150200MC->Draw("same");
  c2->cd(2);
  fHistd0zITSMIoneSPDInAccP500700MC->SetMinimum(1);
  fHistd0zITSMIoneSPDInAccS500700MC->SetMinimum(1);
  fHistd0zITSMIoneSPDInAccP500700MC->SetLineColor(1);
  fHistd0zITSMIoneSPDInAccP500700MC->Draw();
  fHistd0zITSMIoneSPDInAccS500700MC->SetLineColor(6);
  fHistd0zITSMIoneSPDInAccS500700MC->Draw("same");
  c2->cd(3);
  fHistd0zITSMIoneSPDInAccP10001500MC->SetMinimum(1);
  fHistd0zITSMIoneSPDInAccS10001500MC->SetMinimum(1);
  fHistd0zITSMIoneSPDInAccP10001500MC->SetLineColor(1);
  fHistd0zITSMIoneSPDInAccP10001500MC->Draw();
  fHistd0zITSMIoneSPDInAccS10001500MC->SetLineColor(6);
  fHistd0zITSMIoneSPDInAccS10001500MC->Draw("same");


  TCanvas *c3 = new TCanvas("c3","c3");
  c3->Divide(3,1);
  c3_1->SetLogx();
  c3_2->SetLogx();
  c3_3->SetLogx();
  Float_t d0cut[15]={0.05,0.06,0.08,0.10,0.15,0.2,0.25,0.30,0.35,0.40,0.50,0.60,0.80,1.,3.40};

  c3->cd(1);
  Float_t fracP150200[15],fracS150200[15];
  Float_t intPtot150200=fHistd0zITSMIoneSPDInAccP150200MC->Integral(1,fHistd0zITSMIoneSPDInAccP150200MC->GetNbinsX());
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0zITSMIoneSPDInAccP150200MC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0zITSMIoneSPDInAccP150200MC->FindBin(+d0cut[i]);
    Float_t intPcut150200=fHistd0zITSMIoneSPDInAccP150200MC->Integral(bin1,bin2);
    Float_t intScut150200=fHistd0zITSMIoneSPDInAccS150200MC->Integral(bin1,bin2);
    fracP150200[i]=intPcut150200/intPtot150200;
    fracS150200[i]=1.-intScut150200/(intPcut150200+intScut150200);
  }
  TGraph *gfracP150200=new TGraph(15,d0cut,fracP150200);
  gfracP150200->SetMarkerColor(2);
  gfracP150200->SetMarkerStyle(20);
  gfracP150200->Draw("ap");
  TGraph *gfracS150200=new TGraph(15,d0cut,fracS150200);
  gfracS150200->SetMarkerColor(4);
  gfracS150200->SetMarkerStyle(21);
  gfracS150200->Draw("p");

  c3->cd(2);
  Float_t fracP500700[15],fracS500700[15];
  Float_t intPtot500700=fHistd0zITSMIoneSPDInAccP500700MC->Integral(1,fHistd0zITSMIoneSPDInAccP500700MC->GetNbinsX());
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0zITSMIoneSPDInAccP500700MC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0zITSMIoneSPDInAccP500700MC->FindBin(+d0cut[i]);
    Float_t intPcut500700=fHistd0zITSMIoneSPDInAccP500700MC->Integral(bin1,bin2);
    Float_t intScut500700=fHistd0zITSMIoneSPDInAccS500700MC->Integral(bin1,bin2);
    fracP500700[i]=intPcut500700/intPtot500700;
    fracS500700[i]=1.-intScut500700/(intPcut500700+intScut500700);
  }
  TGraph *gfracP500700=new TGraph(15,d0cut,fracP500700);
  gfracP500700->SetMarkerColor(2);
  gfracP500700->SetMarkerStyle(20);
  gfracP500700->Draw("ap");
  TGraph *gfracS500700=new TGraph(15,d0cut,fracS500700);
  gfracS500700->SetMarkerColor(4);
  gfracS500700->SetMarkerStyle(21);
  gfracS500700->Draw("p");

  c3->cd(3);
  Float_t fracP10001500[15],fracS10001500[15];
  Float_t intPtot10001500=fHistd0zITSMIoneSPDInAccP10001500MC->Integral(1,fHistd0zITSMIoneSPDInAccP10001500MC->GetNbinsX());
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0zITSMIoneSPDInAccP10001500MC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0zITSMIoneSPDInAccP10001500MC->FindBin(+d0cut[i]);
    Float_t intPcut10001500=fHistd0zITSMIoneSPDInAccP10001500MC->Integral(bin1,bin2);
    Float_t intScut10001500=fHistd0zITSMIoneSPDInAccS10001500MC->Integral(bin1,bin2);
    fracP10001500[i]=intPcut10001500/intPtot10001500;
    fracS10001500[i]=1.-intScut10001500/(intPcut10001500+intScut10001500);
  }
  TGraph *gfracP10001500=new TGraph(15,d0cut,fracP10001500);
  gfracP10001500->SetMarkerColor(2);
  gfracP10001500->SetMarkerStyle(20);
  gfracP10001500->Draw("ap");
  TGraph *gfracS10001500=new TGraph(15,d0cut,fracS10001500);
  gfracS10001500->SetMarkerColor(4);
  gfracS10001500->SetMarkerStyle(21);
  gfracS10001500->Draw("p");

  return;
}
//----------------------------------------------------------------------------
void Corrections() {

  TFile *f_cutits= new TFile("ITS.Performance_lhc09d10_500k_noz0cut_nodcaTPC.root");

  TList *list=(TList*)f_cutits->Get("cOutputITS");

  TH1F *fHistNEvents = (TH1F*)list->FindObject("fHistNEvents");
  Float_t nEventsSelWithVertex = fHistNEvents->GetBinContent(5);

  TH1F *fHistPtITSTPCsel = (TH1F*)list->FindObject("fHistPtITSTPCsel");
  NormalizePtHist(fHistPtITSTPCsel,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselP = (TH1F*)list->FindObject("fHistPtITSTPCselP");
  NormalizePtHist(fHistPtITSTPCselP,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselS = (TH1F*)list->FindObject("fHistPtITSTPCselS");
  NormalizePtHist(fHistPtITSTPCselS,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselPfromStrange = (TH1F*)list->FindObject("fHistPtITSTPCselPfromStrange");
  NormalizePtHist(fHistPtITSTPCselPfromStrange,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselSfromStrange = (TH1F*)list->FindObject("fHistPtITSTPCselSfromStrange");
  NormalizePtHist(fHistPtITSTPCselSfromStrange,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselSfromMat = (TH1F*)list->FindObject("fHistPtITSTPCselSfromMat");
  NormalizePtHist(fHistPtITSTPCselSfromMat,nEventsSelWithVertex);


  TFile *f_cuttpc= new TFile("ITS.Performance_lhc09d10_500k_noz0cut.root");

  TList *list=(TList*)f_cuttpc->Get("cOutputITS");

  TH1F *fHistNEvents = (TH1F*)list->FindObject("fHistNEvents");
  nEventsSelWithVertex = fHistNEvents->GetBinContent(5);

  TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  NormalizePtHist(fHistPtTPCInAcc,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccP = (TH1F*)list->FindObject("fHistPtTPCInAccP");
  NormalizePtHist(fHistPtTPCInAccP,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccS = (TH1F*)list->FindObject("fHistPtTPCInAccS");
  NormalizePtHist(fHistPtTPCInAccS,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccPfromStrange = (TH1F*)list->FindObject("fHistPtTPCInAccPfromStrange"); 
  NormalizePtHist(fHistPtTPCInAccPfromStrange,nEventsSelWithVertex);
 TH1F *fHistPtTPCInAccSfromStrange = (TH1F*)list->FindObject("fHistPtTPCInAccSfromStrange");
  NormalizePtHist(fHistPtTPCInAccSfromStrange,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccSfromMat = (TH1F*)list->FindObject("fHistPtTPCInAccSfromMat");
  NormalizePtHist(fHistPtTPCInAccSfromMat,nEventsSelWithVertex);


  TFile *f_nocut= new TFile("ITS.Performance_lhc09d10_500k_nod0z0cut.root");

  TList *list=(TList*)f_nocut->Get("cOutputITS");

  TH1F *fHistNEvents = (TH1F*)list->FindObject("fHistNEvents");
  nEventsSelWithVertex = fHistNEvents->GetBinContent(5);

  TH1F *fHistPtITSTPCselnocut = (TH1F*)list->FindObject("fHistPtITSTPCsel");
  NormalizePtHist(fHistPtITSTPCselnocut,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselnocutP = (TH1F*)list->FindObject("fHistPtITSTPCselP");
  NormalizePtHist(fHistPtITSTPCselnocutP,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselnocutS = (TH1F*)list->FindObject("fHistPtITSTPCselS");
  NormalizePtHist(fHistPtITSTPCselnocutS,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselnocutPfromStrange = (TH1F*)list->FindObject("fHistPtITSTPCselPfromStrange");
  NormalizePtHist(fHistPtITSTPCselnocutPfromStrange,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselnocutSfromStrange = (TH1F*)list->FindObject("fHistPtITSTPCselSfromStrange");
  NormalizePtHist(fHistPtITSTPCselnocutSfromStrange,nEventsSelWithVertex);
  TH1F *fHistPtITSTPCselnocutSfromMat = (TH1F*)list->FindObject("fHistPtITSTPCselSfromMat");
  NormalizePtHist(fHistPtITSTPCselnocutSfromMat,nEventsSelWithVertex);

  TH1F *fHistPtTPCInAccnocut = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  NormalizePtHist(fHistPtTPCInAccnocut,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccnocutP = (TH1F*)list->FindObject("fHistPtTPCInAccP");
  NormalizePtHist(fHistPtTPCInAccnocutP,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccnocutS = (TH1F*)list->FindObject("fHistPtTPCInAccS");
  NormalizePtHist(fHistPtTPCInAccnocutS,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccnocutPfromStrange = (TH1F*)list->FindObject("fHistPtTPCInAccPfromStrange"); 
  NormalizePtHist(fHistPtTPCInAccnocutPfromStrange,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccnocutSfromStrange = (TH1F*)list->FindObject("fHistPtTPCInAccSfromStrange");
  NormalizePtHist(fHistPtTPCInAccnocutSfromStrange,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccnocutSfromMat = (TH1F*)list->FindObject("fHistPtTPCInAccSfromMat");
  NormalizePtHist(fHistPtTPCInAccnocutSfromMat,nEventsSelWithVertex);

  // cut efficiencies
  
  TH1F *fHistPtITSTPCeffselP = (TH1F*)fHistPtITSTPCselP->Clone("fHistPtITSTPCeffselP");
  fHistPtITSTPCeffselP->Divide(fHistPtITSTPCselnocutP);
  TH1F *fHistPtITSTPCeffselS = (TH1F*)fHistPtITSTPCselS->Clone("fHistPtITSTPCeffselS");
  fHistPtITSTPCeffselS->Divide(fHistPtITSTPCselnocutS);
  TH1F *fHistPtITSTPCeffselSfromStrange = (TH1F*)fHistPtITSTPCselSfromStrange->Clone("fHistPtITSTPCeffselSfromStrange");
  fHistPtITSTPCeffselSfromStrange->Divide(fHistPtITSTPCselnocutSfromStrange);
  TH1F *fHistPtITSTPCeffselSfromMat = (TH1F*)fHistPtITSTPCselSfromMat->Clone("fHistPtITSTPCeffselSfromMat");
  fHistPtITSTPCeffselSfromMat->Divide(fHistPtITSTPCselnocutSfromMat);

  TH1F *fHistPtTPCeffInAccP = (TH1F*)fHistPtTPCInAccP->Clone("fHistPtTPCeffInAccP");
  fHistPtTPCeffInAccP->Divide(fHistPtTPCInAccnocutP);
  TH1F *fHistPtTPCeffInAccS = (TH1F*)fHistPtTPCInAccS->Clone("fHistPtTPCeffInAccS");
  fHistPtTPCeffInAccS->Divide(fHistPtTPCInAccnocutS);
  TH1F *fHistPtTPCeffInAccSfromStrange = (TH1F*)fHistPtTPCInAccSfromStrange->Clone("fHistPtTPCeffInAccSfromStrange");
  fHistPtTPCeffInAccSfromStrange->Divide(fHistPtTPCInAccnocutSfromStrange);
  TH1F *fHistPtTPCeffInAccSfromMat = (TH1F*)fHistPtTPCInAccSfromMat->Clone("fHistPtTPCeffInAccSfromMat");
  fHistPtTPCeffInAccSfromMat->Divide(fHistPtTPCInAccnocutSfromMat);


  TLegend *l0=new TLegend(0.5,0.5,0.9,0.9);
  TCanvas *c0=new TCanvas("c0","c0");
  c0->Divide(2,1);
  c0->cd(1);
  TH1F *fHistPtITSTPCselFracS = (TH1F*)fHistPtITSTPCselS->Clone("fHistPtITSTPCselFracS");
  fHistPtITSTPCselFracS->Divide(fHistPtITSTPCsel);
  fHistPtITSTPCselFracS->SetLineColor(2);
  l0->AddEntry(fHistPtITSTPCselFracS,"secondaries","l");
  fHistPtITSTPCselFracS->Draw();
  TH1F *fHistPtITSTPCselFracSfromStrange = (TH1F*)fHistPtITSTPCselSfromStrange->Clone("fHistPtITSTPCselFracSfromStrange");
  fHistPtITSTPCselFracSfromStrange->Divide(fHistPtITSTPCsel);
  fHistPtITSTPCselFracSfromStrange->SetLineColor(6);
  l0->AddEntry(fHistPtITSTPCselFracSfromStrange,"sec. from strange","l");
  fHistPtITSTPCselFracSfromStrange->Draw("same");
  TH1F *fHistPtITSTPCselFracSfromMat = (TH1F*)fHistPtITSTPCselSfromMat->Clone("fHistPtITSTPCselFracSfromMat");
  fHistPtITSTPCselFracSfromMat->Divide(fHistPtITSTPCsel);
  fHistPtITSTPCselFracSfromMat->SetLineColor(kOrange+7);
  l0->AddEntry(fHistPtITSTPCselFracSfromMat,"sec. from material","l");
  fHistPtITSTPCselFracSfromMat->Draw("same");
  l0->Draw();
  c0->cd(2);
  TH1F *fHistPtTPCInAccFracS = (TH1F*)fHistPtTPCInAccS->Clone("fHistPtTPCInAccFracS");
  fHistPtTPCInAccFracS->Divide(fHistPtTPCInAcc);
  fHistPtTPCInAccFracS->SetLineColor(2);
  fHistPtTPCInAccFracS->Draw();
  TH1F *fHistPtTPCInAccFracSfromStrange = (TH1F*)fHistPtTPCInAccSfromStrange->Clone("fHistPtTPCInAccFracSfromStrange");
  fHistPtTPCInAccFracSfromStrange->Divide(fHistPtTPCInAcc);
  fHistPtTPCInAccFracSfromStrange->SetLineColor(6);
  fHistPtTPCInAccFracSfromStrange->Draw("same");
  TH1F *fHistPtTPCInAccFracSfromMat = (TH1F*)fHistPtTPCInAccSfromMat->Clone("fHistPtTPCInAccFracSfromMat");
  fHistPtTPCInAccFracSfromMat->Divide(fHistPtTPCInAcc);
  fHistPtTPCInAccFracSfromMat->SetLineColor(kOrange+7);
  fHistPtTPCInAccFracSfromMat->Draw("same");

  TLegend *l1=new TLegend(0.5,0.5,0.9,0.9);
  TCanvas *c1=new TCanvas("c1","c1");
  c1->Divide(2,1);
  c1->cd(1);
  fHistPtITSTPCeffselP->SetYTitle("dca cut efficiency");
  fHistPtITSTPCeffselP->SetLineColor(1);
  l1->AddEntry(fHistPtITSTPCeffselP,"primaries","l");
  fHistPtITSTPCeffselP->Draw();
  fHistPtITSTPCeffselS->SetLineColor(2);
  l1->AddEntry(fHistPtITSTPCeffselS,"secondaries","l");
  fHistPtITSTPCeffselS->Draw("same");
  fHistPtITSTPCeffselSfromStrange->SetLineColor(6);
  l1->AddEntry(fHistPtITSTPCeffselSfromStrange,"sec. from strange","l");
  fHistPtITSTPCeffselSfromStrange->Draw("same");
  fHistPtITSTPCeffselSfromMat->SetLineColor(kOrange+7);
  l1->AddEntry(fHistPtITSTPCeffselSfromMat,"sec. from material","l");
  fHistPtITSTPCeffselSfromMat->Draw("same");
  l1->Draw();
  c1->cd(2);
  fHistPtTPCeffInAccP->SetYTitle("dca cut efficiency");
  fHistPtTPCeffInAccP->SetLineColor(1);
  fHistPtTPCeffInAccP->Draw();
  fHistPtTPCeffInAccS->SetLineColor(2);
  fHistPtTPCeffInAccS->Draw("same");
  fHistPtTPCeffInAccSfromStrange->SetLineColor(6);
  fHistPtTPCeffInAccSfromStrange->Draw("same");
  fHistPtTPCeffInAccSfromMat->SetLineColor(kOrange+7);
  fHistPtTPCeffInAccSfromMat->Draw("same");



  Float_t weightStrange=1.;
  Float_t weightMat=1.;

  TH1F *fHistPtITSTPCselCocktail = MakeCocktail(fHistPtITSTPCselP,
						fHistPtITSTPCselPfromStrange,
						fHistPtITSTPCselS,
						fHistPtITSTPCselSfromStrange,
						fHistPtITSTPCselSfromMat,
						weightStrange,
						weightMat,
						"fHistPtITSTPCselCocktail");

  TH1F *fHistPtITSTPCselPCocktail = CorrectToP(fHistPtITSTPCselCocktail,
					       fHistPtITSTPCselSfromStrange,
					       fHistPtITSTPCselSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtITSTPCselPCocktail");

  TH1F *fHistPtITSTPCallPCocktail = CorrectToPnocut(fHistPtITSTPCselPCocktail,
						    fHistPtITSTPCeffselP,
						    "fHistPtITSTPCallPCocktail");


  TH1F *fHistPtTPCInAccCocktail = MakeCocktail(fHistPtTPCInAccP,
						fHistPtTPCInAccPfromStrange,
						fHistPtTPCInAccS,
						fHistPtTPCInAccSfromStrange,
						fHistPtTPCInAccSfromMat,
						weightStrange,
					       weightMat,
					       "fHistPtTPCInAccCocktail");

  TH1F *fHistPtTPCInAccPCocktail = CorrectToP(fHistPtTPCInAccCocktail,
					      fHistPtTPCInAccSfromStrange,
					      fHistPtTPCInAccSfromMat,
					      weightStrange,
					      weightMat,
					      "fHistPtTPCInAccPCocktail");

  TH1F *fHistPtTPCallPCocktail = CorrectToPnocut(fHistPtTPCInAccPCocktail,
						 fHistPtTPCeffInAccP,
						 "fHistPtTPCallPCocktail");



  TLegend *l2=new TLegend(0.5,0.5,0.9,0.9);
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(2,1);
  c2->cd(1);
  fHistPtITSTPCselCocktail->SetLineColor(4);
  l2->AddEntry(fHistPtITSTPCselCocktail,"all, with dca cut","l");
  fHistPtITSTPCselCocktail->Draw();
  fHistPtITSTPCselPCocktail->SetLineColor(7);
  l2->AddEntry(fHistPtITSTPCselPCocktail,"primaries, with dca cut","l");
  fHistPtITSTPCselPCocktail->Draw("same");
  fHistPtITSTPCallPCocktail->SetLineColor(6);
  l2->AddEntry(fHistPtITSTPCallPCocktail,"primaries, w/o dca cut","l");
  fHistPtITSTPCallPCocktail->Draw("same");
  l2->Draw();
  c2->cd(2);
  fHistPtTPCInAccCocktail->SetLineColor(4);
  fHistPtTPCInAccCocktail->Draw();
  fHistPtTPCInAccPCocktail->SetLineColor(7);
  fHistPtTPCInAccPCocktail->Draw("same");
  fHistPtTPCallPCocktail->SetLineColor(6);
  fHistPtTPCallPCocktail->Draw("same");


  // systematic from secondaries
  weightStrange=1.3;
  weightMat=1.0;
  TH1F *fHistPtITSTPCselPCocktail1 = CorrectToP(fHistPtITSTPCselCocktail,
					       fHistPtITSTPCselSfromStrange,
					       fHistPtITSTPCselSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtITSTPCselPCocktail1");
  TH1F *fHistPtITSTPCallPCocktail1 = CorrectToPnocut(fHistPtITSTPCselPCocktail1,
						    fHistPtITSTPCeffselP,
						    "fHistPtITSTPCallPCocktail1");
  weightStrange=0.7;
  weightMat=1.0;
  TH1F *fHistPtITSTPCselPCocktail2 = CorrectToP(fHistPtITSTPCselCocktail,
					       fHistPtITSTPCselSfromStrange,
					       fHistPtITSTPCselSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtITSTPCselPCocktail2");
  TH1F *fHistPtITSTPCallPCocktail2 = CorrectToPnocut(fHistPtITSTPCselPCocktail2,
						    fHistPtITSTPCeffselP,
						    "fHistPtITSTPCallPCocktail2");
  weightStrange=1.0;
  weightMat=1.1;
  TH1F *fHistPtITSTPCselPCocktail3 = CorrectToP(fHistPtITSTPCselCocktail,
					       fHistPtITSTPCselSfromStrange,
					       fHistPtITSTPCselSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtITSTPCselPCocktail3");
  TH1F *fHistPtITSTPCallPCocktail3 = CorrectToPnocut(fHistPtITSTPCselPCocktail3,
						    fHistPtITSTPCeffselP,
						    "fHistPtITSTPCallPCocktail3");
  weightStrange=1.0;
  weightMat=0.9;
  TH1F *fHistPtITSTPCselPCocktail4 = CorrectToP(fHistPtITSTPCselCocktail,
					       fHistPtITSTPCselSfromStrange,
					       fHistPtITSTPCselSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtITSTPCselPCocktail4");
  TH1F *fHistPtITSTPCallPCocktail4 = CorrectToPnocut(fHistPtITSTPCselPCocktail4,
						    fHistPtITSTPCeffselP,
						    "fHistPtITSTPCallPCocktail4");
  weightStrange=1.3;
  weightMat=1.1;
  TH1F *fHistPtITSTPCselPCocktail5 = CorrectToP(fHistPtITSTPCselCocktail,
					       fHistPtITSTPCselSfromStrange,
					       fHistPtITSTPCselSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtITSTPCselPCocktail5");
  TH1F *fHistPtITSTPCallPCocktail5 = CorrectToPnocut(fHistPtITSTPCselPCocktail5,
						    fHistPtITSTPCeffselP,
						    "fHistPtITSTPCallPCocktail5");
  weightStrange=0.7;
  weightMat=0.9;
  TH1F *fHistPtITSTPCselPCocktail6 = CorrectToP(fHistPtITSTPCselCocktail,
					       fHistPtITSTPCselSfromStrange,
					       fHistPtITSTPCselSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtITSTPCselPCocktail6");
  TH1F *fHistPtITSTPCallPCocktail6 = CorrectToPnocut(fHistPtITSTPCselPCocktail6,
						    fHistPtITSTPCeffselP,
						    "fHistPtITSTPCallPCocktail6");

  weightStrange=1.3;
  weightMat=1.0;
  TH1F *fHistPtTPCInAccPCocktail1 = CorrectToP(fHistPtTPCInAccCocktail,
					       fHistPtTPCInAccSfromStrange,
					       fHistPtTPCInAccSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtTPCInAccPCocktail1");
  TH1F *fHistPtTPCallPCocktail1 = CorrectToPnocut(fHistPtTPCInAccPCocktail1,
						    fHistPtTPCeffInAccP,
						    "fHistPtTPCallPCocktail1");
  weightStrange=0.7;
  weightMat=1.0;
  TH1F *fHistPtTPCInAccPCocktail2 = CorrectToP(fHistPtTPCInAccCocktail,
					       fHistPtTPCInAccSfromStrange,
					       fHistPtTPCInAccSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtTPCInAccPCocktail2");
  TH1F *fHistPtTPCallPCocktail2 = CorrectToPnocut(fHistPtTPCInAccPCocktail2,
						    fHistPtTPCeffInAccP,
						    "fHistPtTPCallPCocktail2");
  weightStrange=1.0;
  weightMat=1.1;
  TH1F *fHistPtTPCInAccPCocktail3 = CorrectToP(fHistPtTPCInAccCocktail,
					       fHistPtTPCInAccSfromStrange,
					       fHistPtTPCInAccSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtTPCInAccPCocktail3");
  TH1F *fHistPtTPCallPCocktail3 = CorrectToPnocut(fHistPtTPCInAccPCocktail3,
						    fHistPtTPCeffInAccP,
						    "fHistPtTPCallPCocktail3");
  weightStrange=1.0;
  weightMat=0.9;
  TH1F *fHistPtTPCInAccPCocktail4 = CorrectToP(fHistPtTPCInAccCocktail,
					       fHistPtTPCInAccSfromStrange,
					       fHistPtTPCInAccSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtTPCInAccPCocktail4");
  TH1F *fHistPtTPCallPCocktail4 = CorrectToPnocut(fHistPtTPCInAccPCocktail4,
						    fHistPtTPCeffInAccP,
						    "fHistPtTPCallPCocktail4");
  weightStrange=1.3;
  weightMat=1.1;
  TH1F *fHistPtTPCInAccPCocktail5 = CorrectToP(fHistPtTPCInAccCocktail,
					       fHistPtTPCInAccSfromStrange,
					       fHistPtTPCInAccSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtTPCInAccPCocktail5");
  TH1F *fHistPtTPCallPCocktail5 = CorrectToPnocut(fHistPtTPCInAccPCocktail5,
						    fHistPtTPCeffInAccP,
						    "fHistPtTPCallPCocktail5");
  weightStrange=0.7;
  weightMat=0.9;
  TH1F *fHistPtTPCInAccPCocktail6 = CorrectToP(fHistPtTPCInAccCocktail,
					       fHistPtTPCInAccSfromStrange,
					       fHistPtTPCInAccSfromMat,
					       weightStrange,
					       weightMat,
					       "fHistPtTPCInAccPCocktail6");
  TH1F *fHistPtTPCallPCocktail6 = CorrectToPnocut(fHistPtTPCInAccPCocktail6,
						    fHistPtTPCeffInAccP,
						    "fHistPtTPCallPCocktail6");

  TLegend *l3=new TLegend(0.5,0.5,0.9,0.9);
  TCanvas *c3 = new TCanvas("c3","c3");
  c3->Divide(2,1);
  c3->cd(1);
  l3->AddEntry(fHistPtITSTPCallPCocktail,"central","l");
  fHistPtITSTPCallPCocktail->Draw();
  l3->AddEntry(fHistPtITSTPCallPCocktail1,"strange #pm 30%","l");
  fHistPtITSTPCallPCocktail1->SetLineColor(1);
  fHistPtITSTPCallPCocktail2->SetLineColor(1);
  fHistPtITSTPCallPCocktail1->Draw("same");
  fHistPtITSTPCallPCocktail2->Draw("same");
  l3->AddEntry(fHistPtITSTPCallPCocktail3,"material #pm 10%","l");
  fHistPtITSTPCallPCocktail3->SetLineColor(8);
  fHistPtITSTPCallPCocktail4->SetLineColor(8);
  fHistPtITSTPCallPCocktail3->Draw("same");
  fHistPtITSTPCallPCocktail4->Draw("same");
  l3->AddEntry(fHistPtITSTPCallPCocktail5,"strange #pm 30% & material #pm 10%","l");
  fHistPtITSTPCallPCocktail5->SetLineColor(4);
  fHistPtITSTPCallPCocktail6->SetLineColor(4);
  fHistPtITSTPCallPCocktail5->Draw("same");
  fHistPtITSTPCallPCocktail6->Draw("same");
  l3->Draw();
  c3->cd(2);
  fHistPtTPCallPCocktail->Draw();
  fHistPtTPCallPCocktail1->SetLineColor(1);
  fHistPtTPCallPCocktail2->SetLineColor(1);
  fHistPtTPCallPCocktail1->Draw("same");
  fHistPtTPCallPCocktail2->Draw("same");
  fHistPtTPCallPCocktail3->SetLineColor(8);
  fHistPtTPCallPCocktail4->SetLineColor(8);
  fHistPtTPCallPCocktail3->Draw("same");
  fHistPtTPCallPCocktail4->Draw("same");
  fHistPtTPCallPCocktail5->SetLineColor(4);
  fHistPtTPCallPCocktail6->SetLineColor(4);
  fHistPtTPCallPCocktail5->Draw("same");
  fHistPtTPCallPCocktail6->Draw("same");

  /*
  TCanvas *c4 = new TCanvas("c4","c4");
  c4->Divide(2,1);
  c4->cd(1);
  fHistPtITSTPCallPCocktail2->Add(fHistPtITSTPCallPCocktail,-1);
  fHistPtITSTPCallPCocktail2->Divide(fHistPtITSTPCallPCocktail);
  fHistPtITSTPCallPCocktail2->Draw();
  fHistPtITSTPCallPCocktail4->Add(fHistPtITSTPCallPCocktail,-1);
  fHistPtITSTPCallPCocktail4->Divide(fHistPtITSTPCallPCocktail);
  fHistPtITSTPCallPCocktail4->Draw("same");
  fHistPtITSTPCallPCocktail6->Add(fHistPtITSTPCallPCocktail,-1);
  fHistPtITSTPCallPCocktail6->Divide(fHistPtITSTPCallPCocktail);
  fHistPtITSTPCallPCocktail6->Draw("same");
  l3->Draw();
  c4->cd(2);
  fHistPtTPCallPCocktail2->Add(fHistPtTPCallPCocktail,-1);
  fHistPtTPCallPCocktail2->Divide(fHistPtTPCallPCocktail);
  fHistPtTPCallPCocktail2->Draw();
  fHistPtTPCallPCocktail4->Add(fHistPtTPCallPCocktail,-1);
  fHistPtTPCallPCocktail4->Divide(fHistPtTPCallPCocktail);
  fHistPtTPCallPCocktail4->Draw("same");
  fHistPtTPCallPCocktail6->Add(fHistPtTPCallPCocktail,-1);
  fHistPtTPCallPCocktail6->Divide(fHistPtTPCallPCocktail);
  fHistPtTPCallPCocktail6->Draw("same");
  */

  return;
}
//-------------------------------------------------------------
void NormalizePtHist(TH1F *h,Float_t n) {

  for(Int_t i=1; i<=h->GetNbinsX(); i++) {
    if(h->GetBinContent(i)) {
      h->SetBinContent(i,h->GetBinContent(i)/h->GetBinWidth(i));
      h->SetBinError(i,h->GetBinError(i)/h->GetBinWidth(i));
    }
  }
  h->Scale(1./n);

  return;
}
//----------------------------------------------------------------------
TH1F *MakeCocktail(TH1F *fHistPtselP,
		   TH1F *fHistPtselPfromStrange,
		   TH1F *fHistPtselS,
		   TH1F *fHistPtselSfromStrange,
		   TH1F *fHistPtselSfromMat,
		   Float_t weightStrange,
		   Float_t weightMat,TString name) {

  // build histogram with yield passing cut
  TH1F *fHistPtselCocktail = (TH1F*)fHistPtselP->Clone(name.Data());
  fHistPtselCocktail->Add(fHistPtselS);
  fHistPtselCocktail->Add(fHistPtselPfromStrange,weightStrange-1.);
  fHistPtselCocktail->Add(fHistPtselSfromStrange,weightStrange-1.);
  fHistPtselCocktail->Add(fHistPtselSfromMat,weightMat-1.);


  return fHistPtselCocktail;
}
//----------------------------------------------------------------------
TH1F *CorrectToP(TH1F *fHistPtselCocktail,
		TH1F *fHistPtselSfromStrange,
		TH1F *fHistPtselSfromMat,
		Float_t weightStrange,
		Float_t weightMat,TString name) {


  // build histogram with yield of primaries passing cut
  TH1F* fHistPtselPCocktail = (TH1F*)fHistPtselCocktail->Clone(name.Data());
  fHistPtselPCocktail->Add(fHistPtselSfromStrange,-weightStrange);
  fHistPtselPCocktail->Add(fHistPtselSfromMat,-weightMat);

  return fHistPtselPCocktail;
}
//----------------------------------------------------------------------
TH1F *CorrectToPnocut(TH1F *fHistPtselPCocktail,
		     TH1F *fHistPteffselP,TString name) {

  // build histogram with yield of primaries without cut
  TH1F *fHistPtallPCocktail = (TH1F*)fHistPtselPCocktail->Clone(name.Data());
  //printf("%f\n",fHistPtallPCocktail->GetBinContent(6));
  fHistPtallPCocktail->Divide(fHistPteffselP);
  //printf("%f\n",fHistPtallPCocktail->GetBinContent(6));

  return fHistPtallPCocktail;
}
//----------------------------------------------------------------------
void CompareTPCPt() {

  TFile *f_tpcrefit= new TFile("ITS.Performance_104892pass2_etale05.root");

  TList *list=(TList*)f_tpcrefit->Get("cOutputITS");

  TH1F *fHistNEvents = (TH1F*)list->FindObject("fHistNEvents");
  Float_t nEventsSelWithVertex = fHistNEvents->GetBinContent(5);

  TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  NormalizePtHist(fHistPtTPCInAcc,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccNoTRDout = (TH1F*)list->FindObject("fHistPtTPCInAccNoTRDout");
  NormalizePtHist(fHistPtTPCInAccNoTRDout,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccNoTOFout = (TH1F*)list->FindObject("fHistPtTPCInAccNoTOFout");
  NormalizePtHist(fHistPtTPCInAccNoTOFout,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccWithPtTPCAtVtx = (TH1F*)list->FindObject("fHistPtTPCInAccWithPtTPCAtVtx");
  NormalizePtHist(fHistPtTPCInAccWithPtTPCAtVtx,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccWithPtTPCAtInnerWall = (TH1F*)list->FindObject("fHistPtTPCInAccWithPtTPCAtInnerWall");
  NormalizePtHist(fHistPtTPCInAccWithPtTPCAtInnerWall,nEventsSelWithVertex);


  TFile *f_notpcrefit= new TFile("ITS.Performance_104892pass3_etale05.root");

  TList *list=(TList*)f_notpcrefit->Get("cOutputITS");

  TH1F *fHistNEventsNoRefit = (TH1F*)list->FindObject("fHistNEvents");
  nEventsSelWithVertex = fHistNEventsNoRefit->GetBinContent(5);

  TH1F *fHistPtTPCInAccNoRefit = (TH1F*)list->FindObject("fHistPtTPCInAcc");
  NormalizePtHist(fHistPtTPCInAccNoRefit,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccNoTRDoutNoRefit = (TH1F*)list->FindObject("fHistPtTPCInAccNoTRDout");
  NormalizePtHist(fHistPtTPCInAccNoTRDoutNoRefit,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccNoTOFoutNoRefit = (TH1F*)list->FindObject("fHistPtTPCInAccNoTOFout");
  NormalizePtHist(fHistPtTPCInAccNoTOFoutNoRefit,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccWithPtTPCAtVtxNoRefit = (TH1F*)list->FindObject("fHistPtTPCInAccWithPtTPCAtVtx");
  NormalizePtHist(fHistPtTPCInAccWithPtTPCAtVtxNoRefit,nEventsSelWithVertex);
  TH1F *fHistPtTPCInAccWithPtTPCAtInnerWallNoRefit = (TH1F*)list->FindObject("fHistPtTPCInAccWithPtTPCAtInnerWall");
  NormalizePtHist(fHistPtTPCInAccWithPtTPCAtInnerWallNoRefit,nEventsSelWithVertex);



  TLegend *l1 = new TLegend(0.5,0.5,0.9,0.9);
  TCanvas *c1=new TCanvas("c1","c1",0,0,1000,600);
  c1->Divide(2,1);
  c1->cd(1);
  fHistPtTPCInAcc->SetMarkerStyle(20);
  l1->AddEntry(fHistPtTPCInAcc,"p_{t} measurement:","");
  l1->AddEntry(fHistPtTPCInAcc,"at vtx, from main params","p");
  fHistPtTPCInAcc->Draw();
  fHistPtTPCInAccWithPtTPCAtVtx->SetMarkerStyle(21);
  fHistPtTPCInAccWithPtTPCAtVtx->SetMarkerColor(2);
  l1->AddEntry(fHistPtTPCInAccWithPtTPCAtVtx,"at vtx, from TPC params","p");
  fHistPtTPCInAccWithPtTPCAtVtx->Draw("same");
  fHistPtTPCInAccWithPtTPCAtInnerWall->SetMarkerStyle(25);
  fHistPtTPCInAccWithPtTPCAtInnerWall->SetMarkerColor(2);
  l1->AddEntry(fHistPtTPCInAccWithPtTPCAtInnerWall,"at r=80cm, from TPC params","p");
  fHistPtTPCInAccWithPtTPCAtInnerWall->Draw("same");
  fHistPtTPCInAccNoTOFout->SetMarkerStyle(23);
  fHistPtTPCInAccNoTOFout->SetMarkerColor(3);
  fHistPtTPCInAccNoTOFout->Draw("same");
  fHistPtTPCInAccNoTRDout->SetMarkerStyle(24);
  fHistPtTPCInAccNoTRDout->SetMarkerColor(4);
  fHistPtTPCInAccNoTRDout->Draw("same");
  l1->Draw();
  c1->cd(2);
  fHistPtTPCInAccNoRefit->SetMarkerStyle(20);
  fHistPtTPCInAccNoRefit->Draw();
  fHistPtTPCInAccWithPtTPCAtVtxNoRefit->SetMarkerStyle(21);
  fHistPtTPCInAccWithPtTPCAtVtxNoRefit->SetMarkerColor(2);
  fHistPtTPCInAccWithPtTPCAtVtxNoRefit->Draw("same");
  fHistPtTPCInAccWithPtTPCAtInnerWallNoRefit->SetMarkerStyle(25);
  fHistPtTPCInAccWithPtTPCAtInnerWallNoRefit->SetMarkerColor(2);
  fHistPtTPCInAccWithPtTPCAtInnerWallNoRefit->Draw("same");
  fHistPtTPCInAccNoTOFoutNoRefit->SetMarkerStyle(23);
  fHistPtTPCInAccNoTOFoutNoRefit->SetMarkerColor(3);
  fHistPtTPCInAccNoTOFoutNoRefit->Draw("same");
  fHistPtTPCInAccNoTRDoutNoRefit->SetMarkerStyle(24);
  fHistPtTPCInAccNoTRDoutNoRefit->SetMarkerColor(4);
  fHistPtTPCInAccNoTRDoutNoRefit->Draw("same");
  //l1->Draw();

  return;
}

TH1F*  gHistTemplateStrange;
TH1F*  gHistTemplatePandMat;
Float_t gd0min=0.15;
Float_t gd0max=1.5;


Double_t FitMCTemplate(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   if(TMath::Abs(xx)<gd0min) return 0;
   Int_t bin = gHistTemplateStrange->GetXaxis()->FindBin(xx);
   Double_t pandmat = par[0]*gHistTemplatePandMat->GetBinContent(bin);
   Double_t strange = par[1]*gHistTemplateStrange->GetBinContent(bin);
   return strange + pandmat;
}

void FitImpPar_rphi() {

  Int_t rebin=1;

  TFile *fMC= new TFile("ITS.Performance_lhc10a8a12_trkvtx_noTPCdca_120310.root");

  TList *list=(TList*)fMC->Get("cOutputITS");
  TH1F *fHistd0rphiITSMIoneSPDInAccPMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP150200");
  TH1F *fHistd0rphiITSMIoneSPDInAccSMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS150200");
  TH1F *fHistd0rphiITSMIoneSPDInAccSfromMatMC = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS150200fromMat");
  fHistd0rphiITSMIoneSPDInAccPMC->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAccSMC->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAccSfromMatMC->Rebin(rebin);

  TH1F *fHistd0rphiITSMIoneSPDInAccSfromStrangeMC=(TH1F*)fHistd0rphiITSMIoneSPDInAccSMC->Clone("fHistd0rphiITSMIoneSPDInAccSfromStrangeMC");
  fHistd0rphiITSMIoneSPDInAccSfromStrangeMC->Add(fHistd0rphiITSMIoneSPDInAccSfromMatMC,-1.);
  TH1F *fHistd0rphiITSMIoneSPDInAccMC=(TH1F*)fHistd0rphiITSMIoneSPDInAccPMC->Clone("fHistd0rphiITSMIoneSPDInAccMC");
  fHistd0rphiITSMIoneSPDInAccMC->Add(fHistd0rphiITSMIoneSPDInAccSMC);
  //fHistd0rphiITSMIoneSPDInAccMC->Scale(1./fHistd0rphiITSMIoneSPDInAccMC->GetEntries());
  TH1F *fHistd0rphiITSMIoneSPDInAccPandMatMC=(TH1F*)fHistd0rphiITSMIoneSPDInAccPMC->Clone("fHistd0rphiITSMIoneSPDInAccPandMatMC");
  fHistd0rphiITSMIoneSPDInAccPandMatMC->Add(fHistd0rphiITSMIoneSPDInAccSfromMatMC);

  // templates
  gHistTemplateStrange = (TH1F*)fHistd0rphiITSMIoneSPDInAccSfromStrangeMC->Clone("gHistTemplateStrange");
  gHistTemplatePandMat = (TH1F*)fHistd0rphiITSMIoneSPDInAccPandMatMC->Clone("gHistTemplatePandMat");

  


  TFile *f= new TFile("ITS.Performance_104065_104892_pass4_trkvtx_noTPCdca.root");
  //TFile *f= new TFile("ITS.Performance_lhc10a8_300k_trkvtx_noTPCdca.root");

  TList *list=(TList*)f->Get("cOutputITS");
  TH1F *fHistd0rphiITSMIoneSPDInAccP = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccP150200");
  TH1F *fHistd0rphiITSMIoneSPDInAccS = (TH1F*)list->FindObject("fHistd0rphiITSMIoneSPDInAccS150200");
  fHistd0rphiITSMIoneSPDInAccP->Rebin(rebin);
  fHistd0rphiITSMIoneSPDInAccS->Rebin(rebin);

  TH1F *fHistd0rphiITSMIoneSPDInAcc=(TH1F*)fHistd0rphiITSMIoneSPDInAccP->Clone("fHistd0rphiITSMIoneSPDInAcc");
  fHistd0rphiITSMIoneSPDInAcc->Add(fHistd0rphiITSMIoneSPDInAccS);


  TF1 *fitMCtemplate = new TF1("fitMCtemplate",FitMCTemplate,-gd0max,+gd0max,2);

  TH1F* fHistd0rphiITSMIoneSPDInAccCopy = (TH1F*)fHistd0rphiITSMIoneSPDInAcc->Clone("fHistd0rphiITSMIoneSPDInAccCopy");
  for(Int_t ib=1;ib<=fHistd0rphiITSMIoneSPDInAccCopy->GetNbinsX();ib++) {
    if(TMath::Abs(fHistd0rphiITSMIoneSPDInAccCopy->GetBinCenter(ib))<gd0min) {
      fHistd0rphiITSMIoneSPDInAccCopy->SetBinContent(ib,0);
      fHistd0rphiITSMIoneSPDInAccCopy->SetBinError(ib,0);
    }
  }
  fHistd0rphiITSMIoneSPDInAccCopy->Fit("fitMCtemplate","RL");


  fHistd0rphiITSMIoneSPDInAccPandMatMC->Scale(fitMCtemplate->GetParameter(0));
  fHistd0rphiITSMIoneSPDInAccSfromStrangeMC->Scale(fitMCtemplate->GetParameter(1));
  TH1F *fHistd0rphiITSMIoneSPDInAccMCFit=(TH1F*)fHistd0rphiITSMIoneSPDInAccPandMatMC->Clone("fHistd0rphiITSMIoneSPDInAccMCFit");
  fHistd0rphiITSMIoneSPDInAccMCFit->Add(fHistd0rphiITSMIoneSPDInAccSfromStrangeMC);

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->SetLogy();
  fHistd0rphiITSMIoneSPDInAccCopy->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccCopy->SetMaximum(2*fHistd0rphiITSMIoneSPDInAcc->GetMaximum());
  fHistd0rphiITSMIoneSPDInAcc->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccSfromStrangeMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccPandMatMC->SetMinimum(1);
  fHistd0rphiITSMIoneSPDInAccCopy->SetLineColor(4);
  fHistd0rphiITSMIoneSPDInAccCopy->Draw();
  fHistd0rphiITSMIoneSPDInAcc->SetLineColor(4);
  fHistd0rphiITSMIoneSPDInAcc->Draw("same");
  fHistd0rphiITSMIoneSPDInAccMCFit->SetLineColor(2);
  fHistd0rphiITSMIoneSPDInAccMCFit->Draw("same");
  fHistd0rphiITSMIoneSPDInAccSfromStrangeMC->SetLineColor(3);
  fHistd0rphiITSMIoneSPDInAccSfromStrangeMC->Draw("same");
  fHistd0rphiITSMIoneSPDInAccPandMatMC->SetLineColor(15);
  fHistd0rphiITSMIoneSPDInAccPandMatMC->Draw("same");



  TCanvas *c3 = new TCanvas("c3","c3");
  c3->SetLogx();
  c3->SetGridx();
  c3->SetGridy();
  Float_t d0cut[15]={0.0301,0.0401,0.0501,0.0601,0.0801,0.101,0.151,0.201,0.251,0.301,0.401,0.601,0.801,1.001,1.401};

  Float_t fracSfromStrange[15],fracSfromStrangeMC[15],fracSfromMatMC[15],fracS[15],fracSMC[15],fracSRatio[15];
  for(Int_t i=0;i<15;i++) {
    Int_t bin1=fHistd0rphiITSMIoneSPDInAccPandMatMC->FindBin(-d0cut[i]);
    Int_t bin2=fHistd0rphiITSMIoneSPDInAccPandMatMC->FindBin(+d0cut[i]);
    Float_t intDatacut=fHistd0rphiITSMIoneSPDInAcc->Integral(bin1,bin2);
    Float_t intMCcut=fHistd0rphiITSMIoneSPDInAccMC->Integral(bin1,bin2);
    Float_t intSfromStrangecut=fHistd0rphiITSMIoneSPDInAccSfromStrangeMC->Integral(bin1,bin2);
    Float_t intSfromMatcut=fHistd0rphiITSMIoneSPDInAccSfromMatMC->Integral(bin1,bin2);
    fracSfromStrange[i]=intSfromStrangecut/intDatacut;
    fracSfromStrangeMC[i]=intSfromStrangecut/fitMCtemplate->GetParameter(1)/intMCcut;
    fracSfromMatMC[i]=intSfromMatcut/intMCcut;

    fracS[i]=fracSfromStrange[i]+fracSfromMatMC[i];
    fracSMC[i]=fracSfromStrangeMC[i]+fracSfromMatMC[i];

    fracSRatio[i]=fracS[i]/fracSMC[i];
  }
  TGraph *gfracSfromStrange=new TGraph(15,d0cut,fracSfromStrange);
  gfracSfromStrange->SetMarkerColor(2);
  gfracSfromStrange->SetMarkerStyle(20);
  gfracSfromStrange->Draw("ap");

  TGraph *gfracSfromStrangeMC=new TGraph(15,d0cut,fracSfromStrangeMC);
  gfracSfromStrangeMC->SetMarkerColor(1);
  gfracSfromStrangeMC->SetMarkerStyle(22);
  gfracSfromStrangeMC->Draw("p");

  TGraph *gfracSfromMatMC=new TGraph(15,d0cut,fracSfromMatMC);
  gfracSfromMatMC->SetMarkerColor(4);
  gfracSfromMatMC->SetMarkerStyle(21);
  gfracSfromMatMC->Draw("p");

  TGraph *gfracS=new TGraph(15,d0cut,fracS);
  gfracS->SetMarkerColor(2);
  gfracS->SetMarkerStyle(24);
  gfracS->Draw("p");

  TGraph *gfracSMC=new TGraph(15,d0cut,fracSMC);
  gfracSMC->SetMarkerColor(1);
  gfracSMC->SetMarkerStyle(25);
  gfracSMC->Draw("p");


  TCanvas *c4 = new TCanvas("c4","c4");
  c4->SetLogx();
  c4->SetGridx();
  c4->SetGridy();

  TGraph *gfracSRatio=new TGraph(15,d0cut,fracSRatio);
  gfracSRatio->SetMarkerColor(9);
  gfracSRatio->SetMarkerStyle(21);
  gfracSRatio->Draw("ap");


  return;
}
//-----------------------------------------------------------------------------
void PlotEffVSErrors() {

  TFile *f1= new TFile("eff_104892_pd_testnew1.root");

  TH1F *fHistPtITSMI6InAcctestnew1 = (TH1F*)f1->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMIge2InAcctestnew1 = (TH1F*)f1->Get("fHistPtITSMIge2InAcc");
  TH1F *fHistPtITSMISPDInAcctestnew1 = (TH1F*)f1->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcctestnew1 = (TH1F*)f1->Get("fHistPtITSMIoneSPDInAcc");

  fHistPtITSMI6InAcctestnew1->SetMarkerStyle(20);
  fHistPtITSMIge2InAcctestnew1->SetMarkerStyle(20);
  fHistPtITSMISPDInAcctestnew1->SetMarkerStyle(20);
  fHistPtITSMIoneSPDInAcctestnew1->SetMarkerStyle(20);
  fHistPtITSMI6InAcctestnew1->SetMarkerColor(2);
  fHistPtITSMIge2InAcctestnew1->SetMarkerColor(1);
  fHistPtITSMISPDInAcctestnew1->SetMarkerColor(9);
  fHistPtITSMIoneSPDInAcctestnew1->SetMarkerColor(15);


  TFile *f2= new TFile("eff_104892_pd_testnew2.root");

  TH1F *fHistPtITSMI6InAcctestnew2 = (TH1F*)f2->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMIge2InAcctestnew2 = (TH1F*)f2->Get("fHistPtITSMIge2InAcc");
  TH1F *fHistPtITSMISPDInAcctestnew2 = (TH1F*)f2->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcctestnew2 = (TH1F*)f2->Get("fHistPtITSMIoneSPDInAcc");

  fHistPtITSMI6InAcctestnew2->SetMarkerStyle(21);
  fHistPtITSMIge2InAcctestnew2->SetMarkerStyle(21);
  fHistPtITSMISPDInAcctestnew2->SetMarkerStyle(21);
  fHistPtITSMIoneSPDInAcctestnew2->SetMarkerStyle(21);
  fHistPtITSMI6InAcctestnew2->SetMarkerColor(2);
  fHistPtITSMIge2InAcctestnew2->SetMarkerColor(1);
  fHistPtITSMISPDInAcctestnew2->SetMarkerColor(9);
  fHistPtITSMIoneSPDInAcctestnew2->SetMarkerColor(15);


  TFile *f3= new TFile("eff_104892_pd_testnew3.root");

  TH1F *fHistPtITSMI6InAcctestnew3 = (TH1F*)f3->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMIge2InAcctestnew3 = (TH1F*)f3->Get("fHistPtITSMIge2InAcc");
  TH1F *fHistPtITSMISPDInAcctestnew3 = (TH1F*)f3->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcctestnew3 = (TH1F*)f3->Get("fHistPtITSMIoneSPDInAcc");

  fHistPtITSMI6InAcctestnew3->SetMarkerStyle(22);
  fHistPtITSMIge2InAcctestnew3->SetMarkerStyle(22);
  fHistPtITSMISPDInAcctestnew3->SetMarkerStyle(22);
  fHistPtITSMIoneSPDInAcctestnew3->SetMarkerStyle(22);
  fHistPtITSMI6InAcctestnew3->SetMarkerColor(2);
  fHistPtITSMIge2InAcctestnew3->SetMarkerColor(1);
  fHistPtITSMISPDInAcctestnew3->SetMarkerColor(9);
  fHistPtITSMIoneSPDInAcctestnew3->SetMarkerColor(15);


  TFile *f5= new TFile("eff_104892_pd_testnew5.root");

  TH1F *fHistPtITSMI6InAcctestnew5 = (TH1F*)f5->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMIge2InAcctestnew5 = (TH1F*)f5->Get("fHistPtITSMIge2InAcc");
  TH1F *fHistPtITSMISPDInAcctestnew5 = (TH1F*)f5->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcctestnew5 = (TH1F*)f5->Get("fHistPtITSMIoneSPDInAcc");

  fHistPtITSMI6InAcctestnew5->SetMarkerStyle(24);
  fHistPtITSMIge2InAcctestnew5->SetMarkerStyle(24);
  fHistPtITSMISPDInAcctestnew5->SetMarkerStyle(24);
  fHistPtITSMIoneSPDInAcctestnew5->SetMarkerStyle(24);
  fHistPtITSMI6InAcctestnew5->SetMarkerColor(2);
  fHistPtITSMIge2InAcctestnew5->SetMarkerColor(1);
  fHistPtITSMISPDInAcctestnew5->SetMarkerColor(9);
  fHistPtITSMIoneSPDInAcctestnew5->SetMarkerColor(15);


  TFile *f6= new TFile("eff_104892_pd_testnew6.root");

  TH1F *fHistPtITSMI6InAcctestnew6 = (TH1F*)f6->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMIge2InAcctestnew6 = (TH1F*)f6->Get("fHistPtITSMIge2InAcc");
  TH1F *fHistPtITSMISPDInAcctestnew6 = (TH1F*)f6->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcctestnew6 = (TH1F*)f6->Get("fHistPtITSMIoneSPDInAcc");

  fHistPtITSMI6InAcctestnew6->SetMarkerStyle(25);
  fHistPtITSMIge2InAcctestnew6->SetMarkerStyle(25);
  fHistPtITSMISPDInAcctestnew6->SetMarkerStyle(25);
  fHistPtITSMIoneSPDInAcctestnew6->SetMarkerStyle(25);
  fHistPtITSMI6InAcctestnew6->SetMarkerColor(2);
  fHistPtITSMIge2InAcctestnew6->SetMarkerColor(1);
  fHistPtITSMISPDInAcctestnew6->SetMarkerColor(9);
  fHistPtITSMIoneSPDInAcctestnew6->SetMarkerColor(15);


  TFile *f7= new TFile("eff_104892_pd_testnew7.root");

  TH1F *fHistPtITSMI6InAcctestnew7 = (TH1F*)f7->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMIge2InAcctestnew7 = (TH1F*)f7->Get("fHistPtITSMIge2InAcc");
  TH1F *fHistPtITSMISPDInAcctestnew7 = (TH1F*)f7->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcctestnew7 = (TH1F*)f7->Get("fHistPtITSMIoneSPDInAcc");

  fHistPtITSMI6InAcctestnew7->SetMarkerStyle(23);
  fHistPtITSMIge2InAcctestnew7->SetMarkerStyle(23);
  fHistPtITSMISPDInAcctestnew7->SetMarkerStyle(23);
  fHistPtITSMIoneSPDInAcctestnew7->SetMarkerStyle(23);
  fHistPtITSMI6InAcctestnew7->SetMarkerColor(2);
  fHistPtITSMIge2InAcctestnew7->SetMarkerColor(1);
  fHistPtITSMISPDInAcctestnew7->SetMarkerColor(9);
  fHistPtITSMIoneSPDInAcctestnew7->SetMarkerColor(15);


  TFile *f8= new TFile("eff_104892_pd_testnew8.root");

  TH1F *fHistPtITSMI6InAcctestnew8 = (TH1F*)f8->Get("fHistPtITSMI6InAcc");
  TH1F *fHistPtITSMIge2InAcctestnew8 = (TH1F*)f8->Get("fHistPtITSMIge2InAcc");
  TH1F *fHistPtITSMISPDInAcctestnew8 = (TH1F*)f8->Get("fHistPtITSMISPDInAcc");
  TH1F *fHistPtITSMIoneSPDInAcctestnew8 = (TH1F*)f8->Get("fHistPtITSMIoneSPDInAcc");

  fHistPtITSMI6InAcctestnew8->SetMarkerStyle(28);
  fHistPtITSMIge2InAcctestnew8->SetMarkerStyle(28);
  fHistPtITSMISPDInAcctestnew8->SetMarkerStyle(28);
  fHistPtITSMIoneSPDInAcctestnew8->SetMarkerStyle(28);
  fHistPtITSMI6InAcctestnew8->SetMarkerColor(2);
  fHistPtITSMIge2InAcctestnew8->SetMarkerColor(1);
  fHistPtITSMISPDInAcctestnew8->SetMarkerColor(9);
  fHistPtITSMIoneSPDInAcctestnew8->SetMarkerColor(15);


  TCanvas *c= new TCanvas("c","c");
  fHistPtITSMI6InAcctestnew1->Draw();
  fHistPtITSMIge2InAcctestnew1->Draw("same");
  fHistPtITSMISPDInAcctestnew1->Draw("same");
  fHistPtITSMIoneSPDInAcctestnew1->Draw("same");
  fHistPtITSMI6InAcctestnew2->Draw("same");
  fHistPtITSMIge2InAcctestnew2->Draw("same");
  fHistPtITSMISPDInAcctestnew2->Draw("same");
  fHistPtITSMIoneSPDInAcctestnew2->Draw("same");
  fHistPtITSMI6InAcctestnew3->Draw("same");
  fHistPtITSMIge2InAcctestnew3->Draw("same");
  fHistPtITSMISPDInAcctestnew3->Draw("same");
  fHistPtITSMIoneSPDInAcctestnew3->Draw("same");
  fHistPtITSMI6InAcctestnew5->Draw("same");
  fHistPtITSMIge2InAcctestnew5->Draw("same");
  fHistPtITSMISPDInAcctestnew5->Draw("same");
  fHistPtITSMIoneSPDInAcctestnew5->Draw("same");
  fHistPtITSMI6InAcctestnew6->Draw("same");
  fHistPtITSMIge2InAcctestnew6->Draw("same");
  fHistPtITSMISPDInAcctestnew6->Draw("same");
  fHistPtITSMIoneSPDInAcctestnew6->Draw("same");
  fHistPtITSMI6InAcctestnew7->Draw("same");
  fHistPtITSMIge2InAcctestnew7->Draw("same");
  fHistPtITSMISPDInAcctestnew7->Draw("same");
  fHistPtITSMIoneSPDInAcctestnew7->Draw("same");
  fHistPtITSMI6InAcctestnew8->Draw("same");
  fHistPtITSMIge2InAcctestnew8->Draw("same");
  fHistPtITSMISPDInAcctestnew8->Draw("same");
  fHistPtITSMIoneSPDInAcctestnew8->Draw("same");


  return;
}
//--------------------------------------------------------------------------
void ReweightStrange(TH1F *hPt,TH1F* hPtPfromStrange,TH1F* hPtSfromStrange) {

  for(Int_t i=1;i<=hPt->GetNbinsX();i++) {
    Double_t pt=hPt->GetBinCenter(i);
    Double_t weight=1.;
    Double_t weightP=1.;
    if(pt<.25) {
      weight=1.3;//1.;
    } else if(pt>=0.25 && pt<0.5) {
      weight=1.5;//1.3;
    } else if(pt>=0.5 && pt<1.1) {
      weight=1.7;//1.44;
    } else if(pt>=1.1) {
      weight=2;//1.7;
    }
    
    weightP*=1.3;
    hPt->SetBinContent(i,hPt->GetBinContent(i)+(weight-1.)*hPtSfromStrange->GetBinContent(i));
    if(hPtPfromStrange) hPt->SetBinContent(i,hPt->GetBinContent(i)+(weightP-1.)*hPtPfromStrange->GetBinContent(i));
  }

  return;
}
//--------------------------------------------------------------------------
void ITSTrackingTrending(Int_t firstrun=158508,Int_t lastrun=158511,
			 //TString pathBeforeRun="/alice/sim/LHC11a10a_bis/",
			 //TString pathAfterRun="/QA67/Stage_4/004/QAresults.root",
			 TString pathBeforeRun="/alice/data/2011/LHC11d/000",
			 TString pathAfterRun="/ESDs/pass1/QAresults.root",
			 //TString pathBeforeRun="/alice/data/2010/LHC10b/000",
			 //TString pathAfterRun="/ESDs/pass2/QA9/QAresults.root",
			 TString pathAfterRun2="") 
{
  //
  // Make ITS efficiency trending plots from QAresults.root files
  //
  gStyle->SetOptStat(0);

  TGrid::Connect("alien://");

  Float_t ioValues[100],ioErrors[100];
  Int_t index;
  Int_t nruns=lastrun-firstrun+1;


  TH1F *hFracSPD1 = new TH1F("hFracSPD1","SPD inner; run number; Fraction of HSs",nruns,firstrun-0.5,lastrun+0.5);
  hFracSPD1->SetLineColor(3);
  hFracSPD1->SetMarkerColor(3);
  hFracSPD1->SetMarkerStyle(20);
  TH1F *hFracSPD2 = new TH1F("hFracSPD2","SPD outer; run number; Fraction of HSs",nruns,firstrun-0.5,lastrun+0.5);
  hFracSPD2->SetLineColor(8);
  hFracSPD2->SetMarkerColor(8);


  TH1F *hEffge2Pt02 = new TH1F("hEffge2Pt02","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEffge2Pt02->SetLineWidth(2);
  hEffge2Pt02->SetLineColor(1);
  hEffge2Pt02->SetMarkerColor(1);
  hEffge2Pt02->SetMarkerStyle(20);
  TH1F *hEffge2Pt1 = new TH1F("hEffge2Pt1","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEffge2Pt1->SetLineWidth(2);
  hEffge2Pt1->SetLineColor(1);
  hEffge2Pt1->SetMarkerColor(1);
  hEffge2Pt1->SetMarkerStyle(20);
  TH1F *hEffge2Pt10 = new TH1F("hEffge2Pt10","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEffge2Pt10->SetLineWidth(2);
  hEffge2Pt10->SetLineColor(1);
  hEffge2Pt10->SetMarkerColor(1);
  hEffge2Pt10->SetMarkerStyle(20);

  TH1F *hEff6Pt02 = new TH1F("hEff6Pt02","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEff6Pt02->SetLineWidth(2);
  hEff6Pt02->SetLineColor(2);
  hEff6Pt02->SetMarkerColor(2);
  hEff6Pt02->SetMarkerStyle(20);
  TH1F *hEff6Pt1 = new TH1F("hEff6Pt1","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEff6Pt1->SetLineWidth(2);
  hEff6Pt1->SetLineColor(2);
  hEff6Pt1->SetMarkerColor(2);
  hEff6Pt1->SetMarkerStyle(20);
  TH1F *hEff6Pt10 = new TH1F("hEff6Pt10","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEff6Pt10->SetLineWidth(2);
  hEff6Pt10->SetLineColor(2);
  hEff6Pt10->SetMarkerColor(2);
  hEff6Pt10->SetMarkerStyle(20);

  TH1F *hEffSPDPt02 = new TH1F("hEffSPDPt02","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEffSPDPt02->SetLineWidth(2);
  hEffSPDPt02->SetLineColor(9);
  hEffSPDPt02->SetMarkerColor(9);
  hEffSPDPt02->SetMarkerStyle(20);
  TH1F *hEffSPDPt1 = new TH1F("hEffSPDPt1","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEffSPDPt1->SetLineWidth(2);
  hEffSPDPt1->SetLineColor(9);
  hEffSPDPt1->SetMarkerColor(9);
  hEffSPDPt1->SetMarkerStyle(20);
  TH1F *hEffSPDPt10 = new TH1F("hEffSPDPt10","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEffSPDPt10->SetLineWidth(2);
  hEffSPDPt10->SetLineColor(9);
  hEffSPDPt10->SetMarkerColor(9);
  hEffSPDPt10->SetMarkerStyle(20);

  TH1F *hEffoneSPDPt02 = new TH1F("hEffoneSPDPt02","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEffoneSPDPt02->SetLineWidth(2);
  hEffoneSPDPt02->SetLineColor(15);
  hEffoneSPDPt02->SetMarkerColor(15);
  hEffoneSPDPt02->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt1 = new TH1F("hEffoneSPDPt1","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEffoneSPDPt1->SetLineWidth(2);
  hEffoneSPDPt1->SetLineColor(15);
  hEffoneSPDPt1->SetMarkerColor(15);
  hEffoneSPDPt1->SetMarkerStyle(20);
  TH1F *hEffoneSPDPt10 = new TH1F("hEffoneSPDPt10","Efficiency; run number; TPC+ITS / TPC",nruns,firstrun-0.5,lastrun+0.5);
  hEffoneSPDPt10->SetLineWidth(2);
  hEffoneSPDPt10->SetLineColor(15);
  hEffoneSPDPt10->SetMarkerColor(15);
  hEffoneSPDPt10->SetMarkerStyle(20);

  // loop on runs
  for(Int_t irun=firstrun; irun<=lastrun; irun++) {
    if(!KeepRun(irun)) continue;

    TString path=pathBeforeRun;
    path+=irun;
    path.Append(pathAfterRun.Data());
    path.Prepend("alien://");

    printf("%s\n",path.Data());

    if(!PlotITSTrackingHists(path.Data(),ioValues,ioErrors)) {
      if(pathAfterRun2.Data()=="") continue;
      TString path2=pathBeforeRun;
      path2+=irun;
      path2.Append(pathAfterRun2.Data());
      path2.Prepend("alien://");
      if(!PlotITSTrackingHists(path2.Data(),ioValues,ioErrors)) continue;
    }

    Int_t bin=hEffge2Pt1->FindBin(irun);


    // fill histos
    index=0;
    hFracSPD1->SetBinContent(bin,ioValues[index]);
    hFracSPD1->SetBinError(bin,.01);
    index=1;
    hFracSPD2->SetBinContent(bin,ioValues[index]);
    hFracSPD2->SetBinError(bin,.01);

    index=2;
    hEffge2Pt02->SetBinContent(bin,ioValues[index]);
    hEffge2Pt02->SetBinError(bin,ioErrors[index]);
    index=3;
    hEffge2Pt1->SetBinContent(bin,ioValues[index]);
    hEffge2Pt1->SetBinError(bin,ioErrors[index]);
    index=4;
    hEffge2Pt10->SetBinContent(bin,ioValues[index]);
    hEffge2Pt10->SetBinError(bin,ioErrors[index]);

    index=5;
    hEff6Pt02->SetBinContent(bin,ioValues[index]);
    hEff6Pt02->SetBinError(bin,ioErrors[index]);
    index=6;
    hEff6Pt1->SetBinContent(bin,ioValues[index]);
    hEff6Pt1->SetBinError(bin,ioErrors[index]);
    index=7;
    hEff6Pt10->SetBinContent(bin,ioValues[index]);
    hEff6Pt10->SetBinError(bin,ioErrors[index]);

    index=8;
    hEffSPDPt02->SetBinContent(bin,ioValues[index]);
    hEffSPDPt02->SetBinError(bin,ioErrors[index]);
    index=9;
    hEffSPDPt1->SetBinContent(bin,ioValues[index]);
    hEffSPDPt1->SetBinError(bin,ioErrors[index]);
    index=10;
    hEffSPDPt10->SetBinContent(bin,ioValues[index]);
    hEffSPDPt10->SetBinError(bin,ioErrors[index]);

    index=11;
    hEffoneSPDPt02->SetBinContent(bin,ioValues[index]);
    hEffoneSPDPt02->SetBinError(bin,ioErrors[index]);
    index=12;
    hEffoneSPDPt1->SetBinContent(bin,ioValues[index]);
    hEffoneSPDPt1->SetBinError(bin,ioErrors[index]);
    index=13;
    hEffoneSPDPt10->SetBinContent(bin,ioValues[index]);
    hEffoneSPDPt10->SetBinError(bin,ioErrors[index]);

  }

  TCanvas *cSPD = new TCanvas("cSPD","cSPD",0,0,1000,300);
  cSPD->SetGridy();
  hFracSPD1->SetMaximum(1.2);
  hFracSPD1->SetMaximum(0);
  hFracSPD1->Draw("p");
  hFracSPD2->Draw("same,p");
  TLegend* lSPD=new TLegend(0.9,0.8,1,1);
  lSPD->AddEntry(hFracSPD1,"Frac. SPD1 ON","l");
  lSPD->AddEntry(hFracSPD2,"Frac. SPD2 ON","l");
  lSPD->Draw();

  TCanvas *cPt02 = new TCanvas("cPt02","cPt02",0,0,1000,300);
  cPt02->SetGridy();
  hEffge2Pt02->SetMaximum(1.2);
  hEffge2Pt02->SetMaximum(0);
  hEffge2Pt02->Draw();
  hEff6Pt02->Draw("same");
  hEffSPDPt02->Draw("same");
  hEffoneSPDPt02->Draw("same");
  TLegend* lPt02=new TLegend(0.9,0.8,1,1);
  lPt02->AddEntry(hEffge2Pt02,">=2 cls","l");
  lPt02->AddEntry(hEffoneSPDPt02,">=1 SPD + any","l");
  lPt02->AddEntry(hEffSPDPt02,"2 SPD + any","l");
  lPt02->AddEntry(hEff6Pt02,"6 cls","l");
  lPt02->Draw();

  TCanvas *cPt1 = new TCanvas("cPt1","cPt1",0,0,1000,300);
  cPt1->SetGridy();
  hEffge2Pt1->SetMaximum(1.2);
  hEffge2Pt1->SetMaximum(0);
  hEffge2Pt1->Draw();
  hEff6Pt1->Draw("same");
  hEffSPDPt1->Draw("same");
  hEffoneSPDPt1->Draw("same");
  lPt02->Draw();

  TCanvas *cPt10 = new TCanvas("cPt10","cPt10",0,0,1000,300);
  cPt10->SetGridy();
  hEffge2Pt10->SetMaximum(1.2);
  hEffge2Pt10->SetMaximum(0);
  hEffge2Pt10->Draw("p");
  hEff6Pt10->Draw("same,p");
  hEffSPDPt10->Draw("same,p");
  hEffoneSPDPt10->Draw("same,p");
  lPt02->Draw();

  return;
}
//------------------------------------------------------------------------
Bool_t KeepRun(Int_t irun) {



  // LHC10c good runs
  Int_t nruns10b=33;
  Int_t goodruns10b[33]={117222, 117220, 117116, 117112, 117109, 117099, 117092, 117086, 117077, 117065, 117063, 117060, 117059, 117054, 117053, 117052, 117050, 117048, 116645, 116643, 116574, 116571, 116562, 116403, 116402, 116288, 116102, 115414, 115401, 115393, 115193, 115186, 114931};

  // LHC10c good runs
  Int_t nruns10c=46;
  Int_t goodruns10c[46]={121040, 121039, 120829, 120825, 120824, 120823, 120822, 120821, 120820, 120758, 120750, 120741, 120671, 120617, 120616, 120505, 120504, 120503, 120244, 120079, 120076, 120073, 120072, 120069, 120067, 119862, 119859, 119856, 119853, 119849, 119846, 119845, 119844, 119842, 119841, 119163, 119161, 119159, 118561, 118560, 118558, 118556, 118518, 118512, 118507, 118506};
  //

  // LHC10d good runs
  Int_t nruns10d=80;
  Int_t goodruns10d[80]={126437, 126432, 126425, 126424, 126422, 126409, 126408, 126407, 126406, 126405, 126404, 126403, 126359, 126352, 126351, 126350, 126285, 126284, 126283, 126168, 126167, 126160, 126158, 126097, 126090, 126088, 126082, 126081, 126078, 126073, 126008, 126007, 126004, 125855, 125851, 125850, 125849, 125848, 125847, 125844, 125843, 125842, 125633, 125632, 125630, 125628, 125296, 125295, 125186, 125156, 125140, 125139, 125134, 125133, 125101, 125100, 125097, 125085, 125083, 125023, 124751, 124750, 124746, 124702, 124608, 124607, 124606, 124605, 124604, 124381, 124380, 124378, 124367, 124362, 124358, 124355, 124191, 124187, 122375, 122374};
  //

  // LHC10e good runs
  Int_t nruns10e=158;
  Int_t goodruns10e[158]={130369, 130365, 130360, 130358, 130356, 130354, 130353, 130348, 130343, 130342, 130179, 130178, 130172, 130170, 130168, 130158, 130157, 130156, 130151, 130149, 130148, 129983, 129966, 129962, 129961, 129960, 129959, 129763, 129760, 129750, 129748, 129747, 129745, 129744, 129742, 129738, 129736, 129735, 129734, 129731, 129729, 129726, 129725, 129723, 129667, 129666, 129659, 129655, 129654, 129653, 129652, 129651, 129650, 129649, 129648, 129647, 129641, 129639, 129599, 129598, 129597, 129587, 129586, 129541, 129540, 129536, 129529, 129528, 129527, 129526, 129525, 129524, 129523, 129522, 129521, 129520, 129519, 129516, 129515, 129514, 129513, 129512, 129510, 129508, 129042, 129041, 128913, 128912, 128911, 128910, 128855, 128853, 128850, 128849, 128843, 128836, 128835, 128833, 128824, 128823, 128820, 128819, 128817, 128813, 128777, 128776, 128678, 128677, 128621, 128615, 128611, 128609, 128605, 128596, 128594, 128592, 128590, 128582, 128581, 128507, 128506, 128505, 128504, 128503, 128498, 128495, 128494, 128486, 128483, 128452, 128366, 128260, 128257, 128192, 128191, 128189, 128186, 128185, 128182, 128180, 128175, 128053, 127942, 127941, 127937, 127936, 127935, 127933, 127932, 127931, 127822, 127817, 127815, 127814, 127730, 127729, 127724, 127719};



  // LHC10h good runs
  Int_t nruns10h=49;
  Int_t goodruns10h[49]={138624, 138637, 138666, 138730, 138731, 138732, 138737, 138740, 138742, 138830, 138836, 138871, 138872, 138924, 138965, 138972, 138977, 138978, 138980, 138982, 138983, 139029, 139036, 139037, 139038, 139042, 139105, 139107, 139309, 139310, 139314, 139360, 139438, 139439, 139440, 139441, 139465, 139466, 139467, 139470, 139471, 139503, 139505, 139507, 139510, 139511, 139513, 139514, 139517};

  /*
  Int_t nruns10h=91;
  Int_t goodruns10h[91]={139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137162, 137161};
  */

  // LHC11b good runs
  Int_t nruns11b=83;
  Int_t goodruns11b[83]={148531, 148541, 148544, 148648, 148800, 148847, 148857, 149102, 149115, 149122, 149123, 149367, 149382, 149389, 149395, 149398, 149432, 149434, 149442, 149449, 149452, 149453, 149456, 149461, 149462, 149472, 149484, 149487, 149499, 149529, 149530, 149531, 149533, 149534, 149540, 149543, 149544, 149545, 149547, 149549, 149582, 149584, 149585, 149588, 149590, 149592, 149594, 149596, 149598, 149600, 149604, 149657, 149669, 149680, 149682, 149712, 149723, 149726, 149727, 149734, 149735, 149740, 149746, 149749, 149757, 149761, 149769, 149771, 149772, 149773, 149774, 149775, 149776, 149781, 149785, 149786, 149789, 149790, 150375, 150421, 150586, 150589, 150705};

  // LHC11c good runs
  Int_t nruns11c=244;
  Int_t goodruns11c[244]={154570, 154495, 154485, 154483, 154482, 154480, 154478, 154448, 154385, 154383, 154382, 154296, 154293, 154289, 154286, 154283, 154281, 154273, 154270, 154269, 154266, 154264, 154261, 154257, 154252, 154222, 154221, 154220, 154219, 154211, 154207, 154163, 154158, 154151, 154145, 154143, 154141, 154138, 154136, 154132, 154131, 154130, 154129, 154126, 154125, 154091, 154083, 154081, 154039, 154031, 154030, 154026, 154024, 154018, 153954, 153946, 153944, 153939, 153938, 153935, 153929, 153924, 153916, 153911, 153909, 153906, 153876, 153875, 153873, 153812, 153808, 153807, 153805, 153798, 153796, 153794, 153784, 153781, 153779, 153778, 153777, 153776, 153738, 153733, 153728, 153727, 153726, 153725, 153718, 153709, 153702, 153594, 153591, 153589, 153588, 153587, 153583, 153578, 153571, 153570, 153566, 153560, 153558, 153552, 153548, 153544, 153542, 153541, 153539, 153536, 153533, 153513, 153465, 153415, 153413, 153373, 153371, 153369, 153363, 153362, 153360, 153296, 153234, 153232, 153227, 153223, 153127, 153121, 153120, 153117, 153116, 153115, 153059, 153056, 152935, 152934, 152907, 152823, 152822, 152821, 152820, 152819, 152817, 152780, 152773, 152772, 152751, 152750, 152718, 152717, 152716, 152715, 152708, 152702, 152701, 152698, 152697, 152696, 152695, 152658, 152581, 152513, 152512, 152488, 152455, 152377, 152371, 152369, 152368, 152367, 152334, 152332, 152323, 152322, 152321, 152320, 152319, 152318, 152314, 152313, 152312, 152311, 152310, 152309, 152307, 152306, 152304, 152285, 152284, 152257, 152256, 152214, 152209, 152208, 152207, 152206, 152146, 152138, 152137, 152136, 152091, 152090, 152089, 152087, 152086, 152083, 152082, 152081, 152079, 152078, 152046, 152015, 152011, 152008, 152007, 152005, 152003, 152002, 151852, 151851, 151850, 151849, 151810, 151809, 151752, 151751, 151732, 151724, 151689, 151681, 151680, 151678, 151674, 151672, 151671, 151669, 151666, 151665, 151664, 151661, 151660, 151655, 151638, 151636};

  // LHC11d good runs
  Int_t nruns11d=81;
  Int_t goodruns11d[81]={158508, 158509, 158511, 158518, 158520, 158521, 158526, 158528, 158531, 158592, 158602, 158604, 158611, 158615, 158622, 158626, 158706, 158717, 158722, 158745, 158777, 158780, 158784, 158788, 158791, 158793, 158794, 158848, 158868, 158876, 158877, 158879, 159028, 159042, 159076, 159090, 159120, 159121, 159128, 159147, 159154, 159162, 159168, 159173, 159185, 159191, 159193, 159199, 159201, 159205, 159207, 159214, 159216, 159218, 159221, 159223, 159255, 159258, 159260, 159285, 159286, 159318, 159379, 159451, 159503, 159517, 159521, 159535, 159538, 159571, 159577, 159580, 159581, 159582, 159586, 159593, 159595, 159599, 159602, 159606, 159635};


  Bool_t found=kFALSE;
  Int_t i=0;

  for(i=0; i<nruns10b; i++) {
    if(irun==goodruns10b[i]) {found=kTRUE; break;}
  }  
  if(found) return kTRUE;

  for(i=0; i<nruns10c; i++) {
    if(irun==goodruns10c[i]) {found=kTRUE; break;}
  }  
  if(found) return kTRUE;

  for(i=0; i<nruns10d; i++) {
    if(irun==goodruns10d[i]) {found=kTRUE; break;}
  }  
  if(found) return kTRUE;

  for(i=0; i<nruns10e; i++) {
    if(irun==goodruns10e[i]) {found=kTRUE; break;}
  }  
  if(found) return kTRUE;

  for(i=0; i<nruns10h; i++) {
    if(irun==goodruns10h[i]) {found=kTRUE; break;}
  }  
  if(found) return kTRUE;

  for(i=0; i<nruns11b; i++) {
    if(irun==goodruns11b[i]) {found=kTRUE; break;}
  }  
  if(found) return kTRUE;

  for(i=0; i<nruns11c; i++) {
    if(irun==goodruns11c[i]) {found=kTRUE; break;}
  }  
  if(found) return kTRUE;

  for(i=0; i<nruns11d; i++) {
    if(irun==goodruns11d[i]) {found=kTRUE; break;}
  }  
  if(found) return kTRUE;


  return kFALSE;
}
//-------------------------------------------------------------------------
void FakesWithChi2Cut(TString fname) {

  gStyle->SetOptStat(0);

  TFile *f= TFile::Open(fname.Data());
  if(!f) return kFALSE;

  TList *list=(TList*)f->Get("cOutputITS");
  TDirectoryFile *dir=0;
  if(!list) {
    dir=(TDirectoryFile*)f->GetDirectory("ITS_Performance");
    if(dir) list = (TList*)dir->Get("cOutputITS");
  }

  if(!list) return kFALSE;




  TH1F *fHistRedChi2AllPt02 = (TH1F*)list->FindObject("fHistITSRedChi2NonFakePt02");
  TH1F *fHistRedChi2FakesPt02 = (TH1F*)list->FindObject("fHistITSRedChi2FakePt02");
  TH1F *fHistRedChi2NonFakesPt02 = (TH1F*)fHistRedChi2AllPt02->Clone("fHistITSRedChi2NonFakePt02");
  fHistRedChi2NonFakesPt02->Add(fHistRedChi2FakesPt02,-1.);

  TH1F *fHistFakeFracVSChi2CutPt02 = (TH1F*)fHistRedChi2AllPt02->Clone("fHistFakeFracVSChi2CutPt02");
  fHistFakeFracVSChi2CutPt02->SetLineColor(2);
  fHistFakeFracVSChi2CutPt02->SetLineWidth(2);
  TH1F *fHistNonFakeEffVSChi2CutPt02 = (TH1F*)fHistRedChi2AllPt02->Clone("fHistFakeEffVSChi2CutPt02");
  fHistNonFakeEffVSChi2CutPt02->SetLineColor(4);
  fHistNonFakeEffVSChi2CutPt02->SetLineWidth(2);

  for(Int_t bin=1;bin<=fHistFakeFracVSChi2CutPt02->GetNbinsX();bin++) {
    Float_t fakesBelowCut=fHistRedChi2FakesPt02->Integral(1,bin);
    Float_t nonfakesBelowCut=fHistRedChi2NonFakesPt02->Integral(1,bin);
    Float_t allBelowCut=fHistRedChi2AllPt02->Integral(1,bin);
    Float_t nonfakesNoCut=fHistRedChi2NonFakesPt02->Integral();
    fHistFakeFracVSChi2CutPt02->SetBinContent(bin,fakesBelowCut/allBelowCut);
    fHistNonFakeEffVSChi2CutPt02->SetBinContent(bin,nonfakesBelowCut/nonfakesNoCut);
  }

  TH1F *fHistRedChi2AllPt05 = (TH1F*)list->FindObject("fHistITSRedChi2NonFakePt05");
  TH1F *fHistRedChi2FakesPt05 = (TH1F*)list->FindObject("fHistITSRedChi2FakePt05");
  TH1F *fHistRedChi2NonFakesPt05 = (TH1F*)fHistRedChi2AllPt05->Clone("fHistITSRedChi2NonFakePt05");
  fHistRedChi2NonFakesPt05->Add(fHistRedChi2FakesPt05,-1.);

  TH1F *fHistFakeFracVSChi2CutPt05 = (TH1F*)fHistRedChi2AllPt05->Clone("fHistFakeFracVSChi2CutPt05");
  fHistFakeFracVSChi2CutPt05->SetLineColor(2);
  fHistFakeFracVSChi2CutPt05->SetLineWidth(2);
  TH1F *fHistNonFakeEffVSChi2CutPt05 = (TH1F*)fHistRedChi2AllPt05->Clone("fHistFakeEffVSChi2CutPt05");
  fHistNonFakeEffVSChi2CutPt05->SetLineColor(4);
  fHistNonFakeEffVSChi2CutPt05->SetLineWidth(2);

  for(Int_t bin=1;bin<=fHistFakeFracVSChi2CutPt05->GetNbinsX();bin++) {
    Float_t fakesBelowCut=fHistRedChi2FakesPt05->Integral(1,bin);
    Float_t nonfakesBelowCut=fHistRedChi2NonFakesPt05->Integral(1,bin);
    Float_t allBelowCut=fHistRedChi2AllPt05->Integral(1,bin);
    Float_t nonfakesNoCut=fHistRedChi2NonFakesPt05->Integral();
    fHistFakeFracVSChi2CutPt05->SetBinContent(bin,fakesBelowCut/allBelowCut);
    fHistNonFakeEffVSChi2CutPt05->SetBinContent(bin,nonfakesBelowCut/nonfakesNoCut);
  }

  TH1F *fHistRedChi2AllPt1 = (TH1F*)list->FindObject("fHistITSRedChi2NonFakePt1");
  TH1F *fHistRedChi2FakesPt1 = (TH1F*)list->FindObject("fHistITSRedChi2FakePt1");
  TH1F *fHistRedChi2NonFakesPt1 = (TH1F*)fHistRedChi2AllPt1->Clone("fHistITSRedChi2NonFakePt1");
  fHistRedChi2NonFakesPt1->Add(fHistRedChi2FakesPt1,-1.);

  TH1F *fHistFakeFracVSChi2CutPt1 = (TH1F*)fHistRedChi2AllPt1->Clone("fHistFakeFracVSChi2CutPt1");
  fHistFakeFracVSChi2CutPt1->SetLineColor(2);
  fHistFakeFracVSChi2CutPt1->SetLineWidth(2);
  TH1F *fHistNonFakeEffVSChi2CutPt1 = (TH1F*)fHistRedChi2AllPt1->Clone("fHistFakeEffVSChi2CutPt1");
  fHistNonFakeEffVSChi2CutPt1->SetLineColor(4);
  fHistNonFakeEffVSChi2CutPt1->SetLineWidth(2);

  for(Int_t bin=1;bin<=fHistFakeFracVSChi2CutPt1->GetNbinsX();bin++) {
    Float_t fakesBelowCut=fHistRedChi2FakesPt1->Integral(1,bin);
    Float_t nonfakesBelowCut=fHistRedChi2NonFakesPt1->Integral(1,bin);
    Float_t allBelowCut=fHistRedChi2AllPt1->Integral(1,bin);
    Float_t nonfakesNoCut=fHistRedChi2NonFakesPt1->Integral();
    fHistFakeFracVSChi2CutPt1->SetBinContent(bin,fakesBelowCut/allBelowCut);
    fHistNonFakeEffVSChi2CutPt1->SetBinContent(bin,nonfakesBelowCut/nonfakesNoCut);
  }


  TCanvas *c=new TCanvas("c","c",0,0,1000,500);
  c->Divide(3,1);
  c->cd(1);
  fHistRedChi2NonFakesPt02->Scale(1./fHistRedChi2NonFakesPt02->Integral());
  fHistRedChi2NonFakesPt02->SetLineColor(4);
  fHistRedChi2NonFakesPt02->Draw();
  fHistRedChi2FakesPt02->Scale(1./fHistRedChi2FakesPt02->Integral());
  fHistRedChi2FakesPt02->SetLineColor(2);
  fHistRedChi2FakesPt02->Draw("same");
  c->cd(2);
  fHistRedChi2NonFakesPt05->Scale(1./fHistRedChi2NonFakesPt05->Integral());
  fHistRedChi2NonFakesPt05->SetLineColor(4);
  fHistRedChi2NonFakesPt05->Draw();
  fHistRedChi2FakesPt05->Scale(1./fHistRedChi2FakesPt05->Integral());
  fHistRedChi2FakesPt05->SetLineColor(2);
  fHistRedChi2FakesPt05->Draw("same");
  c->cd(3);
  fHistRedChi2NonFakesPt1->Scale(1./fHistRedChi2NonFakesPt1->Integral());
  fHistRedChi2NonFakesPt1->SetLineColor(4);
  fHistRedChi2NonFakesPt1->Draw();
  fHistRedChi2FakesPt1->Scale(1./fHistRedChi2FakesPt1->Integral());
  fHistRedChi2FakesPt1->SetLineColor(2);
  fHistRedChi2FakesPt1->Draw("same");

  TCanvas *cc=new TCanvas("cc","cc",0,0,1000,500);
  cc->Divide(3,1);
  cc_1->SetLogx();
  cc_2->SetLogx();
  cc_3->SetLogx();
  cc->cd(1);
  fHistFakeFracVSChi2CutPt02->SetXTitle("maximum ITS #chi^{2}/nclusters");
  fHistFakeFracVSChi2CutPt02->SetMinimum(0);
  fHistFakeFracVSChi2CutPt02->SetMaximum(1);
  fHistFakeFracVSChi2CutPt02->Draw();
  fHistNonFakeEffVSChi2CutPt02->Draw("same");
  cc->cd(2);
  fHistFakeFracVSChi2CutPt05->SetXTitle("maximum ITS #chi^{2}/nclusters");
  fHistFakeFracVSChi2CutPt05->SetMinimum(0);
  fHistFakeFracVSChi2CutPt05->SetMaximum(1);
  fHistFakeFracVSChi2CutPt05->Draw();
  fHistNonFakeEffVSChi2CutPt05->Draw("same");
  cc->cd(3);
  fHistFakeFracVSChi2CutPt1->SetXTitle("maximum ITS #chi^{2}/nclusters");
  fHistFakeFracVSChi2CutPt1->SetMinimum(0);
  fHistFakeFracVSChi2CutPt1->SetMaximum(1);
  fHistFakeFracVSChi2CutPt1->Draw();
  fHistNonFakeEffVSChi2CutPt1->Draw("same");



  return;
}
