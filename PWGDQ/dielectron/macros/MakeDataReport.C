void SetupStyle();

void MakeDataReport(const char* outputFile="JpsiDataReport.pdf",
                          const char* histos="jpsi_HistosSE.root",
                          const char* cf="jpsi_CF.root")
{
  //
  // Make a pdf file with the efficiency report
  //

  SetupStyle();
  
  AliDielectronCFdraw d(cf);
  d.SetRangeUser("PairType",1,1);
  d.SetRangeUser("Y",-.89,.9,"0");
  
  
  TFile f("jpsi_HistosSE.root");
  
  AliDielectronHistos h;
  TIter nextHists((TList*)f.Get("Dielectron_Histos"));
  
  TPaveText pt(.02,.6,.98,.8);
  TText *t1=pt.AddText("");
  TText *t2=pt.AddText("");
  
  TCanvas *c1=new TCanvas;
  
  TPDF p(outputFile);

  //
  // Invariant mass plots
  //
  
  
  //
  // Make QA info
  //
  
  t1->SetTitle("QA summary plots for");
  THashList *list=0x0;
  while ( (list=(THashList*)nextHists()) ){
    h.SetHistogramList(*list);
    t2->SetTitle(list->GetName());
    pt.Draw();
    c1->Update();
    h.Draw();
    c1->Clear();
  }
  p.Close();
  delete c1;
}

void SetupStyle()
{
  const Int_t NCont=255;
  
  TStyle *st = new TStyle("mystyle","mystyle");
  gROOT->GetStyle("Plain")->Copy((*st));
  st->SetTitleX(0.1);
  st->SetTitleW(0.8);
  st->SetTitleH(0.08);
  st->SetStatX(.9);
  st->SetStatY(.9);
  st->SetNumberContours(NCont);
  st->SetPalette(1,0);
  st->SetOptStat("erm");
  st->SetOptFit(0);
  st->SetGridColor(kGray+1);
  st->SetPadGridX(kTRUE);
  st->SetPadGridY(kTRUE);
  st->SetPadTickX(kTRUE);
  st->SetPadTickY(kTRUE);
  st->cd();
  
  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  
}

void DrawUnbinned(){
  TFile f("jpsi_debug.root");
//   if (!f.IsOpen()) return;

  TTree *t=(TTree*)f.Get("Pair");
//   if (!t) return;

  TCanvas c1;
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  
  TLegend *leg=new TLegend(0.59,.79,.99,.99);
  TLine l;
  
  l.SetLineColor(kGreen-5);
  l.SetLineWidth(2);
  l.SetLineStyle(2);
  leg->SetFillColor(10);

  leg->Clear();
  

  t->SetLineColor(kBlack);
  t->Draw("M>>hAll(200,-.01,3.99)","","histe");
  TH1 *hAll=(TH1*)gROOT->FindObject("hAll");
  hAll->SetMinimum(0.1);
  hAll->SetTitle(";M [GeV]; yield");
  leg->AddEntry(hAll,"|n#sigma e|<2 + pt>0.3 GeV","l");

  l.DrawLine(3.097,1,3.097,1e4);
  
  t->SetLineColor(kOrange-5);
  t->Draw("M>>hC11(200,-.01,3.99)","abs(Leg1_ImpactParXY)<.004&&abs(Leg2_ImpactParXY)<.004","histesame");
  hAll=(TH1*)gROOT->FindObject("hC11");
  leg->AddEntry(hAll,"|n#sigma e|<2 + pt>0.3 GeV + |dXY|<40#mum","l");
  
  TCut d1_1="abs(Leg1_TPC_nSigma_Electrons)<1";
  TCut d2_1="abs(Leg2_TPC_nSigma_Electrons)<1";
  TCut d_1=d1_1+d2_1;
  
  t->SetLineColor(kRed);
  t->Draw("M>>hC1(200,-.01,3.99)","Leg2_Pt>1&&Leg1_Pt>1","histesame");
  hAll=(TH1*)gROOT->FindObject("hC1");
  leg->AddEntry(hAll,"|n#sigma e|<2 + pt>1 GeV","l");

  
  t->SetLineColor(kGreen);
  t->Draw("M>>hC2(200,-.01,3.99)",d_1+"Leg2_Pt>1&&Leg1_Pt>1","histesame");
  hAll=(TH1*)gROOT->FindObject("hC2");
  leg->AddEntry(hAll,"|n#sigma e|<1 + pt>1 GeV","l");

  t->SetLineColor(kMagenta);
  t->Draw("M>>hC3(200,-.01,3.99)","Leg1_Pt>2&&Leg2_Pt>2","histesame");
  hAll=(TH1*)gROOT->FindObject("hC3");
  leg->AddEntry(hAll,"|n#sigma e|<2 + pt>2 GeV","l");
  
  t->SetLineColor(kBlue);
  t->Draw("M>>hC4(200,-.01,3.99)","Leg1_Pt>3&&Leg2_Pt>3","histesame");
  hAll=(TH1*)gROOT->FindObject("hC4");
  leg->AddEntry(hAll,"|n#sigma e|<2 + pt>3 GeV","l");
  
  leg->Draw();
}

/*
  Double_t alephParameters[5];
  // simulation
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;
  Color_t color=kRed;

  TF1 *foProton = new TF1("foProton", "50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foPion = new TF1("foPion", "50*AliExternalTrackParam::BetheBlochAleph(x/0.13957,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foElec = new TF1("foElec", "50*AliExternalTrackParam::BetheBlochAleph(x/0.000511,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foKaon = new TF1("foKaon", "50*AliExternalTrackParam::BetheBlochAleph(x/0.493677,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foMuon = new TF1("foMuon", "50*AliExternalTrackParam::BetheBlochAleph(x/0.105658,[0],[1],[2],[3],[4])",0.05,20);
  //
  foProton->SetParameters(alephParameters);
  foPion->SetParameters(alephParameters);
  foElec->SetParameters(alephParameters);
  foKaon->SetParameters(alephParameters);
  foMuon->SetParameters(alephParameters);
  //
  foProton->SetLineColor(color);
  foPion->SetLineColor(color);
  foElec->SetLineColor(color);
  foKaon->SetLineColor(color);
  foMuon->SetLineColor(color);
  //
  Int_t lineWidth=1;
  foProton->SetLineWidth(lineWidth);
  foPion->SetLineWidth(lineWidth);
  foElec->SetLineWidth(lineWidth);
  foKaon->SetLineWidth(lineWidth);
  foMuon->SetLineWidth(lineWidth);

  //
  foProton->SetNpx(200);
  foPion->SetNpx(200);
  foElec->SetNpx(200);
  foKaon->SetNpx(200);
  foMuon->SetNpx(200);
  //
  foProton->Draw("same");
  foPion->Draw("same");
  foElec->Draw("same");
  foKaon->Draw("same");
  foMuon->Draw("same");





  // data
  Double_t res=5.2e-2;
  alephParameters[0] = 0.0283086;
  alephParameters[1] = 2.63394e+01;
  alephParameters[2] = 5.04114e-11;
  alephParameters[3] = 2.12543e+00;
  alephParameters[4] = 4.88663e+00;
  Color_t color=kRed;

    alephParameters[0] = 0.0283086/0.97;
    //alephParameters[0] = 0.0283086;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;


  TF1 *foDataProton = new TF1("foDataProton", "50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foDataProtonP = new TF1("foDataProtonP",Form( "50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1+%f)",3*res),0.05,20);
  TF1 *foDataProtonM = new TF1("foDataProtonM", Form("50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1-%f)",res),0.05,20);

  TF1 *foDataPion = new TF1("foDataPion", "50*AliExternalTrackParam::BetheBlochAleph(x/0.13957,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foDataPionP = new TF1("foDataPionP",Form( "50*AliExternalTrackParam::BetheBlochAleph(x/0.13957,[0],[1],[2],[3],[4])*(1+%f)",res),0.05,20);
  TF1 *foDataPionM = new TF1("foDataPionM", Form("50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1-%f)",res),0.05,20);

  TF1 *foDataElec = new TF1("foDataElec", "50*AliExternalTrackParam::BetheBlochAleph(x/0.000511,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foDataElecP = new TF1("foDataElecP",Form( "50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1+%f)",res),0.05,20);
  TF1 *foDataElecM = new TF1("foDataElecM", Form("50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1-%f)",res),0.05,20);

  TF1 *foDataKaon = new TF1("foDataKaon", "50*AliExternalTrackParam::BetheBlochAleph(x/0.493677,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foDataKaonP = new TF1("foDataKaonP",Form( "50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1+%f)",res),0.05,20);
  TF1 *foDataKaonM = new TF1("foDataKaonM", Form("50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1-%f)",res),0.05,20);

  TF1 *foDataMuon = new TF1("foDataMuon", "50*AliExternalTrackParam::BetheBlochAleph(x/0.105658,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foDataMuonP = new TF1("foDataMuonP",Form( "50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1+%f)",res),0.05,20);
  TF1 *foDataMuonM = new TF1("foDataMuonM", Form("50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1-%f)",res),0.05,20);

  //
  foDataProton->SetParameters(alephParameters);
  foDataProtonP->SetParameters(alephParameters);
  foDataProtonM->SetParameters(alephParameters);
  foDataPion->SetParameters(alephParameters);
  foDataPionP->SetParameters(alephParameters);
  foDataPionM->SetParameters(alephParameters);
  foDataElec->SetParameters(alephParameters);
  foDataElecP->SetParameters(alephParameters);
  foDataElecM->SetParameters(alephParameters);
  foDataKaon->SetParameters(alephParameters);
  foDataKaonP->SetParameters(alephParameters);
  foDataKaonM->SetParameters(alephParameters);
  foDataMuon->SetParameters(alephParameters);
  foDataMuonP->SetParameters(alephParameters);
  foDataMuonM->SetParameters(alephParameters);
  //
  foDataProton->SetLineColor(color);
  foDataProtonP->SetLineColor(color-4);
  foDataProtonM->SetLineColor(color-4);
  foDataPion->SetLineColor(color);
  foDataPionP->SetLineColor(color-4);
  foDataPionM->SetLineColor(color-4);
  foDataElec->SetLineColor(color);
  foDataElecP->SetLineColor(color-4);
  foDataElecM->SetLineColor(color-4);
  foDataKaon->SetLineColor(color);
  foDataKaonP->SetLineColor(color-4);
  foDataKaonM->SetLineColor(color-4);
  foDataMuon->SetLineColor(color);
  foDataMuonP->SetLineColor(color-4);
  foDataMuonM->SetLineColor(color-4);
  //
  Int_t lineWidth=1;
  foDataProton->SetLineWidth(lineWidth);
  foDataProtonP->SetLineWidth(lineWidth);
  foDataProtonM->SetLineWidth(lineWidth);
  foDataPion->SetLineWidth(lineWidth);
  foDataPionP->SetLineWidth(lineWidth);
  foDataPionM->SetLineWidth(lineWidth);
  foDataElec->SetLineWidth(lineWidth);
  foDataElecP->SetLineWidth(lineWidth);
  foDataElecM->SetLineWidth(lineWidth);
  foDataKaon->SetLineWidth(lineWidth);
  foDataKaonP->SetLineWidth(lineWidth);
  foDataKaonM->SetLineWidth(lineWidth);
  foDataMuon->SetLineWidth(lineWidth);
  foDataMuonP->SetLineWidth(lineWidth);
  foDataMuonM->SetLineWidth(lineWidth);

  //
  foDataProtonP->SetLineStyle(2);
  foDataProtonM->SetLineStyle(2);
  foDataPionP->SetLineStyle(2);
  foDataPionM->SetLineStyle(2);
  foDataElecP->SetLineStyle(2);
  foDataElecM->SetLineStyle(2);
  foDataKaonP->SetLineStyle(2);
  foDataKaonM->SetLineStyle(2);
  foDataMuonP->SetLineStyle(2);
  foDataMuonM->SetLineStyle(2);

  //
  foDataProton->SetNpx(200);
  foDataProtonP->SetNpx(200);
  foDataProtonM->SetNpx(200);
  foDataPion->SetNpx(200);
  foDataPionP->SetNpx(200);
  foDataPionM->SetNpx(200);
  foDataElec->SetNpx(200);
  foDataKaon->SetNpx(200);
  foDataMuon->SetNpx(200);
  //
  foDataProton->Draw("same");
  foDataProtonP->Draw("same");
  foDataProtonM->Draw("same");
  foDataPion->Draw("same");
  foDataElec->Draw("same");
  foDataKaon->Draw("same");
  foDataMuon->Draw("same");





{

  Int_t baseColors[5]={kRed, kGreen+1, kAzure-4, kMagenta, kCyan+1};
  Int_t sigmaColorOffset=1;
  
Int_t baseColors[5]={kRed, kGreen+1, kAzure-4, kMagenta, kCyan+1};
  Int_t sigmaColorOffset=0;
Int_t baseColors[5]={kRed, kRed, kRed, kRed, kRed};

  Double_t sigmas[5]={3,3,3,3,3};
  Double_t masses[5];

  for (Int_t i=0; i<AliPID::kSPECIES; ++i) masses[i]=AliPID::ParticleMass(i);
  
  Double_t res=7.e-2;
  Double_t alephParameters[5];

  alephParameters[0] = 0.0283086;
  alephParameters[1] = 2.63394e+01;
  alephParameters[2] = 5.04114e-11;
  alephParameters[3] = 2.12543e+00;
  alephParameters[4] = 4.88663e+00;
      alephParameters[0] = 1.25202/50.;   //was 1.79571/55.;
      alephParameters[1] = 2.74992e+01;   //was 22.0028;
      alephParameters[2] = TMath::Exp(-3.31517e+01);  //was1.55354e-11;
      alephParameters[3] = 2.46246;       //was 2.39804;
      alephParameters[4] = 6.78938;       //was 5.1209;

  Double_t mip=50;

  Color_t color=kRed;
  Int_t lineWidth=2;

TF1 *fBethe[5];
TF1 *fBetheP[5];
TF1 *fBetheM[5];

for (Int_t i=0; i<5; ++i){
  fBethe[i] = new TF1(Form("fBethe%d",i), Form("%f*AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])",mip,masses[i]),0.05,20);
  fBetheP[i] = new TF1(Form("fBethe%d",i), Form("%f*AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])*(1+%f*%f)",mip,masses[i],res,sigmas[i]),0.05,20);
  fBetheM[i] = new TF1(Form("fBethe%d",i), Form("%f*AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])*(1-%f*%f)",mip,masses[i],res,sigmas[i]),0.05,20);

  fBethe[i]->SetParameters(alephParameters);
  fBetheP[i]->SetParameters(alephParameters);
  fBetheM[i]->SetParameters(alephParameters);

  fBethe[i]->SetLineColor(baseColors[i]);
  fBetheP[i]->SetLineColor(baseColors[i]-sigmaColorOffset);
  fBetheM[i]->SetLineColor(baseColors[i]-sigmaColorOffset);

  fBethe[i]->SetLineWidth(lineWidth);
  fBetheP[i]->SetLineWidth(lineWidth);
  fBetheM[i]->SetLineWidth(lineWidth);
  
  fBetheP[i]->SetLineStyle(2);
  fBetheM[i]->SetLineStyle(2);

  fBethe[i]->SetNpx(200);
  fBetheP[i]->SetNpx(200);
  fBetheM[i]->SetNpx(200);
}

for (Int_t i=0; i<5; ++i){
  fBethe[i]->Draw("same");
//   fBetheP[i]->Draw("same");
//   fBetheM[i]->Draw("same");
}
}

*/


/*

    Double_t resolution=0.052;
    Double_t nSigma=3.;
    TF1 *ffPio=new TF1(Form("fBethe%d",AliPID::kPion), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigma,resolution, AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);
    ffPio->SetParameters(0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00);
fPPio->Draw("same");

  TF1 f("fP","-8*exp(-0.6*x)",0,40);
f.Draw("same")

//Unbinned fit

RooWorkspace w("w",kTRUE);
w.factory("CBShape::cb(M[0,5],x0[0],sigma[.04],alpha[1],n[1])")
RooRealVar cut("cut","cut")
RooDataSet data("data","data",c,w::M,"cut")



abs(Leg1_TPC_nSigma_Electrons)<3&&abs(Leg2_TPC_nSigma_Electrons)<3&&(Leg1_TPC_nSigma_Pions)>=3.5&&(Leg2_TPC_nSigma_Pions)>=3.5&&(Leg1_TPC_nSigma_Protons)>=3.0&&(Leg2_TPC_nSigma_Protons)>=3.0&&abs(Leg1_Eta)<0.9&&abs(Leg2_Eta)<0.9&&Leg1_NclsTPC>=70&&Leg2_NclsTPC>=70&&Leg1_Pt>=1&&Leg2_Pt>=1&&abs(Y)<0.9&&Leg1_ITSLayerFirstCls<2&&Leg2_ITSLayerFirstCls<2
*/

/*
c=Pair
c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<3&&abs(Leg2_TPC_nSigma_Electrons)<3");
c->SetAlias("cutPi","(Leg1_TPC_nSigma_Pions)>=3.5&&(Leg2_TPC_nSigma_Pions)>=3.5");
c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>=1.5&&(Leg2_TPC_nSigma_Protons)>=1.5");
// c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>=3.&&(Leg2_TPC_nSigma_Protons)>=3.");
c->SetAlias("pidSig","cutE&&cutPi&&cutP");

c->SetAlias("LegEta","abs(Leg1_Eta)<0.9&&abs(Leg2_Eta)<0.9");
c->SetAlias("LegNcl","Leg1_NclsTPC>=70&&Leg2_NclsTPC>=70");
c->SetAlias("LegPt","Leg1_Pt>=.6&&Leg2_Pt>=.6");
c->SetAlias("Rap","abs(Y)<0.9");
c->SetAlias("QA","LegNcl&&LegEta&&Rap&&LegPt&&Rap");

c->SetAlias("ITS1","Leg1_ITSLayerFirstCls<1&&Leg2_ITSLayerFirstCls<1")
c->SetAlias("ITS2","Leg1_ITSLayerFirstCls<2&&Leg2_ITSLayerFirstCls<2")
c->SetAlias("ITS3","Leg1_ITSLayerFirstCls<3&&Leg2_ITSLayerFirstCls<3")
c->SetAlias("ITS4","Leg1_ITSLayerFirstCls<4&&Leg2_ITSLayerFirstCls<4")
c->SetAlias("ITS5","Leg1_ITSLayerFirstCls<5&&Leg2_ITSLayerFirstCls<5")
c->SetAlias("ITS6","Leg1_ITSLayerFirstCls<6&&Leg2_ITSLayerFirstCls<6")

c->SetAlias("cut","PairType==1 && QA && pidSig && ITS5");

TFile f("/data/Work/data/pp2.76/m.root","recreate")
TTree *tnew=c->CopyTree("cut")
tnew->Write("Pair");
f.Save()
f.Close()

c->SetMarkerStyle(20);
c->SetMarkerSize(1);
c->Draw("M>>hM(250,0.,5.)","cut","e");
hM=hM
hM->SetMarkerColor(kRed)
hM->SetLineColor(kRed)

c->SetAlias("cutLS","(PairType==0||PairType==2) && QA && pidSig && ITS5");

c->Draw("M>>hMLS(250,0.,5.)","cutLS","esame");
hMLS=hMLS
hMLS->SetLineColor(kBlue);
hMLS->SetMarkerStyle(27)
Double_t factor=hM->Integral(hM->FindBin(3.201),hM->FindBin(4.9))/hMLS->Integral(hMLS->FindBin(3.201),hMLS->FindBin(4.9));
printf("scale factor: %.2f\n",factor)
hMLS->Scale(factor);

hSub=(TH1*)hM->Clone("hSub");
hSub->Add(hMLS,-1);
hSub->SetLineColor(kBlack);
hSub->SetMarkerColor(kBlack);
hSub->Draw("same")

hSub->Integral(hSub->FindBin(2.961),hSub->FindBin(3.159))

*/



/*
//
// Default
//
c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<3&&abs(Leg2_TPC_nSigma_Electrons)<3");
c->SetAlias("cutPi","(Leg1_TPC_nSigma_Pions)>=3.5&&(Leg2_TPC_nSigma_Pions)>=3.5");
c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>=3&&(Leg2_TPC_nSigma_Protons)>=3");
c->SetAlias("pidSig","cutE&&cutPi&&cutP");

c->SetAlias("LegEta","abs(Leg1_Eta)<0.9&&abs(Leg2_Eta)<0.9");
c->SetAlias("LegNcl","Leg1_NclsTPC>=70&&Leg2_NclsTPC>=70");
c->SetAlias("LegPt","Leg1_Pt>=.6&&Leg2_Pt>=0.6");
c->SetAlias("Rap","abs(Y)<0.9");
c->SetAlias("QA","LegNcl&&LegEta&&Rap&&LegPt&&Rap");

c->SetAlias("cutSigN","Leg2_TPCsignalN>=Leg2_NclsTPC&&Leg1_TPCsignalN>=Leg1_NclsTPC");
c->SetAlias("etaAdd","Leg1_Eta<0&&Leg2_Eta<0")
c->SetAlias("ITS1","Leg1_ITSLayerFirstCls<1&&Leg2_ITSLayerFirstCls<1")
c->SetAlias("ITS2","Leg1_ITSLayerFirstCls<2&&Leg2_ITSLayerFirstCls<2")
c->SetAlias("ITS3","Leg1_ITSLayerFirstCls<3&&Leg2_ITSLayerFirstCls<3")
c->SetAlias("ITS4","Leg1_ITSLayerFirstCls<4&&Leg2_ITSLayerFirstCls<4")
c->SetAlias("ITS5","Leg1_ITSLayerFirstCls<5&&Leg2_ITSLayerFirstCls<5")
c->SetAlias("ITS6","Leg1_ITSLayerFirstCls<6&&Leg2_ITSLayerFirstCls<6")

c->SetAlias("cut","PairType==1 && QA && pidSig");

c->SetMarkerStyle(20);
c->SetMarkerSize(.5);
c->Draw("M>>hM(125,0.,5.)","cut","e");
c->Draw("M>>hM1(125,0.,5.)","cut&&ITS1","esame"); hM1->SetLineColor(kRed); hM1->SetMarkerColor(kRed);
c->Draw("M>>hM2(125,0.,5.)","cut&&ITS2","esame"); hM2->SetLineColor(kBlue); hM2->SetMarkerColor(kBlue);
c->Draw("M>>hM3(125,0.,5.)","cut&&ITS3","esame"); hM3->SetLineColor(kGreen);  hM3->SetMarkerColor(kGreen);
c->Draw("M>>hM4(125,0.,5.)","cut&&ITS4","esame"); hM4->SetLineColor(kMagenta);  hM4->SetMarkerColor(kMagenta);
c->Draw("M>>hM5(125,0.,5.)","cut&&ITS5","esame"); hM5->SetLineColor(kYellow);  hM5->SetMarkerColor(kYellow);
c->Draw("M>>hM6(125,0.,5.)","cut&&ITS6","esame"); hM6->SetLineColor(kAzure);  hM6->SetMarkerColor(kAzure);

c->SetAlias("cut","PairType==1&&QA&&pidSig&&Leg1_ITSLayerFirstCls<1&&Leg2_ITSLayerFirstCls<1");
c->Draw("M>>hM1(125,0.,5.)","cut","esame");
hM1->SetMarkerColor(kBlue); hM1->SetLineColor(kBlue)

c->SetAlias("cut","PairType==1&&QA&&pidSig&&Leg1_ITSLayerFirstCls<2&&Leg2_ITSLayerFirstCls<2&&Run!=146805");
c->Draw("M>>hM2(125,0.,5.)","cut","esame");
hM2->SetMarkerColor(kMagenta); hM2->SetLineColor(kMagenta); hM2->SetMarkerStyle(20);

c->SetAlias("cut","PairType==1&&QA&&pidSig&&Leg1_ITSLayerFirstCls<5&&Leg2_ITSLayerFirstCls<5");
c->Draw("M>>hM3(125,0.,5.)","cut","esame");
hM3->SetMarkerColor(kGreen); hM3->SetLineColor(kGreen)

//
c->SetMarkerColor(kBlack); c->SetLineColor(kBlack); c->Draw("M>>hM(125,0.,5.)","cut","e");
new TCanvas; c->SetMarkerColor(kBlue); c->SetLineColor(kBlue); c->Draw("M>>hM2(125,0.01,5.01)","cut","e");
new TCanvas; c->SetMarkerColor(kRed); c->SetLineColor(kRed); c->Draw("M>>hM3(125,0.02,5.02)","cut","e");
new TCanvas; c->SetMarkerColor(kGreen); c->SetLineColor(kGreen); c->Draw("M>>hM4(125,0.03,5.03)","cut","e");



// trd cut

c->SetAlias("stp","1./9.*TMath::Pi()")
c->SetAlias("trdcut","((Leg1_Phi>1.5*stp&&Leg1_Phi<5.5*stp || Leg1_Phi>9.5*stp&&Leg1_Phi<16.5*stp) && (Leg2_Phi>1.5*stp&&Leg2_Phi<5.5*stp || Leg2_Phi>9.5*stp&&Leg2_Phi<16
.5*stp))")

c->SetAlias("trdcutn","(!(Leg1_Phi>1.5*stp&&Leg1_Phi<5.5*stp || Leg1_Phi>9.5*stp&&Leg1_Phi<16.5*stp) && !(Leg2_Phi>1.5*stp&&Leg2_Phi<5.5*stp || Leg2_Phi>9.5*stp&&Leg2_Phi
<16.5*stp))")

c->SetAlias("trdcutn2","((Leg1_Phi>1.5*stp&&Leg1_Phi<5.5*stp || Leg1_Phi>9.5*stp&&Leg1_Phi<16.5*stp) && !(Leg2_Phi>1.5*stp&&Leg2_Phi<5.5*stp || Leg2_Phi>9.5*stp&&Leg2_Phi<16.5*stp))")

c->SetAlias("trdcutn3","(!(Leg1_Phi>1.5*stp&&Leg1_Phi<5.5*stp || Leg1_Phi>9.5*stp&&Leg1_Phi<16.5*stp) && (Leg2_Phi>1.5*stp&&Leg2_Phi<5.5*stp || Leg2_Phi>9.5*stp&&Leg2_Phi<16.5*stp))")

c1->Divide(2,2)
c1->cd(1)
c->Draw("M>>hM1(125,0.,5.)","cut&&trdcut","e");
c1->cd(2)
c->Draw("M>>hM2(125,0.,5.)","cut&&trdcutn","e");
c1->cd(3)
c->Draw("M>>hM3(125,0.,5.)","cut&&trdcutn2","e");
c1->cd(4)
c->Draw("M>>hM4(125,0.,5.)","cut&&trdcutn3","e");


c->SetAlias("notrd","((Leg1_TrackStatus&512)==0 && (Leg2_TrackStatus&512)==0)")
c->SetAlias("withtrd","((Leg1_TrackStatus&512)==512 && (Leg2_TrackStatus&512)==512)")

c->SetAlias("notrd2","((Leg1_TrackStatus&512)==0 && (Leg2_TrackStatus&512)==512)")
c->SetAlias("notrd3","((Leg1_TrackStatus&512)==512 && (Leg2_TrackStatus&512)==0)")

c3->Divide(2,2)
c3->cd(1)
c->Draw("M>>hM5(125,0.,5.)","cut&&notrd","e");
c3->cd(2)
c->Draw("M>>hM6(125,0.,5.)","cut&&withtrd","e");
c3->cd(3)
c->Draw("M>>hM7(125,0.,5.)","cut&&notrd2","e");
c3->cd(4)
c->Draw("M>>hM8(125,0.,5.)","cut&&notrd3","e");

//pt ele only pid
c->SetAlias("cutE","Leg1_TPC_nSigma_Electrons>-2.&&Leg1_TPC_nSigma_Electrons<3&&Leg2_TPC_nSigma_Electrons>-2.&&Leg2_TPC_nSigma_Electrons<3");
c->SetAlias("pidSig","cutE");
c->SetAlias("LegPt","Leg1_Pt>1.2&&Leg2_Pt>1.2");

c->SetAlias("LegEta","abs(Leg1_Eta)<0.9&&abs(Leg2_Eta)<0.9");
c->SetAlias("LegNcl","Leg1_NclsTPC>70&&Leg2_NclsTPC>70");
c->Draw("M>>hM(125,0.,5.)","cut","e");

//---------------------

c->SetAlias("cut","PairType==1&&QA")


c->SetAlias("nCls","Leg1_NclsTPC>90&&Leg2_NclsTPC>90");


c->SetAlias("cutE","Leg1_TPC_nSigma_Electrons>-1 && Leg2_TPC_nSigma_Electrons>-1");
c->SetAlias("LegPt","Leg1_Pt>1.5&&Leg2_Pt>1.5");

//distance to chamber edge
c->SetAlias("tgs","pi/180.*20");
c->SetAlias("ed1","Leg1_Phi-int(Leg1_Phi/tgs+.5)*tgs");
c->SetAlias("ed2","Leg2_Phi-int(Leg2_Phi/tgs+.5)*tgs");

c->Draw("M:abs(ed1)>>hMd(15,0.,.18,65,0.,5.2)","cut","colz");


//--------PID
//-Param dEdx
c->SetAlias("cutPipardEdx","Leg1_TPC_signal>75-20*exp(-.7*Leg1_P_InnerParam)&&Leg2_TPC_signal>75-20*exp(-.7*Leg2_P_InnerParam)")

//require TOF PID
  c->SetAlias("TOFe1r","(((Leg1_TrackStatus&32768)==32768)&&abs(Leg1_TOF_nSigma_Electrons)<3)");
  c->SetAlias("TOFe2r","(((Leg2_TrackStatus&32768)==32768)&&abs(Leg2_TOF_nSigma_Electrons)<3)");

//if avail TOF PID
  c->SetAlias("TOFe1","((Leg1_TrackStatus&32768)==0)||(((Leg1_TrackStatus&32768)==32768)&&abs(Leg1_TOF_nSigma_Electrons)<3)");
  c->SetAlias("TOFe2","((Leg2_TrackStatus&32768)==0)||(((Leg2_TrackStatus&32768)==32768)&&abs(Leg2_TOF_nSigma_Electrons)<3)");
  c->SetAlias("TOFe","TOFe1&&TOFe2");


c->SetAlias("cutPspecial","(abs(Leg1_TPC_nSigma_Protons)>3||(abs(Leg1_TPC_nSigma_Protons)<=3&&TOFe1r))&&(abs(Leg2_TPC_nSigma_Protons)>3||(abs(Leg2_TPC_nSigma_Protons)<=3&&TOFe2r))")

c->SetAlias("cutP","(Leg1_Pt>1.2||TOFe1r)) && (Leg2_Pt>1.2 || TOFe2r)")

//-- nsigma
c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<3&&abs(Leg2_TPC_nSigma_Electrons)<3");
// c->SetAlias("cutE","Leg1_TPC_nSigma_Electrons>-1 && Leg2_TPC_nSigma_Electrons>-1");
c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>3&&abs(Leg2_TPC_nSigma_Pions)>3");
c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>3.&&(Leg2_TPC_nSigma_Protons)>3.3");
c->SetAlias("pidSig","cutE&&cutPi&&cutP");
//-- Pi param
// c->SetAlias("eleParam","Leg1_TPC_nSigma_Electrons<5&&Leg2_TPC_nSigma_Electrons<5&&Leg1_TPC_nSigma_Electrons>-2.65*exp(-0.9*Leg1_P_InnerParam)&&Leg2_TPC_nSigma_Electrons>-8*exp(-0.6*Leg2_P_InnerParam)");
c->SetAlias("eleParam","Leg1_TPC_nSigma_Electrons<5&&Leg2_TPC_nSigma_Electrons<5&&Leg1_TPC_nSigma_Electrons>-3.7*exp(-0.9*Leg1_P_InnerParam)-0.1&&Leg2_TPC_nSigma_Electrons>-3.7*exp(-0.6*Leg2_P_InnerParam)-0.1");
c->SetAlias("pidParam","eleParam&&cutP");



c->SetAlias("LegEta","abs(Leg1_Eta)<0.9&&abs(Leg2_Eta<0.9)");
c->SetAlias("LegNcl","Leg1_NclsTPC>90&&Leg2_NclsTPC>90");
c->SetAlias("LegPt","Leg1_Pt>1&&Leg2_Pt>1");
c->SetAlias("Rap","abs(Y)<0.9");
c->SetAlias("QA","LegNcl&&LegEta&&Rap");
c->SetAlias("spdFirst","(Leg1_ITS_clusterMap&1)==1 && (Leg2_ITS_clusterMap&1)==1");
c->SetAlias("LegNclDiffIter1","abs(Leg1_NclsTPC-Leg1_NclsTPCiter1)<10&&abs(Leg2_NclsTPC-Leg2_NclsTPCiter1)<10")
c->SetAlias("LegNclPID","(Leg1_NclsTPC-Leg1_TPCsignalN)<20&&(Leg2_NclsTPC-Leg2_TPCsignalN)<20")

c->SetAlias("LegNcl","Leg1_NFclsTPCrobust&&Leg2_NFclsTPCrobust");
c->SetAlias("cut","PairType==1&&QA&&pidSig&&LegPt");

c->SetAlias("cut","PairType==1&&QA&&cutPipardEdx&&cutPspecial")
c->SetAlias("cut","PairType==1&&QA&&pidSig&&LegPt")

c->SetMarkerStyle(22);
c->SetMarkerSize(.8);

//-------- nsigma

c->SetMarkerColor(kBlack);
c->SetLineColor(kBlack);

c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>3&&abs(Leg2_TPC_nSigma_Pions)>3");
c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>3&&(Leg2_TPC_nSigma_Protons)>3");

c->Draw("M>>hM(301,-.01,6.01)","cut","e");
c->GetHistogram()->GetXaxis()->SetRangeUser(2,4);

//--------

c->SetMarkerColor(kBlue);
c->SetLineColor(kBlue);

c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>3&&abs(Leg2_TPC_nSigma_Pions)>3");
c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>3&&(Leg2_TPC_nSigma_Protons)>3");

c->Draw("M>>hM2(301,-.01,6.01)","cut","esame");

//--------

c->SetMarkerColor(kGreen);
c->SetLineColor(kGreen);

c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>4&&abs(Leg2_TPC_nSigma_Pions)>4");
c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>3.5&&(Leg2_TPC_nSigma_Protons)>3.5");

c->Draw("M>>hM3(301,-.01,6.01)","cut","esame");

//--------

c->SetMarkerColor(kMagenta);
c->SetLineColor(kMagenta);

c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>3.5&&abs(Leg2_TPC_nSigma_Pions)>3.5");
c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>4&&(Leg2_TPC_nSigma_Protons)>4");

c->Draw("M>>hM4(301,-.01,6.01)","cut","esame");





//
//-------- rapidity
//

c->SetMarkerColor(kBlack);
c->SetLineColor(kBlack);

c->Draw("M>>hM(301,-.01,6.01)","cut&&Y<=0","e");

//--------

c->SetMarkerColor(kBlue);
c->SetLineColor(kBlue);

c->Draw("M>>hM2(301,-.01,6.01)","cut&&Y>0","esame");





c->SetAlias("cutE","Leg1_TPC_nSigma_Electrons>-1 && Leg2_TPC_nSigma_Electrons>-1");
c->SetAlias("LegPt","Leg1_Pt>1.2&&Leg2_Pt>1.2");

c->SetAlias("cut","PairType==1&&QA&&pidSig");


//-------- binning

c->SetMarkerColor(kBlack);
c->SetLineColor(kBlack);

// c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>3&&abs(Leg2_TPC_nSigma_Pions)>3");
// c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>3&&(Leg2_TPC_nSigma_Protons)>3");

c->Draw("M>>hM(601,-.015,6.005)","cut","e");

TGraphErrors gr0;
TH1 *h=c->GetHistogram();

gr0->SetLineColor(h->GetLineColor());
gr0->SetMarkerColor(h->GetMarkerColor());
for (Int_t i=0;i<h->GetNbinsX();++i){
 gr0.SetPoint(i,h->GetXaxis()->GetBinCenter(i+1),h->GetBinContent(i+1));
//  gr0.SetPointError(i,h->GetXaxis()->GetBinWidth(i+1)/2,h->GetBinError(i+1));
 gr0.SetPointError(i,0,h->GetBinError(i+1));
}

//--------

c->SetMarkerColor(kBlue);
c->SetLineColor(kBlue);

c->Draw("M>>hM2(601,-.01,6.01)","cut","egoff");

TGraphErrors gr1;
TH1 *h=c->GetHistogram();

gr1->SetLineColor(h->GetLineColor());
gr1->SetMarkerColor(h->GetMarkerColor());
for (Int_t i=0;i<h->GetNbinsX();++i){
 gr1.SetPoint(i,h->GetXaxis()->GetBinCenter(i+1),h->GetBinContent(i+1));
//  gr1.SetPointError(i,h->GetXaxis()->GetBinWidth(i+1)/2,h->GetBinError(i+1));
 gr1.SetPointError(i,0,h->GetBinError(i+1));
}

//--------

c->SetMarkerColor(kGreen);
c->SetLineColor(kGreen);

c->Draw("M>>hM3(601,-.005,6.015)","cut","egoff");
c->GetHistogram()->GetXaxis()->SetRangeUser(2.,4);

TGraphErrors gr2;
TH1 *h=c->GetHistogram();

gr2->SetLineColor(h->GetLineColor());
gr2->SetMarkerColor(h->GetMarkerColor());
for (Int_t i=0;i<h->GetNbinsX();++i){
 gr2.SetPoint(i,h->GetXaxis()->GetBinCenter(i+1),h->GetBinContent(i+1));
//  gr2.SetPointError(i,h->GetXaxis()->GetBinWidth(i+1)/2,h->GetBinError(i+1));
 gr2.SetPointError(i,0,h->GetBinError(i+1));
}

//--------

c->SetMarkerColor(kMagenta);
c->SetLineColor(kMagenta);

c->Draw("M>>hM4(601,-.0,6.02)","cut","egoff");

TGraphErrors gr3;
TH1 *h=c->GetHistogram();

gr3->SetLineColor(h->GetLineColor());
gr3->SetMarkerColor(h->GetMarkerColor());
for (Int_t i=0;i<h->GetNbinsX();++i){
 gr3.SetPoint(i,h->GetXaxis()->GetBinCenter(i+1),h->GetBinContent(i+1));
//  gr3.SetPointError(i,h->GetXaxis()->GetBinWidth(i+1)/2,h->GetBinError(i+1));
 gr3.SetPointError(i,0,h->GetBinError(i+1));
}

gr0->Draw("ap");
gr0->GetHistogram()->GetXaxis()->SetRangeUser(2.,4);
gr1->Draw("p");
gr2->Draw("p");
gr3->Draw("p");




c->SetAlias("LegNclDiffIter1","(Leg1_NclsTPC-Leg1_NclsTPCiter1)>-1&&(Leg2_NclsTPC-Leg2_NclsTPCiter1)>-1")
c->SetAlias("cut","PairType==1&&LegNclDiffIter1")
// histos
AliDielectronHistos h("h","h");
h.AddClass("TPCsignal");
h.UserHistogram("TPCsignal","sigTPC","TPC signal;P [GeV];TPC signal [arb. Units]",400,.3,40,400,0.,200.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","sigTPC")->SetDirectory(gDirectory)

h.UserHistogram("TPCsignal","nSigE","TPC n #sigma Electrons;P [GeV];TPC n #sigma Electrons",200,.3,40.,100,-10.,10.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigE")->SetDirectory(gDirectory)
h.UserHistogram("TPCsignal","nSigMu","TPC n #sigma Muons;P [GeV];TPC n #sigma Muons",400,.3,40.,500,-10.,10.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigMu")->SetDirectory(gDirectory)
h.UserHistogram("TPCsignal","nSigPi","TPC n #sigma Pions;P [GeV];TPC n #sigma Pions",400,.3,40.,500,-10.,10.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigPi")->SetDirectory(gDirectory)
h.UserHistogram("TPCsignal","nSigK","TPC n #sigma Kaons;P [GeV];TPC n #sigma Kaons",400,.3,40,500,-10,10,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigK")->SetDirectory(gDirectory)
h.UserHistogram("TPCsignal","nSigP","TPC n #sigma Protons;P [GeV];TPC n #sigma Protons",400,.3,40.,500,-10,10.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigP")->SetDirectory(gDirectory)

h.UserHistogram("TPCsignal","nSigDiffP","ncls-nclsXX;P [GeV];ncls-nclsXX",400,.3,40.,200,-40,160.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigDiffP")->SetDirectory(gDirectory)


c->Draw("Leg1_TPC_signal:Leg1_P_InnerParam>>sigTPC","cut","colz")
c->Draw("Leg2_TPC_signal:Leg2_P_InnerParam>>+sigTPC","cut","colz")


c->Draw("Leg1_TPC_nSigma_Electrons:Leg1_P_InnerParam>>nSigE","cut","colz")
c->Draw("Leg2_TPC_nSigma_Electrons:Leg2_P_InnerParam>>+nSigE","cut","colzsame")

c->Draw("Leg1_TPC_nSigma_Muos:Leg1_P_InnerParam>>nSigMu","cut","goff")
c->Draw("Leg2_TPC_nSigma_Muons:Leg1_P_InnerParam>>nSigMu","cut","goff")
c->Draw("Leg2_TPC_nSigma_Muons:Leg2_P_InnerParam>>+nSigMu","cut","colz")

c->Draw("Leg1_TPC_nSigma_Pions:Leg1_P_InnerParam>>nSigPi","cut","goff")
c->Draw("Leg2_TPC_nSigma_Pions:Leg2_P_InnerParam>>+nSigPi","cut","colz")

c->Draw("Leg1_TPC_nSigma_Kaons:Leg1_P_InnerParam>>nSigK","cut","goff")
c->Draw("Leg2_TPC_nSigma_Kaons:Leg2_P_InnerParam>>+nSigK","cut","colz")

c->Draw("Leg1_TPC_nSigma_Protons:Leg1_P_InnerParam>>nSigP","cut","goff")
c->Draw("Leg2_TPC_nSigma_Protons:Leg2_P_InnerParam>>+nSigP","cut","colz")



c->Draw("Leg1_TOF_nSigma_Electrons:Leg1_P_InnerParam>>nSigE","cut","colz")
c->Draw("Leg2_TOF_nSigma_Electrons:Leg2_P_InnerParam>>+nSigE","cut","colz")


c->Draw("Leg1_TOF_nSigma_Protons:Leg1_P_InnerParam>>nSigP","cut","goff")
c->Draw("Leg2_TOF_nSigma_Protons:Leg2_P_InnerParam>>+nSigP","cut","colz")

c->SetScanField(0)
c->Scan("M:EventInFile:Leg1_ID:Leg2_ID:File.GetString()","","colsize=1 col=.4f:6.d:5.d:5.d:130.s")

c->Scan("1:Run:EventInFile:File.GetString()","","colsize=1 col=.2f:8.d:8.d:130.s");



AliDielectronSignalFunc sig;
sig.SetDefaults(1);


//----------
c->Draw("Leg1_NclsTPC>>hNcls(160,-0.5,159.5)","cut","goff");
c->Draw("Leg2_NclsTPC>>+hNcls","cut","goff");
hNclsPID->SetLineColor(kBlack)

c->Draw("Leg1_TPCsignalN>>hNclsPID(160,-0.5,159.5)","cut","goff");
c->Draw("Leg2_TPCsignalN>>+hNclsPID","cut","goff");
hNclsPID->SetLineColor(kBlue)

c->Draw("Leg1_NclsTPCiter1>>hNclsIter1(160,-0.5,159.5)","cut","goff");
c->Draw("Leg2_NclsTPCiter1>>+hNclsIter1","cut","goff");
hNclsIter1->SetLineColor(kGreen)

hNcls->Draw();
hNclsPID->Draw("same");
hNclsIter1->Draw("same");
//-----------



c->Draw("Leg1_TPCsignalN:Leg1_NclsTPC>>hNclsPIDNcls(160,-0.5,159.5,160,-0.5,159.5)","cut","colz");
c->Draw("Leg2_TPCsignalN:Leg2_NclsTPC>>+hNclsPIDNcls","cut","colz");


c->Draw("Leg1_NclsTPC-Leg1_TPCsignalN:Leg1_P_InnerParam>>nSigDiffP","cut","colz");
c->Draw("Leg2_NclsTPC-Leg2_TPCsignalN:Leg1_P_InnerParam>>+nSigDiffP","cut","colz");


c->Draw("Leg1_NclsTPC-Leg1_NclsTPCiter1:Leg1_P_InnerParam>>nSigDiffP","cut","colz");
c->Draw("Leg2_NclsTPC-Leg2_NclsTPCiter1:Leg1_P_InnerParam>>+nSigDiffP","cut","colz");









c->Draw("Leg1_NclsTPC-Leg1_TPCsignalN:Leg1_TPC_nSigma_Electrons:Leg1_P_InnerParam>>hXX(100,0,10,20,-4,4)","cut","profcolz")
c->Draw("Leg1_NclsTPC-Leg1_TPCsignalN:Leg2_TPC_nSigma_Electrons:Leg2_P_InnerParam>>+hXX","cut","profcolz")









//WooJins cuts:
c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<3&&abs(Leg2_TPC_nSigma_Electrons)<3");
c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>5&&abs(Leg2_TPC_nSigma_Pions)>5");
// c->SetAlias("cutPi","Leg1_TPC_signal>65&&Leg2_TPC_signal>65");
c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>5&&(Leg2_TPC_nSigma_Protons)>5");
c->SetAlias("cutK","abs(Leg1_TPC_nSigma_Kaons)>2&&abs(Leg2_TPC_nSigma_Kaons)>2");
c->SetAlias("pid","cutE&&cutPi&&cutP");
c->SetAlias("etaLeg","abs(Leg1_Eta)<.9&&abs(Leg2_Eta)<.9");
c->SetAlias("rap","abs(Y)<.9");


c->SetAlias("cutAdd","PairType==1&&abs(Leg1_ImpactParXY)<.02&&abs(Leg2_ImpactParXY)<.02&&Leg2_Pt>1.")
c->Draw("M>>hM(50,2,4)","PairType==1&&pid","e");
h.Rebin();
h.Rebin();
h.Rebin();
sig.Process(hM);
sig.Draw("samestat");




//test
c->SetAlias("cutAdd","PairType==1&&abs(Leg1_ImpactParXY)<.03&&abs(Leg2_ImpactParXY)<.03&&Leg2_Pt>1")
c->SetAlias("cut","cutAdd&&pid&etaLeg&&rap")

c->Draw("M>>hM(50,1.99,3.99)","cut","e");
sig.Process(hM);
sig.Draw("samestat");


c->SetAlias("cut","PairType==1")
c->Draw("Leg1_TPC_signal:Leg1_P_InnerParam>>sigTPC","cut","colz")
c->Draw("Leg2_TPC_signal:Leg2_P_InnerParam>>+sigTPC","cut","colz")

////
c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<2&&abs(Leg2_TPC_nSigma_Electrons)<2");
c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>2&&abs(Leg2_TPC_nSigma_Pions)>2");
c->SetAlias("cutPi2","Leg1_TPC_signal>65&&Leg2_TPC_signal>65");
c->SetAlias("cutPi3","abs(Leg1_TPC_nSigma_Pions)>2.5&&abs(Leg2_TPC_nSigma_Pions)>2.5");
c->SetAlias("cutP","abs(Leg1_TPC_nSigma_Protons)>2&&abs(Leg2_TPC_nSigma_Protons)>2");
c->SetAlias("cutP2","abs(Leg1_TPC_nSigma_Protons)>1.5&&abs(Leg2_TPC_nSigma_Protons)>1.5");
c->SetAlias("cutK","abs(Leg1_TPC_nSigma_Kaons)>2&&abs(Leg2_TPC_nSigma_Kaons)>2");
c->SetAlias("cutdXY","abs(Leg1_ImpactParXY)<.02&&abs(Leg2_ImpactParXY)<.02")
c->SetAlias("cutPt","Leg2_Pt>1.")
c->SetAlias("cutPt2","Leg2_Pt>1.2")

//----
c->SetMarkerSize(0.7);

c->SetMarkerStyle(20);
c->SetLineColor(kRed);
c->SetMarkerColor(kRed);
c->Draw("M>>hM(100,2,4)","cutPi&&cutPt","e");

c->SetMarkerStyle(21);
c->SetLineColor(kRed-1);
c->SetMarkerColor(kRed-2);
c->Draw("M>>hM2(100,2,4)","cutPi2&&cutPt","esame");

c->SetMarkerStyle(22);
c->SetLineColor(kRed-2);
c->SetMarkerColor(kRed-2);
c->Draw("M>>hM3(100,2,4)","cutPi3&&cutPt","esame");

//----
c->SetMarkerStyle(20);
c->SetLineColor(kBlue);
c->SetMarkerColor(kBlue);
c->Draw("M>>hM4(100,2,4)","cutPi&&cutPt&&cutP","esame");

c->SetMarkerStyle(21);
c->SetLineColor(kBlue-1);
c->SetMarkerColor(kBlue-1);
c->Draw("M>>hM5(100,2,4)","cutPi2&&cutPt&&cutP","esame");

c->SetMarkerStyle(22);
c->SetLineColor(kBlue-2);
c->SetMarkerColor(kBlue-2);
c->Draw("M>>hM6(100,2,4)","cutPi3&&cutPt&&cutP","esame");

//----

c->SetMarkerStyle(20);
c->SetLineColor(kGreen);
c->SetMarkerColor(kGreen);
c->Draw("M>>hM7(100,2,4)","cutPi&&cutPt&&cutP2","esame");

c->SetMarkerStyle(21);
c->SetLineColor(kGreen-1);
c->SetMarkerColor(kGreen-1);
c->Draw("M>>hM8(100,2,4)","cutPi2&&cutPt&&cutP2","esame");

c->SetMarkerStyle(22);
c->SetLineColor(kGreen-2);
c->SetMarkerColor(kGreen-2);
c->Draw("M>>hM9(100,2,4)","cutPi3&&cutPt&&cutP2","esame");

//----


c->SetMarkerStyle(20);
c->SetLineColor(kMagentha);
c->SetMarkerColor(kMagentha);
c->Draw("M>>hM7(100,2,4)","cutPi&&cutPt&&cutP2","esame");

c->SetMarkerStyle(21);
c->SetLineColor(kMagentha-1);
c->SetMarkerColor(kMagentha-1);
c->Draw("M>>hM8(100,2,4)","cutPi2&&cutPt&&cutP2","esame");

c->SetMarkerStyle(22);
c->SetLineColor(kMagentha-2);
c->SetMarkerColor(kMagentha-2);
c->Draw("M>>hM9(100,2,4)","cutPi3&&cutPt&&cutP2","esame");


c->SetLineColor(kBlack);
c->SetMarkerColor(kBlue);
c->Draw("M>>hM4(100,2,4)","cutE&&cutPi&&cutK&&cutP&&cutdXY&&cutPt","esame");




//
c->SetAlias("cutE","Leg1_TPC_nSigma_Electrons>-1.5&&Leg2_TPC_nSigma_Electrons>-1.5");
// c->SetAlias("cutE","Leg1_TPC_signal>60&&Leg2_TPC_signal>60");
c->SetAlias("cutP","abs(Leg1_TPC_nSigma_Protons)>3&&abs(Leg2_TPC_nSigma_Protons)>3")

c->SetAlias("cutAdd","PairType==1&&abs(Leg1_ImpactParXY)<.03&&abs(Leg2_ImpactParXY)<.03&&Leg2_Pt>0")
c->SetAlias("cut","Leg2_Pt>1&&cutE&&cutP")

c->Draw("M>>hM(50,1.99,3.99)","cut","e");

c->SetAlias("cutAdd","PairType==1&&abs(Leg1_ImpactParXY)<.03&&abs(Leg2_ImpactParXY)<.03&&Leg2_Pt>.8")
c->Draw("M>>hM2(50,1.99,3.99)","cut","esame");

c->SetAlias("cutAdd","PairType==1&&abs(Leg1_ImpactParXY)<.03&&abs(Leg2_ImpactParXY)<.03&&Leg2_Pt>1")
c->Draw("M>>hM3(50,1.99,3.99)","cut","esame");

c->SetAlias("cutAdd","PairType==1&&abs(Leg1_ImpactParXY)<.03&&abs(Leg2_ImpactParXY)<.03&&Leg2_Pt>1.2")
c->Draw("M>>hM4(50,1.99,3.99)","cut","esame");


c->Draw("Leg1_TPC_signal:Leg1_P_InnerParam>>sigTPC","cut","goff")
c->Draw("Leg2_TPC_signal:Leg2_P_InnerParam>>+sigTPC","cut","goff")
c1->Modified();c1->Update()


c->SetAlias("cutE","Leg1_TPC_nSigma_Electrons>-1.5&&Leg2_TPC_nSigma_Electrons>-1.5");
c->SetAlias("cut","Leg2_P_InnerParam>1.5&&cutE")
c->SetMarkerStyle(21);
c->Draw("M>>hM(50,1.99,3.99)","cut","e");

c->SetAlias("cutE","Leg1_TPC_nSigma_Electrons>-1.5+.8&&Leg2_TPC_nSigma_Electrons>-1.5+.8");
c->SetAlias("cut","Leg2_P_InnerParam>1.5&&cutE")
c->SetMarkerStyle(20);
c->SetMarkerColor(kRed);
c->Draw("M>>hM3(50,1.99,3.99)","cut","esame");




//=================== pass 1 pass 2 comparison =======================


//------ cuts ------------


c->SetAlias("cutPro","Leg1_TPC_signal>50*exp(-.5*(Leg1_P_InnerParam-2))&&Leg2_TPC_signal>50*exp(-.5*(Leg2_P_InnerParam-2))&&Leg1_TPC_signal<85&&Leg2_TPC_signal<85");

c->SetAlias("Pcut","Leg1_P_InnerParam>1.3&&Leg2_P_InnerParam>1.3")
c->SetAlias("Ptcut","Leg1_Pt>1&&Leg2_Pt>1")
c->SetAlias("TOFcut","abs(Leg1_TOF_nSigma_Electrons)<3&&abs(Leg2_TOF_nSigma_Electrons)<3");
c->SetAlias("TOFcut2","(Leg1_P_InnerParam<1.3&&abs(Leg1_TOF_nSigma_Electrons)<3||Leg1_P_InnerParam>=1.3)&&(Leg2_P_InnerParam<1.3&&abs(Leg2_TOF_nSigma_Electrons)<3||Leg2_P_InnerParam>=1.3)");
c->SetAlias("TPCcut","(Leg1_TPC_signal>70*(1-exp(-1*(Leg1_P_InnerParam+2))))&&(Leg2_TPC_signal>70*(1-exp(-1*(Leg2_P_InnerParam+2))))")
c->SetAlias("NClcut","Leg1_NclsTPC>90&&Leg2_NclsTPC>90");

c->SetAlias("eleParam","Leg1_TPC_nSigma_Electrons<5&&Leg2_TPC_nSigma_Electrons<5&&Leg1_TPC_nSigma_Electrons>-2.65*exp(-0.6757*Leg1_P_InnerParam)&&Leg2_TPC_nSigma_Electrons>-2.65*exp(-0.6757*Leg2_P_InnerParam)")
c->SetAlias("cut","PairType==1&&eleParam&&Run<127719")
c->SetAlias("cut","1==1")
c->SetAlias("cut","NClcut")

c->SetAlias("cut","TOFcut&&TPCcut&&NClcut")

c->SetAlias("cut","TOFcut2&&TPCcut&&NClcut")

c->SetAlias("cut","cutPro&&TPCcut&&NClcut")

c->SetAlias("cut","Pcut&&TPCcut&&NClcut")

c->SetAlias("cut","Ptcut&&TPCcut&&NClcut")



//------------ plots --------------

//no cut
c->SetAlias("cut","1==1")
c1->SetLogx(0)
c1->SetLogz(0)
c->SetLineColor(kBlack);
c->Draw("M>>hM(50,1.99,3.99)","cut","e");
hM->SetTitle(";Inv. Mass [GeV]; Pair (e+e-) / 40GeV")
c1->Modified()
c1->Update();
c1->Print("pics/M_noCut.png");

//=============
//ncl: No cut
//=============
c->SetAlias("cut","1==1")
c1->SetLogx(1)
c1->SetLogz(0)
gStyle->SetOptStat(0);
c->Draw("Leg1_NclsTPC:Leg1_TPC_signal:Leg1_P_InnerParam>>test(1000,0,40,200,60,100)","cut","profcolz")
c->Draw("Leg2_NclsTPC:Leg2_TPC_signal:Leg2_P_InnerParam>>+test","cut","profcolz")
test->SetMinimum(80)
test->SetTitle("mean TPC number of clusters;P_{TPC} [GeV]; TPC signal [arb. units]")
c1->Modified()
c1->Update();
c1->Print("pics/TPCnCl_P.png");

//=============
//TPCsignal: ncl cut
//=============
c->SetAlias("cut","NClcut")
c->Draw("Leg1_NclsTPC:Leg1_TPC_signal:Leg1_P_InnerParam>>test(1000,0,40,200,60,100)","cut","profcolz")
c->Draw("Leg2_NclsTPC:Leg2_TPC_signal:Leg2_P_InnerParam>>+test","cut","profcolz")
test->SetMinimum(80)
test->SetTitle("mean TPC number of clusters;P_{TPC} [GeV]; TPC signal [arb. units]")
c1->Modified()
c1->Update();
c1->Print("pics/TPCnCl_P_cutNcl.png");


//=============
//tpc signal + signal cut
//=============
c1->SetLogx(1)
c1->SetLogy(0)
c1->SetLogz(1)
h.GetHistogram("TPCsignal","sigTPC")->GetYaxis()->SetRangeUser(60,100);
c->SetAlias("cut","NClcut")
c->Draw("Leg1_TPC_signal:Leg1_P_InnerParam>>sigTPC","cut","colz")
c->Draw("Leg2_TPC_signal:Leg2_P_InnerParam>>+sigTPC","cut","colz")
TF1 f("f1","[0]*(1-exp(-[1]*(x-[2])))",0.3,40);
f.SetParameters(70,1,-2);
f.Draw("same");
c1->Modified();
c1->Update();
c1->Print("pics/TPCsignal_P_cutNcl.png");

//------- Mass

c1->SetLogx(0)
c1->SetLogy(1)
c1->SetLogz(0)
c->SetAlias("cut","1==1")
c->SetLineColor(kBlack);
c->SetMarkerColor(kBlack);
c->Draw("M>>hM(50,1.99,3.99)","cut","e");
hM->SetTitle(";Inv. Mass [GeV]; Pair (e+e-) / 40GeV")
hM->SetMinimum(5e2);
c->SetAlias("cut","NClcut")
c->SetLineColor(kRed);
c->SetMarkerColor(kRed);
c->Draw("M>>hM2(50,1.99,3.99)","cut","esame");
c1->Modified();
c1->Update();
c1->Print("pics/M_nClCut.png");


//==========
//tpc signal: ncl + tpc cut
//==========
c1->SetLogx(1)
c1->SetLogy(0)
c1->SetLogz(1)
c->SetAlias("cut","TPCcut&&NClcut")
c->Draw("Leg1_TPC_signal:Leg1_P_InnerParam>>sigTPC","cut","colz")
c->Draw("Leg2_TPC_signal:Leg2_P_InnerParam>>+sigTPC","cut","colz")
c1->Modified();
c1->Update();
c1->Print("pics/TPCsignal_P_cutNcl_tpc.png");

/--- Mass

c1->SetLogx(0)
c1->SetLogy(1)
c1->SetLogz(0)
c->SetAlias("cut","1==1")
c->SetLineColor(kBlack);
c->SetMarkerColor(kBlack);
c->Draw("M>>hM(50,1.99,3.99)","cut","e");
hM->SetTitle(";Inv. Mass [GeV]; Pair (e+e-) / 40GeV")
hM->SetMinimum(5);
c->SetAlias("cut","NClcut")
c->SetLineColor(kRed);
c->SetMarkerColor(kRed);
c->Draw("M>>hM2(50,1.99,3.99)","cut","esame");
c->SetAlias("cut","TPCcut&&NClcut")
c->SetLineColor(kGreen);
c->SetMarkerColor(kGreen);
c->Draw("M>>hM3(50,1.99,3.99)","cut","esame");
c1->Modified();
c1->Update();
c1->Print("pics/M_nClCut_tpc.png");


//========
//TPC signal: ncl + tpc +  tof cut
//=======
c1->SetLogx(1)
c1->SetLogy(0)
c1->SetLogz(1)
c->SetAlias("cut","TOFcut2&&TPCcut&&NClcut")
c->Draw("Leg1_TPC_signal:Leg1_P_InnerParam>>sigTPC","cut","colz")
c->Draw("Leg2_TPC_signal:Leg2_P_InnerParam>>+sigTPC","cut","colz")
c1->Modified();
c1->Update();
c1->Print("pics/TPCsignal_P_cutNcl_tpc.png");

//--- Mass

c1->SetLogx(0)
c1->SetLogy(0)
c1->SetLogz(0)
c->SetAlias("cut","1==1")
c->SetAlias("cut","TPCcut&&NClcut")
c->SetLineColor(kGreen);
c->Draw("M>>hM(50,1.99,3.99)","cut","e");
hM->SetTitle(";Inv. Mass [GeV]; Pair (e+e-) / 40GeV")
hM->SetMinimum(.1);
c->SetAlias("cut","TOFcut2&&TPCcut&&NClcut")
c->SetLineColor(kBlue);
c->Draw("M>>hM2(50,1.99,3.99)","cut","esame");
c1->Modified();
c1->Update();
c1->Print("pics/M_nClCut_tpc.png");

//========
//Inv Mass: different cuts
//=======

c1->SetLogx(0)
c1->SetLogy(0)
c1->SetLogz(0)
c->SetAlias("cut","Ptcut&&TPCcut&&NClcut")
c->SetLineColor(kMagenta);
c->SetMarkerColor(kMagenta);
c->SetMarkerStyle(22);
c->Draw("M>>hM(50,1.99,3.99)","cut","e");

c->SetAlias("cut","Pcut&&TPCcut&&NClcut")
c->SetLineColor(kCyan+1);
c->SetMarkerColor(kCyan+1);
c->SetMarkerStyle(21);
c->Draw("M>>hM2(50,1.99,3.99)","cut","esame");

c->SetAlias("cut","TOFcut2&&TPCcut&&NClcut")
c->SetMarkerStyle(20);
c->SetLineColor(kBlue);
c->SetMarkerColor(kBlue);
c->Draw("M>>hM3(50,1.99,3.99)","cut","esame");

c1->Modified();
c1->Update();
c1->Print("pics/M_nClCut_tpc_tof.png");




c->SetAlias("NClcut","Leg1_NclsTPC>90&&Leg2_NclsTPC>90");
c->SetAlias("Ptcut","Leg1_Pt>1&&Leg2_Pt>1")
c->SetAlias("PairT","PairType==1");
c->SetAlias("cut","NClcut&&Ptcut&&PairT")


c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<3&&abs(Leg2_TPC_nSigma_Electrons)<3");
c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>3&&abs(Leg2_TPC_nSigma_Pions)>3");
c->SetAlias("cutP","(Leg1_TPC_nSigma_Protons)>3&&(Leg2_TPC_nSigma_Protons)>3");
c->SetAlias("cut","NClcut&&Ptcut&&PairT&&cutE&&cutPi&&cutP")

c->SetAlias("eta","abs(Leg1_Eta)<0.88&&abs(Leg1_Eta)<0.88");
c->SetAlias("rap","abs(Y)<0.88");
c->SetAlias("spdFirst","(Leg1_ITS_clusterMap&1)==1&&(Leg2_ITS_clusterMap&1)==1")
c->SetAlias("cut","NClcut&&Ptcut&&PairT&&cutE&&cutPi&&cutP&&eta&&rap&&spdFirst")
 c->Draw("M>>hM(50,1.98,1.98+50*.04)","cut","e")


c->SetAlias("cutProL1",Form("Leg1_TPC_nSigma_Electrons>(%f*%f+(AliExternalTrackParam::BetheBlochAleph(Leg1_P_InnerParam/%f,0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00)-AliExternalTrackParam::BetheBlochAleph(Leg1_P_InnerParam/%f,0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00)))/%f",nSigma,resolution, AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kElectron), resolution))

c->SetAlias("cutProL2",Form("Leg2_TPC_nSigma_Electrons>(%f*%f+(AliExternalTrackParam::BetheBlochAleph(Leg2_P_InnerParam/%f,0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00)-AliExternalTrackParam::BetheBlochAleph(Leg2_P_InnerParam/%f,0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00)))/%f",nSigma,resolution, AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kElectron), resolution))

c->SetAlias("cutPioL1",Form("Leg1_TPC_nSigma_Electrons>(%f*%f+(AliExternalTrackParam::BetheBlochAleph(Leg1_P_InnerParam/%f,0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00)-AliExternalTrackParam::BetheBlochAleph(Leg1_P_InnerParam/%f,0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00)))/%f",nSigma,resolution, AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kElectron), resolution))

c->SetAlias("cutPioL2",Form("Leg2_TPC_nSigma_Electrons>(%f*%f+(AliExternalTrackParam::BetheBlochAleph(Leg2_P_InnerParam/%f,0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00)-AliExternalTrackParam::BetheBlochAleph(Leg2_P_InnerParam/%f,0.0283086/0.97,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00)))/%f",nSigma,resolution, AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kElectron), resolution))

c->SetAlias("cutPro","cutProL1&&cutProL2");
c->SetAlias("cutPio","cutPioL1&&cutPioL2");


//
// track in or out of TRD
//
c->SetAlias("TRD11","(Leg1_TrackStatus&512)==512 && (Leg2_TrackStatus&512)==512");
c->SetAlias("TRD00","(Leg1_TrackStatus&512)==  0 && (Leg2_TrackStatus&512)==  0");
c->SetAlias("TRD01","(Leg1_TrackStatus&512)==  0 && (Leg2_TrackStatus&512)==512");
c->SetAlias("TRD10","(Leg1_TrackStatus&512)==512 && (Leg2_TrackStatus&512)==  0");

TCanvas *c2=(TCanvas)gROOT->FindObject("c2");
if (!c2) c2=new TCanvas("c2","c2");
c2->Clear()
c2->Divide(2,2);

c->SetAlias("cut","PairType==1 && QA && pidSig && ITS2");

TLine l;
l.SetLineColor(kGray+2);
// l.SetLineWidth(2);
TBox b;
b.SetFillColor(kGray);
b.SetFillStyle(3004);

c2->cd(1);
c->Draw("M>>hM11(125,0.,5.)","cut&TRD11","e");
b.DrawBox(2.92,0,3.16,gPad->GetUymax())
l.DrawLine(3.095,0,3.095,gPad->GetUymax());
c2->cd(2);
c->Draw("M>>hM10(125,0.,5.)","cut&TRD10","e");
b.DrawBox(2.92,0,3.16,gPad->GetUymax())
l.DrawLine(3.095,0,3.095,gPad->GetUymax());
c2->cd(3);
c->Draw("M>>hM01(125,0.,5.)","cut&TRD01","e");
b.DrawBox(2.92,0,3.16,gPad->GetUymax())
l.DrawLine(3.095,0,3.095,gPad->GetUymax());
c2->cd(4);
c->Draw("M>>hM00(125,0.,5.)","cut&TRD00","e");
b.DrawBox(2.92,0,3.16,gPad->GetUymax())
l.DrawLine(3.095,0,3.095,gPad->GetUymax());


//
// number of cluster in ITS
//

TH1F h("hITS","hITS",10,0,10);
for (Int_t i=0; i<6; ++i){
  c->Draw(Form("((Leg2_ITS_clusterMap&%d)==%d)*(%d+1)>>+hITS",i,i,i),"cut");
}

c->SetAlias("nclsITS1","((Leg2_ITS_clusterMap&1)==1) + ((Leg2_ITS_clusterMap&2)==2) + ((Leg2_ITS_clusterMap&4)==4) + ((Leg2_ITS_clusterMap&8)==8) + ((Leg2_ITS_clusterMap&16)==16) + ((Leg2_ITS_clusterMap&32)==32)")

c->SetAlias("nclsITS2","((Leg1_ITS_clusterMap&1)==1) + ((Leg1_ITS_clusterMap&2)==2) + ((Leg1_ITS_clusterMap&4)==4) + ((Leg1_ITS_clusterMap&8)==8) + ((Leg1_ITS_clusterMap&16)==16) + ((Leg1_ITS_clusterMap&32)==32)")




*/

