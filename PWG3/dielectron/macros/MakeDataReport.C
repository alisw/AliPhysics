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
  Double_t res=5.7e-2;
  alephParameters[0] = 0.0283086;
  alephParameters[1] = 2.63394e+01;
  alephParameters[2] = 5.04114e-11;
  alephParameters[3] = 2.12543e+00;
  alephParameters[4] = 4.88663e+00;
  Color_t color=kRed;



  TF1 *foDataProton = new TF1("foDataProton", "50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foDataProtonP = new TF1("foDataProtonP",Form( "50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1+%f)",res),0.05,20);
  TF1 *foDataProtonM = new TF1("foDataProtonM", Form("50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1-%f)",res),0.05,20);

  TF1 *foDataPion = new TF1("foDataPion", "50*AliExternalTrackParam::BetheBlochAleph(x/0.13957,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foDataPionP = new TF1("foDataPionP",Form( "50*AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])*(1+%f)",res),0.05,20);
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

  Double_t sigmas[5]={2,2,2,2,2};
  Double_t masses[5];

  for (Int_t i=0; i<AliPID::kSPECIES; ++i) masses[i]=AliPID::ParticleMass(i);
  
  Double_t res=7.e-2;
  Double_t alephParameters[5];

  alephParameters[0] = 0.0283086;
  alephParameters[1] = 2.63394e+01;
  alephParameters[2] = 5.04114e-11;
  alephParameters[3] = 2.12543e+00;
  alephParameters[4] = 4.88663e+00;
  Double_t mip=49.2;

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
  fBetheP[i]->Draw("same");
  fBetheM[i]->Draw("same");
}
}

*/

/*

c->SetLineColor(kRed); c->Draw("M","Leg1_Pt>1.2","same");
c->SetLineColor(kGreen); c->Draw("M","Leg1_Pt>1.5","same");
c->SetLineColor(kBlue); c->Draw("M","Leg1_Pt>2","same");
c->SetLineColor(kMagenta); c->Draw("M","Leg1_Pt>2.5","same");
c->SetLineColor(kGreen-4); c->Draw("M","Leg1_Pt>3","same");

c->SetLineColor(kRed-2); c->Draw("M","Leg1_Pt>1.2&&Leg2_Pt>.5","same");
c->SetLineColor(kRed-4); c->Draw("M","Leg1_Pt>1.2&&Leg2_Pt>.7","same");
c->SetLineColor(kRed-6); c->Draw("M","Leg1_Pt>1.2&&Leg2_Pt>1","same");
c->SetLineColor(kRed-8); c->Draw("M","Leg1_Pt>1.2&&Leg2_Pt>1.2","same");


c=Pair


c->SetAlias("cutE","Leg1_TPC_nSigma_Electrons>-2&&(Leg2_TPC_nSigma_Electrons)>-2");
c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<2&&abs(Leg2_TPC_nSigma_Electrons)<2");

c->SetAlias("cutPi","Leg1_TPC_nSigma_Pions>2&&Leg2_TPC_nSigma_Pions>2");
c->SetAlias("cutP","abs(Leg1_TPC_nSigma_Protons)>2&&abs(Leg2_TPC_nSigma_Protons)>2");
c->SetAlias("cutK","abs(Leg1_TPC_nSigma_Kaons)>2&&abs(Leg2_TPC_nSigma_Kaons)>2");
c->SetAlias("cutMu","abs(Leg1_TPC_nSigma_Muons)>1.5&&abs(Leg2_TPC_nSigma_Muons)>1.5");

c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<3&&abs(Leg2_TPC_nSigma_Electrons)<3");
c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>2&&abs(Leg2_TPC_nSigma_Pions)>2");
c->SetAlias("cutP","abs(Leg1_TPC_nSigma_Protons+.5)>2&&abs(Leg2_TPC_nSigma_Protons+.5)>2");
c->SetAlias("pid","cutPi&&cutE");



c->SetAlias("pid","cutPi&&cutP&&cutK&&cutE");
c->SetAlias("pid","cutPi&&cutP&&cutK&&cutMu");
c->Draw("M>>h(50,1.99,3.99)","pid&&PairType==1&&Leg2_Pt>1","e");
hist=h
TObjArray arr;
arr.AddAt(hist,1);

AliDielectronSignalFunc fun
fun.SetDefaults(1);
fun.Process(&arr);
fun.Draw("samestat");


c->SetLineColor(kBlack);
c->Draw("M>>h(50,1.99,3.99)","cutPi&&PairType==1","e");

c->SetLineColor(kRed);
c->Draw("M>>h2(50,1.99,3.99)","cutPi&&cutK&&PairType==1","esame");

c->SetLineColor(kBlue);
c->Draw("M>>h3(50,1.99,3.99)","cutPi&&cutP&&PairType==1","esame");

c->SetLineColor(kGreen);
c->Draw("M>>h4(50,1.99,3.99)","cutPi&&cutMu&&PairType==1","esame");







c->SetLineColor(kBlack);
c->Draw("Leg1_Pt>>hPt(50,1.99,3.99)","cutPi&&PairType==1","goff");
c->Draw("Leg2_Pt>>hPt+","cutPi&&PairType==1","e");

c->SetLineColor(kRed);
c->Draw("Leg1_Pt>>hPt2(50,1.99,3.99)","cutPi&&cutK&&PairType==1","goff");
c->Draw("Leg2_Pt>>hPt2+","cutPi&&cutK&&PairType==1","esame");

c->SetLineColor(kBlue);
c->Draw("Leg1_Pt>>hPt3(50,1.99,3.99)","cutPi&&cutP&&PairType==1","goff");
c->Draw("Leg2_Pt>>hPt3+","cutPi&&cutP&&PairType==1","esame");

c->SetLineColor(kGreen);
c->Draw("Leg1_Pt>>hPt4(50,1.99,3.99)","cutPi&&cutLeg1_Ptu&&PairType==1","goff");
c->Draw("Leg1_Pt>>hPt4+","cutPi&&cutMu&&PairType==1","esame");




c->Draw("M>>hM5(100,1.99,3.99)","Leg1_TPC_signal>65&&Leg2_TPC_signal>65&&Leg2_Pt>1.3&&PairType==1","e");


c=Pair
c->Draw("M>>hM5(100,1.99,3.99)","Leg1_TPC_signal>65&&Leg2_TPC_signal>65&&Leg2_Pt>1.&&PairType==1","e")
hist=hM5
AliDielectronSignalFunc fun
fun.SetDefaults(1);
fun.SetFitRange(2,4)
fun.Process(hM)
fun.Draw("samestat");

c->Draw("M>>hM5(50,1.99,3.99)","Leg1_TPC_signal>65&&Leg2_TPC_signal>65&&Leg2_Pt>1.1&&PairType==1","e")




c->SetAlias("cut","Leg1_TPC_signal>65&&Leg2_TPC_signal>65&&Leg2_Pt>1.&&PairType==1")
c->SetAlias("cut","Leg2_Pt>1&&PairType==1")
c->SetAlias("cut","Leg1_TPC_signal>65&&Leg2_TPC_signal>65&&cutP&&cutK&&PairType==1")
c->SetAlias("cut","cutP&&PairType==1")
c->SetAlias("cut","PairType==1&&abs(Leg1_TOF_nSigma_Protons)>2&&abs(Leg2_TOF_nSigma_Protons)>2")

c->SetAlias("cut","abs(Leg2_TOF_nSigma_Protons)>3&&abs(Leg2_TOF_nSigma_Pions)>")
c->SetAlias("cut","abs(Leg2_TOF_nSigma_Protons)>7-5.5/3*Leg2_P_InnerParam")
// histos
AliDielectronHistos h("h","h");
h.AddClass("TPCsignal");
h.UserHistogram("TPCsignal","sigTPC","TPC signal;P [GeV];TPC signal [arb. Units]",400,.3,40,400,0.,200.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","sigTPC")->SetDirectory(gDirectory)

h.UserHistogram("TPCsignal","nSigE","TPC n #sigma Electrons;P [GeV];TPC n #sigma Electrons",400,.3,40.,400,-4.,4.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigE")->SetDirectory(gDirectory)
h.UserHistogram("TPCsignal","nSigMu","TPC n #sigma Muons;P [GeV];TPC n #sigma Muons",400,.3,40.,400,-4.,4.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigMu")->SetDirectory(gDirectory)
h.UserHistogram("TPCsignal","nSigPi","TPC n #sigma Pions;P [GeV];TPC n #sigma Pions",400,.3,40.,400,-4.,4.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigPi")->SetDirectory(gDirectory)
h.UserHistogram("TPCsignal","nSigK","TPC n #sigma Kaons;P [GeV];TPC n #sigma Kaons",400,.3,40,400,-4,4,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigK")->SetDirectory(gDirectory)
h.UserHistogram("TPCsignal","nSigP","TPC n #sigma Protons;P [GeV];TPC n #sigma Protons",400,.3,40.,400,-4,4.,0,0,kTRUE,kFALSE)
h.GetHistogram("TPCsignal","nSigP")->SetDirectory(gDirectory)


c->Draw("Leg1_TPC_signal:Leg1_P_InnerParam>>sigTPC","cut","colz")
c->Draw("Leg2_TPC_signal:Leg2_P_InnerParam>>+sigTPC","cut","colz")


c->Draw("Leg1_TPC_nSigma_Electrons:Leg1_P_InnerParam>>nSigE","cut","colz")
c->Draw("Leg2_TPC_nSigma_Electrons:Leg2_P_InnerParam>>+nSigE","cut","colz")

c->Draw("Leg1_TPC_nSigma_Muos:Leg1_P_InnerParam>>nSigMu","cut","goff")
c->Draw("Leg2_TPC_nSigma_Muons:Leg1_P_InnerParam>>nSigMu","cut","goff")
c->Draw("Leg2_TPC_nSigma_Muons:Leg2_P_InnerParam>>+nSigMu","cut","colz")

c->Draw("Leg1_TPC_nSigma_Pions:Leg1_P_InnerParam>>nSigPi","cut","goff")
c->Draw("Leg2_TPC_nSigma_Pions:Leg2_P_InnerParam>>+nSigPi","cut","colz")

c->Draw("Leg1_TPC_nSigma_Kaons:Leg1_P_InnerParam>>nSigK","cut","goff")
c->Draw("Leg2_TPC_nSigma_Kaons:Leg2_P_InnerParam>>+nSigK","cut","colz")

c->Draw("Leg1_TPC_nSigma_Protons+.5:Leg1_P_InnerParam>>nSigP","cut","goff")
c->Draw("Leg2_TPC_nSigma_Protons+.5:Leg2_P_InnerParam>>+nSigP","cut","colz")




c->Draw("Leg1_TOF_nSigma_Protons:Leg1_P_InnerParam>>nSigP","cut","goff")
c->Draw("Leg2_TOF_nSigma_Protons:Leg2_P_InnerParam>>+nSigP","cut","colz")


Pair->SetScanField(0)
Pair->Scan("EventInFile:File.GetString()","","colsize=1 col=3.d:100.s")


AliDielectronSignalFunc sig;
sig.SetDefaults(1);

//WooJins cuts:
c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<3&&abs(Leg2_TPC_nSigma_Electrons)<3");
c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>2&&abs(Leg2_TPC_nSigma_Pions)>2");
// c->SetAlias("cutPi","Leg1_TPC_signal>65&&Leg2_TPC_signal>65");
c->SetAlias("cutP","abs(Leg1_TPC_nSigma_Protons)>2&&abs(Leg2_TPC_nSigma_Protons)>2");
c->SetAlias("cutK","abs(Leg1_TPC_nSigma_Kaons)>2&&abs(Leg2_TPC_nSigma_Kaons)>2");
c->SetAlias("pid","cutE&&cutPi&&cutP&&cutK");


c->SetAlias("cutAdd","PairType==1&&abs(Leg1_ImpactParXY)<.02&&abs(Leg2_ImpactParXY)<.02&&Leg2_Pt>1.")
c->Draw("M>>hM(100,2,4)","cutAdd&&pid","e");
h.Rebin();
h.Rebin();
h.Rebin();
sig.Process(hM);
sig.Draw("samestat");




//test
c->SetAlias("cutE","abs(Leg1_TPC_nSigma_Electrons)<3&&abs(Leg2_TPC_nSigma_Electrons)<2");
c->SetAlias("cutPi","abs(Leg1_TPC_nSigma_Pions)>1&&abs(Leg2_TPC_nSigma_Pions)>2");
c->SetAlias("cutPi","Leg1_TPC_signal>67&&Leg2_TPC_signal>67");
c->SetAlias("cutP","abs(Leg1_TPC_nSigma_Protons)>2&&abs(Leg2_TPC_nSigma_Protons)>2");
c->SetAlias("cutK","abs(Leg1_TPC_nSigma_Kaons)>2&&abs(Leg2_TPC_nSigma_Kaons)>2");


c->SetAlias("pid","cutE&&cutPi&&cutP&&cutK");
c->SetAlias("cutAdd","PairType==1&&abs(Leg1_ImpactParXY)<.03&&abs(Leg2_ImpactParXY)<.03&&Leg2_Pt>1")
c->SetAlias("cut","cutAdd&&pid")

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
c->Draw("M>>hM2(50,1.99,3.99)","cut","esame");

*/

