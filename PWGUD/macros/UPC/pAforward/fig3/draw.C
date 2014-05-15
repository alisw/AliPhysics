#include "TLatex.h"

void SetGraphStyle(TGraph* g, Color_t mcolor, Style_t mstyle, Size_t msize, Color_t lcolor, Style_t lstyle, Width_t lwidth);
TGraphErrors* read_jmrt(const char* fileName, const Int_t nPoints = 147);
TGraph* read_bsat(const char* fileName, const Int_t nPoints=15, Bool_t skipLowEnergy=1);
TGraphErrors* read_h1(const char* fileName);
TGraphErrors* read_zeus_ee();
TGraphErrors* read_zeus_mm();
TGraphAsymmErrors* read_alice(Bool_t statOnly=0,Int_t option=0);
TGraphErrors* read_lhcb(const char* fileName, const int nPoints=20);
TGraphErrors* read_clark();
TGraphErrors* read_slac();
TGraphErrors* read_gittelman();
TGraphErrors* read_e814();
TGraphErrors* read_e516();

void draw(){
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02,"X");  
  gStyle->SetTickLength(0.02,"Y"); 
  gStyle->SetLineScalePS(2);
  TGraphErrors* gJMRT_LO  = read_jmrt("JMRT-LO-JPSI.txt");
  TGraphErrors* gJMRT_NLO = read_jmrt("JMRT-NLO-JPSI.txt");
  TGraph* gBsat_eik = read_bsat("crossy_jpsi.dat.bSat_7TeV");
  TGraph* gBsat_pom = read_bsat("crossy_jpsi.dat.bSat_7TeVnoeik");
  TF1* fStarlight = new TF1("fStarlight","((1.0-( (4.035*4.035)/(x*x) ))**2.0)*4.1*exp(0.65*log(x))",4.1,2000);
  TGraphErrors* gH1 = read_h1("d13-058.table_w_allH1_v1.txt");
  TGraphErrors* gZeus_ee = read_zeus_ee();
  TGraphErrors* gZeus_mm = read_zeus_mm();
  TGraphAsymmErrors* gALICEstat = read_alice(1);
  TGraphAsymmErrors* gALICEfull = read_alice(0);
  TGraphAsymmErrors* gALICEpPb  = read_alice(0,1);
  TGraphAsymmErrors* gALICEPbp  = read_alice(0,2);
  TGraphErrors* gLHCb = read_lhcb("lhcb_data2.txt");
  TGraphErrors* gClark     = read_clark();
  TGraphErrors* gSLAC      = read_slac();
  TGraphErrors* gGittelman = read_gittelman();
  TGraphErrors* gE814      = read_e814();
  TGraphErrors* gE516      = read_e516();
  fStarlight->SetLineWidth(2);
  fStarlight->SetLineColor(kBlack);
  SetGraphStyle(gJMRT_LO  ,     1,kDot             ,  1,kRed   , 7,2);
  SetGraphStyle(gJMRT_NLO ,     1,kDot             ,  1,kBlue  , 5,2);
  SetGraphStyle(gBsat_eik ,     1,kDot             ,  1,kGreen+2, 9,2);
  SetGraphStyle(gBsat_pom ,     1,kDot             ,  1,kGreen+2,10,2);
  SetGraphStyle(gH1       ,kBlack,kFullCircle      ,1.4,kBlack , 1,1);
  SetGraphStyle(gZeus_ee  ,kBlack,kOpenSquare      ,1.4,kBlack , 1,1);
  SetGraphStyle(gZeus_mm  ,kBlack,kOpenSquare      ,1.4,kBlack , 1,1);
  SetGraphStyle(gLHCb     ,kBlack,kOpenCircle      ,1.4,kBlack , 1,1);
  SetGraphStyle(gSLAC     ,kBlack,kFullTriangleDown,1.2,kBlack , 1,1);
  SetGraphStyle(gClark    ,kBlack,kFullTriangleDown,1.2,kBlack , 1,1);
  SetGraphStyle(gGittelman,kBlack,kFullTriangleDown,1.2,kBlack , 1,1);
  SetGraphStyle(gE814     ,kBlack,kFullTriangleDown,1.2,kBlack , 1,1);
  SetGraphStyle(gALICEfull,kRed  ,kFullSquare      ,1.4,kRed   , 1,2);
  SetGraphStyle(gALICEstat,kRed  ,kFullSquare      ,1.4,kRed   , 1,2);
  SetGraphStyle(gALICEpPb ,kRed  ,kFullSquare      ,1.4,kRed   , 1,2);
  SetGraphStyle(gALICEPbp ,kRed  ,kFullDiamond     ,2.4,kRed   , 1,2);
  
  TCanvas* c1 = new TCanvas("c1","c1",1000,700);
  gPad->SetLeftMargin(0.09);
  gPad->SetRightMargin(0.003);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.12);
  gPad->SetLogx();
  gPad->SetLogy();
  TH1F* frame1 = gPad->DrawFrame(19.9999,9.99,1400,1000);
  frame1->SetTitle(";W_{#gammap} [GeV];#sigma(#gamma+p #rightarrow J/#psi+p) [nb]");
  frame1->GetXaxis()->SetTitleOffset(1.35);
  frame1->GetXaxis()->SetLabelSize(0.045);
  frame1->GetYaxis()->SetLabelSize(0.045);
  frame1->GetXaxis()->SetTitleSize(0.045);
  frame1->GetYaxis()->SetTitleSize(0.045);
  fStarlight->Draw("same");
  gJMRT_LO->Draw("cxsame");
  gJMRT_NLO->Draw("cxsame");
  gBsat_eik->Draw("lsame");
  gBsat_pom->Draw("lsame");
  gH1->Draw("pzsame");
  gZeus_ee->Draw("pzsame");
  gZeus_mm->Draw("pzsame");
  gLHCb->Draw("pzsame");
//  gE814->Draw("pzsame");
  gClark->Draw("pzsame");
  gALICEpPb->Draw("pzsame");
  gALICEPbp->Draw("pzsame");
//  gALICEfull->Draw("pzsame");
  TF1* fAliceFit = new TF1("fAliceFit","[0]*x^[1]",0,1000);
  gALICEfull->Fit(fAliceFit,"0");
  
  TArrow* arrow = new TArrow(100,10,1000,10,0.013,"<|>");
  arrow->SetFillColor(kRed);
  arrow->SetLineColor(kRed);
  arrow->SetLineWidth(2);
  arrow->DrawArrow(21,11.2,45,11.2);
  arrow->DrawArrow(577,11.2,952,11.2);
  arrow->DrawArrow(25,700,32,700);
  
  TLegend* leg1 = new TLegend(0.10,0.60,0.6,0.93);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry(arrow,"W_{#gammap} interval probed by ALICE","<|>");
  leg1->AddEntry(gALICEpPb,"ALICE (p-Pb)","p");
  leg1->AddEntry(gALICEPbp,"ALICE (Pb-p)","p");
  leg1->AddEntry(gLHCb,"LHCb solutions (pp)","p");
  leg1->AddEntry(gH1,"H1","p");
  leg1->AddEntry(gZeus_ee,"ZEUS","p");
//  leg1->AddEntry(gE814,"Fixed target experiments","p");
  leg1->Draw("same");

  TLegend* leg2 = new TLegend(0.45,0.18,1.00,0.45);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(gJMRT_LO,"JMRT LO","l");
  leg2->AddEntry(gJMRT_NLO,"JMRT NLO","l");
  leg2->AddEntry(gBsat_eik,"b-Sat (eikonalized)","l");
  leg2->AddEntry(gBsat_pom,"b-Sat (1-Pomeron)","l");
  leg2->AddEntry(fStarlight,"STARLIGHT parameterization","l");
  leg2->Draw("same");
  
//  TLatex* l = new TLatex();
//  l->SetTextAlign(22);
//  l->SetTextFont(42);
//  l->DrawLatex(160,7,"W_{#gammap} interval probed by ALICE");
  

  c1->Print("fig3.eps");
  c1->Print("fig3.png");
  
}

void SetGraphStyle(TGraph* g, Color_t mcolor, Style_t mstyle, Size_t msize, Color_t lcolor, Style_t lstyle, Width_t lwidth){
  g->SetMarkerColor(mcolor);
  g->SetMarkerSize(msize);
  g->SetMarkerStyle(mstyle);
  g->SetLineColor(lcolor);
  g->SetLineStyle(lstyle);
  g->SetLineWidth(lwidth);
}

TGraphErrors* read_jmrt(const char* fileName, const Int_t nPoints){
  ifstream f;
  f.open(fileName);
  char b[1000];
  getline(f,b);
  getline(f,b);
  Double_t q,w[nPoints],cs[nPoints],dcs[nPoints];
  for (Int_t i=0;i<nPoints;i++)  f >> q >> w[i] >> cs[i] >> dcs[i];
  f.close();
  return new TGraphErrors(nPoints,w,cs,NULL,dcs);
}

TGraph* read_bsat(const char* fileName, const Int_t nPoints, Bool_t skipLowEnergy){
  ifstream f;
  f.open(fileName);
  char b[1000];
  Double_t y,w[nPoints],q,cs[nPoints];
  if (skipLowEnergy) for (Int_t i=0;i<6;i++) getline(f,b);
  for (Int_t i=0;i<nPoints;i++) f >> y >> w[i] >> q >> cs[i];
  f.close();
  return new TGraph(nPoints,w,cs);
}

TGraphErrors* read_h1(const char* fileName){
  // read H1 data from http://www-h1.desy.de/psfiles/figures/d13-058.table_w_allH1_v1.txt
  ifstream f;
  f.open(fileName);
  char b[1000];
  Double_t sigm[29];
  Double_t dsig[29];
  Double_t wavg[29];
  Double_t n,pd,wmin,wmax,PhiT;
  // header
  for (Int_t i=0;i<65;i++) getline(f,b);
  Int_t nExclusive=0;
  for (Int_t i=0;i<40;i++){
    f >> n >> pd; 
    if (pd==1) { getline(f,b); continue; }
    f >> wmin >> wmax >> wavg[nExclusive] >> PhiT >> sigm[nExclusive] >> dsig[nExclusive];
    getline(f,b);
    nExclusive++;
  }
  return new TGraphErrors(nExclusive,wavg,sigm,NULL,dsig);
}


TGraphErrors* read_zeus_mm(){
  // ZEUS points
  const int nZEUSmm=8;
  Double_t wavgZEUSmm[]    = { 25.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0 };
  Double_t sigmZEUSmm[]    = { 32.6, 41.5, 55.8, 66.6,  73.4,  86.7, 104.0, 110.0 };
  Double_t statZEUSmm[]    = { 5.4, 1.1, 1.5, 2.0, 2.3, 3.2, 5.0, 11.0 };
  Double_t systZEUSmm[]    = { 5.2, 3.3, 4.6, 7.0, 6.0, 6.5,11.0, 12.0 };
  Double_t dsigZEUSmm[nZEUSmm];
  for (Int_t i=0;i<nZEUSmm;i++) dsigZEUSmm[i]=sqrt(pow(statZEUSmm[i],2.)+pow(systZEUSmm[i],2.));
  return new TGraphErrors(nZEUSmm,wavgZEUSmm,sigmZEUSmm,NULL,dsigZEUSmm);
}

TGraphErrors* read_zeus_ee(){
  const int nZEUSee=14;
  Double_t wavgZEUSee[]    = { 27.5, 42.5, 55.0, 65.0, 75.0, 85.0, 100.0, 117.5, 132.5, 155.0, 185.0, 215.0, 245.0, 275.0 };
  Double_t sigmZEUSee[]    = { 33.6, 43.8, 57.2, 62.5, 68.9, 72.1, 81.9, 95.7, 103.9, 115.0, 129.1, 141.7, 140.3, 189.0 };
  Double_t statZEUSee[]    = { 1.6, 2.0, 1.8, 2.3, 2.6, 2.9, 2.3, 3.2, 3.6, 3.3, 4.7, 6.1, 7.4, 13.0 };
  Double_t systZEUSee[]    = { 2.4, 2.8, 3.5, 3.9, 4.5, 4.5, 4.8, 5.4, 5.8, 6.7, 7.7, 8.7, 9.9, 26.0 };
  Double_t dsigZEUSee[nZEUSee];
  for (Int_t i=0;i<nZEUSee;i++) dsigZEUSee[i]=sqrt(pow(statZEUSee[i],2.)+pow(systZEUSee[i],2.));
  return new TGraphErrors(nZEUSee,wavgZEUSee,sigmZEUSee,NULL,dsigZEUSee);
}

TGraphAsymmErrors* read_alice(Bool_t statOnly, Int_t option){
  const int nALICE=4;
  Double_t wavgALICE[]     = { 24.1, 30.9, 39.6,   706};
  Double_t sigmALICE[]     = { 26.6, 33.6, 36.9,   275};
  Double_t statALICE[]     = {  3.5,  3.0,  5.3,    35};
  Double_t systALICEl[]    = {  2.7,  2.7,  3.9,    31};
  Double_t systALICEh[]    = {  2.7,  2.7,  4.0,    26};
  Double_t fluxALICE[]     = {  0.5,  0.7,  0.7,    25};
  if (statOnly) {
    if (option==0) return new TGraphAsymmErrors(nALICE,wavgALICE,sigmALICE,NULL,NULL,statALICE,statALICE);
    if (option==1) return new TGraphAsymmErrors(3,wavgALICE,sigmALICE,NULL,NULL,statALICE,statALICE);
    if (option==2) return new TGraphAsymmErrors(1,&(wavgALICE[3]),&(sigmALICE[3]),NULL,NULL,&(statALICE[3]),&(statALICE[3]));
  }
  Double_t dsigALICEl[nALICE];
  Double_t dsigALICEh[nALICE];
  for (Int_t i=0;i<nALICE;i++) dsigALICEl[i] = sqrt(pow(statALICE[i],2)+pow(systALICEl[i],2));
  for (Int_t i=0;i<nALICE;i++) dsigALICEh[i] = sqrt(pow(statALICE[i],2)+pow(systALICEh[i],2));
  
  if (option==1) return new TGraphAsymmErrors(3,wavgALICE,sigmALICE,NULL,NULL,dsigALICEl,dsigALICEh);
  if (option==2) return new TGraphAsymmErrors(1,&(wavgALICE[3]),&(sigmALICE[3]),NULL,NULL,&(dsigALICEl[3]),&(dsigALICEh[3]));
  return new TGraphAsymmErrors(nALICE,wavgALICE,sigmALICE,NULL,NULL,dsigALICEl,dsigALICEh);
}

TGraphErrors* read_lhcb(const char* fileName, const int nPoints){
  ifstream f;
  f.open(fileName);
  Double_t w[nPoints],cs[nPoints],dcs[nPoints];
  for (Int_t i=0;i<nPoints;i++)  f >> w[i] >> cs[i] >> dcs[i];
  f.close();
  return new TGraphErrors(nPoints,w,cs,NULL,dcs);
}

TGraphErrors* read_clark(){
  const int n=5;
  Double_t w[]   = {8.04, 10.29, 12.06, 14.3, 16.23};
  Double_t cs[]  = {11.4,  14.4,  17.8, 19.6,  22.0};
  Double_t dcs[] = { 2.2,   3.6,   3.4,  3.7,   6.3};
  return new TGraphErrors(n,w,cs,NULL,dcs);
}

TGraphErrors* read_slac(){
  const int n=5;
  Double_t w[]   = {5.03, 5.56, 5.72, 6.04, 6.35};
  Double_t cs[]  = {1.31, 2.83, 3.72, 4.14, 5.03};
  Double_t dcs[] = {0.31, 0.48, 0.52, 0.57, 0.66};
  return new TGraphErrors(n,w,cs,NULL,dcs);
}

TGraphErrors* read_gittelman(){
  const int n=1;
  Double_t w[]   = {4.64};
  Double_t cs[]  = {0.484};
  Double_t dcs[] = {0.123};
  return new TGraphErrors(n,w,cs,NULL,dcs);
}

TGraphErrors* read_e814(){
  const int n=4;
  Double_t w[]   = {11.79, 15.36, 18.11, 20.48};
  Double_t cs[]  = { 17.7,  20.9,  26.6,  34.2};
  Double_t dcs[] = {  2.6,   2.6,   4.8,   3.0};
  return new TGraphErrors(n,w,cs,NULL,dcs);
}

TGraphErrors* read_e516(){
  const int n=1;
  Double_t w[]   = {14.1};
  Double_t cs[]  = { 9.8};
  Double_t dcs[] = {2.05};
  return new TGraphErrors(n,w,cs,NULL,dcs);
}
