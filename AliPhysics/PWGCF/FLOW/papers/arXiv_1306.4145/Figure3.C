// contact: Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
//
// Macro to reproduce Fig. 3 of the directed flow paper http://arxiv.org/abs/1306.4145
//
#include "STAR_v1_PRL_v101_i25_e252301_y2008.C"

float myMarkerSize=1.5;
const int xRes=600, yRes=800;
float mainFont = 20;
const int nSystOpt = 4;
float rangeYa = 0.99e-3;
float scale62=0.12;
float scale200=0.37;

TGraphErrors *v1_even_cen_stat_05_80;
TGraphErrors *v1_even_cen_syst_05_80;
TGraphErrors *v1_odd_cen_stat_05_80;
TGraphErrors *v1_odd_cen_syst_05_80;

void Figure3(bool rWrite = false,TString dataFileName="ALICE_v1_arxiv_1306_4145.root")
{
  TGaxis::SetMaxDigits(3);
  myOptions();
  gROOT->ForceStyle();
  directedFlow_2008_STARdataCentrality(0.01*scale200, 0.01*scale62);
  
  float minPad[3]={0,0.36,0.665};
  TCanvas *myCan = new TCanvas("myCan","Figure 3",xRes,yRes);
  myCan->cd();
  TPad *myPad1 = new TPad("myPad1","myPad1",0.,minPad[2],1,1);
  myPadSetUp(myPad1,0.08, 0.09, 0.01, 0.00);
  TPad *myPad2 = new TPad("myPad2","myPad2",0.,minPad[1],1,minPad[2]);
  myPadSetUp(myPad2,0.08, 0.0, 0.01, 0.0);
  TPad *myPad3 = new TPad("myPad3","myPad3",0.,minPad[0],1,minPad[1]);
  myPadSetUp(myPad3,0.08, 0.00, 0.01, 0.16);
  myPad1->Draw();
  myPad2->Draw();
  myPad3->Draw();
  
  TH1F *myBlankHisto1 = new TH1F("myBlankHisto1","Blank Histogram",10,0,82);
  myHistoSetUp(myBlankHisto1,"centrality (%)","v_{1}",-rangeYa,rangeYa*0.25, 1, 2, 115, 305);
  TH1F *myBlankHisto2 = new TH1F("myBlankHisto2","Blank Histogram",10,0,82);
  myHistoSetUp(myBlankHisto2,"centrality (%)","#LTp_{x}#GT/#LTp_{T}#GT",-rangeYa,rangeYa*0.25, 1, 2, 115, 305);
  TH1F *myBlankHisto3 = new TH1F("myBlankHisto3","Blank Histogram",0,0,82);
  myHistoSetUp(myBlankHisto3,"centrality (%)","v_{1}",-rangeYa,rangeYa*0.25, 1, 2, 115, 305);
  
  // get graphs
  TFile *dataFile = TFile::Open(dataFileName,"READ");
  
  v1_even_cen_stat_05_80 = (TGraphErrors*)dataFile->Get("v1_even_cen_stat_05_80");
  v1_even_cen_syst_05_80 = (TGraphErrors*)dataFile->Get("v1_even_cen_syst_05_80");
  v1_odd_cen_stat_05_80 = (TGraphErrors*)dataFile->Get("v1_odd_cen_stat_05_80");
  v1_odd_cen_syst_05_80 = (TGraphErrors*)dataFile->Get("v1_odd_cen_syst_05_80");
  px_even_cen_stat_05_80 = (TGraphErrors*)dataFile->Get("px_even_cen_stat_05_80");
  px_even_cen_syst_05_80 = (TGraphErrors*)dataFile->Get("px_even_cen_syst_05_80");
  px_odd_cen_stat_05_80 = (TGraphErrors*)dataFile->Get("px_odd_cen_stat_05_80");
  px_odd_cen_syst_05_80 = (TGraphErrors*)dataFile->Get("px_odd_cen_syst_05_80");
  
  // setup graphs
  int SystlineWidth =6;
  int colorStat = kBlue; int colorSyst = colorStat-10;
  myTGraphSetUp(v1_even_cen_syst_05_80,kOpenCircle,colorSyst,0,1,colorSyst,SystlineWidth,1001,colorSyst);
  myTGraphSetUp(v1_even_cen_stat_05_80,kOpenCircle,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  colorStat = kMagenta+2; colorSyst = colorStat-10;
  myTGraphSetUp(v1_odd_cen_syst_05_80,kFullDiamond,colorSyst,0,1,colorSyst,SystlineWidth,1001,colorSyst);
  myTGraphSetUp(v1_odd_cen_stat_05_80,kFullDiamond,colorStat,1.7*myMarkerSize,1,colorStat,2,1001,colorStat);
  
  colorStat = kBlue; colorSyst = colorStat-10;
  myTGraphSetUp(px_even_cen_syst_05_80,kOpenSquare,colorSyst,0,1,colorSyst,SystlineWidth,1001,colorSyst);
  myTGraphSetUp(px_even_cen_stat_05_80,kOpenSquare,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  colorStat = kMagenta+2; colorSyst = colorStat-10;
  myTGraphSetUp(px_odd_cen_syst_05_80,kFullSquare,colorSyst,0,1,colorSyst,SystlineWidth,1001,colorSyst);
  myTGraphSetUp(px_odd_cen_stat_05_80,kFullSquare,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  myTGraphSetUp(v1_star_AuAu200_cent,kOpenCross,kGreen+2,myMarkerSize,1,kGreen+2,2,1001,kGreen+2);
  myTGraphSetUp(v1_star_AuAu62_cent,kOpenStar,kBlue+2,myMarkerSize*1.3,1,kBlue+2,2,1001,kBlue+2);  
  
  
  TGraphErrors *v1_odd_cen_syst_05_80_clone=(TGraphErrors*)v1_odd_cen_syst_05_80->Clone("v1_odd_cen_syst_05_80_clone"); v1_odd_cen_syst_05_80_clone->SetLineWidth(1);
  TGraphErrors *v1_even_cen_syst_05_80_clone=(TGraphErrors*)v1_even_cen_syst_05_80->Clone("v1_even_cen_syst_05_80_clone"); v1_even_cen_syst_05_80_clone->SetLineWidth(1);
  
  TLegend *myLegend1SysA = new TLegend(0.11,0.09,0.47,0.2);
  myLegendSetUp(myLegend1SysA,mainFont);
  myLegend1SysA->AddEntry(v1_odd_cen_syst_05_80_clone," ","F");
  TLegend *myLegend1a = new TLegend(0.11,0.09,0.47,0.2);
  myLegendSetUp(myLegend1a,mainFont);
  myLegend1a->AddEntry(v1_odd_cen_stat_05_80," ","P");
  
  TLegend *myLegend1b = new TLegend(0.215,0.09,0.575,0.2);
  myLegendSetUp(myLegend1b,mainFont);
  myLegend1b->AddEntry(v1_even_cen_stat_05_80," ","P");
  TLegend *myLegend1SysB = new TLegend(0.215,0.09,0.575,0.2);
  myLegendSetUp(myLegend1SysB,mainFont);
  myLegend1SysB->AddEntry(v1_even_cen_syst_05_80_clone," v_{1}","F");
  
  TGraphErrors *px_odd_cen_syst_05_80_clone=(TGraphErrors*)px_odd_cen_syst_05_80->Clone("px_odd_cen_syst_05_80_clone"); px_odd_cen_syst_05_80_clone->SetLineWidth(1);
  TGraphErrors *px_even_cen_syst_05_80_clone=(TGraphErrors*)px_even_cen_syst_05_80->Clone("px_even_cen_syst_05_80_clone"); px_even_cen_syst_05_80_clone->SetLineWidth(1);
  
  TLegend *myLegend1pTSysA = new TLegend(0.11,0.09,0.47,0.2);
  myLegendSetUp(myLegend1pTSysA,mainFont);
  myLegend1pTSysA->AddEntry(px_odd_cen_syst_05_80_clone," ","F");
  TLegend *myLegend1pTa = new TLegend(0.11,0.09,0.47,0.2);
  myLegendSetUp(myLegend1pTa,mainFont);
  myLegend1pTa->AddEntry(px_odd_cen_stat_05_80," ","P");
  
  TLegend *myLegend1pTb = new TLegend(0.215,0.09,0.575,0.2);
  myLegendSetUp(myLegend1pTb,mainFont);
  myLegend1pTb->AddEntry(px_even_cen_stat_05_80," ","P");
  TLegend *myLegend1pTSysB = new TLegend(0.215,0.09,0.575,0.2);
  myLegendSetUp(myLegend1pTSysB,mainFont);
  myLegend1pTSysB->AddEntry(px_even_cen_syst_05_80_clone," #LTp_{x}#GT/#LTp_{T}#GT","F");
  
  TLegend *myLegend3 = new TLegend(0.11,0.2,0.47,0.38);
  myLegendSetUp(myLegend3,mainFont);
  myLegend3->AddEntry(v1_star_AuAu200_cent," #times 0.37  Au-Au@200GeV","P");
  myLegend3->AddEntry(v1_star_AuAu62_cent," #times 0.12  Au-Au@62GeV","P");
  
  // pad1
  myPad1->cd();  
  myBlankHisto1->Draw();
  
  TLatex *myText = new TLatex();
  myText->SetNDC();
  myText->SetTextSize(mainFont);
  myText->SetTextColor(1);
  myText->DrawLatex(0.1,0.81,"ALICE Pb-Pb@2.76TeV  |#eta|<0.8 p_{T}>0.15 GeV/c");
  myText->DrawLatex(0.13,0.21,"odd     even");
  myText->SetTextSize(mainFont*1.5);
  myText->DrawLatex(0.93,0.81,"(a)");
  
  myLegend1SysA->Draw();
  myLegend1SysB->Draw();
  myLegend1a->Draw();
  myLegend1b->Draw();
  
  float shift =-1.2;
  ShiftAlongXaxis(v1_even_cen_syst_05_80, shift);
  ShiftAlongXaxis(v1_even_cen_stat_05_80, shift);
  ShiftAlongXaxis(px_even_cen_syst_05_80, shift);
  ShiftAlongXaxis(px_even_cen_stat_05_80, shift);
  
  v1_even_cen_syst_05_80->Draw("P,eZ");  
  v1_odd_cen_syst_05_80->Draw("eZ");  
  v1_even_cen_stat_05_80->Draw("P,eZ");
  v1_odd_cen_stat_05_80->Draw("P,eZ");
  
  // pad2
  myPad2->cd();
  myBlankHisto2->Draw();
  
  TLatex *myText2 = new TLatex();
  myText2->SetNDC();
  myText2->SetTextColor(1);
  myText2->SetTextSize(mainFont);
  myText2->DrawLatex(0.13,0.21,"odd     even");
  myText2->SetTextSize(mainFont*1.5);
  myText2->DrawLatex(0.93,0.88,"(b)");
  
  myLegend1pTSysA->Draw();
  myLegend1pTSysB->Draw();
  myLegend1pTa->Draw();
  myLegend1pTb->Draw();
  
  px_even_cen_syst_05_80->Draw("P,eZ");
  px_odd_cen_syst_05_80->Draw("P,eZ");
  px_even_cen_stat_05_80->Draw("P,eZ");
  px_odd_cen_stat_05_80->Draw("P,eZ");
  
  // pad3
  myPad3->cd();
  myBlankHisto3->Draw();
  
  TLatex *myText3 = new TLatex();
  myText3->SetNDC();
  myText3->SetTextColor(1);
  myText3->SetTextSize(mainFont*1.5);
  myText3->DrawLatex(0.93,0.9,"(c)");
  myText->SetTextSize(mainFont);
  myText->DrawLatex(0.11,0.4,"STAR (scaled) |#eta|<1.3");
  
  myLegend3->Draw();
  
  ShiftAlongXaxis(v1_star_AuAu62_cent, -1.5);
  ShiftAlongXaxis(v1_star_AuAu200_cent, 1.5);
  
  v1_star_AuAu62_cent->GetXaxis()->SetLimits(10,80);
  v1_odd_cen_syst_05_80->Draw("eZ");
  v1_star_AuAu62_cent->Draw("P,eZ");
  v1_star_AuAu200_cent->Draw("P,eZ");
  v1_odd_cen_stat_05_80->Draw("P,eZ");
  
  // save output if option rWrite=true
  TString fileName="Figure3";
  myCan->Update();
  if (rWrite)  
  {
    myCan->SaveAs(fileName+".png");
    myCan->SaveAs(fileName+".eps");
    myCan->SaveAs(fileName+".pdf");
  }
}

// only helper functions below

void myHistoSetUp
(
  TH1F *hist=0, TString xTitle="xTitle", TString yTitle="yTitle", float minY=-1, float maxY=1, int lineColor=1, int lineStyle=2,
 int nDivisionsX=305,int nDivisionsY=305
)
{
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle(yTitle);
  hist->SetMinimum(minY);
  hist->SetMaximum(maxY);
  hist->SetLineColor(lineColor);
  hist->SetLineStyle(lineStyle);
  hist->SetNdivisions(nDivisionsX,"x");
  hist->SetNdivisions(nDivisionsY,"y");
  return;
}


void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}

void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.0,
		  int currentMarkerStyle=20, int currentMarkerColor=0,
		  int currentLineStyle=1, int currentLineColor=0)
{
  currentGraph->SetMarkerSize(currentMarkerSize);
  currentGraph->SetMarkerStyle(currentMarkerStyle);
  currentGraph->SetMarkerColor(currentMarkerColor);
  currentGraph->SetLineStyle(currentLineStyle);
  currentGraph->SetLineColor(currentLineColor);
  return;
}

void myOptions(Int_t lStat=0){
  int font = 43;
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetLegendFont(font);
  gStyle->SetStatFontSize(20);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(20,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");
  gStyle->SetTitleOffset(1.2,"y");  
  gStyle->SetTitleOffset(3,"xz");  
  gStyle->SetTitleSize(21,"x");  
  gStyle->SetTitleSize(24,"y");  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0); 
  
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
  }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}

TGraphErrors* makeGraphH1(TH1* hist, TString name="")
{
  name.ReplaceAll(" ","");
  Int_t nbins = hist->GetNbinsX();
  Double_t* x = new Double_t[nbins];
  Double_t* y = new Double_t[nbins];
  Double_t* xerr = new Double_t[nbins]; 
  Double_t* yerr = new Double_t[nbins];
  Int_t n=0;
  for (Int_t i=0; i<nbins; i++)
  {
    x[n] = hist->GetXaxis()->GetBinCenter(i+1);
    y[n] = hist->GetBinContent(i+1);
    xerr[n] = 0.0;
    yerr[n] = hist->GetBinError(i+1);
    n++;
  }
  TGraphErrors* gr = new TGraphErrors(n,x,y,xerr,yerr);
  delete [] x;
  delete [] y;
  delete [] xerr;
  delete [] yerr;
  return gr;
}

TGraphErrors* makeGraphPr(TProfile* hist, TString name="")
{
  name.ReplaceAll(" ","");
  Int_t nbins = hist->GetNbinsX();
  Double_t* x = new Double_t[nbins];
  Double_t* y = new Double_t[nbins];
  Double_t* xerr = new Double_t[nbins]; 
  Double_t* yerr = new Double_t[nbins];
  Int_t n=0;
  for (Int_t i=0; i<nbins; i++)
  {
    x[n] = hist->GetXaxis()->GetBinCenter(i+1);
    y[n] = hist->GetBinContent(i+1);
    xerr[n] = 0.0;
    yerr[n] = hist->GetBinError(i+1);
    n++;
  }
  TGraphErrors* gr = new TGraphErrors(n,x,y,xerr,yerr);
  delete [] x;
  delete [] y;
  delete [] xerr;
  delete [] yerr;
  return gr;
}

void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)
{
  if(!ge)
  {
    printf("\n WARNING: ge is NULL in ShiftAlongXaxis() !!!! \n\n");
    return;
  }
  Int_t nPoints = ge->GetN();
  Double_t x = 0.;
  Double_t y = 0.;
  for(Int_t p=0;p<nPoints;p++)
  {
    ge->GetPoint(p,x,y);
    x+=shift;
    ge->SetPoint(p,x,y);
  }
}

void myTGraphSetUp
(
  TGraphErrors *currentGraph=0,
 int myMarkerStyle=8,
 int myMarkerColor=1,
 float myMarkerSize=1,
 int myLineStyle=1,
 int myLineColor=1,
 float myLineWidth=1,
 int myFillStyle =1001,
 int myFillColor =1 
)
{
  currentGraph->SetMarkerStyle(myMarkerStyle);
  currentGraph->SetMarkerColor(myMarkerColor);
  currentGraph->SetMarkerSize(myMarkerSize);
  currentGraph->SetLineColor(myLineColor);
  currentGraph->SetLineStyle(myLineStyle);
  currentGraph->SetLineWidth(myLineWidth);
  currentGraph->SetFillStyle(myFillStyle);
  currentGraph->SetFillColor(myFillColor);
}
