// contact: Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
//
// Macro to reproduce Fig. 4 of the directed flow paper http://arxiv.org/abs/1306.4145
//
#include "STAR_v1_PRL_v101_i25_e252301_y2008.C"

float rangeYaMax = 3.79e-3;
float rangeYaMin = 1.29e-3;
float rangeYbMax = 3.2e-3;
float rangeYbMin = 3.49e-3*0.45;
float maxPt =4.7;
float maxPtFit =4.5;
float minPtFit =0.;
float shiftX =0.04;
float scale62=0.12;
float scale200=0.37;

float myMarkerSize=1.5;
const int xRes=600, yRes=600;
int SystlineWidth =6;

TGraphErrors *v1_even_pT_stat_05_80;
TGraphErrors *v1_even_pT_syst_05_80;
TGraphErrors *v1_odd_pT_stat_05_80;
TGraphErrors *v1_odd_pT_syst_05_80;

TGraphErrors *v1_odd_pT_stat_05_40;
TGraphErrors *v1_odd_pT_syst_05_40;
TGraphErrors *v1_odd_pT_stat_40_80;
TGraphErrors *v1_odd_pT_syst_40_80;

float mainFont = 20;
TString PtFitFunction = "[0]*x+[1]*x^2+[2]*x^3";

const int nSystOpt = 4;

void Figure4(bool rWrite = false,TString dataFileName="ALICE_v1_arxiv_1306_4145.root")
{
  TGaxis::SetMaxDigits(3);
  myOptions();
  gROOT->ForceStyle();
  directedFlow_2008_STARdataPt(0.01*scale200, 0.01*scale200);
  
  TCanvas *myCan = new TCanvas("myCan","Figure 4",xRes,yRes);
  myCan->cd();
  TPad *myPad1 = new TPad("myPad1","myPad1",0.,0.52,1,1);
  myPadSetUp(myPad1,0.08,0.078,0.01,0.00);
  TPad *myPad2 = new TPad("myPad2","myPad2",0.,0.,1,0.52);
  myPadSetUp(myPad2,0.08,0.00,0.01,0.16);
  myPad1->Draw();
  myPad2->Draw();
  
  TH1F *myBlankHisto1 = new TH1F("myBlankHisto1","Blank Histogram",10,0,maxPt);
  myHistoSetUp(myBlankHisto1,"p_{T}, GeV/c","v_{1}",-rangeYaMin*0.8,rangeYaMax*0.91, 1, 2, 310, 503);
  TH1F *myBlankHisto2 = new TH1F("myBlankHisto2","Blank Histogram",10,0,maxPt);
  myHistoSetUp(myBlankHisto2,"p_{T} (GeV/c)","v_{1}",-rangeYaMin*1.1,rangeYaMax*0.75, 1, 2, 310, 503);
  
  // get graphs
  TFile *dataFile = TFile::Open(dataFileName,"READ");
  
  v1_even_pT_stat_05_80 = (TGraphErrors*)dataFile->Get("v1_even_pT_stat_05_80");
  v1_even_pT_syst_05_80 = (TGraphErrors*)dataFile->Get("v1_even_pT_syst_05_80");
  v1_odd_pT_stat_05_80 = (TGraphErrors*)dataFile->Get("v1_odd_pT_stat_05_80");
  v1_odd_pT_syst_05_80 = (TGraphErrors*)dataFile->Get("v1_odd_pT_syst_05_80");
  
  v1_odd_pT_stat_05_40 = (TGraphErrors*)dataFile->Get("v1_odd_pT_stat_05_40");
  v1_odd_pT_syst_05_40 = (TGraphErrors*)dataFile->Get("v1_odd_pT_syst_05_40");
  v1_odd_pT_stat_40_80 = (TGraphErrors*)dataFile->Get("v1_odd_pT_stat_40_80");
  v1_odd_pT_syst_40_80 = (TGraphErrors*)dataFile->Get("v1_odd_pT_syst_40_80");
  
  // setup graphs
  int colorStat = kBlue; int colorSyst = kBlue-10;
  myTGraphSetUp(v1_even_pT_syst_05_80,kOpenCircle,colorSyst,0,1,colorSyst, SystlineWidth,1001,colorSyst);
  myTGraphSetUp(v1_even_pT_stat_05_80,kOpenCircle,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  myTGraphSetUp(v1_odd_pT_syst_05_40,kFullTriangleDown,colorSyst,0,1,colorSyst, SystlineWidth,1001,colorSyst);
  myTGraphSetUp(v1_odd_pT_stat_05_40,kFullTriangleDown,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  colorStat = kMagenta+2; colorSyst = kMagenta-8;
  myTGraphSetUp(v1_odd_pT_syst_05_80,kFullSquare,colorSyst,0,1,colorSyst, SystlineWidth,1001,colorSyst);
  myTGraphSetUp(v1_odd_pT_stat_05_80,kFullSquare,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  myTGraphSetUp(v1_odd_pT_syst_40_80,kFullTriangleUp,colorSyst,0,1,colorSyst, SystlineWidth,1001,colorSyst);
  myTGraphSetUp(v1_odd_pT_stat_40_80,kFullTriangleUp,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  colorStat = kGreen+2; colorSyst = kGreen-8;
  myTGraphSetUp(v1_star_AuAu200_40_80_pT,kOpenCross,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  colorStat = kBlue+2; colorSyst = kBlue-8;
  myTGraphSetUp(v1_star_AuAu200_5_40_pT,kOpenStar,colorStat,myMarkerSize*1.3,1,colorStat,2,1001,colorStat);  
  
  // fitting
  TF1*fit_v1_even_pT_stat_05_80 = new TF1("fit_v1_even_pT_stat_05_80", PtFitFunction, minPtFit, maxPtFit);
  v1_even_pT_stat_05_80->Fit(fit_v1_even_pT_stat_05_80,"0","",minPtFit, maxPtFit);
  fit_v1_even_pT_stat_05_80 = (TF1*)v1_even_pT_stat_05_80->GetListOfFunctions()->At(0)->Clone();
  fit_v1_even_pT_stat_05_80->SetLineColor(v1_even_pT_stat_05_80->GetLineColor()-10);
  
  TF1*fit_v1_odd_pT_stat_05_80 = new TF1("fit_v1_odd_pT_stat_05_80", PtFitFunction, minPtFit, maxPtFit);
  v1_odd_pT_stat_05_80->Fit(fit_v1_odd_pT_stat_05_80,"0","",minPtFit, maxPtFit);
  fit_v1_odd_pT_stat_05_80 = (TF1*)v1_odd_pT_stat_05_80->GetListOfFunctions()->At(0)->Clone();
  fit_v1_odd_pT_stat_05_80->SetLineColor(v1_odd_pT_stat_05_80->GetLineColor()-10);
  
  TF1*fit_v1_odd_pT_stat_05_40 = new TF1("fit_v1_odd_pT_stat_05_40", PtFitFunction, 0,2.2);
  v1_odd_pT_stat_05_40->Fit(fit_v1_odd_pT_stat_05_40,"0","",minPtFit, maxPtFit);
  fit_v1_odd_pT_stat_05_40 = (TF1*)v1_odd_pT_stat_05_40->GetListOfFunctions()->At(0)->Clone();
  fit_v1_odd_pT_stat_05_40->SetLineColor(v1_odd_pT_stat_05_40->GetLineColor()-10);
  
  TF1*fit_v1_odd_pT_stat_40_80 = new TF1("fit_v1_odd_pT_stat_40_80", PtFitFunction, minPtFit, maxPtFit);
  v1_odd_pT_stat_40_80->Fit(fit_v1_odd_pT_stat_40_80,"0","",minPtFit, maxPtFit);
  fit_v1_odd_pT_stat_40_80 = (TF1*)v1_odd_pT_stat_40_80->GetListOfFunctions()->At(0)->Clone();
  fit_v1_odd_pT_stat_40_80->SetLineColor(v1_odd_pT_stat_40_80->GetLineColor()-10);
  
  // close graphs to fill legends
  TGraphErrors *v1_odd_pT_syst_05_80_clone=(TGraphErrors*)v1_odd_pT_syst_05_80->Clone("v1_odd_pT_syst_05_80_clone"); v1_odd_pT_syst_05_80_clone->SetLineWidth(1);
  TGraphErrors *v1_even_pT_syst_05_80_clone=(TGraphErrors*)v1_even_pT_syst_05_80->Clone("v1_even_pT_syst_05_80_clone"); v1_even_pT_syst_05_80_clone->SetLineWidth(1);
  TGraphErrors *v1_odd_pT_syst_05_40_clone=(TGraphErrors*)v1_odd_pT_syst_05_40->Clone("v1_odd_pT_syst_05_40_clone"); v1_odd_pT_syst_05_40_clone->SetLineWidth(1);
  TGraphErrors *v1_odd_pT_syst_40_80_clone=(TGraphErrors*)v1_odd_pT_syst_40_80->Clone("v1_odd_pT_syst_40_80_clone"); v1_odd_pT_syst_40_80_clone->SetLineWidth(1);
  
  //create and fill legends
  TLegend *myLegend1SysA = new TLegend(0.1,0.64,0.42,0.74);
  myLegendSetUp(myLegend1SysA,mainFont);
  myLegend1SysA->AddEntry(v1_odd_pT_syst_05_80_clone," ","F");
  TLegend *myLegend1a = new TLegend(0.1,0.64,0.42,0.74);
  myLegendSetUp(myLegend1a,mainFont);
  myLegend1a->AddEntry(v1_odd_pT_stat_05_80," ","P");
  
  TLegend *myLegend1b = new TLegend(0.2,0.64,0.52,0.74);
  myLegendSetUp(myLegend1b,mainFont);
  myLegend1b->AddEntry(v1_even_pT_stat_05_80," ","P");
  TLegend *myLegend1SysB = new TLegend(0.2,0.64,0.52,0.74);
  myLegendSetUp(myLegend1SysB,mainFont);
  myLegend1SysB->AddEntry(v1_even_pT_syst_05_80_clone,"  5-80%","F");
  
  TLegend *myLegend2Sys = new TLegend(0.1,0.45,0.42,0.61);
  myLegendSetUp(myLegend2Sys,mainFont);
  myLegend2Sys->AddEntry(v1_odd_pT_syst_05_40_clone," ","F");
  myLegend2Sys->AddEntry(v1_odd_pT_syst_40_80_clone," ","F");
  
  TLegend *myLegend2 = new TLegend(0.1,0.45,0.42,0.61);
  myLegendSetUp(myLegend2,mainFont);
  myLegend2->AddEntry(v1_odd_pT_stat_05_40,"   5-40%","P");
  myLegend2->AddEntry(v1_odd_pT_stat_40_80,"   40-80%","P");
  
  TLegend *myLegend3 = new TLegend(0.1,0.73,0.42,0.9);
  myLegendSetUp(myLegend3,mainFont);
  myLegend3->AddEntry(v1_star_AuAu200_5_40_pT,"  #times 0.37  5-40%","P");
  myLegend3->AddEntry(v1_star_AuAu200_40_80_pT,"  #times 0.37  40-80%","P");
  
  TLegend *myLegend4 = new TLegend(0.1,0.53,0.42,0.57);
  myLegendSetUp(myLegend4,mainFont);
  myLegend4->AddEntry(v1_odd_pT_stat_05_80," ","L");
  myLegend4->AddEntry("NULL"," polynomial fits","");
  myLegend4->AddEntry(v1_even_pT_stat_05_80," ","L");
  
  // pad1
  myPad1->cd();
  myBlankHisto1->Draw();
  
  TLatex *myText = new TLatex();
  myText->SetNDC();
  myText->SetTextSize(mainFont);
  myText->SetTextColor(1);
  myText->DrawLatex(0.11,0.85,"ALICE Pb-Pb@2.76TeV |#eta|<0.8");
  myText->DrawLatex(0.115,0.76,"odd     even  v_{1}");
  myText->SetTextSize(1.5*mainFont);
  myText->DrawLatex(0.93,0.83,"(a)");
  
  myLegend1SysA->Draw();
  myLegend1SysB->Draw();
  myLegend1a->Draw();
  myLegend1b->Draw();
  myLegend4->Draw();
  
  ShiftAlongXaxis(v1_odd_pT_syst_05_80, shiftX);
  ShiftAlongXaxis(v1_odd_pT_stat_05_80, shiftX);
  ShiftAlongXaxis(v1_even_pT_syst_05_80, -shiftX);
  ShiftAlongXaxis(v1_even_pT_stat_05_80, -shiftX);
  
  fit_v1_even_pT_stat_05_80->Draw("same");
  fit_v1_odd_pT_stat_05_80->Draw("same");
  
  v1_even_pT_syst_05_80->Draw("eZ");
  v1_odd_pT_syst_05_80->Draw("eZ");
  v1_even_pT_stat_05_80->Draw("P,eZ");
  v1_odd_pT_stat_05_80->Draw("P,eZ");
  
  // pad2
  myPad2->cd();
  myBlankHisto2->Draw();
  
  myText->SetTextSize(mainFont);
  myText->DrawLatex(0.095,0.63,"ALICE  odd v_{1}");
  myText->SetTextSize(1.5*mainFont);
  myText->DrawLatex(0.93,0.91,"(b)");
  myText->SetTextSize(mainFont);
  myText->DrawLatex(0.1,0.92,"STAR (scaled) Au-Au@200GeV |#eta|<1.3");
  
  myLegend2Sys->Draw();
  myLegend2->Draw();
  myLegend3->Draw();
  
  ShiftAlongXaxis(v1_odd_pT_stat_40_80, shiftX);
  ShiftAlongXaxis(v1_odd_pT_syst_40_80, shiftX);
  ShiftAlongXaxis(v1_odd_pT_stat_05_40, -shiftX);
  ShiftAlongXaxis(v1_odd_pT_syst_05_40, -shiftX);
  
  fit_v1_odd_pT_stat_05_40->Draw("same");
  fit_v1_odd_pT_stat_40_80->Draw("same");
  
  v1_odd_pT_syst_05_40->Draw("eZ");
  v1_odd_pT_syst_40_80->Draw("eZ");
  v1_star_AuAu200_5_40_pT->Draw("P,eZ");
  v1_star_AuAu200_40_80_pT->Draw("P,eZ");
  v1_odd_pT_stat_05_40->Draw("P,eZ");
  v1_odd_pT_stat_40_80->Draw("P,eZ");
  
  
  // save output if option rWrite=true
  TString fileName="Figure4";
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
  gStyle->SetTitleOffset(0.8,"y");  
  gStyle->SetTitleOffset(2,"xz");  
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
