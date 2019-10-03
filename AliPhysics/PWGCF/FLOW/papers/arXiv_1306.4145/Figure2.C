// contact: Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
//
// Macro to reproduce Fig. 2 of the directed flow paper http://arxiv.org/abs/1306.4145
//
#include "STAR_v1_PRL_v101_i25_e252301_y2008.C"

float myMarkerSize=1.8;
const int xRes=600, yRes=800;
int fillStyle = 1001;
float mainFont = 20;
const int nSystOpt = 4;
float rangeYa = 0.99e-3;
float rangeYb = 0.99e-3;
float scale62=0.12;
float scale200=0.37;
int SystlineWidth =6;

TGraphErrors *v1_odd_eta_syst_10_20;
TGraphErrors *v1_odd_eta_syst_30_40;
TGraphErrors *v1_odd_eta_stat_10_20;
TGraphErrors *v1_odd_eta_stat_30_40;

TGraphErrors *px_even_eta_syst_10_60;
TGraphErrors *px_even_eta_stat_10_60;
TGraphErrors *px_odd_eta_syst_10_60;
TGraphErrors *px_odd_eta_stat_10_60;

TGraphErrors *v1_even_eta_syst_10_60;
TGraphErrors *v1_even_eta_stat_10_60;
TGraphErrors *v1_odd_eta_syst_10_60;
TGraphErrors *v1_odd_eta_stat_10_60;

void Figure2(bool rWrite = false,TString dataFileName="ALICE_v1_arxiv_1306_4145.root")
{
  TGaxis::SetMaxDigits(3);
  myOptions();
  gROOT->ForceStyle();
  directedFlow_2008_STARdataEta(0.01*scale200, 0.01*scale62);
  
  float minPad[3]={0,0.36,0.665};
  TCanvas *myCan = new TCanvas("myCan","Figure 3",xRes,yRes);
  myCan->cd();
  TPad *myPad1 = new TPad("myPad1","myPad1",0.,minPad[2],1,1);
  myPadSetUp(myPad1,0.11, 0.09, 0.01, 0.00);
  TPad *myPad2 = new TPad("myPad2","myPad2",0.,minPad[1],1,minPad[2]);
  myPadSetUp(myPad2,0.11, 0.0, 0.01, 0.0);
  TPad *myPad3 = new TPad("myPad3","myPad3",0.,minPad[0],1,minPad[1]);
  myPadSetUp(myPad3,0.11, 0.00, 0.01, 0.16);
  myPad1->Draw();
  myPad2->Draw();
  myPad3->Draw();
  
  TH1F *myBlankHisto1 = new TH1F("myBlankHisto1","Blank Histogram",10,-0.85,0.85);
  myHistoSetUp(myBlankHisto1,"#eta","v_{1}",-rangeYb,rangeYb, 1, 2, 505, 305);
  TH1F *myBlankHisto2 = new TH1F("myBlankHisto2","Blank Histogram",10,-0.85,0.85);
  myHistoSetUp(myBlankHisto2,"#eta","#LTp_{x}#GT/#LTp_{T}#GT",-rangeYa,rangeYa, 1, 2, 505, 305);
  TH1F *myBlankHisto3 = new TH1F("myBlankHisto3","Blank Histogram",100,-0.85,0.85);
  myHistoSetUp(myBlankHisto3,"#eta","v_{1}",-rangeYb,rangeYb, 1, 2, 505, 305);
  
  // get graphs
  TFile *dataFile = TFile::Open(dataFileName,"READ");
  
  v1_odd_eta_stat_10_20 = (TGraphErrors*)dataFile->Get("v1_odd_eta_stat_10_20");
  v1_odd_eta_syst_10_20 = (TGraphErrors*)dataFile->Get("v1_odd_eta_syst_10_20");
  v1_odd_eta_stat_30_40 = (TGraphErrors*)dataFile->Get("v1_odd_eta_stat_30_40");
  v1_odd_eta_syst_30_40 = (TGraphErrors*)dataFile->Get("v1_odd_eta_syst_30_40");
  v1_odd_eta_stat_30_60 = (TGraphErrors*)dataFile->Get("v1_odd_eta_stat_30_60");
  v1_odd_eta_syst_30_60 = (TGraphErrors*)dataFile->Get("v1_odd_eta_syst_30_60");
  
  v1_even_eta_stat_10_20 = (TGraphErrors*)dataFile->Get("v1_even_eta_stat_10_20");
  v1_even_eta_syst_10_20 = (TGraphErrors*)dataFile->Get("v1_even_eta_syst_10_20");
  v1_even_eta_stat_30_40 = (TGraphErrors*)dataFile->Get("v1_even_eta_stat_30_40");
  v1_even_eta_syst_30_40 = (TGraphErrors*)dataFile->Get("v1_even_eta_syst_30_40");
  
  v1_odd_eta_stat_10_60 = (TGraphErrors*)dataFile->Get("v1_odd_eta_stat_10_60");
  v1_odd_eta_syst_10_60 = (TGraphErrors*)dataFile->Get("v1_odd_eta_syst_10_60");
  v1_even_eta_stat_10_60 = (TGraphErrors*)dataFile->Get("v1_even_eta_stat_10_60");
  v1_even_eta_syst_10_60 = (TGraphErrors*)dataFile->Get("v1_even_eta_syst_10_60");
  
  px_odd_eta_stat_10_60 = (TGraphErrors*)dataFile->Get("px_odd_eta_stat_10_60");
  px_odd_eta_syst_10_60 = (TGraphErrors*)dataFile->Get("px_odd_eta_syst_10_60");
  px_even_eta_stat_10_60 = (TGraphErrors*)dataFile->Get("px_even_eta_stat_10_60");
  px_even_eta_syst_10_60 = (TGraphErrors*)dataFile->Get("px_even_eta_syst_10_60");
  
  // setup graphs
  int colorStat = kBlue; int colorSyst = kBlue-10;
  
  myTGraphSetUp(v1_odd_eta_syst_30_40,kFullTriangleDown,colorSyst,0,1,colorSyst,SystlineWidth,fillStyle,colorSyst);
  myTGraphSetUp(v1_odd_eta_stat_30_40,kFullTriangleDown,colorStat,myMarkerSize,1,colorStat,2,fillStyle,colorStat);
  
  myTGraphSetUp(v1_even_eta_syst_30_40,kOpenTriangleDown,colorSyst,0,1,colorSyst,SystlineWidth,fillStyle,colorSyst);
  myTGraphSetUp(v1_even_eta_stat_30_40,kOpenTriangleDown,colorStat,myMarkerSize,1,colorStat,2,fillStyle,colorStat);
  
  myTGraphSetUp(px_even_eta_syst_10_60,kOpenSquare,colorSyst,0,1,colorSyst,SystlineWidth,1001,colorSyst);
  myTGraphSetUp(px_even_eta_stat_10_60,kOpenSquare,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  colorStat = kGreen+2; colorSyst = kGreen-8;
  
  myTGraphSetUp(v1_odd_eta_syst_10_20,kFullTriangleUp,colorSyst,0,1,colorSyst,SystlineWidth,fillStyle,colorSyst);
  myTGraphSetUp(v1_odd_eta_stat_10_20,kFullTriangleUp,colorStat,myMarkerSize,1,colorStat,2,fillStyle,colorStat);
  
  myTGraphSetUp(v1_even_eta_syst_10_20,kOpenTriangleUp,colorSyst,0,1,colorSyst,SystlineWidth,fillStyle,colorSyst);
  myTGraphSetUp(v1_even_eta_stat_10_20,kOpenTriangleUp,colorStat,myMarkerSize,1,colorStat,2,fillStyle,colorStat);
  
  colorStat = kMagenta+2; colorSyst = kMagenta-8;
  
  myTGraphSetUp(v1_even_eta_syst_10_60,kOpenDiamond,colorSyst,0,1,colorSyst,SystlineWidth,1001,colorSyst);
  myTGraphSetUp(v1_even_eta_stat_10_60,kOpenDiamond,colorStat,myMarkerSize*1.3,1,colorStat,2,1001,colorStat);
  
  myTGraphSetUp(px_odd_eta_syst_10_60,kFullSquare,colorSyst,0,1,colorSyst,SystlineWidth,1001,colorSyst);
  myTGraphSetUp(px_odd_eta_stat_10_60,kFullSquare,colorStat,myMarkerSize,1,colorStat,2,1001,colorStat);
  
  myTGraphSetUp(v1_odd_eta_syst_10_60,kFullDiamond,colorSyst,0,1,colorSyst,SystlineWidth,1001,colorSyst);
  myTGraphSetUp(v1_odd_eta_stat_10_60,kFullDiamond,colorStat,myMarkerSize*1.4,1,colorStat,2,1001,colorStat);
  
  myTGraphSetUp(v1_odd_eta_syst_30_60,kFullCircle,colorSyst,0,1,colorSyst,SystlineWidth,fillStyle,colorSyst);
  myTGraphSetUp(v1_odd_eta_stat_30_60,kFullCircle,colorStat,myMarkerSize,1,colorStat,2,fillStyle,colorStat);
  
  myTGraphSetUp(v1_star_AuAu200_30_60_eta,kOpenCross,kGreen+2,myMarkerSize,1,kGreen+2,2,fillStyle,kGreen+2);
  myTGraphSetUp(v1_star_AuAu62_30_60_eta,kOpenStar,kBlue+2,myMarkerSize*1.3,1,kBlue+2,2,fillStyle,kBlue+2);  
  
  TLatex *myText = new TLatex();
  myText->SetNDC();
  
  // fitting
  TF1*fit_px_odd_eta_stat_10_60 = new TF1("fit_v1eta_pTweight_odd_Stat_10_60", "[0]*x", -0.8,0.8);
  px_odd_eta_stat_10_60->Fit(fit_px_odd_eta_stat_10_60,"0","",-0.78,0.78);
  fit_px_odd_eta_stat_10_60 = (TF1*)px_odd_eta_stat_10_60->GetListOfFunctions()->At(0)->Clone();
  fit_px_odd_eta_stat_10_60->SetLineColor(px_odd_eta_stat_10_60->GetLineColor());
  fit_px_odd_eta_stat_10_60->SetLineStyle(1);
  fit_px_odd_eta_stat_10_60->SetLineWidth(2);
  
  TF1*fit_px_even_eta_stat_10_60 = new TF1("fit_v1eta_pTweight_even_Stat_10_60", "[0]", -0.8,0.8);
  px_even_eta_stat_10_60->Fit(fit_px_even_eta_stat_10_60,"0","",-0.78,0.78);
  fit_px_even_eta_stat_10_60 = (TF1*)px_even_eta_stat_10_60->GetListOfFunctions()->At(0)->Clone();
  fit_px_even_eta_stat_10_60->SetLineColor(px_even_eta_stat_10_60->GetLineColor());
  fit_px_even_eta_stat_10_60->SetLineStyle(2);
  fit_px_even_eta_stat_10_60->SetLineWidth(2);
  
  TF1*fit_v1_odd_eta_stat_10_60 = new TF1("fit_v1eta_odd_Stat_10_60", "[0]*x", -0.8,0.8);
  v1_odd_eta_stat_10_60->Fit(fit_v1_odd_eta_stat_10_60,"0","",-0.78,0.78);
  fit_v1_odd_eta_stat_10_60 = (TF1*)v1_odd_eta_stat_10_60->GetListOfFunctions()->At(0)->Clone();
  fit_v1_odd_eta_stat_10_60->SetLineColor(v1_odd_eta_stat_10_60->GetLineColor());
  fit_v1_odd_eta_stat_10_60->SetLineStyle(1);
  fit_v1_odd_eta_stat_10_60->SetLineWidth(2);
  
  TF1*fit_v1_even_eta_stat_10_60 = new TF1("fit_v1eta_even_Stat_10_60", "[0]", -0.8,0.8);
  v1_even_eta_stat_10_60->Fit(fit_v1_even_eta_stat_10_60,"0","",-0.78,0.78);
  fit_v1_even_eta_stat_10_60 = (TF1*)v1_even_eta_stat_10_60->GetListOfFunctions()->At(0)->Clone();
  fit_v1_even_eta_stat_10_60->SetLineColor(v1_even_eta_stat_10_60->GetLineColor());
  fit_v1_even_eta_stat_10_60->SetLineStyle(2);
  fit_v1_even_eta_stat_10_60->SetLineWidth(2);
  
  TF1*fit_v1_odd_eta_stat_30_60 = new TF1("fit_v1_odd_eta_stat_30_60", "[0]*x", -0.78,0.78);
  v1_odd_eta_stat_30_60->Fit(fit_v1_odd_eta_stat_30_60,"0","",-0.78,0.78);
  fit_v1_odd_eta_stat_30_60 = (TF1*)v1_odd_eta_stat_30_60->GetListOfFunctions()->At(0)->Clone();
  fit_v1_odd_eta_stat_30_60->SetLineColor(v1_odd_eta_stat_30_60->GetLineColor());
  fit_v1_odd_eta_stat_30_60->SetLineStyle(1);
  fit_v1_odd_eta_stat_30_60->SetLineWidth(2);
  
  // clone histos for legend drawing
  TGraphErrors *px_odd_eta_syst_10_60_clone=(TGraphErrors*)px_odd_eta_syst_10_60->Clone("px_odd_eta_syst_10_60_clone"); px_odd_eta_syst_10_60_clone->SetLineWidth(1);
  TGraphErrors *px_even_eta_syst_10_60_clone=(TGraphErrors*)px_even_eta_syst_10_60->Clone("px_even_eta_syst_10_60_clone"); px_even_eta_syst_10_60_clone->SetLineWidth(1);
  TGraphErrors *v1_odd_eta_syst_10_60_clone=(TGraphErrors*)v1_odd_eta_syst_10_60->Clone("v1_odd_eta_syst_10_60_clone"); v1_odd_eta_syst_10_60_clone->SetLineWidth(1);
  TGraphErrors *v1_odd_eta_syst_10_20_clone=(TGraphErrors*)v1_odd_eta_syst_10_20->Clone("v1_odd_eta_syst_10_20_clone"); v1_odd_eta_syst_10_20_clone->SetLineWidth(1);
  TGraphErrors *v1_odd_eta_syst_30_40_clone=(TGraphErrors*)v1_odd_eta_syst_30_40->Clone("v1_odd_eta_syst_30_40_clone"); v1_odd_eta_syst_30_40_clone->SetLineWidth(1);
  TGraphErrors *v1_even_eta_syst_10_60_clone=(TGraphErrors*)v1_even_eta_syst_10_60->Clone("v1_even_eta_syst_10_60_clone"); v1_even_eta_syst_10_60_clone->SetLineWidth(1);
  TGraphErrors *v1_even_eta_syst_10_20_clone=(TGraphErrors*)v1_even_eta_syst_10_20->Clone("v1_even_eta_syst_10_20_clone"); v1_even_eta_syst_10_20_clone->SetLineWidth(1);
  TGraphErrors *v1_even_eta_syst_30_40_clone=(TGraphErrors*)v1_even_eta_syst_30_40->Clone("v1_even_eta_syst_30_40_clone"); v1_even_eta_syst_30_40_clone->SetLineWidth(1);
  TGraphErrors *v1_odd_eta_syst_30_60_clone=(TGraphErrors*)v1_odd_eta_syst_30_60->Clone("v1_odd_eta_syst_30_60_clone"); v1_odd_eta_syst_30_60_clone->SetLineWidth(1);  
  
  // create and fill legends
  TLegend *myLegend1SysA = new TLegend(0.59,0.52,0.86,0.82);
  myLegendSetUp(myLegend1SysA,mainFont);
  myLegend1SysA->AddEntry(v1_odd_eta_syst_10_20_clone," ","F");
  myLegend1SysA->AddEntry(v1_odd_eta_syst_30_40_clone," ","F");
  myLegend1SysA->AddEntry(v1_odd_eta_syst_10_60_clone," ","F");
  
  TLegend *myLegend1a = new TLegend(0.59,0.52,0.86,0.82);
  myLegendSetUp(myLegend1a,mainFont);
  myLegend1a->AddEntry(v1_odd_eta_stat_10_20," ","P");
  myLegend1a->AddEntry(v1_odd_eta_stat_30_40," ","P");
  myLegend1a->AddEntry(v1_odd_eta_stat_10_60," ","P");
  
  TLegend *myLegend1Sys1C = new TLegend(0.59,0.52,0.86,0.82);
  myLegendSetUp(myLegend1Sys1C,mainFont);
  myLegend1Sys1C->AddEntry("NULL"," ","");
  myLegend1Sys1C->AddEntry("NULL"," ","");
  myLegend1Sys1C->AddEntry(fit_v1_odd_eta_stat_10_60," ","L");
  
  TLegend *myLegend1b = new TLegend(0.69,0.52,0.95,0.82);
  myLegendSetUp(myLegend1b,mainFont);
  myLegend1b->AddEntry(v1_even_eta_syst_10_20_clone," ","F");
  myLegend1b->AddEntry(v1_even_eta_syst_30_40_clone," ","F");
  myLegend1b->AddEntry(v1_even_eta_syst_10_60_clone," ","F");
  
  TLegend *myLegend1SysB = new TLegend(0.69,0.52,0.95,0.82);
  myLegendSetUp(myLegend1SysB,mainFont);
  myLegend1SysB->AddEntry(v1_even_eta_stat_10_20,"  10-20%","P");
  myLegend1SysB->AddEntry(v1_even_eta_stat_30_40,"  30-40%","P");
  myLegend1SysB->AddEntry(v1_even_eta_stat_10_60,"  10-60% with fit","P");
  
  TLegend *myLegend1SysC = new TLegend(0.69,0.52,0.95,0.82);
  myLegendSetUp(myLegend1SysC,mainFont);
  myLegend1SysC->AddEntry("NULL"," ","");
  myLegend1SysC->AddEntry("NULL"," ","");
  myLegend1SysC->AddEntry(fit_v1_even_eta_stat_10_60," ","L");
  
  TLegend *myLegend1SysA2 = new TLegend(0.55,0.73,0.82,0.85);
  myLegendSetUp(myLegend1SysA2,mainFont);
  myLegend1SysA2->AddEntry(px_odd_eta_syst_10_60_clone," ","F");
  
  TLegend *myLegend1a2 = new TLegend(0.55,0.73,0.82,0.85);
  myLegendSetUp(myLegend1a2,mainFont);
  myLegend1a2->AddEntry(px_odd_eta_stat_10_60," ","P");
  
  TLegend *myLegend1c2 = new TLegend(0.55,0.73,0.82,0.85);
  myLegendSetUp(myLegend1c2,mainFont);
  myLegend1c2->AddEntry(fit_v1_odd_eta_stat_10_60," ","L");
  
  TLegend *myLegend1b2 = new TLegend(0.64,0.73,0.91,0.85);
  myLegendSetUp(myLegend1b2,mainFont);
  myLegend1b2->AddEntry(px_even_eta_syst_10_60_clone," ","F");
  
  TLegend *myLegend1SysB2 = new TLegend(0.64,0.73,0.91,0.85);
  myLegendSetUp(myLegend1SysB2,mainFont);
  myLegend1SysB2->AddEntry(px_even_eta_stat_10_60,"   10-60% with fit","P");
  
  TLegend *myLegend1SysC2 = new TLegend(0.64,0.73,0.91,0.85);
  myLegendSetUp(myLegend1SysC2,mainFont);
  myLegend1SysC2->AddEntry(fit_px_even_eta_stat_10_60," ","L");
  
  TLegend *myLegend1SysA3 = new TLegend(0.65,0.73,0.91,0.85);
  myLegendSetUp(myLegend1SysA3,mainFont);
  myLegend1SysA3->AddEntry(v1_odd_eta_syst_30_60_clone," ","F");
  
  TLegend *myLegend1a3 = new TLegend(0.65,0.73,0.91,0.85);
  myLegendSetUp(myLegend1a3,mainFont);
  myLegend1a3->AddEntry(v1_odd_eta_stat_30_60," ","P");
  
  TLegend *myLegend1c3 = new TLegend(0.65,0.73,0.91,0.85);
  myLegendSetUp(myLegend1c3,mainFont);
  myLegend1c3->AddEntry(fit_v1_odd_eta_stat_30_60," 30-60% with fit","L");
  
  TLegend *myLegend1Sys = new TLegend(0.78,0.5,1.06,0.8);
  myLegendSetUp(myLegend1Sys,mainFont);
  myLegend1Sys->AddEntry(v1_odd_eta_syst_10_20_clone," ","F");
  myLegend1Sys->AddEntry(v1_odd_eta_syst_30_40_clone," ","F");
  myLegend1Sys->AddEntry(v1_odd_eta_syst_30_60_clone," ","F");
  
  TLegend *myLegend1 = new TLegend(0.78,0.5,1.06,0.8);
  myLegendSetUp(myLegend1,mainFont);
  myLegend1->AddEntry(v1_odd_eta_stat_10_20,"  10-20%","P");
  myLegend1->AddEntry(v1_odd_eta_stat_30_40,"  30-40%","P");
  myLegend1->AddEntry(v1_odd_eta_stat_30_60,"  30-60%","P");
  
  TLegend *myLegend3 = new TLegend(0.135,0.3,0.43,0.5);
  myLegendSetUp(myLegend3,0.95*mainFont);
  myLegend3->AddEntry(v1_star_AuAu200_30_60_eta,"   #times 0.37  200GeV","P");
  myLegend3->AddEntry(v1_star_AuAu62_30_60_eta,"   #times 0.12  62.4GeV","P");
  
  // shift points along x-axis for visibility
  float shift =0.01;
  
  ShiftAlongXaxis(v1_odd_eta_stat_10_20, 4*shift);
  ShiftAlongXaxis(v1_odd_eta_syst_10_20, 4*shift);
  ShiftAlongXaxis(v1_odd_eta_stat_30_40, -4*shift);
  ShiftAlongXaxis(v1_odd_eta_syst_30_40, -4*shift);
  
  ShiftAlongXaxis(v1_even_eta_stat_10_20, 4*shift);
  ShiftAlongXaxis(v1_even_eta_syst_10_20, 4*shift);
  ShiftAlongXaxis(v1_even_eta_stat_30_40, -4*shift);
  ShiftAlongXaxis(v1_even_eta_syst_30_40, -4*shift);  
  
  //pad1
  myPad1->cd();
  myBlankHisto1->Draw();
  myText->SetTextSize(1.5*mainFont);
  myText->DrawLatex(0.93,0.81,"(a)");
  myText->SetTextSize(mainFont);
  myText->DrawLatex(0.6,0.84,"odd    even     v_{1}");
  
  myLegend1SysA->Draw();
  myLegend1a->Draw();
  myLegend1Sys1C->Draw();
  myLegend1b->Draw();
  myLegend1SysB->Draw();
  myLegend1SysC->Draw();
  
  fit_v1_odd_eta_stat_10_60->Draw("same");
  fit_v1_even_eta_stat_10_60->Draw("same");  
  v1_odd_eta_syst_10_20->Draw("eZ");
  v1_odd_eta_syst_30_40->Draw("eZ");
  v1_odd_eta_syst_10_60->Draw("eZ");
  v1_even_eta_syst_10_20->Draw("eZ");
  v1_even_eta_syst_30_40->Draw("eZ");
  v1_even_eta_syst_10_60->Draw("eZ");
  v1_odd_eta_stat_10_20->Draw("P,eZ");
  v1_odd_eta_stat_30_40->Draw("P,eZ");
  v1_odd_eta_stat_10_60->Draw("P,eZ");
  v1_even_eta_stat_10_20->Draw("P,eZ");
  v1_even_eta_stat_30_40->Draw("P,eZ");
  v1_even_eta_stat_10_60->Draw("P,eZ");
  
  //pad2
  myPad2->cd();
  myBlankHisto2->Draw();
  myText->SetTextSize(1.5*mainFont);
  myText->DrawLatex(0.93,0.88,"(b)");
  myText->SetTextSize(mainFont);
  myText->DrawLatex(0.15,0.1,"ALICE Pb-Pb@2.76TeV  p_{T}>0.15 GeV/c");
  myText->SetTextSize(mainFont);
  myText->DrawLatex(0.555,0.89,"odd    even   #LTp_{x}#GT/#LTp_{T}#GT");
  
  myLegend1SysA2->Draw();
  myLegend1a2->Draw();
  myLegend1b2->Draw();
  myLegend1SysB2->Draw();
  myLegend1SysC2->Draw();
  myLegend1c2->Draw();
  
  fit_px_odd_eta_stat_10_60->Draw("same");
  fit_px_even_eta_stat_10_60->Draw("same");  
  px_even_eta_syst_10_60->Draw("eZ");
  px_even_eta_stat_10_60->Draw("P,eZ");
  px_odd_eta_syst_10_60->Draw("eZ");
  px_odd_eta_stat_10_60->Draw("P,eZ");
  
  //pad3
  myPad3->cd();
  myBlankHisto3->Draw();
  myText->SetTextSize(1.5*mainFont);
  myText->DrawLatex(0.93,0.9,"(c)");  
  myText->SetTextSize(mainFont);
  myText->DrawLatex(0.655,0.89,"odd   v_{1}");
  myText->SetTextSize(0.95*mainFont);
  myText->DrawLatex(0.13,0.51,"STAR  (scaled)");
  myText->SetTextSize(0.95*mainFont);
  myText->DrawLatex(0.15,0.23,"Au-Au 30-60% p_{T}>0.15 GeV/c");
  
  myLegend1SysA3->Draw();
  myLegend1a3->Draw();
  myLegend1c3->Draw();
  myLegend3->Draw();
  
  ShiftAlongXaxis(v1_star_AuAu200_30_60_eta, 2*shift);
  ShiftAlongXaxis(v1_star_AuAu62_30_60_eta, -2*shift);
  
  fit_v1_odd_eta_stat_30_60->Draw("same");  
  v1_star_AuAu62_30_60_eta->Draw("P,eZ");
  v1_star_AuAu200_30_60_eta->Draw("P,eZ");
  v1_odd_eta_syst_30_60->Draw("eZ");
  v1_odd_eta_stat_30_60->Draw("P,eZ");
  
  // save output if option rWrite=true
  TString fileName="Figure2";
  myCan->Update();
  if (rWrite)  
  {
    myCan->SaveAs(fileName+".png");
    myCan->SaveAs(fileName+".eps");
    myCan->SaveAs(fileName+".pdf");
  }
}

// only helper functions below

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
  gStyle->SetTitleOffset(1.7,"y");  
  gStyle->SetTitleOffset(2,"xz");  
  gStyle->SetTitleSize(24,"x");  
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
  TString title = name;
  title.ReplaceAll("_"," ");
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
  gr->SetName(name);
  gr->SetTitle(title);
  delete [] x;
  delete [] y;
  delete [] xerr;
  delete [] yerr;
  return gr;
}

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

void ScaleXaxis(TGraphErrors *ge, Double_t scale)
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
    x*=scale;
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
 int myFillStyle =fillStyle,
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
