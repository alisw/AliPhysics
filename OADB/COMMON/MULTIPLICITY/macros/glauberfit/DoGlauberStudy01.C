#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "AliMultGlauberNBDFitter.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLegend.h"

//________________________________________________________________
Double_t FastIntegrate(TF1 *f1, Double_t a, Double_t b, Int_t n = 5){
  //Do fast integration with N sampling points
  const Int_t nc = n;
  Double_t x[nc], y[nc];
  Double_t lWidth = (b-a)/((double)(n-1));
  for(Int_t ii=0; ii<n; ii++){
    x[ii] = a + ((double)(ii))*lWidth;
    y[ii] = f1->Eval( x[ii] );
  }
  //Now go via trapezoids, please (this probably has a name)
  Double_t lIntegral = 0;
  for(Int_t ii=0; ii<n-1; ii++){
    lIntegral += 0.5*lWidth*(y[ii]+y[ii+1]);
  }
  return lIntegral/(b-a);
}

void DoGlauberStudy01(Double_t lFitRange = 132, Int_t lMode = 2, Bool_t lFreek = 1){
  gStyle->SetLineScalePS(1);
  gStyle->SetOptStat(0);
  cout<<"This macro does a full glauber fit of data"<<endl;
  
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // Acquire data to start
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  TFile *file = new TFile("AnalysisResults_244918.root", "READ");
  TTree *fTree = (TTree*) file->Get("MultSelection/fTreeEvent");
  fTree->Print();
  
  TProfile *fProVtxZ = new TProfile("fProVtxZ","", 10,-10,10);
  fTree->Draw("(fAmplitude_V0A+fAmplitude_V0C):fEvSel_VtxZ>>fProVtxZ", "fnContributors>0&&fEvSel_Triggered&&TMath::Abs(fEvSel_VtxZ)<10", "goff");
  
  TCanvas *c0 = new TCanvas("c0", "", 800,600);
  fProVtxZ->Draw();
  
  TF1 *f1 = new TF1("f1", "[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x)", -10,10);
  fProVtxZ->Fit("f1", "REM0");
  f1->Draw("same");
  
  Int_t lRebin = 25;
  
  TH1D *hV0M = new TH1D("hV0M", "", 38000/lRebin, 0, 38000);
  TString lExpression = Form("(fAmplitude_V0A+fAmplitude_V0C)/(1+(%.10f)*fEvSel_VtxZ+(%.10f)*fEvSel_VtxZ*fEvSel_VtxZ+(%.10f)*fEvSel_VtxZ*fEvSel_VtxZ*fEvSel_VtxZ)",
                             f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3));
  //lExpression = "(fAmplitude_V0A+fAmplitude_V0C)";
  fTree->Draw(Form("%s>>hV0M",lExpression.Data()), "fnContributors>0&&fEvSel_Triggered&&TMath::Abs(fEvSel_VtxZ)<10", "goff");
  
  TCanvas *c1 = new TCanvas("c1", "", 1300,900);
  c1->Divide(1,2);
  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetTicks(1,1);
  c1->cd(1)->SetPad(0,0.5,1,1);
  c1->cd(2)->SetPad(0,0.0,1,.5);
  
  c1->cd(1)->SetBottomMargin(0.001);
  c1->cd(1)->SetRightMargin(0.25);
  c1->cd(1)->SetTopMargin(0.02);
  c1->cd(1)->SetLeftMargin(0.07);
  
  c1->cd(2)->SetBottomMargin(0.14);
  c1->cd(2)->SetRightMargin(0.25);
  c1->cd(2)->SetTopMargin(0.001);
  c1->cd(2)->SetLeftMargin(0.07);
  c1->cd(2)->SetTicks(1,1);
  c1->cd(1);
  
  hV0M->GetXaxis()->SetRangeUser(0,47000);
  hV0M->GetYaxis()->SetRangeUser(3e-2,3990);
  hV0M->SetLineColor(kBlack);
  hV0M->SetMarkerStyle(20);
  hV0M->SetMarkerColor(kBlack);
  hV0M->SetMarkerSize(0.5);
  hV0M->GetYaxis()->SetTitleSize(0.07);
  hV0M->GetYaxis()->SetLabelSize(0.05);
  hV0M->GetYaxis()->SetTitle("Count");
  hV0M->GetYaxis()->SetTitleOffset(0.5);
  hV0M->GetXaxis()->SetLabelSize(0.05);
  hV0M->GetXaxis()->SetTitleSize(0.06);
  hV0M->GetXaxis()->SetTitle("V0M Amplitude");
  hV0M->GetYaxis()->SetTickLength(0.015);
  hV0M->Draw("E");
  
  TFile *file2 = new TFile("V0M_fit_244918_2_0_toia.root", "READ");
  TH1F *hToia = (TH1F*) file2->Get("fHOutMultV0M");
  TH1F *hToiaGlau = (TH1F*) file2->Get("fHOutMultV0M_GLAU");
  cout<<"hToia Nbins native : "<<hToia->GetNbinsX()<<endl;
  cout<<"hToia max : "<<hToia->GetBinLowEdge(hToia->GetNbinsX()+2)<<endl;
  hToia->SetLineColor(kRed);
  hToiaGlau->SetLineColor(kMagenta);
  
//  hToia->Rebin(lRebin);
//  hToiaGlau->Rebin(lRebin);
//  hToia->Draw("same");
//  hToiaGlau->Draw("same");
//  return;
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  
  //Proof of concept fitter that makes use of the AliMultGlauberNBDFitter class
  //Step 0: create object of class
  AliMultGlauberNBDFitter *g = new AliMultGlauberNBDFitter(lMode,"fBazooka");
  
  //Step 1: open the (Npart, Ncoll) pair information, provide
  TFile *fbasefile = new TFile("basehistos.root","READ");
  TH2D *hNpNc = (TH2D*) fbasefile->Get("hNpNc");
  g->SetNpartNcollCorrelation(hNpNc);
  
  
  Double_t lFitRangeMax = 43000;
  
  g->SetInputV0M(hV0M);
  g->SetFitRange(lFitRange,lFitRangeMax);
  //Step 3: go for it ...
  g->SetNorm(4000000);
  TString lString = "REM0";
  g->SetFitOptions(lString.Data());
  g->SetFitNpx(100000);
  
  TF1 *fitfunc = g->GetGlauberNBD();
  
  fitfunc->SetParLimits(2,0.65,95);
  fitfunc->SetParLimits(0,10,100);
  fitfunc->SetParameter(0,45);
  //fitfunc->SetParameter(1,0.2);
  if(!lFreek){
    fitfunc->FixParameter(1,1.5);
  }else{
    fitfunc->SetParLimits(1,0.1,1.9);
    fitfunc->SetParameter(1,0.45);
  }
  fitfunc->FixParameter(2,0.8);
  fitfunc->SetParameter(3,7.11231e+06);
  
  //g->InitializeNpNc();
  g->DoFit();
  
  fitfunc->SetLineColor(kRed);
  fitfunc->SetLineWidth(2);
  fitfunc->Draw("same");
  
//  return;
  //Do a ratio plot
  TH1D *hRatio = (TH1D*) hV0M->Clone("hRatio");
  TH1D *hRatioWide = (TH1D*) hV0M->Clone("hRatioWide");
  hRatioWide->Rebin(20);
  hRatioWide->Scale(1./20.);
  //hRatio->Scale(1.,"width") ;
  
  for(Int_t ii=1; ii<hRatio->GetNbinsX()+1; ii++){
    Double_t lRatio = hRatio->GetBinContent(ii);
    //Double_t lFuncVal = fitfunc->Eval( hRatio->GetBinCenter(ii) );
    Double_t lFuncVal = FastIntegrate( fitfunc, hRatio->GetBinLowEdge(ii), hRatio->GetBinLowEdge(ii+1), 10);
    if ( lRatio < 2 ){
      hRatio->SetBinContent(ii,-1);
      hRatio->SetBinError(ii,1e-9);
      continue;
    }
    lRatio /= lFuncVal;
    hRatio->SetBinContent(ii, lRatio);
    hRatio->SetBinError(ii, hRatio->GetBinError(ii)/lFuncVal);
  }
  
  for(Int_t ii=1; ii<hRatioWide->GetNbinsX()+1; ii++){
    Double_t lRatio = hRatioWide->GetBinContent(ii);
    Double_t lFuncVal = FastIntegrate( fitfunc, hRatioWide->GetBinLowEdge(ii), hRatioWide->GetBinLowEdge(ii+1), 20);
    if ( lRatio < 2 ){
      hRatioWide->SetBinContent(ii,-1);
      hRatioWide->SetBinError(ii,1e-9);
      continue;
    }
    lRatio /= lFuncVal;
    hRatioWide->SetBinContent(ii, lRatio);
    hRatioWide->SetBinError(ii, hRatioWide->GetBinError(ii)/lFuncVal);
  }
  
  c1->cd(2);
  hRatio->GetYaxis()->SetTitle("Data/Fit");
  hRatio->GetYaxis()->SetRangeUser(0.45,1.55);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerColor(kGray+2);
  hRatio->SetLineColor(kGray+2);
  //hRatio->SetMarkerSize(1.0);
  hRatio->SetMarkerSize(.7);
  
  hRatio->Draw("hist");
  hRatioWide->SetLineColor(kRed);
  hRatioWide->SetMarkerColor(kRed);
  hRatioWide->SetMarkerStyle(20);
  hRatioWide->SetMarkerSize(1.2);
  hRatioWide->SetLineWidth(2);
  
  TLine *line = new TLine(0,1,38000,1);
  line->SetLineStyle(7);
  line->SetLineColor(kGray+1);
  line->Draw();
  
  TLine *lFitRangeLine = new TLine( lFitRange, 0.35, lFitRange, 1.65) ;
  lFitRangeLine->SetLineColor(kBlue);
  lFitRangeLine->SetLineWidth(2);
  lFitRangeLine->SetLineStyle(2);
  lFitRangeLine->Draw();
  
  TH1D *hRatioGrayed = (TH1D*) hRatio->Clone("hRatioGrayed");
  hRatioGrayed->SetMarkerColor(kGray+2);
  hRatioGrayed->SetLineColor(kGray+2);
  hRatioGrayed->Draw("same");
  
  hRatio->SetLineWidth(1);
  hRatio->Draw("same hist");
  hRatioWide->Draw("same");
  
  c1->cd(1);
  TLatex *lat = new TLatex();
  lat->SetNDC();
  Float_t lPosText = 0.76;
  Float_t lYShift = 0.25;
  lat->SetTextSize(0.042);
  lat->DrawLatex(lPosText,0.67+lYShift, "Pb-Pb 5.02 TeV");
  lat->DrawLatex(lPosText,0.61+lYShift, "Glauber + NBD fit");
  lat->DrawLatex(lPosText,0.55+lYShift, "2015 data, Run 244918");
  lat->SetTextFont(42);
  lat->DrawLatex(lPosText,0.49+lYShift, Form("Fit range: %.1f-%.1f", lFitRange, lFitRangeMax) );
  lat->DrawLatex(lPosText,0.43+lYShift, Form("#Chi^{2}/ndf: %.1f / %i = %.3f", fitfunc->GetChisquare(), fitfunc->GetNDF(), fitfunc->GetChisquare() / ((Double_t)(fitfunc->GetNDF() ) ) ) );
  lat->DrawLatex(lPosText,0.37+lYShift, Form("Fit options: %s", lString.Data() ) );
  lat->DrawLatex(lPosText,0.31+lYShift, Form("Glauber f: %.3f", fitfunc->GetParameter(2) ) );
  lat->DrawLatex(lPosText,0.25+lYShift, Form("Glauber #mu: %.3f", fitfunc->GetParameter(0) ) );
  lat->DrawLatex(lPosText,0.19+lYShift, Form("Glauber k: %.3f", fitfunc->GetParameter(1) ) );
  
  //Now extract hyper-fine cumulative function from fit function
  const Long_t lSamplePoints = 1e+5; //because, because. Just because.

  Double_t lDelta = 50000./((Double_t)(lSamplePoints));

  Double_t lX[lSamplePoints], lY[lSamplePoints];
  cout<<"Calculating cumulative function..."<<endl;
  lX[0] = 0; lY[0] = 0;
  for(Long_t ii=1; ii<lSamplePoints; ii++){
    if(ii%5000==0) cout<<"At sample #"<<ii<<endl;
    lX[ii] = ((Double_t) ii) * lDelta ;
    lY[ii] = lY[ii-1] + fitfunc->Eval(lX[ii]);
  }
  Int_t lFirstPointAbove = -1;
  for(Long_t ii=1; ii<lSamplePoints; ii++){
    lX[ii] = ((Double_t) ii) * lDelta ;
    lY[ii] = lY[ii] / lY[lSamplePoints-1] ; //Normalize
    if( lFirstPointAbove<0 && lY[ii]>0.1 ) lFirstPointAbove=ii;
  }
  cout<<"Integrating, please wait..."<<endl;
  Double_t lIntegralTotal = fitfunc->Integral(0.,38000.);
  cout<<"integral = "<<lIntegralTotal<<endl;
  cout<<"At ll = "<<132<<", evaluate fractional integral: "<<fitfunc->Integral(0., 132)/lIntegralTotal<<endl;

  TGraph *gr = new TGraph(lSamplePoints, lX, lY);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.3);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kBlack);

  TCanvas *cCumu = new TCanvas("cCumu", "", 800,600);
  cCumu->SetTicks(1,1);
  gr->Draw("ALP") ;

  //Determine anchor point based on hyperfinely binned TGraph, please
  Double_t lAnchorPoint;
  Double_t lFrac = (0.1 - gr->GetY()[lFirstPointAbove-1])/(gr->GetY()[lFirstPointAbove] - gr->GetY()[lFirstPointAbove-1]);
  lAnchorPoint = gr->GetX()[lFirstPointAbove-1] + lFrac*(gr->GetX()[lFirstPointAbove] - gr->GetX()[lFirstPointAbove-1]);
  cout<<"Anchor point determined to be: "<<lAnchorPoint<<endl;
  
  TGraph *grInverse = new TGraph(lSamplePoints, lY, lX);
  grInverse->SetMarkerStyle(20);
  grInverse->SetMarkerSize(0.3);
  grInverse->SetMarkerColor(kPink+3);
  grInverse->SetLineColor(kPink+3);
  
  TCanvas *cCumuInv = new TCanvas("cCumuInv", "", 800,600);
  cCumuInv->SetTicks(1,1);
  grInverse->Draw("ALP") ;
  
  Double_t lCentrality[] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0};
  Double_t lCentralityShiftedPlu[] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0};
  Double_t lCentralityShiftedMin[] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0};
  Double_t lCentralityValue[20];
  Double_t lCentralityValuePlu[20];
  Double_t lCentralityValueMin[20];
  
  Int_t lNCentrality = sizeof(lCentrality)/sizeof(Double_t);
  for(Int_t ii=1; ii<lNCentrality; ii++) lCentralityShiftedPlu[ii] *= 90.5/90.0;
  for(Int_t ii=1; ii<lNCentrality; ii++) lCentralityShiftedMin[ii] *= 89.5/90.0;
  
  lCentralityValue[0] = 0.0;
  lCentralityValue[lNCentrality-1] = 5e+4;
  lCentralityValuePlu[0] = 0.0;
  lCentralityValuePlu[lNCentrality-1] = 5e+4;
  lCentralityValueMin[0] = 0.0;
  lCentralityValueMin[lNCentrality-1] = 5e+4;
  for(Int_t ii=1; ii<lNCentrality-1; ii++){
    lCentralityValue[ii] = grInverse->Eval((100-lCentrality[ii])/100.0);
    lCentralityValuePlu[ii] = grInverse->Eval((100-lCentralityShiftedPlu[ii])/100.0);
    lCentralityValueMin[ii] = grInverse->Eval((100-lCentralityShiftedMin[ii])/100.0);
  }
  for(Int_t ii=0; ii<lNCentrality; ii++)
  cout<<"Centrality boundary for "<<lCentrality[ii]<<" is "<<lCentralityValue[ii]<<endl;
  
  TProfile *gRefMult = new TProfile("gRefMult", "", lNCentrality-1, lCentralityValue);
  TProfile *gRefMultPlu = new TProfile("gRefMultPlu", "", lNCentrality-1, lCentralityValuePlu);
  TProfile *gRefMultMin = new TProfile("gRefMultMin", "", lNCentrality-1, lCentralityValueMin);
  
  fTree->Draw(Form("fnTracklets08/1.6:%s>>gRefMult",lExpression.Data()), "fnContributors>0&&fEvSel_Triggered&&TMath::Abs(fEvSel_VtxZ)<10", "goff");
  fTree->Draw(Form("fnTracklets08/1.6:%s>>gRefMultPlu",lExpression.Data()), "fnContributors>0&&fEvSel_Triggered&&TMath::Abs(fEvSel_VtxZ)<10", "goff");
  fTree->Draw(Form("fnTracklets08/1.6:%s>>gRefMultMin",lExpression.Data()), "fnContributors>0&&fEvSel_Triggered&&TMath::Abs(fEvSel_VtxZ)<10", "goff");
  
  TCanvas *cNch = new TCanvas("cNch", "", 800,600);
  cNch->SetTicks(1,1);
  gRefMult->Draw("") ;
  
  /*
   TH1D *hCumulative = new TH1D("hCumulative", "", lSamplePoints, 0, 50000);
   cout<<"Calculating cumulative..."<<endl;
   hCumulative->SetBinContent(1,0);
   for(Long_t ii=1; ii<lSamplePoints; ii++){
   if(ii%50000==0) cout<<"At sample #"<<ii<<endl;
   hCumulative->SetBinContent(ii+1, hCumulative->GetBinContent(ii) +
   lDelta*0.5*(fitfunc->Eval(hCumulative->GetBinLowEdge(ii)) ) );
   }
   
   TCanvas *cCumu = new TCanvas("cCumu", "", 800,600);
   cCumu->SetTicks(1,1);
   hCumulative->Scale(1./hCumulative->GetBinContent( lSamplePoints-1 ));
   hCumulative->Draw();
   
   Double_t lAnchorPoint = 0;
   //Find precise 90% anchor point
   for(Long_t ii=1; ii<lSamplePoints; ii++){
   if( hCumulative->GetBinContent(ii) > 0.1){
   cout<<"Located anchor point for 90% anchoring!"<<endl;
   cout<<"It's at approximately: "<<hCumulative->GetBinCenter(ii)<<endl;
   cout<<"Will now smoothen it slightly for you, hang on"<<endl;
   Double_t lVal2 = hCumulative->GetBinLowEdge(ii);
   Double_t lVal1 = hCumulative->GetBinLowEdge(ii-1);
   
   //Proportion:
   Double_t lFrac = (10.-hCumulative->GetBinContent(ii-1))/
   (hCumulative->GetBinContent(11)-hCumulative->GetBinContent(ii-1));
   lAnchorPoint = lVal1 + lFrac*(lVal2-lVal1);
   cout<<"Better value: "<<lAnchorPoint<<endl;
   
   break;
   }
   }
   */
  
  Double_t lFracAnchoredOut = gr->Eval(lFitRange);
  Double_t lFracAnchoredOut400 = gr->Eval(400);
  Double_t lFracAnchoredOut600 = gr->Eval(600);
  c1->cd(1);
  //lat->DrawLatex(lPosText,0.17-0.05, Form("Glauber k: %.3f", fitfunc->GetParameter(1) ) );
  lat->DrawLatex(lPosText,0.38-0.00, Form("90%% anchor point: %.3f", lAnchorPoint ) );
  lat->DrawLatex(lPosText,0.32-0.00, Form("Perc. above 132: %.2f%%", 100.*(1-lFracAnchoredOut) ) );
  
  //Sample Nchs
  lat->SetTextSize(0.042);
  lat->DrawLatex(lPosText,0.26-0.00, Form("d#it{N}_{ch}^{raw}/d#eta(0-5%%): %.2f", gRefMult->GetBinContent(lNCentrality-1) ) );
  lat->DrawLatex(lPosText,0.20-0.00, Form("d#it{N}_{ch}^{raw}/d#eta(60-70%%): %.2f", gRefMult->GetBinContent(4) ) );
  lat->DrawLatex(lPosText,0.14-0.00, Form("d#it{N}_{ch}^{raw}/d#eta(70-80%%): %.2f", gRefMult->GetBinContent(3) ) );
  lat->DrawLatex(lPosText,0.08-0.00, Form("d#it{N}_{ch}^{raw}/d#eta(80-90%%): %.2f", gRefMult->GetBinContent(2) ) );
  lat->DrawLatex(lPosText,0.02-0.00, Form("d#it{N}_{ch}^{raw}/d#eta(90-100%%): %.2f", gRefMult->GetBinContent(1) ) );
  
  
  /* DEPRECATED
   for(Long_t ii=1; ii<lSamplePoints; ii++){
   if( hCumulative->GetBinCenter(ii) > lFitRange){
   cout<<"Located equivalent percentage of hadronic cross section"<<endl;
   //Proportion:
   Double_t lVal2 = hCumulative->GetBinContent(ii);
   Double_t lVal1 = hCumulative->GetBinContent(ii-1);
   Double_t lFrac = (lFitRange-hCumulative->GetBinCenter(ii-1))/
   (hCumulative->GetBinCenter(11)-hCumulative->GetBinCenter(ii-1));
   lFracAnchoredOut = lVal1 + lFrac*(lVal2-lVal1);
   cout<<"Fraction Anchored out: "<<lFracAnchoredOut<<endl;
   cout<<"Anchor percentile for fit: "<<1-lFracAnchoredOut<<endl;
   break;
   }
   }
   */
  
  
  
  

  
  c1->cd(2);
  TLegend *leg = new TLegend(0.15, 0.185, 0.400, 0.384);
  leg->SetBorderSize(0);
  leg->AddEntry(hRatio, "Data-to-fit (fine bins)", "l");
  leg->AddEntry(hRatioWide, "Data-to-fit (wide bins)", "lp");
  TH1D *hDummy = new TH1D("hDummy", "", 10,0,10); hDummy->SetLineColor(kBlue);
  hDummy->SetLineWidth(2);
  hDummy->SetLineStyle(2);
  leg->AddEntry(hDummy, "ALICE anchor point (132.5)");
  leg->Draw();
  
  c1->SaveAs(Form("glauberfit_%.0f.pdf",lFitRange));
  
  hRatio->GetYaxis()->SetRangeUser(.3,1.7);
  c1->SaveAs(Form("glauberfit-mode-%i-freek-%i.pdf",lMode, lFreek));
  hV0M->GetXaxis()->SetRangeUser(0,2500);
  hV0M->GetYaxis()->UnZoom();
  hRatio->GetXaxis()->SetRangeUser(0,2500);
  hRatio->GetYaxis()->SetRangeUser(.3,1.7);
  c1->cd(1)->SetLogy(kFALSE);
  line->SetX2(2500);
  c1->SaveAs(Form("glauberfitzoom-mode-%i-freek-%i.pdf",lMode, lFreek));
  
  cout<<"mode: "<<g->GetAncestorMode()<<endl;


  
  TCanvas *cdeb = new TCanvas("cdeb", "", 800, 600);
  TH1D *hAnc = g->GetAncestorHistogram();
  hAnc->Draw();

  //Estimate Npart, Ncoll in multiplicity bins based on the glauber ntuple
  TProfile *gNpart = new TProfile("gNpart", "", lNCentrality-1, lCentralityValue);
  TProfile *gNcoll = new TProfile("gNcoll", "", lNCentrality-1, lCentralityValue);
  
  g->CalculateAvNpNc(gNpart, gNcoll);

  
  
  
  
  TFile *fNch = new TFile(Form("nchtest-mode-%i-freek-%i.root", lMode, lFreek), "RECREATE");
  gRefMult->Write();
  gRefMultPlu->Write();
  gRefMultMin->Write();
  gNpart->Write();
  gNcoll->Write();
  fNch->Write();
  fNch->Close();
}
