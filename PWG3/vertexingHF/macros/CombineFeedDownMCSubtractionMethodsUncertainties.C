#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"

#include "AliHFSystErr.h"

#include <Riostream.h>

//_________________________________________________________________________________________
//
//  Macro to combine the the MonteCarlo B feed-down subtraction uncertainties
//
//   Take as input the output files from the HFPtSpectrum class 
//    from both fc & Nb subtraction methods and combine the uncertainties. 
//   The final central value is set as the one from the Nb-method. 
//   The final uncertainties are defined as the envelope of both fc & Nb
//      uncertainties with respect to the new central-value.
//   The final global uncertainties are also defined and a preliminary drawing done. 
//
//
//   Usage parameters:
//      1. HFPtSpectrum fc subtraction file 
//      2. HFPtSpectrum Nb subtraction file 
//      3. Output file name
//      4. FONLL theoretical predictions file to draw on top
//      5. Decay channel as defined in the AliHFSystErr class
//
//_________________________________________________________________________________________

enum centrality{ kpp7, kpp276, k010, k020, k2040, k4060, k6080, k4080, k80100 };

void CombineFeedDownMCSubtractionMethodsUncertainties(const char *fcfilename="HFPtSpectrum_D0Kpi_method1_221110_newnorm.root",
						      const char *nbfilename="HFPtSpectrum_D0Kpi_method2_221110_newnorm.root",
						      const char *outfilename="HFPtSpectrum_D0Kpi_combinedFD.root",
						      const char *thfilename="D0DplusDstarPredictions_y05.root",
						      Int_t decay=1, Int_t centrality=kpp7)
{
  
  // 
  // Get fc file inputs
  TFile * fcfile = new TFile(fcfilename,"read");
  TH1D * histoSigmaCorrFc = (TH1D*)fcfile->Get("histoSigmaCorr");
  histoSigmaCorrFc->SetNameTitle("histoSigmaCorrFc","histoSigmaCorrFc");
  TGraphAsymmErrors * gSigmaCorrFc = (TGraphAsymmErrors*)fcfile->Get("gSigmaCorr");
  gSigmaCorrFc->SetNameTitle("gSigmaCorrFc","gSigmaCorrFc");
  TGraphAsymmErrors * gSigmaCorrConservativeFc = (TGraphAsymmErrors*)fcfile->Get("gSigmaCorrConservative");
  gSigmaCorrConservativeFc->SetNameTitle("gSigmaCorrConservativeFc","gSigmaCorrConservativeFc");

  // 
  // Get Nb file inputs
  TFile * nbfile = new TFile(nbfilename,"read");
  TH1D * histoSigmaCorrNb = (TH1D*)nbfile->Get("histoSigmaCorr");
  histoSigmaCorrNb->SetNameTitle("histoSigmaCorrNb","histoSigmaCorrNb");
  TGraphAsymmErrors * gSigmaCorrNb = (TGraphAsymmErrors*)nbfile->Get("gSigmaCorr");
  gSigmaCorrNb->SetNameTitle("gSigmaCorrNb","gSigmaCorrNb");
  TGraphAsymmErrors * gSigmaCorrConservativeNb = (TGraphAsymmErrors*)nbfile->Get("gSigmaCorrConservative");
  gSigmaCorrConservativeNb->SetNameTitle("gSigmaCorrConservativeNb","gSigmaCorrConservativeNb");

  //
  // Get the predictions input
  TFile *thfile = new TFile(thfilename,"read");
  TGraphAsymmErrors * thD0KpifromBprediction = (TGraphAsymmErrors*)thfile->Get("D0Kpiprediction");
  TGraphAsymmErrors * thDpluskpipiprediction = (TGraphAsymmErrors*)thfile->Get("Dpluskpipiprediction");
  TGraphAsymmErrors * thDstarD0piprediction = (TGraphAsymmErrors*)thfile->Get("DstarD0piprediction");
  thD0KpifromBprediction->SetLineColor(4);
  thD0KpifromBprediction->SetFillColor(kAzure+9);
  thDpluskpipiprediction->SetLineColor(4);
  thDpluskpipiprediction->SetFillColor(kAzure+9);
  thDstarD0piprediction->SetLineColor(4);
  thDstarD0piprediction->SetFillColor(kAzure+9);

  //
  // Get the spectra bins & limits
  Int_t nbins = histoSigmaCorrFc->GetNbinsX();
  Double_t *limits = new Double_t[nbins+1];
  Double_t xlow=0., binwidth=0.;
  for (Int_t i=1; i<=nbins; i++) {
    binwidth = histoSigmaCorrFc->GetBinWidth(i);
    xlow = histoSigmaCorrFc->GetBinLowEdge(i);
    limits[i-1] = xlow;
  }
  limits[nbins] = xlow + binwidth;


  //
  // Define a new histogram with the real-data reconstructed spectrum binning 
  //   they will be filled with central value equal to the Nb result
  //   and uncertainties taken from the envelope of the result uncertainties
  //   The systematical unc. (but FD) will also be re-calculated
  //
  TH1D * histoSigmaCorr = new TH1D("histoSigmaCorr","corrected cross-section (combined fc and Nb MC feed-down subtraction)",nbins,limits);
  TGraphAsymmErrors * gSigmaCorr = new TGraphAsymmErrors(nbins+1);
  gSigmaCorr->SetNameTitle("gSigmaCorr","gSigmaCorr (combined fc and Nb MC FD)");
  TGraphAsymmErrors * gSigmaCorrConservative = new TGraphAsymmErrors(nbins+1);
  gSigmaCorrConservative->SetNameTitle("gSigmaCorrConservative","Conservative gSigmaCorr (combined fc and Nb MC FD)");
  TGraphAsymmErrors * gSigmaCorrConservativePC = new TGraphAsymmErrors(nbins+1);
  gSigmaCorrConservativePC->SetNameTitle("gSigmaCorrConservativePC","Conservative gSigmaCorr (combined fc and Nb MC FD) in percentages [for drawing with AliHFSystErr]");

  //
  // Call the systematics uncertainty class for a given decay
  //   will help to compute the systematical unc. (but FD) 
  AliHFSystErr systematics;
  if( centrality==kpp276 ) {
    systematics.SetIsLowEnergy(true);
  } else if( centrality!=kpp7 )  {
    systematics.SetCollisionType(1);
    if ( centrality == k020 ) {
      systematics.SetCentrality("020");
    }
    else if ( centrality == k4080 ) {
      systematics.SetCentrality("4080");
    }
    else { 
      cout << " Systematics not yet implemented " << endl;
      return;
    }
  }
  else { systematics.SetCollisionType(0); }
  systematics.Init(decay);

  // 
  // Loop on all the bins to do the calculations
  //
  Double_t pt=0., average = 0., averageStatUnc=0., avErrx=0., avErryl=0., avErryh=0., avErryfdl=0., avErryfdh=0.;
  Double_t avErrylPC=0., avErryhPC=0., avErryfdlPC=0., avErryfdhPC=0.;
  Double_t valFc = 0., valFcErrstat=0., valFcErrx=0., valFcErryl=0., valFcErryh=0., valFcErryfdl=0., valFcErryfdh=0.;
  Double_t valNb = 0., valNbErrstat=0., valNbErrx=0., valNbErryl=0., valNbErryh=0., valNbErryfdl=0., valNbErryfdh=0.;
  //
  for(Int_t ibin=1; ibin<=nbins; ibin++){

    // Get input values from fc method
    valFc = histoSigmaCorrFc->GetBinContent(ibin);
    pt = histoSigmaCorrFc->GetBinCenter(ibin);
    valFcErrstat = histoSigmaCorrFc->GetBinError(ibin);
    Double_t value =0., ptt=0.;
    gSigmaCorrConservativeFc->GetPoint(ibin,ptt,value);
    if (value<=0.) continue;
    if ( TMath::Abs(valFc-value)>0.1 || TMath::Abs(pt-ptt)>0.1 ) 
      cout << "Hey you ! There might be a problem with the fc input file, please, have a look !" << endl;
    valFcErrx = gSigmaCorrFc->GetErrorXlow(ibin);
    valFcErryl = gSigmaCorrFc->GetErrorYlow(ibin);
    valFcErryh = gSigmaCorrFc->GetErrorYhigh(ibin);
    valFcErryfdl = TMath::Abs( gSigmaCorrConservativeFc->GetErrorYlow(ibin) );
    valFcErryfdh = TMath::Abs( gSigmaCorrConservativeFc->GetErrorYhigh(ibin) );

    // Get input values from Nb method
    valNb = histoSigmaCorrNb->GetBinContent(ibin);
    pt = histoSigmaCorrNb->GetBinCenter(ibin);
    valNbErrstat = histoSigmaCorrNb->GetBinError(ibin);
    gSigmaCorrConservativeNb->GetPoint(ibin,ptt,value);
    if ( TMath::Abs(valNb-value)>0.1 || TMath::Abs(pt-ptt)>0.1 ) 
      cout << "Hey you ! There might be a problem with the Nb input file, please, have a look !" << endl;
    valNbErrx = gSigmaCorrNb->GetErrorXlow(ibin);
    valNbErryl = gSigmaCorrNb->GetErrorYlow(ibin);
    valNbErryh = gSigmaCorrNb->GetErrorYhigh(ibin);
    valNbErryfdl = gSigmaCorrConservativeNb->GetErrorYlow(ibin);
    valNbErryfdh = gSigmaCorrConservativeNb->GetErrorYhigh(ibin);
    

    // Compute the FD combined value
    //    average = valNb
    average = valNb ;
    avErrx = valFcErrx;
    if ( TMath::Abs( valFcErrx - valNbErrx ) > 0.1 ) 
      cout << "Hey you ! There might be consistency problem with the fc & Nb input files, please, have a look !" << endl;
    averageStatUnc = valNbErrstat ;
//     cout << " pt=" << pt << ", average="<<average<<endl;
//     cout << "   stat unc (pc)=" << averageStatUnc/average << ", stat-fc (pc)="<<(valFcErrstat/valFc) << ", stat-Nb (pc)="<<(valNbErrstat/valNb)<<endl;
    
    // now estimate the new feed-down combined uncertainties
    Double_t minimum[2] = { (valFc - valFcErryfdl), (valNb - valNbErryfdl) };
    Double_t maximum[2] = { (valFc + valFcErryfdh), (valNb + valNbErryfdh) };
    avErryfdl = average - TMath::MinElement(2,minimum);
    avErryfdh = TMath::MaxElement(2,maximum) - average;
    avErryfdlPC = avErryfdl / average ; // in percentage
    avErryfdhPC = avErryfdh / average ; // in percentage
//     cout << " fc : val " << valFc << " + " << valFcErryfdh <<" - " << valFcErryfdl <<endl;
//     cout << " Nb : val " << valNb << " + " << valNbErryfdh <<" - " << valNbErryfdl <<endl;
//     cout << " fc  & Nb: val " << average << " + " << avErryfdh <<" - " << avErryfdl <<endl;


    // compute the global systematics
    avErrylPC = systematics.GetTotalSystErr(pt,avErryfdlPC); // in percentage
    avErryhPC = systematics.GetTotalSystErr(pt,avErryfdhPC); // in percentage
    avErryl = avErrylPC * average ;
    avErryh = avErryhPC * average ;
//     cout << "   syst av-l="<<avErryl<<", av-h="<<avErryh<<endl;
//     cout << "     fd-l-pc="<<avErryfdlPC<<", fd-h-pc="<<avErryfdhPC<<", syst err(no fd)-pc="<<systematics.GetTotalSystErr(pt)<<", av-l-pc="<<avErrylPC<<", av-h-pc="<<avErryhPC<<endl;

    // fill in the histos and TGraphs
    //   fill them only when for non empty bins
    if ( average > 0.1 ) {
      histoSigmaCorr->SetBinContent(ibin,average);
      histoSigmaCorr->SetBinError(ibin,averageStatUnc);
      gSigmaCorr->SetPoint(ibin,pt,average);
      gSigmaCorr->SetPointError(ibin,valFcErrx,valFcErrx,avErryl,avErryh);
      gSigmaCorrConservative->SetPoint(ibin,pt,average);
      gSigmaCorrConservative->SetPointError(ibin,valFcErrx,valFcErrx,avErryfdl,avErryfdh);
      gSigmaCorrConservativePC->SetPoint(ibin,pt,0.);
      gSigmaCorrConservativePC->SetPointError(ibin,valFcErrx,valFcErrx,avErryfdlPC,avErryfdhPC);
    }

  }


  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);

  //
  // Plot the results
  TH2F *histo2Draw = new TH2F("histo2Draw","histo2 (for drawing)",100,0,20.,100,1e4,1e10);
  histo2Draw->SetStats(0);
  histo2Draw->GetXaxis()->SetTitle("p_{T}  [GeV]");
  histo2Draw->GetXaxis()->SetTitleSize(0.05);
  histo2Draw->GetXaxis()->SetTitleOffset(0.95);
  histo2Draw->GetYaxis()->SetTitle("#frac{1}{BR} #times #frac{d#sigma}{dp_{T}} |_{|y|<0.5}");
  histo2Draw->GetYaxis()->SetTitleSize(0.05);
  //
  TCanvas *combinefdunc = new TCanvas("combinefdunc","show the FD results combination");
  //
  histo2Draw->Draw();
  //
  histoSigmaCorrFc->SetMarkerStyle(20);
  histoSigmaCorrFc->SetMarkerColor(kGreen+2);
  histoSigmaCorrFc->SetLineColor(kGreen+2);
  histoSigmaCorrFc->Draw("esame");
  gSigmaCorrConservativeFc->SetMarkerStyle(20);
  gSigmaCorrConservativeFc->SetMarkerColor(kGreen+2);
  gSigmaCorrConservativeFc->SetLineColor(kGreen+2);
  gSigmaCorrConservativeFc->SetFillStyle(3002);
  gSigmaCorrConservativeFc->SetFillColor(kGreen);
  gSigmaCorrConservativeFc->Draw("2[]same");
  //
  histoSigmaCorrNb->SetMarkerStyle(25);
  histoSigmaCorrNb->SetMarkerColor(kViolet+5);
  histoSigmaCorrNb->SetLineColor(kViolet+5);
  histoSigmaCorrNb->Draw("esame");
  gSigmaCorrConservativeNb->SetMarkerStyle(25);
  gSigmaCorrConservativeNb->SetMarkerColor(kViolet+5);
  gSigmaCorrConservativeNb->SetLineColor(kViolet+5);
  gSigmaCorrConservativeNb->SetFillStyle(3002);
  gSigmaCorrConservativeNb->SetFillColor(kMagenta);
  gSigmaCorrConservativeNb->Draw("2[]same");
  //
  gSigmaCorrConservative->SetLineColor(kRed);
  gSigmaCorrConservative->SetLineWidth(1);
  gSigmaCorrConservative->SetFillColor(kRed);
  gSigmaCorrConservative->SetFillStyle(0);
  gSigmaCorrConservative->Draw("2");
  histoSigmaCorr->SetMarkerColor(kRed);
  histoSigmaCorr->Draw("esame");
  //
  combinefdunc->SetLogy();
  combinefdunc->Update();

  //
  // Plot the results
  TCanvas *finalresults = new TCanvas("finalresults","show all combined results");
  //
  if ( decay==1 ) {
    thD0KpifromBprediction->SetLineColor(kGreen+2);
    thD0KpifromBprediction->SetLineWidth(3);
    thD0KpifromBprediction->SetFillColor(kGreen-6);
    thD0KpifromBprediction->Draw("3CA");
    thD0KpifromBprediction->Draw("CX");
  }
  else if ( decay==2 ) {
    thDpluskpipiprediction->SetLineColor(kGreen+2);
    thDpluskpipiprediction->SetLineWidth(3);
    thDpluskpipiprediction->SetFillColor(kGreen-6);
    thDpluskpipiprediction->Draw("3CA");
    thDpluskpipiprediction->Draw("CX");
  }
  else if ( decay==3 ) {
    thDstarD0piprediction->SetLineColor(kGreen+2);
    thDstarD0piprediction->SetLineWidth(3);
    thDstarD0piprediction->SetFillColor(kGreen-6);
    thDstarD0piprediction->Draw("3CA");
    thDstarD0piprediction->Draw("CX");
  }
  //
  gSigmaCorr->SetLineColor(kRed);
  gSigmaCorr->SetLineWidth(1);
  gSigmaCorr->SetFillColor(kRed);
  gSigmaCorr->SetFillStyle(0);
  gSigmaCorr->Draw("2");
  histoSigmaCorr->SetMarkerStyle(21);
  histoSigmaCorr->SetMarkerColor(kRed);
  histoSigmaCorr->Draw("esame");
  //
  TLegend * leg = new TLegend(0.7,0.75,0.87,0.5);
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  if ( decay==1 ) leg->AddEntry(thD0KpifromBprediction,"FONLL ","fl");
  else if ( decay==2 ) leg->AddEntry(thDpluskpipiprediction,"FONLL ","fl");
  else if ( decay==3 ) leg->AddEntry(thDstarD0piprediction,"FONLL ","fl");
  leg->AddEntry(histoSigmaCorr,"data stat. unc.","pl");
  leg->AddEntry(gSigmaCorr,"data syst. unc.","f");
  leg->Draw();
  //
  finalresults->SetLogy();
  finalresults->Update();


  //
  // Draw all the systematics independently
  systematics.DrawErrors(gSigmaCorrConservativePC);


  // Write the output to a file
  TFile * out = new TFile(outfilename,"recreate");
  histoSigmaCorr->Write();
  gSigmaCorr->Write();
  gSigmaCorrConservative->Write();
  gSigmaCorrConservativePC->Write();
  out->Write();

}
