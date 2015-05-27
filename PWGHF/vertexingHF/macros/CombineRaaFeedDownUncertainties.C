#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLine.h"

#include <Riostream.h>

//_________________________________________________________________________________________
//
//  Macro to combine the the MonteCarlo B feed-down subtraction uncertainties on RAA
//
//   Take as input the output files from the HFPtSpectrumRaa macro
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
//      4. Flag to switch between Raa/Rcp files
//
//
//   Questions and complains to Z. Conesa del Valle
//_________________________________________________________________________________________

void CombineRaaFeedDownUncertainties(const char *raafileMethod1="",
				     const char *raafileMethod2="",
				     const char *raaOutFilename="", Bool_t isRaa=true)
{

  const char *namesRaa[7] = { "hRABvsPt","gRAB_Norm","gRAB_DataSystematics","gRAB_FeedDownSystematics","gRAB_ElossHypothesis","gRAB_FeedDownSystematicsElossHypothesis","gRAB_GlobalSystematics" };
  const char *namesRcp[7] = { "hRCPvsPt","gRCP_Norm","gRCP_DataSystematics","gRCP_FeedDownSystematics","gRCP_ElossSystematics","gRCP_FeedDownElossSystematics","gRCP_GlobalSystematics" };
  const char *namesRead[7];
  for(int i=0; i<7; i++) {
    if(isRaa) namesRead[i]=namesRaa[i];
    else namesRead[i]=namesRcp[i];
  }
  const char *suffix="AB"; if(!isRaa) suffix="CP";
  //
  // Nb FD method Raa
  //
  TFile * raaf2 = new TFile(raafileMethod2,"read");
  TH1D * hRaa2 = (TH1D*)raaf2->Get(namesRead[0]);
  TGraphAsymmErrors * gRAB_Norm = (TGraphAsymmErrors*)raaf2->Get(namesRead[1]);
  TGraphAsymmErrors * gRAB_Data2 = (TGraphAsymmErrors*)raaf2->Get(namesRead[2]);
  TGraphAsymmErrors * gRAB_Data2PP;
  TGraphAsymmErrors * gRAB_Data2AB;
  if(isRaa){
    gRAB_Data2PP = (TGraphAsymmErrors*)raaf2->Get("gRAB_DataSystematicsPP");
    gRAB_Data2AB = (TGraphAsymmErrors*)raaf2->Get("gRAB_DataSystematicsAB");
  }
  TGraphAsymmErrors * gRAB_FD2 = (TGraphAsymmErrors*)raaf2->Get(namesRead[3]);
  TGraphAsymmErrors * gRAB_Eloss2 = (TGraphAsymmErrors*)raaf2->Get(namesRead[4]);
  TGraphAsymmErrors * gRAB_FDEloss2 = (TGraphAsymmErrors*)raaf2->Get(namesRead[5]);
  TGraphAsymmErrors * gRAB_Global2 =  (TGraphAsymmErrors*)raaf2->Get(namesRead[6]);
  hRaa2->SetTitle( Form("Nb : %s",hRaa2->GetTitle()) );
  gRAB_Data2->SetTitle( Form("Nb : %s",gRAB_Data2->GetTitle()) );
//   gRAB_Data2PP->SetTitle( Form("Nb : %s",gRAB_Data2PP->GetTitle()) );
//   gRAB_Data2AB->SetTitle( Form("Nb : %s",gRAB_Data2AB->GetTitle()) );
  gRAB_FDEloss2->SetTitle( Form("Nb : %s",gRAB_FDEloss2->GetTitle()) );
  gRAB_Global2->SetTitle( Form("Nb : %s",gRAB_Global2->GetTitle()) );
  cout << " Read Raa with Nb method file "<< endl;
  //
  // fc FD method Raa
  //
  TFile * raaf1 = new TFile(raafileMethod1,"read");
  TH1D * hRaa1 = (TH1D*)raaf1->Get(namesRead[0]);
  TGraphAsymmErrors * gRAB_Data1 = (TGraphAsymmErrors*)raaf1->Get(namesRead[2]);
  TGraphAsymmErrors * gRAB_FD1 = (TGraphAsymmErrors*)raaf1->Get(namesRead[3]);
  TGraphAsymmErrors * gRAB_Eloss1 = (TGraphAsymmErrors*)raaf1->Get(namesRead[4]);
  TGraphAsymmErrors * gRAB_FDEloss1 = (TGraphAsymmErrors*)raaf1->Get(namesRead[5]);
  TGraphAsymmErrors * gRAB_Global1 =  (TGraphAsymmErrors*)raaf1->Get(namesRead[6]);
  hRaa1->SetTitle( Form("fc : %s",hRaa1->GetTitle()) );
  gRAB_Data1->SetTitle( Form("fc : %s",gRAB_Data1->GetTitle()) );
  gRAB_FDEloss1->SetTitle( Form("fc : %s",gRAB_FDEloss1->GetTitle()) );
  gRAB_Global1->SetTitle( Form("fc : %s",gRAB_Global1->GetTitle()) );
  cout << " Read Raa with fc method file "<< endl;
  //
  // Output combined Raa
  //
  TFile * out = new TFile(raaOutFilename,"recreate");
  TH1D * hRaa = (TH1D*)hRaa2->Clone(Form("hR%svsPt",suffix));
  hRaa->SetNameTitle(Form("%s",namesRead[0]),Form(" R_{%s}(c) vs p_{T} for both Nb & fc",suffix));
  TGraphAsymmErrors * gRAB_DataSystematics =  new TGraphAsymmErrors();
  gRAB_DataSystematics->SetNameTitle(Form("gR%s_DataSystematics",suffix),Form("Data only systematics on R_{%s}",suffix));
  TGraphAsymmErrors * gRAB_FeedDownSystematicsElossHypothesis = new TGraphAsymmErrors();
  gRAB_FeedDownSystematicsElossHypothesis->SetNameTitle(Form("gR%s_FeedDownSystematicsElossHypothesis",suffix),"Feed-down + Eloss systematics for both Nb & fc");
  TGraphAsymmErrors * gRAB_FeedDownSystematics = new TGraphAsymmErrors();
  gRAB_FeedDownSystematics->SetNameTitle(Form("gR%s_FeedDownSystematics",suffix),"Feed-down systematics for both Nb & fc");
  TGraphAsymmErrors * gRAB_ElossHypothesis = new TGraphAsymmErrors();
  gRAB_ElossHypothesis->SetNameTitle(Form("gR%s_ElossHypothesis",suffix),"Eloss systematics for both Nb & fc");
  TGraphAsymmErrors * gRAB_GlobalSystematics = new TGraphAsymmErrors();
  gRAB_GlobalSystematics->SetNameTitle(Form("gR%s_GlobalSystematics",suffix),"Global systematics (data + FD + Eloss) for both Nb & fc");
  cout << " Created output file : "<< raaOutFilename << endl;

  //
  // Loop on the pt bins
  //
  Double_t pt=0., Raa=0., ptwidth=0.;
  Int_t istartDS=0, istartFDEl1=0,istartFDEl2=0, istartFD1=0,istartFD2=0, istartEl1=0,istartEl2=0;
  Double_t dataHigh=0., dataLow=0., gfdElossHigh=0, gfdElossLow=0., gfdHigh=0, gfdLow=0., gElossHigh=0, gElossLow=0., globalHigh=0, globalLow=0.;
  Double_t ElossHigh1=0., ElossLow1=0., ElossHigh2=0., ElossLow2=0., fdHigh1=0., fdLow1=0., fdHigh2=0., fdLow2=0.;
  Double_t vfdHigh1=0., vfdLow1=0., vfdHigh2=0., vfdLow2=0.;
  Double_t fdElossHigh1=0., fdElossLow1=0., fdElossHigh2=0., fdElossLow2=0.;
  Double_t vfdElossHigh1=0., vfdElossLow1=0., vfdElossHigh2=0., vfdElossLow2=0.;
  Double_t x=0., y=0.;
  Int_t nbins = hRaa->GetNbinsX();
  Int_t jbins = 0.;
  for(Int_t ibin=1; ibin<=nbins; ibin++){

    cout <<" loop, i="<<ibin<<endl;
    pt = hRaa->GetBinCenter(ibin);
    Raa = hRaa->GetBinContent(ibin);
    if( Raa<=0. ) continue;
    cout <<"  pt="<<pt<<", Raa="<<Raa<<endl;

    // Data syst reading
    Int_t n = gRAB_Data2->GetN();
    for(Int_t j=0; j<=n; j++){
      gRAB_Data2->GetPoint(j,x,y);
      if ( TMath::Abs ( x -pt ) < 0.4 ) { 
	istartDS = j; 
	break;
      }
    }
    dataHigh = gRAB_Data2->GetErrorYhigh(istartDS);
    dataLow = gRAB_Data2->GetErrorYlow(istartDS);
    cout <<" Data syst : + "<< dataHigh << " - "<< dataLow << " = +"<< dataHigh/Raa <<" -"<< dataLow/Raa << " %"<<endl;

    // Feed-Down systematics reading (Nb method)
    n = gRAB_FD2->GetN();
    for(Int_t j=0; j<=n; j++){
      gRAB_FD2->GetPoint(j,x,y);
      if ( TMath::Abs ( x -pt ) < 0.4 ) { 
	istartFD2 = j; 
	break;
      }
    }
    fdHigh2 = gRAB_FD2->GetErrorYhigh(istartFD2);
    fdLow2 = gRAB_FD2->GetErrorYlow(istartFD2);
    vfdHigh2 = y + fdHigh2;
    vfdLow2 = y - fdLow2;
    cout << " FD syst : Nb "<< Raa<<" +"<< fdHigh2 << " -"<< fdLow2 <<" = ("<<vfdLow2<<","<<vfdHigh2<<") = +" << fdHigh2/Raa << " -" << fdLow2/Raa << " %" <<endl;

    // Feed-Down systematics reading (fc method)
    n = gRAB_FD1->GetN();
    for(Int_t j=0; j<=n; j++){
      gRAB_FD1->GetPoint(j,x,y);
      if ( TMath::Abs ( x -pt ) < 0.4 ) { 
	istartFD1 = j; 
	break;
      }
    }
    fdHigh1 = gRAB_FD1->GetErrorYhigh(istartFD1);
    fdLow1 = gRAB_FD1->GetErrorYlow(istartFD1);
    vfdHigh1 = y + fdHigh1;
    vfdLow1 = y - fdLow1;
    cout << " FD syst : fc "<< y<<" +"<< fdHigh1 << " -"<< fdLow1 <<" = ("<<vfdLow1<<","<<vfdHigh1
	 <<") = +" << fdHigh1/Raa << " -" << fdLow1/Raa << " %" <<endl;


    // Eloss systematics reading (Nb method)
    n = gRAB_Eloss2->GetN();
    for(Int_t j=0; j<=n; j++){
      gRAB_Eloss2->GetPoint(j,x,y);
      if ( TMath::Abs ( x -pt ) < 0.4 ) { 
	istartEl2 = j; 
	break;
      }
    }
    ElossHigh2 = gRAB_Eloss2->GetErrorYhigh(istartEl2);
    ElossLow2 = gRAB_Eloss2->GetErrorYlow(istartEl2);
    cout << " Eloss syst : Nb "<< Raa<<" +"<< ElossHigh2 << " -"<< ElossLow2 <<" = + "<< ElossHigh2/Raa <<" - "<< ElossLow2/Raa << " %" <<endl;

    // Eloss systematics reading (fc method)
    n = gRAB_Eloss1->GetN();
    for(Int_t j=0; j<=n; j++){
      gRAB_Eloss1->GetPoint(j,x,y);
      if ( TMath::Abs ( x -pt ) < 0.4 ) { 
	istartEl1 = j; 
	break;
      }
    }
    ElossHigh1 = gRAB_Eloss1->GetErrorYhigh(istartEl1);
    ElossLow1 = gRAB_Eloss1->GetErrorYlow(istartEl1);
    cout << " Eloss syst : fc "<< Raa<<" +"<< ElossHigh1 << " -"<< ElossLow1 <<" = + "<< ElossHigh1/Raa <<" - "<< ElossLow1/Raa << " %" <<endl;


    // Feed-Down + Eloss systematics reading (Nb method)
    n = gRAB_FDEloss2->GetN();
    for(Int_t j=0; j<=n; j++){
      gRAB_FDEloss2->GetPoint(j,x,y);
      if ( TMath::Abs ( x -pt ) < 0.4 ) { 
	istartFDEl2 = j; 
	break;
      }
    }
    fdElossHigh2 = gRAB_FDEloss2->GetErrorYhigh(istartFDEl2);
    fdElossLow2 = gRAB_FDEloss2->GetErrorYlow(istartFDEl2);
    vfdElossHigh2 = y + fdElossHigh2;
    vfdElossLow2 = y - fdElossLow2;
    cout << " FD&Eloss syst : Nb "<< Raa<<" +"<< fdElossHigh2 << " -"<< fdElossLow2 <<" = ("<<vfdElossLow2<<","<<vfdElossHigh2<<")" <<endl;

    // Feed-Down + Eloss systematics reading (fc method)
    n = gRAB_FDEloss1->GetN();
    for(Int_t j=0; j<=n; j++){
      gRAB_FDEloss1->GetPoint(j,x,y);
      if ( TMath::Abs ( x -pt ) < 0.4 ) { 
	istartFDEl1 = j; 
	break;
      }
    }
    fdElossHigh1 = gRAB_FDEloss1->GetErrorYhigh(istartFDEl1);
    fdElossLow1 = gRAB_FDEloss1->GetErrorYlow(istartFDEl1);
    vfdElossHigh1 = y + fdElossHigh1;
    vfdElossLow1 = y - fdElossLow1;
    cout << " FD&Eloss syst : fc "<< y<<" +"<< fdElossHigh1 << " -"<< fdElossLow1 <<" = ("<<vfdElossLow1<<","<<vfdElossHigh1
	 <<") = +" << fdElossHigh1/Raa << " -" << fdElossLow1/Raa << " %" <<endl;

    //
    // Combine the FD uncertainties
    gfdLow = ( vfdLow2 < vfdLow1 ) ? vfdLow2 : vfdLow1;
    //    cout << " FD low : "<< gfdLow << " - Raa = ";
    gfdLow = TMath::Abs(gfdLow - Raa);
    //    cout << gfdLow<<endl;
    gfdHigh = ( vfdHigh2 > vfdHigh1 ) ? vfdHigh2 : vfdHigh1 ;
    //    cout << " FD high : "<< gfdHigh << " - Raa = ";
    gfdHigh = TMath::Abs(gfdHigh - Raa);
    //    cout << gfdHigh << endl;
    cout << " FD syst tot : +"<< gfdHigh << " -"<< gfdLow <<" = +"<< gfdHigh/Raa << " - "<< gfdLow/Raa <<" %"<<endl;

    //
    // Combine the Eloss uncertainties
    gElossLow = ( ElossLow2/Raa > ElossLow1/Raa ) ? ElossLow2/Raa : ElossLow1/Raa;
    gElossLow *= Raa ;
    //    cout << " low : "<< gElossLow << endl;
    gElossHigh = ( ElossHigh2/Raa > ElossHigh1/Raa ) ? ElossHigh2/Raa : ElossHigh1/Raa ;
    gElossHigh *= Raa;
    //    cout << " high : "<< gElossHigh << endl;
    cout << " Eloss syst tot : +"<< gElossHigh << " -"<< gElossLow << " = +"<< gElossHigh/Raa <<" -"<< gElossLow/Raa <<" %"<<endl;

    //
    // Combine the FD & Eloss uncertainties
    gfdElossLow = ( vfdElossLow2 < vfdElossLow1 ) ? vfdElossLow2 : vfdElossLow1;
    //    cout << " low : "<< gfdElossLow << " - Raa = ";
    gfdElossLow = TMath::Abs(gfdElossLow - Raa);
    //    cout << gfdElossLow<<endl;
    gfdElossHigh = ( vfdElossHigh2 > vfdElossHigh1 ) ? vfdElossHigh2 : vfdElossHigh1 ;
    //    cout << " high : "<< gfdElossHigh << " - Raa = ";
    gfdElossHigh = TMath::Abs(gfdElossHigh - Raa);
    //    cout << gfdElossHigh << endl;
    cout << " FD&Eloss syst tot : +"<< gfdElossHigh << " -"<< gfdElossLow << " = +"<< gfdElossHigh/Raa <<" -"<< gfdElossLow/Raa << " %"<<endl;

    //
    // Global uncertainties
    globalHigh = TMath::Sqrt( dataHigh*dataHigh + gfdElossHigh*gfdElossHigh );
    globalLow = TMath::Sqrt( dataLow*dataLow + gfdElossLow*gfdElossLow );
    cout << " Global syst tot : +"<< globalHigh << " -"<< globalLow << " = + "<< globalHigh/Raa <<" - " << globalLow/Raa << " %"<<endl;

    ptwidth = hRaa->GetBinWidth(ibin)/2.;
    gRAB_DataSystematics->SetPoint(jbins,pt,Raa);
    gRAB_DataSystematics->SetPointError(jbins,ptwidth,ptwidth,dataLow,dataHigh);
    ptwidth = 0.15;
    gRAB_FeedDownSystematicsElossHypothesis->SetPoint(jbins,pt,Raa);
    gRAB_FeedDownSystematicsElossHypothesis->SetPointError(jbins,ptwidth,ptwidth,gfdElossLow,gfdElossHigh);
    ptwidth = 0.3;
    gRAB_GlobalSystematics->SetPoint(jbins,pt,Raa);
    gRAB_GlobalSystematics->SetPointError(jbins,ptwidth,ptwidth,globalLow,globalHigh);
    jbins++;
  }


  TH2D *hRaaCanvas = new TH2D("hRaaCanvas",Form(" R_{%s}(c) vs p_{T} (for both Nb & fc); p_{t} (GeV/c) ; R_{%s} prompt D",suffix,suffix),25,0.,25.,100,0.,2.0);
  hRaaCanvas->GetXaxis()->SetTitleSize(0.05);
  hRaaCanvas->GetXaxis()->SetTitleOffset(0.9);
  hRaaCanvas->GetYaxis()->SetTitleSize(0.05);
  hRaaCanvas->GetYaxis()->SetTitleOffset(0.9);


  TCanvas *cdata = new TCanvas("Compare data syst unc.");
  hRaaCanvas->Draw();
  hRaa->Draw("same");
  gRAB_DataSystematics->SetFillStyle(3017);
  gRAB_DataSystematics->SetLineColor(4);
  gRAB_DataSystematics->SetLineWidth(3);
  gRAB_DataSystematics->SetFillColor(4);
  gRAB_DataSystematics->Draw("2");
  hRaa2->SetMarkerStyle(28);
  hRaa2->SetMarkerColor(kGreen+2);
  hRaa2->Draw("same");
  gRAB_Data2->SetLineColor(kGreen+2);
  gRAB_Data2->Draw("2");
  hRaa1->SetMarkerStyle(24);
  hRaa1->SetMarkerColor(kOrange+7);
  hRaa1->Draw("same");
  gRAB_Data1->SetLineColor(kOrange+7);
  gRAB_Data1->Draw("2");
  TLegend* leg = cdata->BuildLegend();
  leg->SetFillStyle(0);
  cdata->Update();
  TLine *line = new TLine(0.0172415,1.0,25.,1.0);
  line->SetLineStyle(2);
  line->Draw();


  TCanvas *cfd = new TCanvas("Compare FD & Eloss syst unc.");
  hRaaCanvas->Draw();
  hRaa->Draw("same");
  gRAB_FeedDownSystematicsElossHypothesis->SetFillStyle(3017);
  gRAB_FeedDownSystematicsElossHypothesis->SetLineColor(4);
  gRAB_FeedDownSystematicsElossHypothesis->SetLineColor(kBlack);
  gRAB_FeedDownSystematicsElossHypothesis->SetFillStyle(0);
  gRAB_FeedDownSystematicsElossHypothesis->SetFillColor(kViolet+1);
  gRAB_FeedDownSystematicsElossHypothesis->Draw("2");
  hRaa2->Draw("same");
  gRAB_FDEloss2->SetLineColor(kGreen+2);
  gRAB_FDEloss2->Draw("2");
  hRaa1->Draw("same");
  gRAB_FDEloss1->SetLineColor(kOrange+7);
  gRAB_FDEloss1->Draw("2");
  leg = cfd->BuildLegend();
  leg->SetFillStyle(0);
  cfd->Update();
  line->Draw();

  TCanvas *cglobal = new TCanvas("Compare global syst unc.");
  hRaaCanvas->Draw();
  hRaa->Draw("same");
  gRAB_GlobalSystematics->SetFillStyle(3017);
  gRAB_GlobalSystematics->SetLineColor(4);
  gRAB_GlobalSystematics->Draw("2");
  hRaa2->Draw("same");
  gRAB_Global2->SetLineColor(kGreen+2);
  gRAB_Global2->Draw("2");
  hRaa1->Draw("same");
  gRAB_Global1->SetLineColor(kOrange+7);
  gRAB_Global1->Draw("2");
  leg = cglobal->BuildLegend();
  leg->SetFillStyle(0);
  cglobal->Update();
  line->Draw();

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TCanvas *RaaPlot = new TCanvas("RaaPlot",Form("R%s vs pt, plot all",suffix));
  RaaPlot->SetTopMargin(0.085);
  RaaPlot->SetBottomMargin(0.1);
  RaaPlot->SetTickx();
  RaaPlot->SetTicky();
//   TH2D *hRaaCanvas = new TH2D("hRaaCanvas"," R_{AB}(c) vs p_{T} (no Eloss hypothesis); p_{t} [GeV/c] ; R_{AA} prompt D",25,0.,25.,100,0.,2.0);
//   hRaaCanvas->GetXaxis()->SetTitleSize(0.05);
//   hRaaCanvas->GetXaxis()->SetTitleOffset(0.9);
//   hRaaCanvas->GetYaxis()->SetTitleSize(0.05);
//   hRaaCanvas->GetYaxis()->SetTitleOffset(0.9);
  hRaaCanvas->Draw();
  gRAB_Norm->SetFillStyle(1001);
  gRAB_Norm->SetFillColor(kGray+2);
  gRAB_Norm->Draw("2");
  line = new TLine(0.0172415,1.0,25.,1.0);
  line->SetLineStyle(2);
  line->Draw();
  //  hRaa->SetLineColor(kBlack);
  hRaa->Draw("psame");
  //  gRAB_GlobalSystematics->SetLineColor(kBlack);
  gRAB_GlobalSystematics->SetFillStyle(0);
  gRAB_GlobalSystematics->Draw("2");
  RaaPlot->Draw();

  //
  //
  // Write output file
  //
  out->cd();
  hRaa->Write();
  gRAB_Norm->Write();
  gRAB_DataSystematics->Write();
  if(isRaa){
    gRAB_Data2PP->Write();
    gRAB_Data2AB->Write();
  }
  gRAB_FeedDownSystematics->Write();
  gRAB_ElossHypothesis->Write();
  gRAB_FeedDownSystematicsElossHypothesis->Write();
  gRAB_GlobalSystematics->Write();

}
