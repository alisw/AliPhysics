///
/// \file CompareInvMassGraphsDataProd.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Plot InvMassFit.C output from different inputs in same canvas 
///
/// Example file to read output from InvMassFit.C, mass, width counts vs pT, per SM or all SM 
/// executed for different analysis options and data productions and plot them in the same canvas.
/// It only recovers and plots mass vs pT.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TObject.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include "PlotUtils.C"

#endif

TString fileFormat = ".eps"; /// File format: ".eps" ".pdf" ...

//------------------------------------------------------------------------------
/// Main method
///
/// \param particle  : "Pi0", "Eta" 
/// \param calorimeter : "EMCAL", "DCAL","PHOS"
/// \param trigger   : calo trigger type
//------------------------------------------------------------------------------
void CompareInvMassGraphsDataProd
( 
 TString particle    = "Eta",
 TString calorimeter = "DCAL",
 TString trigger     = "default"
)
{    
  // SM
  const Int_t nSM = 20; 
  Int_t ncol = 4;
  Int_t nrow = 4;
  
  Int_t firstSM  = 0;
  Int_t lastSM   = 11;
  
  if(calorimeter=="DCAL")
  {
    ncol = 3;
    nrow = 3;
    
    firstSM  = 12;
    lastSM   = 19;
  }
  //------------------------------------------------
  // Put the path to the files in prod string array
  // also the explanation title on prodLeg string array
  //------------------------------------------------
  
  const Int_t nProd = 5;
  
  TString prod[] = 
  {
    "LHC17_5TeV/CENTwoSDD"
    ,"LHC17_5TeV/FAST"
    ,"LHC17_13TeV"
    ,"LHC16_13TeV"
    ,"LHC15n"
  };
  
  TString prodLeg[] = 
  {    
     "LHC17pq, #sqrt{s} = 5 TeV, CENT"
    ,"LHC17pq, #sqrt{s} = 5 TeV, FAST"
    ,"LHC17*, #sqrt{s} = 13 TeV"
    ,"LHC16*, #sqrt{s} = 13 TeV"
    ,"LHC15n, #sqrt{s} = 5 TeV"
  };

  TString directory = Form("ComparisonGraphs/%s/%s/%s",calorimeter.Data(),trigger.Data(),particle.Data());

//  const Int_t nProd = 8;
//  
//  TString prod[] = 
//  {
//     "LHC17_5TeV/LHC17p/CENTwoSDD"
//    ,"LHC17_13TeV/LHC17i"
//    ,"LHC17_13TeV/LHC17l"
//    ,"LHC17_13TeV/LHC17o"
//    ,"LHC16_13TeV/LHC16i"
//    ,"LHC16_13TeV/LHC16k"
//    ,"LHC16_13TeV/LHC16o"
//    ,"LHC15n"
//  };
//  
//  TString prodLeg[] = 
//  {    
//    "LHC17p, #sqrt{s} = 5 TeV, CENT"
//    ,"LHC17i, #sqrt{s} = 13 TeV"
//    ,"LHC17l, #sqrt{s} = 13 TeV"
//    ,"LHC17o, #sqrt{s} = 13 TeV"
//    ,"LHC16i, #sqrt{s} = 13 TeV"
//    ,"LHC16l, #sqrt{s} = 13 TeV"
//    ,"LHC16o, #sqrt{s} = 13 TeV"
//    ,"LHC15n, #sqrt{s} = 5 TeV"
//  };
//
//  TString directory = Form("ComparisonGraphsPeriods/%s/%s/%s",calorimeter.Data(),trigger.Data(),particle.Data());
  
  // Create output directories
  TString processline = Form(".! mkdir -p %s",directory.Data()) ;
  gROOT->ProcessLine(processline.Data());
  
  // Init  the graph arrays.
  //
  TGraphErrors* gMass[nSM][nProd];
  TGraphErrors* gMassAllSM[nProd];
  TGraphErrors* gMassAll  [nProd];
  TH1F* hAxis[nSM];
  
  TFile* file[nProd];
  
  TGraphErrors* gMassRat[nSM][nProd];
  TGraphErrors* gMassRatAllSM[nProd];
  TGraphErrors* gMassRatAll  [nProd];
  
  // Init to 0 all arrays
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    gMassAllSM   [iprod] = 0;
    gMassAll     [iprod] = 0;
    gMassRatAllSM[iprod] = 0;
    gMassRatAll  [iprod] = 0;
    file         [iprod] = 0;
    
    for(Int_t ism = 0; ism < nSM; ism++)
    {
      gMass   [ism][iprod] = 0;
      gMassRat[ism][iprod] = 0;
      if ( iprod == 0 ) hAxis[ism] = 0;
    }
  }
  
  //===============================
  // Open the input files and recover the graphs
  // apply some style settings
  //===============================

  Float_t colorProd[] = { 1,kBlue,kRed,kViolet,kRed,kOrange-2,kCyan, kBlue};
  Float_t styleProd[] = {20,   24,  24,     24,  25,       25,   25,    20};  
  Float_t markerSize = 8;
    
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    //printf("iprod %d\n",iprod);
    file[iprod] = new TFile(Form("%s/%s/%s/%s/MassWidthPtHistograms.root",
                                 prod[iprod].Data(),calorimeter.Data(),trigger.Data(),particle.Data()),"read");
    
    printf("iprod %d, %s/%s/%s/%s/MassWidthPtHistograms.root %p\n",
           iprod,prod[iprod].Data(),calorimeter.Data(),trigger.Data(),particle.Data(),file[iprod]);
    
    if(!file[iprod]) continue;
    
    gMassAllSM[iprod] = (TGraphErrors*) file[iprod]->Get("gMass_SameSM");
    //printf("\t gMass_SameSM %p\n",gMassAllSM[iprod]);
    
    gMassAllSM[iprod]->SetMarkerColor(colorProd[iprod]);
    gMassAllSM[iprod]->SetLineColor  (colorProd[iprod]);
    gMassAllSM[iprod]->SetMarkerStyle(styleProd[iprod]);
    gMassAllSM[iprod]->SetMarkerSize (markerSize);
 
    
    gMassAll[iprod] = (TGraphErrors*) file[iprod]->Get("gMass_AllSM");
    //printf("\t gMass_SameSM %p\n",gMassAllSM[iprod]);
    
    gMassAll[iprod]->SetMarkerColor(colorProd[iprod]);
    gMassAll[iprod]->SetLineColor  (colorProd[iprod]);
    gMassAll[iprod]->SetMarkerStyle(styleProd[iprod]);
    gMassAll[iprod]->SetMarkerSize (markerSize);
    
    if ( iprod > 0 && gMassAllSM[0] ) 
    {
      gMassRatAllSM[iprod] = DivideGraphs(gMassAllSM[iprod],gMassAllSM[0]);
      gMassRatAllSM[iprod]->SetMarkerColor(colorProd[iprod]);
      gMassRatAllSM[iprod]->SetLineColor  (colorProd[iprod]);
      gMassRatAllSM[iprod]->SetMarkerStyle(styleProd[iprod]);
      gMassRatAllSM[iprod]->SetMarkerSize (markerSize);
      
      gMassRatAll[iprod] = DivideGraphs(gMassAll[iprod],gMassAll[0]);
      gMassRatAll[iprod]->SetMarkerColor(colorProd[iprod]);
      gMassRatAll[iprod]->SetLineColor  (colorProd[iprod]);
      gMassRatAll[iprod]->SetMarkerStyle(styleProd[iprod]);
      gMassRatAll[iprod]->SetMarkerSize (markerSize);
    }
    
    for(Int_t ism = firstSM; ism <= lastSM; ism++)
    {    
      gMass[ism][iprod] = (TGraphErrors*) file[iprod]->Get(Form("gMass_SM%d",ism));
      //printf("\t gMass_SM%d %p\n",ism,gMass[ism][iprod]);

      if(!gMass[ism][iprod] ) continue;
      
      gMass[ism][iprod]->SetMarkerColor(colorProd[iprod]);
      gMass[ism][iprod]->SetLineColor  (colorProd[iprod]);
      gMass[ism][iprod]->SetMarkerStyle(styleProd[iprod]);
      gMass[ism][iprod]->SetMarkerSize (markerSize);
      
      if ( iprod > 0 && gMass[ism][0] ) 
      {
        gMassRat[ism][iprod] = DivideGraphs(gMass[ism][iprod],gMass[ism][0]);
        gMassRat[ism][iprod]->SetMarkerColor(colorProd[iprod]);
        gMassRat[ism][iprod]->SetLineColor  (colorProd[iprod]);
        gMassRat[ism][iprod]->SetMarkerStyle(styleProd[iprod]);
        gMassRat[ism][iprod]->SetMarkerSize (markerSize);
      }
      
    } // ism

  } // iprod
  
  //===============================
  // PLOTS
  //===============================
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleFontSize(0.06);
  //  gStyle->SetStatFontSize(0.06);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetPadTopMargin(0.07);
  
  //-----------------------------------------------------
  // Plot per SM the Mass distributions and their ratios
  //-----------------------------------------------------
  TString fileName;
  
  TLegend *lE = new TLegend(-0.04,0.1,1,1);
  lE->SetFillColor(0);
  lE->SetFillStyle(0);
  lE->SetLineColor(0);
  lE->SetBorderSize(0);
  lE->SetTextSize(0.07);  
  if      ( trigger == "default" ) lE->SetHeader("     kINT7");
  else if ( trigger == "EMCAL_L0") lE->SetHeader("     kEMC7");
  else if ( trigger == "DCAL_L0" ) lE->SetHeader("     kDMC7"); 
  else if ( trigger == "EMCAL_L1") lE->SetHeader("     kEMCEGA: EG1");
  else if ( trigger == "DCAL_L1" ) lE->SetHeader("     kEMCEGA: DG1"); 
  else if ( trigger == "EMCAL_L2") lE->SetHeader("     kEMCEGA: EG2");
  else if ( trigger == "DCAL_L2" ) lE->SetHeader("     kEMCEGA: DG2");
  
  // Axis E range
  Float_t minE = 1.5;
  Float_t maxE = 12;
  if(particle=="Eta")
  {
    minE = 2;
    maxE = 20;    
  }
  
  TCanvas * cGraphSM = new TCanvas("cGraphSM","Per SM",
                                   ncol*2000,nrow*2000);
  
  cGraphSM->Divide(ncol,nrow);

  for(Int_t ism = firstSM; ism <= lastSM; ism++)
  {  
    cGraphSM->cd(ism+1-firstSM);
    
    //gPad->SetGridx();
    gPad->SetGridy();
 
    hAxis[ism] = new TH1F(Form("hAxisBisRat_SM%d",ism),Form("SM %d",ism),1000,minE,maxE);
    
    hAxis[ism]->GetYaxis()->SetTitle("Mass (MeV/#it{c}^{2})");
    hAxis[ism]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hAxis[ism]->SetMaximum(145);
    hAxis[ism]->SetMinimum(125);    
//    hAxis[ism]->SetMaximum(137);
//    hAxis[ism]->SetMinimum(123);  
    
//    hAxis[ism]->SetMaximum(140);
//    hAxis[ism]->SetMinimum(128);  
    
    if(particle=="Eta")
    {
      hAxis[ism]->SetMaximum(600);
      hAxis[ism]->SetMinimum(450);    
    }
    
    hAxis[ism]->Draw("");
    
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      if ( !file[iprod] || !gMass[ism][iprod] ) continue;
      
      if ( ism == firstSM ) 
       lE->AddEntry(gMass[ism][iprod],prodLeg[iprod],"PL");
      
      gMass[ism][iprod]->Draw("P");
    }
    
  }
  cGraphSM->cd(lastSM-firstSM+2);
  lE->Draw();
  
  fileName = Form("%s/Mass_PerSM",
                  directory.Data());
  
  fileName+=fileFormat;
  cGraphSM->Print(fileName);
  
  // Ratio
  //
  TCanvas * cGraphSMRat = new TCanvas("cGraphSMRat","Ratio Per SM",
                                   ncol*2000,nrow*2000);
  
  cGraphSMRat->Divide(ncol,nrow);
  
  for(Int_t ism = firstSM; ism <= lastSM; ism++)
  {  
    cGraphSMRat->cd(ism+1-firstSM);
    
    //gPad->SetGridx();
    gPad->SetGridy();
    
    hAxis[ism] = new TH1F(Form("hAxisBis_SM%d",ism),Form("SM %d",ism),1000,minE,maxE);
    
    hAxis[ism]->GetYaxis()->SetTitle("Ratio to LHC17pq");
    hAxis[ism]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
//    hAxis[ism]->SetMaximum(1.002);
//    hAxis[ism]->SetMinimum(0.932);       
//    hAxis[ism]->SetMaximum(1.05);
//    hAxis[ism]->SetMinimum(0.90);
    hAxis[ism]->SetMaximum(1.1);
    hAxis[ism]->SetMinimum(0.9);     
    hAxis[ism]->Draw("");
    
    for(Int_t iprod = 1; iprod < nProd; iprod++)
    {
      if ( !gMassRat[ism][iprod] ) continue;
     
      //printf("%d %p\n",iprod,gMassRat[ism][iprod]);
      gMassRat[ism][iprod]->Draw("P");
    }
    
  }
  
  cGraphSMRat->cd(lastSM-firstSM+2);
  
  lE->Draw();
  
  fileName = Form("%s/MassRatio_PerSM",
                  directory.Data());
  
  fileName+=fileFormat;
  cGraphSMRat->Print(fileName);

  
  //-----------------------------------------------------
  // Plot sum of all SM or group of SM
  // the Mass distributions and their ratios
  //-----------------------------------------------------
  
  ncol = 2;
  nrow = 2;

  TH1F* hAxisAll = new TH1F("hAxisAll","All pairs in same SM",1000,minE,maxE);
  
  hAxisAll->GetYaxis()->SetTitle("Mass (MeV/#it{c}^{2})");
  hAxisAll->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hAxisAll->SetMaximum(140);
  hAxisAll->SetMinimum(125);    

  if(particle=="Eta")
  {
    hAxisAll->SetMaximum(600);
    hAxisAll->SetMinimum(450);    
  }
  
//  hAxisAll->SetMaximum(140);
//  hAxisAll->SetMinimum(115);    
  
  TH1F*  hAxisAllRat = new TH1F("hAxisAllRat","Ratio to LHC17pq CENT, all SM, pairs in same SM",1000,minE,maxE);
  hAxisAllRat->GetYaxis()->SetTitle("Ratio to LHC17pq CENT");
  hAxisAllRat->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
//  hAxisAllRat->SetMaximum(1.002);
//  hAxisAllRat->SetMinimum(0.932);    

//  hAxisAllRat->SetMaximum(1.02);
//  hAxisAllRat->SetMinimum(0.92);    
  hAxisAllRat->SetMaximum(1.05);
  hAxisAllRat->SetMinimum(0.95);    

  
  TCanvas * cGraph = new TCanvas("cGraphAllSameSM","All same SM",
                                 ncol*2000,nrow*2000);
  
  cGraph->Divide(ncol,nrow);
  
  cGraph->cd(1);
  
  //gPad->SetGridx();
  gPad->SetGridy();
  
  hAxisAll->Draw("");
  
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    if(!file[iprod]) continue;
    gMassAllSM[iprod]->Draw("P");
  }
  
  gMassAllSM[1]->Draw("P");

  cGraph->cd(2);

  //gPad->SetGridx();
  gPad->SetGridy();
  
  TLegend *lA = new TLegend(0.2,0.8,0.9,0.9);
  lA->SetFillColor(0);
  lA->SetFillStyle(0);
  lA->SetLineColor(0);
  lA->SetBorderSize(0);
  lA->SetTextSize(0.05);  
  
  hAxisAllRat->Draw("");
  
  for(Int_t iprod = 1; iprod < nProd; iprod++)
  {
    if(!file[iprod]) continue;
    
    gMassRatAllSM[iprod]->Draw("P");
    //gMassRatAllSM[iprod]->Fit("pol0","R","",2,8);
    //TF1* func = gMassRatAllSM[iprod]->GetFunction("pol0");
    //lA->AddEntry(func,Form("Fit %2.3f #pm %2.3f",func->GetParameter(0),func->GetParError(0)),"L");
  }
  gMassRatAllSM[1]->Draw("P");

  lA->Draw();

  cGraph->cd(3);

  lE->Draw();
  
  //cGraph->cd(ncol*nrow);
  //lE->Draw();
  
  fileName = Form("%s/Mass_AllSameSM",
                  directory.Data());
  
  fileName+=fileFormat;
  cGraph->Print(fileName);
 
  
  TCanvas * cGraphA = new TCanvas("cGraphAllSM","All SM",
                                 ncol*2000,nrow*2000);
  
  cGraphA->Divide(ncol,nrow);
  
  cGraphA->cd(1);
  
  //gPad->SetGridx();
  gPad->SetGridy();
  
  hAxisAll->SetTitle("all SM combinations");
  hAxisAll->Draw("");
  
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    if(!file[iprod]) continue;
    gMassAll[iprod]->Draw("P");
  }
  
  gMassAll[1]->Draw("P");
  
  cGraphA->cd(2);
  
  //gPad->SetGridx();
  gPad->SetGridy();
  
  hAxisAllRat->SetTitle("Ratio to LHC17pq CENT, all SM combinations");
  hAxisAllRat->Draw("");
  
  for(Int_t iprod = 1; iprod < nProd; iprod++)
  {
    if(!file[iprod]) continue;
    
    gMassRatAllSM[iprod]->Draw("P");
    //gMassRatAllSM[iprod]->Fit("pol0","R","",2,8);
    //TF1* func = gMassRatAllSM[iprod]->GetFunction("pol0");
    //lA->AddEntry(func,Form("Fit %2.3f #pm %2.3f",func->GetParameter(0),func->GetParError(0)),"L");
  }
  
  gMassRatAllSM[1]->Draw("P");
  
  lA->Draw();
  
  cGraphA->cd(3);
  
  lE->Draw();
  
  //cGraph->cd(ncol*nrow);
  //lE->Draw();
  
  fileName = Form("%s/Mass_AllSM",
                  directory.Data());
  
  fileName+=fileFormat;
  cGraphA->Print(fileName);
}


