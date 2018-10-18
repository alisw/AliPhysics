///
/// \file CompareInvMassGraphsMCvsData.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Plot InvMassFit.C output from different inputs, data and MC, in same canvas 
///
/// Example file to read output from InvMassFit.C, mass, width counts vs pT, per SM or all SM 
/// or per cathegory (w/ w/out TRD when applicable) 
/// executed for different analysis options and plot them in the same canvas.
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
/// \param titleMC   : Simplified string acronym of input MC
/// \param titleData : Simplified string acronym of input data
//------------------------------------------------------------------------------
void CompareInvMassGraphsMCvsData
( 
 TString particle    = "Pi0",
 TString calorimeter = "EMCAL",
 TString trigger     = "default",
 TString titleMC     = "simu_pp_7TeV_MB", 
 TString titleData   = "LHC11cd_INT7"     
)
{    
  // SM
  const Int_t nSM = 10; // 10 SM + 1 sum all
  Int_t ncol = 4;
  Int_t nrow = 3;
  
  TString daLeg = "data";
  if     (titleData .Contains("LHC11cd_EMC7"))
    daLeg = "pp@7 TeV, LHC11c+d EMC7";  
  else if(titleData .Contains("LHC11cd_INT7"))
    daLeg = "pp@7 TeV, LHC11c+d INT7";
  else if(titleData == "LHC12_EMC7")
    daLeg = "pp@8 TeV, LHC12 EMC7";   
  else if(titleData == "LHC12_EMG1")
    daLeg = "pp@8 TeV, LHC12 EGA";  
  else if(titleData == "LHC12_INT7")
    daLeg = "pp@8 TeV, LHC12 INT7";  
  else if(titleData == "LHC17pq_EMC7")
    daLeg = "pp@5 TeV, LHC11p+q EMC7";    
  else if(titleData == "LHC17pq_EMG1")
    daLeg = "pp@5 TeV, LHC11p+q EGA";  
  else if(titleData == "LHC17pq_INT7")
    daLeg = "pp@5 TeV, LHC17p+q INT7";
  else if(titleData == "LHC17pq_DCAL_DMC7")
    daLeg = "pp@5 TeV, LHC11p+q DMC7";    
  else if(titleData == "LHC17pq_DCAL_DMG1")
    daLeg = "pp@5 TeV, LHC11p+q DGA";  
  else if(titleData == "LHC17pq_DCAL_INT7")
    daLeg = "pp@5 TeV, LHC17p+q INT7";
  
  //------------------------------------------------
  // Put the path to the files in prod string array
  // also the explanation title on prodLeg string array
  //------------------------------------------------
  
  const Int_t nProd = 3;
  
  TString prod[] = 
  {
    Form("data_%s"    ,titleData.Data())
  , Form("%s_Mimic0_Scaled"  ,titleMC  .Data()) 
  , Form("%s_Mimic10c_EcellCut_Scaled",titleMC  .Data()) 
  };
  
  TString prodLeg[] = 
  {
      "Data"
    , "MC default"
    , "MC xTalk + #it{E}_{inc+cell}>100 MeV"  
  };

  // Create output directories
  TString processline = Form(".! mkdir -p ComparisonGraphs/%s/%s/%s",calorimeter.Data(),trigger.Data(),particle.Data()) ;
  gROOT->ProcessLine(processline.Data());
  
  // Init  the graph arrays.
  //
  TGraphErrors* gMass[nSM][nProd];
  TGraphErrors* gMassAllSM[nProd];
  TGraphErrors* gMassTRDNo[nProd];
  TGraphErrors* gMassTRDYe[nProd];
  TH1F* hAxis[nSM];
  
  TFile* file[nProd];
  
  TGraphErrors* gMassRat[nSM][nProd];
  TGraphErrors* gMassRatAllSM[nProd];
  TGraphErrors* gMassRatTRDNo[nProd];
  TGraphErrors* gMassRatTRDYe[nProd];
  
  //===============================
  // Open the input files and recover the graphs
  // apply some style settings
  //===============================

  Float_t colorProd[] = { 1,kBlue,kRed,kViolet,kRed,kOrange-2,kCyan};
  Float_t styleProd[] = {20,   24,  24,     24,  25,  25,  25};  
  Float_t markerSize = 8;
  
  titleData.ReplaceAll("/","_");
  
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
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
    
    if ( iprod > 0 && gMassAllSM[0] ) 
    {
      gMassRatAllSM[iprod] = DivideGraphs(gMassAllSM[iprod],gMassAllSM[0]);
      gMassRatAllSM[iprod]->SetMarkerColor(colorProd[iprod]);
      gMassRatAllSM[iprod]->SetLineColor  (colorProd[iprod]);
      gMassRatAllSM[iprod]->SetMarkerStyle(styleProd[iprod]);
      gMassRatAllSM[iprod]->SetMarkerSize (markerSize);
    }
    
    gMassTRDNo[iprod] = (TGraphErrors*) file[iprod]->Get("gMass_SameSMTRDNot");
    //printf("\t gMass_SameSMTRDNot %p\n",gMassTRDNo[iprod]);
    
    if ( gMassTRDNo[iprod] )
    {
      gMassTRDNo[iprod]->SetMarkerColor(colorProd[iprod]);
      gMassTRDNo[iprod]->SetLineColor  (colorProd[iprod]);
      gMassTRDNo[iprod]->SetMarkerStyle(styleProd[iprod]);
      gMassTRDNo[iprod]->SetMarkerSize (markerSize);
      
      if ( iprod > 0 && gMassTRDNo[0] ) 
      {
        gMassRatTRDNo[iprod] = DivideGraphs(gMassTRDNo[iprod],gMassTRDNo[0]);
        gMassRatTRDNo[iprod]->SetMarkerColor(colorProd[iprod]);
        gMassRatTRDNo[iprod]->SetLineColor  (colorProd[iprod]);
        gMassRatTRDNo[iprod]->SetMarkerStyle(styleProd[iprod]);
        gMassRatTRDNo[iprod]->SetMarkerSize (markerSize);
      }
    }
    
    gMassTRDYe[iprod] = (TGraphErrors*) file[iprod]->Get("gMass_SameSMTRDYes");
    //printf("\t gMass_SameSMTRDYes %p\n",gMassTRDYe[iprod]);
    
    if ( gMassTRDYe[iprod] )
    {
      gMassTRDYe[iprod]->SetMarkerColor(colorProd[iprod]);
      gMassTRDYe[iprod]->SetLineColor  (colorProd[iprod]);
      gMassTRDYe[iprod]->SetMarkerStyle(styleProd[iprod]);
      gMassTRDYe[iprod]->SetMarkerSize (markerSize);
      
      if ( iprod > 0 && gMassTRDYe[0] ) 
      {
        gMassRatTRDYe[iprod] = DivideGraphs(gMassTRDYe[iprod],gMassTRDYe[0]);
        gMassRatTRDYe[iprod]->SetMarkerColor(colorProd[iprod]);
        gMassRatTRDYe[iprod]->SetLineColor  (colorProd[iprod]);
        gMassRatTRDYe[iprod]->SetMarkerStyle(styleProd[iprod]);
        gMassRatTRDYe[iprod]->SetMarkerSize (markerSize);
      }
    }
    
    for(Int_t ism = 0; ism < nSM; ism++)
    {    
      gMass[ism][iprod] = (TGraphErrors*) file[iprod]->Get(Form("gMass_SM%d",ism));
      // printf("\t gMass_SM%d %p\n",ism,gMass[ism][iprod]);

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
  lE->SetTextSize(0.05);  
  //lE->SetHeader(Form("    SM %d, %s",ism, isoTitle[iso].Data()));
  
  TCanvas * cGraphSM = new TCanvas("cGraphSM","Per SM",
                                   ncol*2000,nrow*2000);
  
  cGraphSM->Divide(ncol,nrow);

  for(Int_t ism = 0; ism < nSM; ism++)
  {  
    cGraphSM->cd(ism+1);
    
    //gPad->SetGridx();
    gPad->SetGridy();
    
    hAxis[ism] = new TH1F(Form("hAxisBisRat_SM%d",ism),Form("SM %d",ism),1000,1,12);
    
    hAxis[ism]->GetYaxis()->SetTitle("Mass (MeV/#it{c}^{2})");
    hAxis[ism]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hAxis[ism]->SetMaximum(140);
    hAxis[ism]->SetMinimum(115);    
//    hAxis[ism]->SetMaximum(137);
//    hAxis[ism]->SetMinimum(123);  
    
//    hAxis[ism]->SetMaximum(140);
//    hAxis[ism]->SetMinimum(128);  
    
    hAxis[ism]->Draw("");
    
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      if(!file[iprod]) continue;
      
      if ( ism == 0 ) 
       lE->AddEntry(gMass[ism][iprod],prodLeg[iprod],"PL");
      
      gMass[ism][iprod]->Draw("P");
    }
    
  }
  cGraphSM->cd(nSM+1);
  lE->Draw();
  
  fileName = Form("ComparisonGraphs/%s/%s/%s/Mass_%s_%s_PerSM",
                  calorimeter.Data(),trigger.Data(),particle.Data(),titleData.Data(),titleMC.Data());
  
  fileName+=fileFormat;
  cGraphSM->Print(fileName);
  
  // Ratio
  //
  TCanvas * cGraphSMRat = new TCanvas("cGraphSMRat","Ratio Per SM",
                                   ncol*2000,nrow*2000);
  
  cGraphSMRat->Divide(ncol,nrow);
  
  for(Int_t ism = 0; ism < nSM; ism++)
  {  
    cGraphSMRat->cd(ism+1);
    
    //gPad->SetGridx();
    gPad->SetGridy();
    
    hAxis[ism] = new TH1F(Form("hAxisBis_SM%d",ism),Form("SM %d",ism),1000,1,12);
    
    hAxis[ism]->GetYaxis()->SetTitle("Mass MC / Data");
    hAxis[ism]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
//    hAxis[ism]->SetMaximum(1.002);
//    hAxis[ism]->SetMinimum(0.932);       
//    hAxis[ism]->SetMaximum(1.05);
//    hAxis[ism]->SetMinimum(0.90);
    hAxis[ism]->SetMaximum(0.98);
    hAxis[ism]->SetMinimum(0.85);     
    hAxis[ism]->Draw("");
    
    for(Int_t iprod = 1; iprod < nProd; iprod++)
    {
      if ( !gMassRat[ism][iprod] ) continue;
     
      //printf("%d %p\n",iprod,gMassRat[ism][iprod]);
      gMassRat[ism][iprod]->Draw("P");
    }
    
  }
  
  cGraphSMRat->cd(nSM+1);
  
  lE->Draw();
  
  fileName = Form("ComparisonGraphs/%s/%s/%s/MassRatio_%s_%s_PerSM",
                  calorimeter.Data(),trigger.Data(),particle.Data(),titleData.Data(),titleMC.Data());
  
  fileName+=fileFormat;
  cGraphSMRat->Print(fileName);

  
  //-----------------------------------------------------
  // Plot sum of all SM or group of SM
  // the Mass distributions and their ratios
  //-----------------------------------------------------
  
  ncol = 2;
  nrow = 2;
  
  if(gMassTRDYe[0])
  {
    ncol = 4;
    nrow = 2;
  }

  TH1F* hAxisAll = new TH1F("hAxisAll","All SM",1000,1,12);
  
  hAxisAll->GetYaxis()->SetTitle("Mass (MeV/#it{c}^{2})");
  hAxisAll->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hAxisAll->SetMaximum(140);
  hAxisAll->SetMinimum(125);    

//  hAxisAll->SetMaximum(140);
//  hAxisAll->SetMinimum(115);    
  
  TH1F*  hAxisAllRat = new TH1F("hAxisAllRat","MC/Data ratio, All SM",1000,1,12);
  hAxisAllRat->GetYaxis()->SetTitle("Mass Data/MC");
  hAxisAllRat->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
//  hAxisAllRat->SetMaximum(1.002);
//  hAxisAllRat->SetMinimum(0.932);    

//  hAxisAllRat->SetMaximum(1.02);
//  hAxisAllRat->SetMinimum(0.92);    
  hAxisAllRat->SetMaximum(1.05);
  hAxisAllRat->SetMinimum(0.98);    

  
  TCanvas * cGraph = new TCanvas("cGraphAllSM","All SM",
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

  if ( gMassTRDYe[0] ) cGraph->cd(5);
  else                 cGraph->cd(3);
  
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
    gMassRatAllSM[iprod]->Fit("pol0","R","",2,8);
    TF1* func = gMassRatAllSM[iprod]->GetFunction("pol0");
    lA->AddEntry(func,Form("Fit %2.3f #pm %2.3f",func->GetParameter(0),func->GetParError(0)),"L");
  }
  gMassRatAllSM[1]->Draw("P");

  lA->Draw();
  
  if ( gMassTRDNo[0] )
  {    
    cGraph->cd(2);
    
    //gPad->SetGridx();
    gPad->SetGridy();
    
    TH1F* hAxisNot = (TH1F*) hAxisAll->Clone("hAxisNot");
    hAxisNot->SetTitle("SM without TRD");
    hAxisNot->Draw("");
    
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      if(!file[iprod]) continue;
      gMassTRDNo[iprod]->Draw("P");
    }
    gMassTRDNo[1]->Draw("P");

    cGraph->cd(6);
    
    TLegend *lN = new TLegend(0.2,0.8,0.9,0.9);
    lN->SetFillColor(0);
    lN->SetFillStyle(0);
    lN->SetLineColor(0);
    lN->SetBorderSize(0);
    lN->SetTextSize(0.05);  
    
    //gPad->SetGridx();
    gPad->SetGridy();
    
    TH1F* hAxisNotRat = (TH1F*) hAxisAllRat->Clone("hAxisNotRat");
    hAxisNotRat->SetTitle("MC/Data ratio, SM without TRD");
    hAxisNotRat->Draw("");   
    
    for(Int_t iprod = 1; iprod < nProd; iprod++)
    {
      if(!file[iprod]) continue;
      
      gMassRatTRDNo[iprod]->Draw("P");
      gMassRatTRDNo[iprod]->Fit("pol0","R","",2,8);
      TF1* func = gMassRatTRDNo[iprod]->GetFunction("pol0");
      lN->AddEntry(func,Form("Fit %2.3f #pm %2.3f",func->GetParameter(0),func->GetParError(0)),"L");

    }
    gMassRatTRDNo[1]->Draw("P");
    lN->Draw();

  }
  
  if ( gMassTRDYe[0] )
  {
    cGraph->cd(3);
    
    //gPad->SetGridx();
    gPad->SetGridy();
    
    TLegend *lY = new TLegend(0.2,0.8,0.9,0.9);
    lY->SetFillColor(0);
    lY->SetFillStyle(0);
    lY->SetLineColor(0);
    lY->SetBorderSize(0);
    lY->SetTextSize(0.05);  
    
    TH1F* hAxisYes = (TH1F*) hAxisAll->Clone("hAxisYes");
    hAxisYes->SetTitle("SM with TRD");
    hAxisYes->Draw("");
    
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      if(!file[iprod]) continue;
      gMassTRDYe[iprod]->Draw("P");
    }
    
    gMassTRDYe[1]->Draw("P");
    
    cGraph->cd(7);
    
    //gPad->SetGridx();
    gPad->SetGridy();
    
    TH1F* hAxisYesRat = (TH1F*) hAxisAllRat->Clone("hAxisYesRat");
    hAxisYesRat->SetTitle("MC/Data ratio, SM with TRD");
    hAxisYesRat->Draw("");   
    
    for(Int_t iprod = 1; iprod < nProd; iprod++)
    {
      if(!file[iprod]) continue;
      
      gMassRatTRDYe[iprod]->Draw("P");
      gMassRatTRDYe[iprod]->Fit("pol0","R","",2,8);
      TF1* func = gMassRatTRDYe[iprod]->GetFunction("pol0");
      lY->AddEntry(func,Form("Fit %2.3f #pm %2.3f",func->GetParameter(0),func->GetParError(0)),"L");
    }
    gMassRatTRDYe[1]->Draw("P");
    lY->Draw();
  }
  
  if ( gMassTRDYe[0] ) cGraph->cd(4);
  else                 cGraph->cd(2);
  
  lE->Draw();
  
  //cGraph->cd(ncol*nrow);
  //lE->Draw();
  
  fileName = Form("ComparisonGraphs/%s/%s/%s/Mass_%s_%s",
                  calorimeter.Data(),trigger.Data(),particle.Data(),titleData.Data(),titleMC.Data());
  
  fileName+=fileFormat;
  cGraph->Print(fileName);
  
}


