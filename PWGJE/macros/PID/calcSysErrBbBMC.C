#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TPaletteAxis.h"

#include <iostream>
#include <fstream>

Bool_t  savePlots = kTRUE;
Bool_t  writeSys  = kTRUE;
TString saveDir = "files/";



// --------------------------------------------------

void setHistoStyleColor(TH1* hist, Bool_t isSpectrum = kFALSE, Int_t color=2)
{
  hist->SetLineColor(1);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.7);
  hist->SetMarkerColor(color);
  hist->SetTitleSize(0.08);
  
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06); 
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.0); 
  hist->GetYaxis()->SetTitleOffset(isSpectrum ? 0.8 : 0.7); 

  hist->GetXaxis()->SetTitleFont(62); 
  hist->GetYaxis()->SetTitleFont(62); 
  hist->GetXaxis()->SetLabelFont(62); 
  hist->GetYaxis()->SetLabelFont(62); 
}


// --------------------------------------------------


TLegend* createLegend(Int_t type)
{
  TLegend* leg = 0x0;
  if (type == 0)
    leg = new TLegend(0.58,0.62,0.92,0.82);
  else if (type == 1)
    leg = new TLegend(0.55,0.18,0.95,0.38);
  else
    leg = new TLegend(0.48,0.60,0.88,0.80);
  
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->SetMargin(0.2);
  
  return leg;
}


// --------------------------------------------------

TCanvas* createCanvas(TString name, Bool_t isSpectrum = kFALSE)
{
  TCanvas* c = new TCanvas(name.Data(),"",760,420); 
  c->Divide(2,2);
  
  for (Int_t i = 1; i <= 4; i++) {
    c->GetPad(i)->SetLeftMargin(isSpectrum ? 0.1 : 0.081);
    c->GetPad(i)->SetRightMargin(0.01);
    c->GetPad(i)->SetBottomMargin(0.13);
    c->GetPad(i)->SetTopMargin(0.12);
  }
  
  return c;
}

// ---------------------------------------------------

Int_t getColor(Int_t col){
  
  Int_t color = kBlack;

  switch(col)
    {
    case 0 : color = kBlue;    break;
    case 1 : color = kCyan;    break;
    case 2 : color = kGreen;   break;
    case 3 : color = kOrange;  break;
    case 4 : color = kRed;     break;
    case 5 : color = kPink;    break;
    case 6 : color = kMagenta; break;
    case 7 : color = kViolet;  break;
    case 8 : color = kBlue-2;  break;
    case 9 : color = kGreen-2;  break;
    case 10: color = kRed+2;  break;
    }

  return color;
  
}

// --------------------------------------------------------------------------------------------------------

void calcSysErrBbBMC(){

  Double_t jetPtLim[] = {5,10,15,20,30,80}; // nBins+1 entries
  const Int_t nJetPtBins = 5;

  TString strSp[] = {"","_pi","_K","_p","_e","_mu"}; 
  const Int_t nSpecies   = 6;

  TString strTitSp[] = {"h^{+} + h^{-}","#pi^{+} + #pi^{-}","K^{+} + K^{-}","p + #bar{p}","e^{+} + e^{-}","#mu^{+} + #mu^{-}"}; 
  //TString strTitSp[] = {"all charged","#pi^{+}+#pi^{-}","p+#bar{p}","K^{+}+K^{-}","e^{+}+e^{-}","#mu^{+}+#mu^{-}"}; 

  TString strSp_10f6a[] = {"","_piPlusPiMinus","_kPlusKMinus","_ppbar","_ePlusEMinus","_muPlusMuMinus"}; 

  // --

  TH1F* hSysErrBbBPt[nSpecies][nJetPtBins];
  TH1F* hSysErrBbBZ[nSpecies][nJetPtBins];
  TH1F* hSysErrBbBXi[nSpecies][nJetPtBins];


  TH1F* corrFacPt_10f6a[nSpecies][nJetPtBins];
  TH1F* corrFacZ_10f6a[nSpecies][nJetPtBins];
  TH1F* corrFacXi_10f6a[nSpecies][nJetPtBins];

  // Pythia Perugia0 

  TH1F* fh1FFTrackPtGenPrim_Perugia0[nSpecies][nJetPtBins];
  TH1F* fh1FFZGenPrim_Perugia0[nSpecies][nJetPtBins];  
  TH1F* fh1FFXiGenPrim_Perugia0[nSpecies][nJetPtBins];
    
  TH1F* fh1FFTrackPtRecPrim_recPt_100_100_Perugia0[nSpecies][nJetPtBins];
  TH1F* fh1FFZRecPrim_recPt_100_100_Perugia0[nSpecies][nJetPtBins];
  TH1F* fh1FFXiRecPrim_recPt_100_100_Perugia0[nSpecies][nJetPtBins];
 
  TH1F* corrFacPt_100_100_Perugia0[nSpecies][nJetPtBins];
  TH1F* corrFacZ_100_100_Perugia0[nSpecies][nJetPtBins];
  TH1F* corrFacXi_100_100_Perugia0[nSpecies][nJetPtBins];



  // --- FFmode3 

  TH1F* fh1FFTrackPtGenPrim_FFmode3[nSpecies][nJetPtBins];
  TH1F* fh1FFZGenPrim_FFmode3[nSpecies][nJetPtBins];  
  TH1F* fh1FFXiGenPrim_FFmode3[nSpecies][nJetPtBins];
    
  TH1F* fh1FFTrackPtRecPrim_recPt_100_100_FFmode3[nSpecies][nJetPtBins];
  TH1F* fh1FFZRecPrim_recPt_100_100_FFmode3[nSpecies][nJetPtBins];
  TH1F* fh1FFXiRecPrim_recPt_100_100_FFmode3[nSpecies][nJetPtBins];
 
  TH1F* corrFacPt_100_100_FFmode3[nSpecies][nJetPtBins];
  TH1F* corrFacZ_100_100_FFmode3[nSpecies][nJetPtBins];
  TH1F* corrFacXi_100_100_FFmode3[nSpecies][nJetPtBins];

 
  // --- FFmode7 

  TH1F* fh1FFTrackPtGenPrim_FFmode7[nSpecies][nJetPtBins];
  TH1F* fh1FFZGenPrim_FFmode7[nSpecies][nJetPtBins];  
  TH1F* fh1FFXiGenPrim_FFmode7[nSpecies][nJetPtBins];
    
  TH1F* fh1FFTrackPtRecPrim_recPt_100_100_FFmode7[nSpecies][nJetPtBins];
  TH1F* fh1FFZRecPrim_recPt_100_100_FFmode7[nSpecies][nJetPtBins];
  TH1F* fh1FFXiRecPrim_recPt_100_100_FFmode7[nSpecies][nJetPtBins];
 
  TH1F* corrFacPt_100_100_FFmode7[nSpecies][nJetPtBins];
  TH1F* corrFacZ_100_100_FFmode7[nSpecies][nJetPtBins];
  TH1F* corrFacXi_100_100_FFmode7[nSpecies][nJetPtBins];

  // ----

  //TString strInFile10f6a = "files/outCorrections_10f6a_tpcCut.root";
  TString strInFile10f6a = "files/outCorrections_10f6a.root";
  
  TString strResDir       = "files";
  TString strResDirP0     = "files";
  TString strResDirModFF  = "files";

  TString strInFilePythia_100_100_P0(Form("%s/outCorrections_PythiaFastJet_eff100_res100.root",strResDirP0.Data()));
  TString strInFilePythia_100_100_FFmode3(Form("%s/outCorrections_PythiaFastJet_eff100_res100_FFmode3.root",strResDirModFF.Data()));
  TString strInFilePythia_100_100_FFmode7(Form("%s/outCorrections_PythiaFastJet_eff100_res100_FFmode7.root",strResDirModFF.Data()));

  // ---

  TFile f1(strInFile10f6a,"READ");
  
  //gDirectory->cd("tpcCut1");

  for(Int_t sp=0; sp<nSpecies; sp++){
    for(Int_t i=0; i<nJetPtBins; i++){
    
      TString strTitlePt(Form("hBbBCorrPt_%02d_%02d%s",(int)jetPtLim[i],(int)jetPtLim[i+1],strSp_10f6a[sp].Data()));
      TString strTitleZ(Form("hBbBCorrZ_%02d_%02d%s",(int)jetPtLim[i],(int)jetPtLim[i+1],strSp_10f6a[sp].Data()));
      TString strTitleXi(Form("hBbBCorrXi_%02d_%02d%s",(int)jetPtLim[i],(int)jetPtLim[i+1],strSp_10f6a[sp].Data()));

      corrFacPt_10f6a[sp][i] = (TH1F*) gDirectory->Get(strTitlePt);
      corrFacZ_10f6a[sp][i]  = (TH1F*) gDirectory->Get(strTitleZ);
      corrFacXi_10f6a[sp][i] = (TH1F*) gDirectory->Get(strTitleXi);

      corrFacPt_10f6a[sp][i]->SetDirectory(0);
      corrFacZ_10f6a[sp][i]->SetDirectory(0);
      corrFacXi_10f6a[sp][i]->SetDirectory(0);
    }
  }
  
  f1.Close();  

  // ---

  for(Int_t sp=0; sp<nSpecies; sp++){

    TFile f2Gen(strInFilePythia_100_100_P0,"READ");
    gDirectory->cd(Form("Gen%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
    
      // for genLevel doesn't matter which histos we take - eff == 1 
 
      TString strTitlePtGen(Form("fh1FFTrackPtGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZGen(Form("fh1FFZGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiGen(Form("fh1FFXiGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
  
      fh1FFTrackPtGenPrim_Perugia0[sp][i] = (TH1F*) gDirectory->Get(strTitlePtGen);
      fh1FFZGenPrim_Perugia0[sp][i]       = (TH1F*) gDirectory->Get(strTitleZGen);
      fh1FFXiGenPrim_Perugia0[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiGen);
      
      fh1FFTrackPtGenPrim_Perugia0[sp][i]->SetDirectory(0);
      fh1FFZGenPrim_Perugia0[sp][i]->SetDirectory(0);
      fh1FFXiGenPrim_Perugia0[sp][i]->SetDirectory(0);
    }
    
    f2Gen.Close();
  }

  // --
  
  for(Int_t sp=0; sp<nSpecies; sp++){

    TFile f2Rec(strInFilePythia_100_100_P0,"READ");

    gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
 
      TString strTitlePtRec(Form("fh1FFTrackPtRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZRec(Form("fh1FFZRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiRec(Form("fh1FFXiRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
  
      fh1FFTrackPtRecPrim_recPt_100_100_Perugia0[sp][i] = (TH1F*) gDirectory->Get(strTitlePtRec);
      fh1FFZRecPrim_recPt_100_100_Perugia0[sp][i]       = (TH1F*) gDirectory->Get(strTitleZRec);
      fh1FFXiRecPrim_recPt_100_100_Perugia0[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiRec);
      
      fh1FFTrackPtRecPrim_recPt_100_100_Perugia0[sp][i]->SetDirectory(0);
      fh1FFZRecPrim_recPt_100_100_Perugia0[sp][i]->SetDirectory(0);
      fh1FFXiRecPrim_recPt_100_100_Perugia0[sp][i]->SetDirectory(0);
    }
    
    f2Rec.Close();
  }

  // ----
  
  for(Int_t sp=0; sp<nSpecies; sp++){

    TFile f2Gen(strInFilePythia_100_100_FFmode3,"READ");

    gDirectory->cd(Form("Gen%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
    
      // for genLevel doesn't matter which histos we take - eff == 1 

      TString strTitlePtGen(Form("fh1FFTrackPtGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZGen(Form("fh1FFZGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiGen(Form("fh1FFXiGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));

 
      fh1FFTrackPtGenPrim_FFmode3[sp][i] = (TH1F*) gDirectory->Get(strTitlePtGen);
      fh1FFZGenPrim_FFmode3[sp][i]       = (TH1F*) gDirectory->Get(strTitleZGen);
      fh1FFXiGenPrim_FFmode3[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiGen);
      
      fh1FFTrackPtGenPrim_FFmode3[sp][i]->SetDirectory(0);
      fh1FFZGenPrim_FFmode3[sp][i]->SetDirectory(0);
      fh1FFXiGenPrim_FFmode3[sp][i]->SetDirectory(0);
    }
    
    f2Gen.Close();
  }


  // ---

  for(Int_t sp=0; sp<nSpecies; sp++){

    TFile f2Rec(strInFilePythia_100_100_FFmode3,"READ");

    gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
    
      // for genLevel doesn't matter which histos we take - eff == 1 
 
      TString strTitlePtRec(Form("fh1FFTrackPtRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZRec(Form("fh1FFZRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiRec(Form("fh1FFXiRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
  
      fh1FFTrackPtRecPrim_recPt_100_100_FFmode3[sp][i] = (TH1F*) gDirectory->Get(strTitlePtRec);
      fh1FFZRecPrim_recPt_100_100_FFmode3[sp][i]       = (TH1F*) gDirectory->Get(strTitleZRec);
      fh1FFXiRecPrim_recPt_100_100_FFmode3[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiRec);
      
      fh1FFTrackPtRecPrim_recPt_100_100_FFmode3[sp][i]->SetDirectory(0);
      fh1FFZRecPrim_recPt_100_100_FFmode3[sp][i]->SetDirectory(0);
      fh1FFXiRecPrim_recPt_100_100_FFmode3[sp][i]->SetDirectory(0);
    }
    
    f2Rec.Close();
  }

  // ----------------------
  
  for(Int_t sp=0; sp<nSpecies; sp++){

    TFile f3Gen(strInFilePythia_100_100_FFmode7,"READ");

    gDirectory->cd(Form("Gen%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
    
      // for genLevel doesn't matter which histos we take - eff == 1 

      TString strTitlePtGen(Form("fh1FFTrackPtGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZGen(Form("fh1FFZGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiGen(Form("fh1FFXiGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
 
      fh1FFTrackPtGenPrim_FFmode7[sp][i] = (TH1F*) gDirectory->Get(strTitlePtGen);
      fh1FFZGenPrim_FFmode7[sp][i]       = (TH1F*) gDirectory->Get(strTitleZGen);
      fh1FFXiGenPrim_FFmode7[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiGen);
      
      fh1FFTrackPtGenPrim_FFmode7[sp][i]->SetDirectory(0);
      fh1FFZGenPrim_FFmode7[sp][i]->SetDirectory(0);
      fh1FFXiGenPrim_FFmode7[sp][i]->SetDirectory(0);
    }
    
    f3Gen.Close();
  }


  // ---

  for(Int_t sp=0; sp<nSpecies; sp++){

    TFile f3Rec(strInFilePythia_100_100_FFmode7,"READ");

    gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
    
      // for genLevel doesn't matter which histos we take - eff == 1 
 
      TString strTitlePtRec(Form("fh1FFTrackPtRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZRec(Form("fh1FFZRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiRec(Form("fh1FFXiRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
  
      fh1FFTrackPtRecPrim_recPt_100_100_FFmode7[sp][i] = (TH1F*) gDirectory->Get(strTitlePtRec);
      fh1FFZRecPrim_recPt_100_100_FFmode7[sp][i]       = (TH1F*) gDirectory->Get(strTitleZRec);
      fh1FFXiRecPrim_recPt_100_100_FFmode7[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiRec);
      
      fh1FFTrackPtRecPrim_recPt_100_100_FFmode7[sp][i]->SetDirectory(0);
      fh1FFZRecPrim_recPt_100_100_FFmode7[sp][i]->SetDirectory(0);
      fh1FFXiRecPrim_recPt_100_100_FFmode7[sp][i]->SetDirectory(0);
    }
    
    f3Rec.Close();
  }

  // ---
  // corr factors

  for(Int_t sp=0; sp<nSpecies; sp++){
    for(Int_t i=0; i<nJetPtBins; i++){

      TString strTitPtCorr_100_100(Form("corrFacPt%s_100_100_P0_%d",strSp[sp].Data(),i));
      TString strTitZCorr_100_100(Form("corrFacZ%s_100_100_P0_%d",strSp[sp].Data(),i));
      TString strTitXiCorr_100_100(Form("corrFacXi%s_100_100_P0_%d",strSp[sp].Data(),i));

      corrFacPt_100_100_Perugia0[sp][i] = (TH1F*) fh1FFTrackPtGenPrim_Perugia0[sp][i]->Clone(strTitPtCorr_100_100);
      corrFacPt_100_100_Perugia0[sp][i]->Divide(fh1FFTrackPtGenPrim_Perugia0[sp][i],fh1FFTrackPtRecPrim_recPt_100_100_Perugia0[sp][i],1,1,"B");
      corrFacZ_100_100_Perugia0[sp][i] = (TH1F*) fh1FFZGenPrim_Perugia0[sp][i]->Clone(strTitZCorr_100_100);
      corrFacZ_100_100_Perugia0[sp][i]->Divide(fh1FFZGenPrim_Perugia0[sp][i],fh1FFZRecPrim_recPt_100_100_Perugia0[sp][i],1,1,"B");

      corrFacXi_100_100_Perugia0[sp][i] = (TH1F*) fh1FFXiGenPrim_Perugia0[sp][i]->Clone(strTitXiCorr_100_100);
      corrFacXi_100_100_Perugia0[sp][i]->Divide(fh1FFXiGenPrim_Perugia0[sp][i],fh1FFXiRecPrim_recPt_100_100_Perugia0[sp][i],1,1,"B");

      strTitPtCorr_100_100.Form("corrFacPt%s_100_100_FFmode3_%d",strSp[sp].Data(),i);
      strTitZCorr_100_100.Form("corrFacZ%s_100_100_FFmode3_%d",strSp[sp].Data(),i);
      strTitXiCorr_100_100.Form("corrFacXi%s_100_100_FFmode3_%d",strSp[sp].Data(),i);

      corrFacPt_100_100_FFmode3[sp][i] = (TH1F*) fh1FFTrackPtGenPrim_FFmode3[sp][i]->Clone(strTitPtCorr_100_100);
      corrFacPt_100_100_FFmode3[sp][i]->Divide(fh1FFTrackPtGenPrim_FFmode3[sp][i],fh1FFTrackPtRecPrim_recPt_100_100_FFmode3[sp][i],1,1,"B");

      corrFacZ_100_100_FFmode3[sp][i] = (TH1F*) fh1FFZGenPrim_FFmode3[sp][i]->Clone(strTitZCorr_100_100);
      corrFacZ_100_100_FFmode3[sp][i]->Divide(fh1FFZGenPrim_FFmode3[sp][i],fh1FFZRecPrim_recPt_100_100_FFmode3[sp][i],1,1,"B");

      corrFacXi_100_100_FFmode3[sp][i] = (TH1F*) fh1FFXiGenPrim_FFmode3[sp][i]->Clone(strTitXiCorr_100_100);
      corrFacXi_100_100_FFmode3[sp][i]->Divide(fh1FFXiGenPrim_FFmode3[sp][i],fh1FFXiRecPrim_recPt_100_100_FFmode3[sp][i],1,1,"B");

      strTitPtCorr_100_100.Form("corrFacPt%s_100_100_FFmode7_%d",strSp[sp].Data(),i);
      strTitZCorr_100_100.Form("corrFacZ%s_100_100_FFmode7_%d",strSp[sp].Data(),i);
      strTitXiCorr_100_100.Form("corrFacXi%s_100_100_FFmode7_%d",strSp[sp].Data(),i);

      corrFacPt_100_100_FFmode7[sp][i] = (TH1F*) fh1FFTrackPtGenPrim_FFmode7[sp][i]->Clone(strTitPtCorr_100_100);
      corrFacPt_100_100_FFmode7[sp][i]->Divide(fh1FFTrackPtGenPrim_FFmode7[sp][i],fh1FFTrackPtRecPrim_recPt_100_100_FFmode7[sp][i],1,1,"B");

      corrFacZ_100_100_FFmode7[sp][i] = (TH1F*) fh1FFZGenPrim_FFmode7[sp][i]->Clone(strTitZCorr_100_100);
      corrFacZ_100_100_FFmode7[sp][i]->Divide(fh1FFZGenPrim_FFmode7[sp][i],fh1FFZRecPrim_recPt_100_100_FFmode7[sp][i],1,1,"B");
      
      corrFacXi_100_100_FFmode7[sp][i] = (TH1F*) fh1FFXiGenPrim_FFmode7[sp][i]->Clone(strTitXiCorr_100_100);
      corrFacXi_100_100_FFmode7[sp][i]->Divide(fh1FFXiGenPrim_FFmode7[sp][i],fh1FFXiRecPrim_recPt_100_100_FFmode7[sp][i],1,1,"B");
    }
  }
  
  
  
  // ---------------------------------------------------
   // sys err: difference between corr fac Pythia/Herwig
   TH1F* corrFacPtSysBbB[nSpecies][nJetPtBins];
   TH1F* corrFacZSysBbB[nSpecies][nJetPtBins];
   TH1F* corrFacXiSysBbB[nSpecies][nJetPtBins];
  
   if(writeSys){
  
     for(Int_t sp=0; sp<nSpecies; sp++){
       for(Int_t i=0; i<nJetPtBins; i++){ // jet slices
     
   Int_t jetPtLoLim = (int) jetPtLim[i];
   Int_t jetPtUpLim = (int) jetPtLim[i+1];
   
   TString strNameBbBPt(Form("hSysErrBbBPt_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
   TString strNameBbBZ(Form("hSysErrBbBZ_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
   TString strNameBbBXi(Form("hSysErrBbBXi_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
   
   TString strNameCorrSysErrBbBPt(Form("corrFacPtSysBbB_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
   TString strNameCorrSysErrBbBZ(Form("corrFacZSysBbB_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
   TString strNameCorrSysErrBbBXi(Form("corrFacXiSysBbB_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
      
     
   hSysErrBbBPt[sp][i] = (TH1F*) corrFacPt_100_100_Perugia0[sp][i]->Clone(strNameBbBPt); 
   hSysErrBbBZ[sp][i]  = (TH1F*) corrFacZ_100_100_Perugia0[sp][i]->Clone(strNameBbBZ);
   hSysErrBbBXi[sp][i] = (TH1F*) corrFacXi_100_100_Perugia0[sp][i]->Clone(strNameBbBXi);
   /*OLD: WRONG BINNING!!!!
   hSysErrBbBPt[sp][i] = (TH1F*) corrFacPt_10f6a[sp][i]->Clone(strNameBbBPt); 
   hSysErrBbBZ[sp][i]  = (TH1F*) corrFacZ_10f6a[sp][i]->Clone(strNameBbBZ);
   hSysErrBbBXi[sp][i] = (TH1F*) corrFacXi_10f6a[sp][i]->Clone(strNameBbBXi);
   */
     
   hSysErrBbBPt[sp][i]->Reset();
   hSysErrBbBZ[sp][i]->Reset();
   hSysErrBbBXi[sp][i]->Reset();  
   
   
   corrFacPtSysBbB[sp][i] = (TH1F*) corrFacPt_100_100_Perugia0[sp][i]->Clone(strNameCorrSysErrBbBPt);
   corrFacZSysBbB[sp][i]  = (TH1F*) corrFacZ_100_100_Perugia0[sp][i]->Clone(strNameCorrSysErrBbBZ);
   corrFacXiSysBbB[sp][i] = (TH1F*) corrFacXi_100_100_Perugia0[sp][i]->Clone(strNameCorrSysErrBbBXi);
       }
     }
     
  
     for(Int_t sp=0; sp<nSpecies; sp++){
       for(Int_t i=0; i<nJetPtBins; i++){
   for(Int_t bin=1; bin<=corrFacPt_100_100_Perugia0[sp][i]->GetNbinsX(); bin++){
     corrFacPtSysBbB[sp][i]->SetBinError(bin, 0);
     if(fh1FFTrackPtGenPrim_Perugia0[sp][i]->GetBinContent(bin) == 0) continue; 
     
     Double_t corrPtP0       = corrFacPt_100_100_Perugia0[sp][i]->GetBinContent(bin);
     Double_t corrPtFFmode3  = corrFacPt_100_100_FFmode3[sp][i]->GetBinContent(bin);
     Double_t corrPtFFmode7  = corrFacPt_100_100_FFmode7[sp][i]->GetBinContent(bin);
      
     const Int_t arrSize = 3;
     Double_t arr[] = {corrPtP0,corrPtFFmode3,corrPtFFmode7};
     Double_t maxEl = TMath::MaxElement(arrSize,arr);
     Double_t minEl = TMath::MinElement(arrSize,arr);
     
     //cout<<" track pt bin "<<bin<<" corr fastsim Pythia P0 "<<corrPtP0<<" P2011 "<<corrPtP2011<<" Herwig "<<corrPtH<<endl;
     
     Double_t sysErrBbB = 0; 
     if(minEl) sysErrBbB = 0.5*fabs((maxEl/minEl)-1);
     /*
     if(maxEl == corrPtP0)      cout<<" pt maxEl: P0:        "<<maxEl<<endl;
     if(maxEl == corrPtFFmode3) cout<<" pt maxEl: FFmod enh: "<<maxEl<<endl;
     if(maxEl == corrPtFFmode7) cout<<" pt maxEl: FFmod depl:"<<maxEl<<endl;
     
     if(minEl == corrPtP0)      cout<<" pt minEl: P0:        "<<minEl<<endl;
     if(minEl == corrPtFFmode3) cout<<" pt minEl: FFmod enh: "<<minEl<<endl;
     if(minEl == corrPtFFmode7) cout<<" pt minEl: FFmod depl:"<<minEl<<endl;
     */
     hSysErrBbBPt[sp][i]->SetBinContent(bin,sysErrBbB); 
     hSysErrBbBPt[sp][i]->SetBinError(bin,0); 
     
     corrFacPtSysBbB[sp][i]->SetBinError(bin, sysErrBbB);
   }
       }
     }

     // --

     for(Int_t sp=0; sp<nSpecies; sp++){
       for(Int_t i=0; i<nJetPtBins; i++){
   for(Int_t bin=1; bin<=corrFacZ_100_100_Perugia0[sp][i]->GetNbinsX(); bin++){
     corrFacZSysBbB[sp][i]->SetBinError(bin, 0);
     if(fh1FFZGenPrim_Perugia0[sp][i]->GetBinContent(bin) == 0) continue; 
     
     Double_t corrZP0       = corrFacZ_100_100_Perugia0[sp][i]->GetBinContent(bin);
     Double_t corrZFFmode3  = corrFacZ_100_100_FFmode3[sp][i]->GetBinContent(bin);
     Double_t corrZFFmode7  = corrFacZ_100_100_FFmode7[sp][i]->GetBinContent(bin);
      
     const Int_t arrSize = 3;
     Double_t arr[] = {corrZP0,corrZFFmode3,corrZFFmode7};
     Double_t maxEl = TMath::MaxElement(arrSize,arr);
     Double_t minEl = TMath::MinElement(arrSize,arr);
     
     //cout<<" track pt bin "<<bin<<" corr fastsim Pythia P0 "<<corrZP0<<" P2011 "<<corrZP2011<<" Herwig "<<corrZH<<endl;
     
     Double_t sysErrBbB = 0; 
     if(minEl) sysErrBbB = 0.5*fabs((maxEl/minEl)-1);
     
     /*
     if(maxEl == corrZP0)      cout<<" z maxEl: P0:        "<<maxEl<<endl;
     if(maxEl == corrZFFmode3) cout<<" z maxEl: FFmod enh: "<<maxEl<<endl;
     if(maxEl == corrZFFmode7) cout<<" z maxEl: FFmod depl:"<<maxEl<<endl;
     
     if(minEl == corrZP0)      cout<<" z minEl: P0:        "<<minEl<<endl;
     if(minEl == corrZFFmode3) cout<<" z minEl: FFmod enh: "<<minEl<<endl;
     if(minEl == corrZFFmode7) cout<<" z minEl: FFmod depl:"<<minEl<<endl;
     */
     hSysErrBbBZ[sp][i]->SetBinContent(bin,sysErrBbB); 
     hSysErrBbBZ[sp][i]->SetBinError(bin,0); 
     
     corrFacZSysBbB[sp][i]->SetBinError(bin, sysErrBbB);
   }
       }
     }
   
     // --

     for(Int_t sp=0; sp<nSpecies; sp++){
       for(Int_t i=0; i<nJetPtBins; i++){
   for(Int_t bin=1; bin<=corrFacXi_100_100_Perugia0[sp][i]->GetNbinsX(); bin++){
     corrFacXiSysBbB[sp][i]->SetBinError(bin, 0);
     
     if(fh1FFTrackPtGenPrim_Perugia0[sp][i]->GetBinContent(bin) == 0) continue; 
     
     Double_t corrXiP0       = corrFacXi_100_100_Perugia0[sp][i]->GetBinContent(bin);
     Double_t corrXiFFmode3  = corrFacXi_100_100_FFmode3[sp][i]->GetBinContent(bin);
     Double_t corrXiFFmode7  = corrFacXi_100_100_FFmode7[sp][i]->GetBinContent(bin);
     
     const Int_t arrSize = 3;
     Double_t arr[] = {corrXiP0,corrXiFFmode3,corrXiFFmode7};
     Double_t maxEl = TMath::MaxElement(arrSize,arr);
     Double_t minEl = TMath::MinElement(arrSize,arr);
     
     //cout<<" track pt bin "<<bin<<" corr fastsim Pythia P0 "<<corrXiP0<<" P2011 "<<corrXiP2011<<" Herwig "<<corrXiH<<endl;
     
     Double_t sysErrBbB = 0; 
     if(minEl) sysErrBbB = 0.5*fabs((maxEl/minEl)-1);
     
     /*
     if(maxEl == corrXiP0)      cout<<" xi maxEl: P0:        "<<maxEl<<endl;
     if(maxEl == corrXiFFmode3) cout<<" xi maxEl: FFmod enh: "<<maxEl<<endl;
     if(maxEl == corrXiFFmode7) cout<<" xi maxEl: FFmod depl:"<<maxEl<<endl;
     
     if(minEl == corrXiP0)      cout<<" xi minEl: P0:        "<<minEl<<endl;
     if(minEl == corrXiFFmode3) cout<<" xi minEl: FFmod enh: "<<minEl<<endl;
     if(minEl == corrXiFFmode7) cout<<" xi minEl: FFmod depl:"<<minEl<<endl;
     */
     hSysErrBbBXi[sp][i]->SetBinContent(bin,sysErrBbB);
     hSysErrBbBXi[sp][i]->SetBinError(bin,0);
     
     corrFacXiSysBbB[sp][i]->SetBinError(bin, sysErrBbB);
   }
       }
     }
     
    TString fnameBbB = Form("%soutSysErr_BbB.root", saveDir.Data());
    TFile fOutBbB(fnameBbB,"RECREATE");
    
    cout<<" write to file "<<fnameBbB<<endl;
    
    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){
  
  hSysErrBbBPt[sp][i]->Write();
  hSysErrBbBZ[sp][i]->Write();
  hSysErrBbBXi[sp][i]->Write();
      }
    }
    
    fOutBbB.Close();
    
   }
   
   

  // --------------------------------------------------------

  gStyle->SetOptStat(0);
  
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleH(0.06);
  
  Int_t selectBin = 0;

  TCanvas* c1 = createCanvas("c1", kTRUE);//new TCanvas("c1","",600,500); 
  
  TLegend* leg1 = createLegend(0);
  
  const Int_t nSpeciesPlot = 4;
  
   for(Int_t sp=0; sp<nSpeciesPlot; sp++){
     for(Int_t i=0; i<nJetPtBins; i++){

       if(i != selectBin) continue;
	   
       c1->cd(sp+1);

       TString strPlotTitleXi(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
       
       fh1FFXiGenPrim_Perugia0[sp][i]->SetTitle(strPlotTitleXi);
       fh1FFXiGenPrim_Perugia0[sp][i]->SetXTitle("#it{#xi}");
       fh1FFXiGenPrim_Perugia0[sp][i]->SetYTitle("1/#it{N}_{Jets} d#it{N}/d#it{#xi}");
       setHistoStyleColor(fh1FFXiGenPrim_Perugia0[sp][i],kTRUE,2);
       fh1FFXiGenPrim_Perugia0[sp][i]->GetXaxis()->SetRangeUser(0,5.2);
       fh1FFXiGenPrim_Perugia0[sp][i]->GetYaxis()->SetRangeUser(0,2.5);
       if(sp>=2) fh1FFXiGenPrim_Perugia0[sp][i]->GetYaxis()->SetRangeUser(0,0.4);
       fh1FFXiGenPrim_Perugia0[sp][i]->DrawCopy();  

       setHistoStyleColor(fh1FFXiGenPrim_FFmode3[sp][i],kTRUE,51);
       fh1FFXiGenPrim_FFmode3[sp][i]->DrawCopy("same");

       setHistoStyleColor(fh1FFXiGenPrim_FFmode7[sp][i],kTRUE,8);
       fh1FFXiGenPrim_FFmode7[sp][i]->DrawCopy("same");
       
       fh1FFXiGenPrim_Perugia0[sp][i]->DrawCopy("same"); //redraw on top   
           
       if(sp==0){
	 leg1->AddEntry(fh1FFXiGenPrim_Perugia0[sp][i],"Perugia-0, gen level","P");
	 leg1->AddEntry(fh1FFXiGenPrim_FFmode3[sp][i],"Low-pt enhanced","P");
	 leg1->AddEntry(fh1FFXiGenPrim_FFmode7[sp][i],"Low-pt depleted","P");
	 
	 leg1->Draw();
       }
       
       gPad->RedrawAxis("");
       gPad->RedrawAxis("G");
     }
   }

   // ---

   TCanvas* c2 = createCanvas("c2", kTRUE);//new TCanvas("c2","",600,500); 
   
   TLegend* leg2 = createLegend(0);

   for(Int_t sp=0; sp<nSpeciesPlot; sp++){
     for(Int_t i=0; i<nJetPtBins; i++){

       if(i != selectBin) continue;
	   
       c2->cd(sp+1);

       gPad->SetLogy();
       
       TString strPlotTitleTrackPt(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
       
       fh1FFTrackPtGenPrim_Perugia0[sp][i]->SetTitle(strPlotTitleTrackPt);
       fh1FFTrackPtGenPrim_Perugia0[sp][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
       fh1FFTrackPtGenPrim_Perugia0[sp][i]->SetYTitle("1/#it{N}_{Jets} d#it{N}/#it{p}_{T}");

       setHistoStyleColor(fh1FFTrackPtGenPrim_Perugia0[sp][i],kTRUE,2);
       fh1FFTrackPtGenPrim_Perugia0[sp][i]->GetXaxis()->SetRangeUser(0,jetPtLim[i+1]);
       fh1FFTrackPtGenPrim_Perugia0[sp][i]->GetYaxis()->SetRangeUser(9e-04,9);
       if(sp>=2) fh1FFTrackPtGenPrim_Perugia0[sp][i]->GetYaxis()->SetRangeUser(5e-05,9e-1);
       fh1FFTrackPtGenPrim_Perugia0[sp][i]->DrawCopy();  

       setHistoStyleColor(fh1FFTrackPtGenPrim_FFmode3[sp][i],kTRUE,51);
       fh1FFTrackPtGenPrim_FFmode3[sp][i]->DrawCopy("same");

       setHistoStyleColor(fh1FFTrackPtGenPrim_FFmode7[sp][i],kTRUE,8);
       fh1FFTrackPtGenPrim_FFmode7[sp][i]->DrawCopy("same");
       
       fh1FFTrackPtGenPrim_Perugia0[sp][i]->DrawCopy("same"); //redraw on top   
           
       if(sp==0){
	 leg2->AddEntry(fh1FFTrackPtGenPrim_Perugia0[sp][i],"Perugia-0, gen level","P");
	 leg2->AddEntry(fh1FFTrackPtGenPrim_FFmode3[sp][i],"Low-pt enhanced","P");
	 leg2->AddEntry(fh1FFTrackPtGenPrim_FFmode7[sp][i],"Low-pt depleted","P");
	 
	 leg2->Draw();
       }
       
       gPad->RedrawAxis("");
       gPad->RedrawAxis("G");
     }
   }

   // ---

   TCanvas* c3 = createCanvas("c3", kTRUE);//new TCanvas("c3","",600,500); 
   
   TLegend* leg3 = createLegend(0);

   for(Int_t sp=0; sp<nSpeciesPlot; sp++){
     for(Int_t i=0; i<nJetPtBins; i++){

       if(i != selectBin) continue;
	   
       c3->cd(sp+1);

       gPad->SetLogy();
       
       TString strPlotTitleZ(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
       
       fh1FFZGenPrim_Perugia0[sp][i]->SetTitle(strPlotTitleZ);
       fh1FFZGenPrim_Perugia0[sp][i]->SetXTitle("#it{z}");
       fh1FFZGenPrim_Perugia0[sp][i]->SetYTitle("1/#it{N}_{Jets} d#it{N}/d#it{z}");

       setHistoStyleColor(fh1FFZGenPrim_Perugia0[sp][i],kTRUE,2);
       fh1FFZGenPrim_Perugia0[sp][i]->GetXaxis()->SetRangeUser(0,1.2);
       fh1FFZGenPrim_Perugia0[sp][i]->GetYaxis()->SetRangeUser(1e-01,40);
       if (sp >= 2) fh1FFZGenPrim_Perugia0[sp][i]->GetYaxis()->SetRangeUser(1e-02,3);
       fh1FFZGenPrim_Perugia0[sp][i]->DrawCopy();  

       setHistoStyleColor(fh1FFZGenPrim_FFmode3[sp][i],kTRUE,51);
       fh1FFZGenPrim_FFmode3[sp][i]->DrawCopy("same");

       setHistoStyleColor(fh1FFZGenPrim_FFmode7[sp][i],kTRUE,8);
       fh1FFZGenPrim_FFmode7[sp][i]->DrawCopy("same");
       
       fh1FFZGenPrim_Perugia0[sp][i]->DrawCopy("same"); //redraw on top   
           
       if(sp==0){
	 leg3->AddEntry(fh1FFZGenPrim_Perugia0[sp][i],"Perugia-0, gen level","P");
	 leg3->AddEntry(fh1FFZGenPrim_FFmode3[sp][i],"Low-pt enhanced","P");
	 leg3->AddEntry(fh1FFZGenPrim_FFmode7[sp][i],"Low-pt depleted","P");
	 
	 leg3->Draw();
       }
       
       gPad->RedrawAxis("");
       gPad->RedrawAxis("G");
     }
   }


   
   // ------------------------------
   // comp corr fac

   TCanvas* c5 = createCanvas("c5");//new TCanvas("c5","",760,420); 
   
   TLegend* leg5 = createLegend(1);
   
   for(Int_t sp=0; sp<nSpeciesPlot; sp++){
     for(Int_t i=0; i<nJetPtBins-1; i++){

       if(i != selectBin) continue;

       c5->cd(sp+1);

       TString strPlotTitleXi(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));

       corrFacXi_100_100_Perugia0[sp][i]->SetTitle(strPlotTitleXi);
       corrFacXi_100_100_Perugia0[sp][i]->SetXTitle("#it{#xi}");
       corrFacXi_100_100_Perugia0[sp][i]->SetYTitle("Correction Factor");
       setHistoStyleColor(corrFacXi_100_100_Perugia0[sp][i],kFALSE,4);
       corrFacXi_100_100_Perugia0[sp][i]->GetXaxis()->SetRangeUser(0,6);
       corrFacXi_100_100_Perugia0[sp][i]->GetYaxis()->SetRangeUser(0.61,2.2);
       corrFacXi_100_100_Perugia0[sp][i]->SetLineColor(4); 
       corrFacXi_100_100_Perugia0[sp][i]->DrawCopy();
       
       setHistoStyleColor(corrFacXiSysBbB[sp][i],kFALSE,4);
       corrFacXiSysBbB[sp][i]->SetFillColor(7);
       corrFacXiSysBbB[sp][i]->DrawCopy("same,E2");

       setHistoStyleColor(corrFacXi_10f6a[sp][i],kFALSE,2);
       corrFacXi_10f6a[sp][i]->SetMarkerStyle(24);
       corrFacXi_10f6a[sp][i]->DrawCopy("same");
       
       if(sp==0){
        leg5->AddEntry(corrFacXi_10f6a[sp][i],"Full simulation","P");
        leg5->AddEntry(corrFacXi_100_100_Perugia0[sp][i],"Fast simulation","P");
        leg5->AddEntry(corrFacXiSysBbB[sp][i],"Low-pt enhanced/depleted","F");
        leg5->Draw();
       }
       /*
       setHistoStyleColor(corrFacXi_100_100_FFmode3[sp][i],kFALSE,51);
       corrFacXi_100_100_FFmode3[sp][i]->SetMarkerStyle(21);
       corrFacXi_100_100_FFmode3[sp][i]->DrawCopy("same");

       setHistoStyleColor(corrFacXi_100_100_FFmode7[sp][i],kFALSE,8);
       corrFacXi_100_100_FFmode7[sp][i]->DrawCopy("same");
       
       setHistoStyleColor(corrFacXi_10f6a[sp][i],kFALSE,2);
       corrFacXi_10f6a[sp][i]->SetMarkerStyle(24);
       corrFacXi_10f6a[sp][i]->DrawCopy("same");
       
       corrFacXi_100_100_Perugia0[sp][i]->DrawCopy("same"); //redraw on top 

       if(sp==0){
	 leg5->AddEntry(corrFacXi_10f6a[sp][i],"Full simulation","P");
	 leg5->AddEntry(corrFacXi_100_100_Perugia0[sp][i],"Fast sim Perugia0","P");
	 leg5->AddEntry(fh1FFTrackPtGenPrim_FFmode3[sp][i],"Low-pt enhanced","P");
	 leg5->AddEntry(fh1FFTrackPtGenPrim_FFmode7[sp][i],"Low-pt depleted","P");
	 leg5->Draw();
       }*/
       
       gPad->RedrawAxis("");
       gPad->RedrawAxis("G");
     }
   }


   // ---

   TCanvas* c6 = createCanvas("c6");//new TCanvas("c6","",760,420); 
   
   TLegend* leg6 = createLegend(2);

   for(Int_t sp=0; sp<nSpeciesPlot; sp++){
     for(Int_t i=0; i<nJetPtBins-1; i++){

       if(i != selectBin) continue;

       c6->cd(sp+1);
     
       TString strPlotTitlePt(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));

       corrFacPt_100_100_Perugia0[sp][i]->SetTitle(strPlotTitlePt);
       corrFacPt_100_100_Perugia0[sp][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
       corrFacPt_100_100_Perugia0[sp][i]->SetYTitle("Correction Factor");
       setHistoStyleColor(corrFacPt_100_100_Perugia0[sp][i],kFALSE,4);
       corrFacPt_100_100_Perugia0[sp][i]->GetXaxis()->SetRangeUser(0,jetPtLim[i+1]);
       corrFacPt_100_100_Perugia0[sp][i]->GetYaxis()->SetRangeUser(0.61,2.2);
       corrFacPt_100_100_Perugia0[sp][i]->SetLineColor(4); 
       corrFacPt_100_100_Perugia0[sp][i]->DrawCopy();  
       
       setHistoStyleColor(corrFacPtSysBbB[sp][i],kFALSE,4);
       corrFacPtSysBbB[sp][i]->SetFillColor(7);
       corrFacPtSysBbB[sp][i]->DrawCopy("same,E2");

       setHistoStyleColor(corrFacPt_10f6a[sp][i],kFALSE,2);
       corrFacPt_10f6a[sp][i]->SetMarkerStyle(24);
       corrFacPt_10f6a[sp][i]->DrawCopy("same");
       
       if(sp==0){
        leg6->AddEntry(corrFacPt_10f6a[sp][i],"Full simulation","P");
        leg6->AddEntry(corrFacPt_100_100_Perugia0[sp][i],"Fast simulation","P");
        leg6->AddEntry(corrFacPtSysBbB[sp][i],"Low-pt enhanced/depleted","F");
        leg6->Draw();
       }
       /*
       setHistoStyleColor(corrFacPt_100_100_FFmode3[sp][i],kFALSE,51);
       corrFacPt_100_100_FFmode3[sp][i]->SetMarkerStyle(21);
       corrFacPt_100_100_FFmode3[sp][i]->DrawCopy("same");
     
       setHistoStyleColor(corrFacPt_100_100_FFmode7[sp][i],kFALSE,8);
       corrFacPt_100_100_FFmode7[sp][i]->DrawCopy("same");
       
       setHistoStyleColor(corrFacPt_10f6a[sp][i],kFALSE,2);
       corrFacPt_10f6a[sp][i]->SetMarkerStyle(24);
       corrFacPt_10f6a[sp][i]->DrawCopy("same");
       
       corrFacPt_100_100_Perugia0[sp][i]->DrawCopy("same"); //redraw on top   
     
       if(sp==0){
	 leg6->AddEntry(corrFacPt_10f6a[sp][i],"Full simulation","P");
	 leg6->AddEntry(corrFacPt_100_100_Perugia0[sp][i],"Fast sim Perugia0","P");
	 leg6->AddEntry(corrFacPt_100_100_FFmode3[sp][i],"Low-pt enhanced","P");
	 leg6->AddEntry(corrFacPt_100_100_FFmode7[sp][i],"Low-pt depleted","P");
	 
	 leg6->Draw();
       }*/
       
       gPad->RedrawAxis("");
       gPad->RedrawAxis("G");
     }
   }

   // ---

   TCanvas* c7 = createCanvas("c7");//new TCanvas("c7","",760,420); 
   
   TLegend* leg7 = createLegend(2);
   
   for(Int_t sp=0; sp<nSpeciesPlot; sp++){
     for(Int_t i=0; i<nJetPtBins-1; i++){
       
       if(i != selectBin) continue;
       
       c7->cd(sp+1);

       TString strPlotTitleZ(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));

       corrFacZ_100_100_Perugia0[sp][i]->SetTitle(strPlotTitleZ);
       corrFacZ_100_100_Perugia0[sp][i]->SetXTitle("#it{z}");
       corrFacZ_100_100_Perugia0[sp][i]->SetYTitle("Correction Factor");
       setHistoStyleColor(corrFacZ_100_100_Perugia0[sp][i],kFALSE,4);
       corrFacZ_100_100_Perugia0[sp][i]->GetXaxis()->SetRangeUser(0,1.0);
       corrFacZ_100_100_Perugia0[sp][i]->GetYaxis()->SetRangeUser(0.61,2.2);
       corrFacZ_100_100_Perugia0[sp][i]->SetLineColor(4); 
       corrFacZ_100_100_Perugia0[sp][i]->DrawCopy();
       
       setHistoStyleColor(corrFacZSysBbB[sp][i],kFALSE,4);
       corrFacZSysBbB[sp][i]->SetFillColor(7);
       corrFacZSysBbB[sp][i]->DrawCopy("same,E2");

       setHistoStyleColor(corrFacZ_10f6a[sp][i],kFALSE,2);
       corrFacZ_10f6a[sp][i]->SetMarkerStyle(24);
       corrFacZ_10f6a[sp][i]->DrawCopy("same");
       
       if(sp==0){
        leg7->AddEntry(corrFacPt_10f6a[sp][i],"Full simulation","P");
        leg7->AddEntry(corrFacPt_100_100_Perugia0[sp][i],"Fast simulation","P");
        leg7->AddEntry(corrFacPtSysBbB[sp][i],"Low-pt enhanced/depleted","F");
        leg7->Draw();
       }
       /*
       setHistoStyleColor(corrFacZ_100_100_FFmode3[sp][i],kFALSE,51);
       corrFacZ_100_100_FFmode3[sp][i]->SetMarkerStyle(21);
       corrFacZ_100_100_FFmode3[sp][i]->DrawCopy("same");

       setHistoStyleColor(corrFacZ_100_100_FFmode7[sp][i],kFALSE,8);
       corrFacZ_100_100_FFmode7[sp][i]->DrawCopy("same");
       
       setHistoStyleColor(corrFacZ_10f6a[sp][i],kFALSE,2);
       corrFacZ_10f6a[sp][i]->SetMarkerStyle(24);
       corrFacZ_10f6a[sp][i]->DrawCopy("same");

       corrFacZ_100_100_Perugia0[sp][i]->DrawCopy("same"); // redraw on top 

       if(sp==0){
	 leg7->AddEntry(corrFacZ_10f6a[sp][i],"Full simulation","P");
	 leg7->AddEntry(corrFacZ_100_100_Perugia0[sp][i],"Fast sim Perugia0","P");
	 leg7->AddEntry(fh1FFTrackPtGenPrim_FFmode3[sp][i],"Low-pt enhanced","P");
	 leg7->AddEntry(fh1FFTrackPtGenPrim_FFmode7[sp][i],"Low-pt depleted","P");
	 leg7->Draw();
       }*/
       
       gPad->RedrawAxis("");
       gPad->RedrawAxis("G");
     }
   }


   if(savePlots){
     c1->SaveAs(Form("%sspectraModFF_dNdXi_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
     c2->SaveAs(Form("%sspectraModFF_dNdPt_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
     c3->SaveAs(Form("%sspectraModFF_dNdZ_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
     
     c5->SaveAs(Form("%scorrFacModFF_dNdXi_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
     c6->SaveAs(Form("%scorrFacModFF_dNdPt_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
     c7->SaveAs(Form("%scorrFacModFF_dNdZ_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
   }

}

