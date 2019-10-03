#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TVirtualFitter.h"

#include "AliPID.h"

#include "THnSparseDefinitions.h"

// E.g. a 'extractSingleTrackEfficiencies.C+("finalCuts/MC_pp/7TeV/LHC10f6a/corrected/finalisedSplines/analytical/Jets/nclCut/noMCidForGen", "output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors__2014_01_22_negCharge.root", "_SS_neg", "pdf")'
// or a 'extractSingleTrackEfficiencies.C+("OliversMacros", "output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors__dummy_negCharge.root", "_SS_neg", "pdf)' -l -b -q
Int_t extractSingleTrackEfficiencies(TString path, TString fileName, TString saveAsSuffix = "", TString saveAs = "")
{
  TVirtualFitter::SetMaxIterations(1e5); 
  
  TString filePathName = Form("%s/%s", path.Data(), fileName.Data());
  TFile* f = TFile::Open(filePathName.Data(), "READ");
  if (!f) {
    printf("Failed to open file: %s!\n", filePathName.Data());
    return -1;
  }
  
  TH1D* hSingleTrackEfficiency[AliPID::kSPECIES] = { 0x0, };
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    hSingleTrackEfficiency[species] = (TH1D*)f->Get(Form("hSingleTrackEfficiency_%s", AliPID::ParticleShortName(species)));
    hSingleTrackEfficiency[species]->SetDirectory(0);
  }
  
  
  f->Close();
  
  TString savePathName = Form("%s/output_SingleTrackEfficiency_%s", path.Data(), fileName.Data());
  
  TFile* fSave = TFile::Open(savePathName.Data(), "RECREATE");
  if (!fSave) {
    printf("Failed to open save file: %s!\n", savePathName.Data());
    return -1;
  }
  
  TString saveAsPathName = "";
  
  // --------------
  // ---- plots ---

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    printf("\n\n_______________________*** Fitting species: %s ***_______________________\n", AliPID::ParticleShortName(species));
    TCanvas *canv = new TCanvas(Form("canv_%s", AliPID::ParticleShortName(species)), "", 1362, 671);//600, 400); 
    
    // ---

    canv->SetTickx();
    canv->SetTicky();
    canv->SetLogx();
    canv->SetGrid(0,1);
    
    canv->SetTopMargin(0.02);
    canv->SetRightMargin(0.015);
    canv->SetLeftMargin(0.085);
    canv->SetBottomMargin(0.155);
    
    
    setAxisTitlesItalic(hSingleTrackEfficiency[species]);
    
    hSingleTrackEfficiency[species]->GetXaxis()->SetLabelSize(0.06);
    hSingleTrackEfficiency[species]->GetXaxis()->SetTitleSize(0.07);
    hSingleTrackEfficiency[species]->GetXaxis()->SetTitleOffset(1);
    
    hSingleTrackEfficiency[species]->GetYaxis()->SetLabelSize(0.06);
    hSingleTrackEfficiency[species]->GetYaxis()->SetTitleSize(0.07);
    hSingleTrackEfficiency[species]->GetYaxis()->SetTitleOffset(0.6);
    
    hSingleTrackEfficiency[species]->SetMarkerStyle(20);
    hSingleTrackEfficiency[species]->SetMarkerSize(1.4);
    hSingleTrackEfficiency[species]->DrawCopy("E1 p");
    /*
    TLegend* leg = new TLegend(0.60,0.76,0.96,0.96);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(1);
    
    leg->AddEntry(hSingleTrackEfficiency[species], hSingleTrackEfficiency[species]->GetTitle(), "P");
    leg->Draw();
    */

    // ---

    // fit in subranges, ensure smoothnes by constraint parameters
    Double_t f0_lowlim = 0.15;
    Double_t f0_uplim = 0.325;
    Double_t f1_uplim = 1.3; 
    Double_t f2_uplim = 2.7; 
    Double_t f3_uplim = 3.3;
    Double_t f4_uplim = 50.;
    Double_t f5_uplim = 50.;

    if (species == AliPID::kPion) {
      f0_uplim = 0.275; //TODO NEW added on Feb 12, 2014
      f1_uplim = 1.3;
      f3_uplim = 4.75;
      f4_uplim = 15.5;
    }
    else if (species == AliPID::kKaon) {
      f0_lowlim = 0.125;
      f0_uplim = 0.275;
      f2_uplim = 3.1;
      f3_uplim = 8.5;
      f4_uplim = 15.5;
    }
    else if (species == AliPID::kProton) {
      f1_uplim = 0.675;
      f2_uplim = 1.5;
      f3_uplim = 2.7;
    }
    else if (species == AliPID::kMuon) {
      f2_uplim = 2.5;
    }
    
    TF1* f0 = new TF1(Form("f0_%s", AliPID::ParticleShortName(species)), "[0]+[1]*TMath::Erf([2]*(x - [3]))", f0_lowlim, f0_uplim);
    //TF1* f0 = new TF1(Form("f0_%s", AliPID::ParticleShortName(species)), "TMath::Max(0, [0]+[1]/x/x/x)", f0_lowlim, f0_uplim);
    //TF1* f0 = new TF1(Form("f0_%s", AliPID::ParticleShortName(species)), "TMath::Max(0, [0]+[1]/TMath::Power(x, [2]))",
    //                  f0_lowlim, f0_uplim);
    //TF1* f1 = new TF1(Form("f1_%s", AliPID::ParticleShortName(species)), "[0]+[2]*(x-[1])+[3]*TMath::Power((x-[1]), 2)",
    //                  f0_uplim, f1_uplim);
    TF1* f1 = new TF1(Form("f1_%s", AliPID::ParticleShortName(species)), 
                      "[0]+[2]*(x-[1])+[3]*TMath::Power((x-[1]),2)+[4]*TMath::Power((x-[1]),3)+[5]*TMath::Power((x-[1]),4)+[6]*TMath::Power((x-[1]),5)",
                      f0_uplim, f1_uplim);
    TF1* f2 = new TF1(Form("f2_%s", AliPID::ParticleShortName(species)),
                      "[0]+[2]*(x-[1])+[3]*TMath::Power((x-[1]),2)+[4]*TMath::Power((x-[1]),3)+[5]*TMath::Power((x-[1]),4)",
                      f1_uplim, f2_uplim);
    TF1* f3 = (species != AliPID::kElectron && species != AliPID::kMuon)
                ? new TF1(Form("f3_%s", AliPID::ParticleShortName(species)), "[0]+[2]*(x-[1])+[3]*TMath::Power((x-[1]),2)",
                          f2_uplim, f3_uplim)
                : new TF1(Form("f3_%s", AliPID::ParticleShortName(species)), "[0]+[2]*(x-[1])", f2_uplim, f3_uplim);
    TF1* f4 = new TF1(Form("f4_%s", AliPID::ParticleShortName(species)), "[0]+[2]*(x-[1])", f3_uplim, f4_uplim);
    TF1* f5 = new TF1(Form("f5_%s", AliPID::ParticleShortName(species)), "[0]+[2]*(x-[1])", f4_uplim, f5_uplim);
    
    /*original
    Double_t f1_uplim = 1; 
    Double_t f2_uplim = 6; 

    TF1 *f1     = new TF1("f1","pol3",0,f1_uplim);
    TF1 *f2     = new TF1("f2","[0]+[2]*(x-[1])+[3]*TMath::Power((x-[1]),2)+[4]*TMath::Power((x-[1]),3)+[5]*TMath::Power((x-[1]),4)",f1_uplim,f2_uplim);
    TF1 *f3     = new TF1("f3","[0]+[2]*(x-[1])",f2_uplim,50);
    */
    
    if (species == AliPID::kProton) {
      f0->SetParameter(0, 0.4);
      f0->SetParameter(1, 0.4);
      f0->SetParameter(2, 22);
      f0->SetParameter(3, 0.3);
    }
    else if (species == AliPID::kKaon) {
      f0->SetParameter(0, 0.24);
      f0->SetParameter(1, 0.1);
      f0->SetParameter(2, 17);
      f0->SetParameter(3, 0.24);
    }
    else {
      f0->SetParameter(0, -28);
      f0->SetParameter(1, 28);
      f0->SetParameter(2, 5.3);
      f0->SetParameter(3, -0.18);
      /*
      f0->SetParameter(2, -3);
      f0->SetParLimits(2, 1, 5);
      */
    }
    printf("\n\n____________________________________\nFitting f0:\n");
    hSingleTrackEfficiency[species]->Fit(f0, "R+");

    // f1: constrain to ensure smoothness 
    Double_t val_f1par0 = f0->Eval(f0_uplim);
    f1->SetParameter(0, val_f1par0);
    f1->SetParLimits(0, val_f1par0, val_f1par0);
    
    f1->SetParameter(1, f0_uplim);
    f1->SetParLimits(1, f0_uplim, f0_uplim);
    
    printf("\n\n____________________________________\nFitting f1:\n");
    hSingleTrackEfficiency[species]->Fit(f1, "R+");

    // f2: constrain to ensure smoothness 
    Double_t val_f2par0 = f1->Eval(f1_uplim);
    f2->SetParameter(0, val_f2par0);
    f2->SetParLimits(0, val_f2par0, val_f2par0);
    
    f2->SetParameter(1, f1_uplim);
    f2->SetParLimits(1, f1_uplim, f1_uplim);

    printf("\n\n____________________________________\nFitting f2:\n");
    hSingleTrackEfficiency[species]->Fit(f2, "R+");

    // f3: constrain to ensure smoothness 
    Double_t val_f3par0 = f2->Eval(f2_uplim);
    f3->SetParameter(0, val_f3par0);
    f3->SetParLimits(0, val_f3par0, val_f3par0);
    
    f3->SetParameter(1, f2_uplim);
    f3->SetParLimits(1, f2_uplim, f2_uplim);
    
    printf("\n\n____________________________________\nFitting f3:\n");
    hSingleTrackEfficiency[species]->Fit(f3, "R+");
    
    // f4: constrain to ensure smoothness 
    Double_t val_f4par0 = f3->Eval(f3_uplim);
    f4->SetParameter(0, val_f4par0);
    f4->SetParLimits(0, val_f4par0, val_f4par0);
    
    f4->SetParameter(1, f3_uplim);
    f4->SetParLimits(1, f3_uplim, f3_uplim);
    
    printf("\n\n____________________________________\nFitting f4:\n");
    hSingleTrackEfficiency[species]->Fit(f4, "R+");
    
    
    if (species == AliPID::kPion || species == AliPID::kKaon) {
      // f5: constrain to ensure smoothness 
      Double_t val_f5par0 = f4->Eval(f4_uplim);
      f5->SetParameter(0, val_f5par0);
      f5->SetParLimits(0, val_f5par0, val_f5par0);
      
      f5->SetParameter(1, f4_uplim);
      f5->SetParLimits(1, f4_uplim, f4_uplim);
      
      printf("\n\n____________________________________\nFitting f5:\n");
      hSingleTrackEfficiency[species]->Fit(f5, "R+");
    }
    
    canv->Modified();
    canv->Update();
    
    ClearTitleFromHistoInCanvas(canv);
    
    canv->Update();
    canv->Modified();
    
    fSave->cd();
    hSingleTrackEfficiency[species]->Write();
    
    if (saveAs != "") {
      saveAsPathName = Form("%s/SingleTrackEff%s_%s.%s", path.Data(), saveAsSuffix.Data(), AliPID::ParticleShortName(species), saveAs.Data());
      canv->SaveAs(saveAsPathName.Data());
      saveAsPathName = saveAsPathName.ReplaceAll(Form(".%s", saveAs.Data()), ".root");
      canv->SaveAs(saveAsPathName.Data());
    }
  }
  
  fSave->Close();

  return 0;
}
