#include <TObjString.h>
#include <TString.h>
#include <TProfile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TFile.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
#include <stdlib.h>
#include <TLegend.h>

Color_t colorPeriod[20]     = { kBlack, kRed+1, kBlue+1, kGreen+2, 807, kAzure+2, kViolet+2, kMagenta+2, kSpring-1, kCyan+2,
                                kGray+1, kRed-6, kBlue-6, kGreen-6, kOrange-7, kAzure-4, kViolet+6, kMagenta-6, kSpring+5, kCyan-6 };
Style_t markerPeriod[20]    = { 20, 24, 21, 25, 29, 30, 33, 27, 34, 28,
                                41, 40, 43, 42, 45, 44, 47, 46, 48, 36};
Style_t linePeriod[20]      = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

Bool_t debug                = kFALSE;

//_____________________________________________________________________________________________________
TString GetFormulaString(TF1* f1) {

  if (!f1 ) return "";

  Double_t xmin, xmax;
  f1->GetRange(xmin, xmax);
  Int_t nPar1                         = f1->GetNpar();
  TString formula1                    = f1->GetExpFormula();
  TString formula2                    = "";

  TObjArray *tempArr  = formula1.Tokenize("+");

  Int_t nCurrParam  = nPar1-1;
  // Put them to the correct variables
  for (Int_t i = 0; i< nPar1; i++){
    TString tempValue       = (TString)((TObjString*)tempArr->At(i))->GetString();
    if(debug) cout << tempValue.Data() << "\t";
    if (TMath::Abs(f1->GetParError(nCurrParam))/TMath::Abs(f1->GetParameter(nCurrParam)) < 2 && f1->GetParameter(nCurrParam) != 0){
      if (i > 0)
        formula2  = formula2+"+"+tempValue;
      else
        formula2  = tempValue;
      if(debug) cout << "taken";
    } else {
      if(debug) cout << "not taken";
    }
    if(debug) cout <<  endl;
    nCurrParam--;
  }

  for (Int_t i = 0; i< nPar1; i++){
    #ifndef __CLING__
      formula2.ReplaceAll(Form("[p%d]",i), Form("%1.1e",f1->GetParameter(i)));
    #else
      formula2.ReplaceAll(Form("[%d]",i), Form("%1.1e",f1->GetParameter(i)));
    #endif
  }
  formula2.ReplaceAll("x*x*x*x","x^{4}");
  formula2.ReplaceAll("x*x*x","x^{3}");
  formula2.ReplaceAll("x*x","x^{2}");
  formula2.ReplaceAll("x","Z_{vtx}");
  formula2.ReplaceAll("*","#upoint ");
  formula2.ReplaceAll("+-","-");
  formula2.ReplaceAll("1.0e+00","1");
  formula2.ReplaceAll("e+0","e");
  formula2.ReplaceAll("e-0","e-");
  for (Int_t i = 0; i < 10; i++) formula2.ReplaceAll(Form("e%d",i),Form("e^{%1d}",i));
  for (Int_t i = 0; i < 10; i++) formula2.ReplaceAll(Form("e-%d",i),Form("e^{-%1d}",i));
  return formula2;
}

//_____________________________________________________________________________________________________
void StyleSettings( TString format = ""){
  //gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptDate(0);   //show day and time
  gStyle->SetOptStat(0);  //show statistic
  gStyle->SetPalette(1,0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTextSize(0.5);
  gStyle->SetLabelSize(0.03,"xyz");
  gStyle->SetLabelOffset(0.002,"xyz");
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetTitleOffset(1,"y");
  gStyle->SetTitleOffset(0.9,"x");
  gStyle->SetCanvasColor(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(1);

  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.09);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.13);


  TGaxis::SetMaxDigits(3);
  gErrorIgnoreLevel=kError;

  if (format.CompareTo("eps") == 0 ||format.CompareTo("pdf") == 0  ) gStyle->SetLineScalePS(1);
}

//__________________________________________________________________________________________________________
void SetStyleTLatex( TLatex* text,
                     Size_t textSize,
                     Width_t lineWidth,
                     Color_t textColor = 1,
                     Font_t textFont = 42,
                     Bool_t kNDC = kTRUE,
                     Short_t align = 11
){
  if (kNDC) {text->SetNDC();}
  text->SetTextFont(textFont);
  text->SetTextColor(textColor);
  text->SetTextSize(textSize);
  text->SetLineWidth(lineWidth);
  text->SetTextAlign(align);
}


//__________________________________________________________________________________________________________
void SetStyleProfile( TProfile* histo,
                    Size_t markerSize,
                    Style_t markerStyle,
                    Color_t lineColor) {
  histo->SetLineWidth(2);
  histo->SetMarkerSize(markerSize);
  histo->SetMarkerStyle(markerStyle);
  histo->SetLineColor(lineColor);
  histo->SetMarkerColor(lineColor);
}


//__________________________________________________________________________________________________________
void SetStyleFit(   TF1* fit,
                    Double_t xRangeStart,
                    Double_t xRangeEnd,
                    Width_t lineWidth,
                    Style_t lineStyle,
                    Color_t lineColor) {
  fit->SetRange(xRangeStart,xRangeEnd);
  fit->SetLineWidth(lineWidth);
  fit->SetLineStyle(lineStyle);
  fit->SetLineColor(lineColor);
}

TLegend *GetAndSetLegend2(  Double_t positionX,
                            Double_t positionY,
                            Double_t positionXRight,
                            Double_t positionYUp,
                            Size_t textSize,
                            Int_t columns               = 1,
                            TString header              = "",
                            Font_t textFont             = 43,
                            Double_t margin             = 0
){

  TLegend *legend = new TLegend(positionX,positionY,positionXRight,positionYUp);
  legend->SetNColumns(columns);
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(textFont);
  legend->SetTextSize(textSize);
  if (margin != 0) legend->SetMargin(margin);
  if (header.CompareTo("")!= 0) legend->SetHeader(header);
  return legend;
}

//________________________________________________________________
void DrawTProfileWithCorrespondingFunction(
  Int_t nInputs,
  TProfile** prof,
  TF1** fit,
  TString* llabels,
  TString xAxisName   ="",
  TString yAxisName   ="",
  TString label       = "",
  TString outputName  = ""
){
  TCanvas* canvas = new TCanvas("canvas1", "canvas1", 600, 450);
  canvas->SetLeftMargin(0.07);
  canvas->SetRightMargin(0.01);
  canvas->SetTopMargin(0.01);
  canvas->SetBottomMargin(0.07);
  canvas->SetFillColor(0);

  TLegend* leg = GetAndSetLegend2(0.09,0.95-1.2*0.035*nInputs,0.8,0.95, 450*0.75*0.035,2,"",43, 0.1);
  Bool_t plottedFirst = 0;
  for (Int_t n = 0; n< nInputs; n++){
    if (prof[n] == NULL || fit[n] == NULL) continue;

    SetStyleProfile(prof[n], 0.7, markerPeriod[n], colorPeriod[n]);
    if (!plottedFirst){
      prof[n]->GetYaxis()->SetTitle("#varepsilon_{vtx}");
      prof[n]->GetXaxis()->SetTitle(xAxisName.Data());
      prof[n]->GetYaxis()->SetRangeUser(prof[n]->GetMinimum(0)*0.9, prof[n]->GetMaximum()*1.2);
      prof[n]->SetTitle("");
      prof[n]->Draw("pe");
      leg->AddEntry(prof[n],llabels[n],"p");
      plottedFirst = 1;
    } else {
      prof[n]->Draw("same,pe");
      leg->AddEntry(prof[n],llabels[n],"p");
    }

    SetStyleFit(fit[n], -10, 10, 7, linePeriod[n], colorPeriod[n]);
    fit[n]->Draw("same");

    leg->AddEntry(fit[n],GetFormulaString(fit[n]),"l");
  }

  leg->Draw();
  TLatex* tlabel = new TLatex(0.95, 0.93, "");
  SetStyleTLatex(tlabel, 0.035, 1, kBlack, 42, kTRUE, 31);
  tlabel->DrawLatex(0.95, 0.925, label.Data());
  tlabel->DrawLatex(0.95, 0.925-1.2*0.035, yAxisName.Data());

  //   TString variableName[5] = {"a","b","c","d","e"};
//   if (fit != NULL){
//     for (Int_t i = 0; i < fit->GetNpar() && i < 5; i++){
//       tlabel->DrawLatex(0.95, 0.90-(i+1)*0.035, Form("%s = %3.4f #pm %3.5f",variableName[i].Data(), fit->GetParameter(i), fit->GetParError(i) ));
//     }
//     tlabel->DrawLatex(0.95, 0.90-(fit->GetNpar()+1)*0.035, Form("#chi^{2}/ndf = %3.4f",fit->GetChisquare()/ fit->GetNDF() ));
//   }
  canvas->SaveAs(outputName.Data());
  delete canvas;
}



//________________________________________________________________
void CompareVtxPositionDifferentPeriods(
                                          TString energy              = "",
                                          TString nameInputFileList   = "",
                                          TString addName             = "",
                                          TString outputFileName      = "ZVtzCorrection.root"
                                        ) {

  StyleSettings("pdf");
  // Function meant to generate calibration OADB
  //
  // --- input : nameInputFile, containing a TTree object

  TString collisionSystem = "";
  if ( energy.CompareTo("pPb_5TeV") == 0 )
      collisionSystem = "p-Pb #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
  else if ( energy.CompareTo("pPb_8TeV") == 0 )
      collisionSystem = "p-Pb #sqrt{#it{s}_{_{NN}}} = 8.16 TeV";
  else if ( energy.CompareTo("XeXe_5TeV") == 0 )
    collisionSystem = "Xe-Xe #sqrt{#it{s}_{_{NN}}} = 5.44 TeV";
  else if ( energy.CompareTo("5TeV") == 0 )
    collisionSystem = "pp #sqrt{#it{s}} = 5.02 TeV";
  else if ( energy.CompareTo("13TeVLowB") == 0 )
    collisionSystem = "pp #sqrt{#it{s}} = 13 TeV, B = 0.2 T";


  TString nameOutputDir = Form("CalibrationQA/CompareDifferentPeriods%s_%s",addName.Data(), energy.Data());
  gSystem->Exec("mkdir -p "+nameOutputDir);

  // read folder and name from file
  ifstream in(nameInputFileList.Data());
  cout<<"Available Cuts:"<<endl;
  cout << nameInputFileList.Data() << endl;
  TString fileNames[20];
  TString periodNames[20];
  TString currFileName, currPeriodName;
  Int_t Number = 0;
  while(!in.eof() && Number < 20){
    in >> currFileName >> currPeriodName;
    currPeriodName.ReplaceAll("_"," ");
    fileNames[Number]     = currFileName;
    periodNames[Number]   = currPeriodName;
    cout<< fileNames[Number]<< "\t" << periodNames[Number]<< endl;
    Number++;
  }
  cout<<"=========================="<<endl;
  Number--;

  TFile* filesInput[20]       = {NULL};
  for (Int_t n = 0; n < Number; n++){
    cout << "reading : " << fileNames[n].Data() << endl;
    filesInput[n]           = new TFile(fileNames[n].Data());
    if (filesInput[n]->IsZombie()) filesInput[n] = NULL;
  }


  TString nameCalibrators[17] = { "V0A", "V0C", "V0M", "V0AEq", "V0CEq",
                                  "V0MEq", "SPDCl0", "SPDCl1", "SPDCl", "RefMultEta5",
                                  "RefMultEta8", "NTracklets", "ZNA", "ZNC", "V0AOnline", "V0COnline",
                                  "V0MOnline",};
  TString axislabelPlot[17]   = { "V0A", "V0C", "V0M", "V0AEq", "V0CEq",
                                  "V0MEq", "CL0", "CL1", "SPD cl.", "mult |#eta| < 0.5",
                                  "mult |#eta| < 0.8", "SPD tracklets", "ZNA", "ZNC", "V0A online", "V0C online", "V0M online",};
  TFile* commonOutputFile     = new TFile(outputFileName.Data(),"UPDATE");

  for (Int_t cal = 0; cal< 17; cal++){
    TProfile* hprofCalibrator[20] = {NULL};
    TF1* fitCalibrator[20]        = {NULL};
    TString currentNameProf       = Form("hprofVtxZvs%s", nameCalibrators[cal].Data());
    TString currentNameFit        = Form("fitVtxZvs%s", nameCalibrators[cal].Data());
    for (Int_t n = 0; n < Number; n++){
      if (filesInput[n] == NULL) continue;

      hprofCalibrator[n]          = (TProfile*)filesInput[n]->Get(currentNameProf.Data());
      fitCalibrator[n]            = (TF1*)filesInput[n]->Get(currentNameFit.Data());
      if (hprofCalibrator[n] != NULL &&  fitCalibrator[n] != NULL){
        Double_t scale        = fitCalibrator[n]->Eval(0);
        hprofCalibrator[n]->Scale(1/scale);
        hprofCalibrator[n]->Fit(fitCalibrator[n],"Q0EMRN","",-10,10);
      } else {
        hprofCalibrator[n]        = NULL;
        fitCalibrator[n]          = NULL;
        continue;
      }
      if (commonOutputFile->Get(periodNames[n].Data()) == NULL){
        commonOutputFile->mkdir(periodNames[n].Data());
      }
      commonOutputFile->cd(periodNames[n].Data());
      fitCalibrator[n]->Write(Form("ZVtxCorr_%s",nameCalibrators[cal].Data()), TObject::kOverwrite);
    }

    DrawTProfileWithCorrespondingFunction(Number, hprofCalibrator, fitCalibrator, periodNames, "Z_{vtx} (cm)", axislabelPlot[cal], collisionSystem, Form("%s/%sZVtzDep.pdf", nameOutputDir.Data(),nameCalibrators[cal].Data()));

  }
  commonOutputFile->Write();

}
