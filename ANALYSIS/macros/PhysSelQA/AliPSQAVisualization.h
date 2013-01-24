/////////////////////////////////////
// Created by: Kevin McDermott     //
// email: kmcderm3@nd.edu          //
// CERN Summer Student 2012        //
// University of Notre Dame du Lac //
//                                 // 
// Revision: 1.0                   //
// Created on: August 6, 2012      //
/////////////////////////////////////

#ifndef ALIPSQAVISUALIZATION_H
#define ALIPSQAVISUALIZATION_H

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include <Riostream.h>
#include "TGaxis.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TSystem.h"
#endif

class AliPSQAVisualization : public TObject
{
 public:
  AliPSQAVisualization();
  virtual ~AliPSQAVisualization();

  void InitializeSelectedPlots(const Char_t * listOfPlots);
  void InitializeColorArray(const Char_t * listOfColors); // called only if fUseColorArray == kTRUE
  
  void PostProcessQA(); // First call of macro
  void ScaleMinAndMax(); // called only if fScaleAuto == kTRUE
  void MakeDir(TString dir); // called if anything is to be saved, first check to make sure directory exists.
  void ImportRunAndFillInfo(const Char_t * listOfRunsAndFills); // import info on run and fill number
  void DrawSelected(Int_t iplot);
  Int_t MatchGoodRunToFillNumber(Int_t runnumber); // match run number to fill number
  Int_t MatchGoodRunToStats(Int_t runnumber); // match run number to it's stats
  void ConvertTGraphErrorsToTH1Ds();
  TString SetDrawSelectedYTitle(TString TGEtitle);
  void SavePDFs(Int_t iplot);
  void DrawSameTriggerOnSameCA(Int_t iplot);
  void DrawOverPlotCA(Int_t iplot);
  void DrawOverPlotNoCA(Int_t iplot);
  void SaveToPDFSeparately(Int_t iplot);
  void SaveOverPlotPDF();
  void SaveOverPlotEPS();

  // Getters and setters from the macro

  void     SetInDirectory(TString indir){fInDirectory = indir;};
  TString  GetInDirectory(){return fInDirectory;};

  void     SetROOTInput(TString RI){fROOTInput = RI;};
  TString  GetROOTInput(){return fROOTInput;};
    
  void     SetSavePDFs(Bool_t pdf){fSavePDFs = pdf;};
  Bool_t   GetSavePDFs(){return fSavePDFs;};
  
  void     SetSaveOverPlotPDF(Bool_t pdf){fSaveOverPlotPDF = pdf;};
  Bool_t   GetSaveOverPlotPDF(){return fSaveOverPlotPDF;};

  void     SetOutPDFName(TString name){fOutPDFName = name;};
  TString  GetOutPDFName(){return fOutPDFName;};

  void     SetSaveOverPlotEPS(Bool_t eps){fSaveOverPlotEPS = eps;};
  Bool_t   GetSaveOverPlotEPS(){return fSaveOverPlotEPS;};
  
  void     SetOutEPSName(TString name){fOutEPSName = name;};
  TString  GetOutEPSName(){return fOutEPSName;};

  void     SetOutDirectory(TString out){fOutDirectory = out;};
  TString  GetOutDirectory(){return fOutDirectory;};  

  void     SetDrawSelected(Bool_t select){fDrawSelected = select;};
  Bool_t   GetDrawSelected(){return fDrawSelected;};

  void     SetUseColorArray(Bool_t ca){fUseColorArray = ca;};
  Bool_t   GetUseColorArray(){return fUseColorArray;};

  void     SetDrawOverPlot(Bool_t op){fDrawOverPlot = op;};
  Bool_t   GetDrawOverPlot(){return fDrawOverPlot;};

  void     SetOverPlotTitle(TString opt){fOverPlotTitle = opt;};
  TString  GetOverPlotTitle(){return fOverPlotTitle;};

  void SetPlotOnSameCanvas(Bool_t plot){fSetPlotOnSameCanvas=plot;}

  // Scaling Options

  void     SetScaleAuto(Bool_t sa){fScaleAuto = sa;};
  Bool_t   GetScaleAuto(){return fScaleAuto;};

  void     SetScaleAutoDivMin(Double_t div_min){fDivMin = div_min;};
  Double_t GetScaleAutoDivMin(){return fDivMin;};

  void     SetScaleAutoMultMax(Double_t mult_max){fMultMax = mult_max;};
  Double_t GetScaleAutoMultMax(){return fMultMax;};

  void     SetScaleManMin(Double_t min){fMinimum = min;};
  Double_t GetScaleManMin(){return fMinimum;};

  void     SetScaleManMax(Double_t max){fMaximum = max;};
  Double_t GetScaleManMax(){return fMaximum;};

 private:

  // Input variables
  TString * fSelectedPlots; // List of plots to be turned from tragphs to th1ds
  Int_t     fNSelectedPlots; // number of plots to be turned from tgraphs to th1ds
  TString   fInDirectory; // input directory of root file, most likely the one just analyzed from AliPSQA
  TString   fROOTInput; // name of root file specified by user
  TString   fRunFillFile;// name of file containing info on various fills
  TFile *   fRootFile; // Root file to be used to produce plots

  // Save DrawSelected
  Bool_t  fDrawSelected;
  Bool_t fSetPlotOnSameCanvas; // plot histos on same canvas

  // Save PDF variables

  Bool_t  fSavePDFs; // if kTRUE, save pdfs of th1ds
  TString fOutDirectory; // save pdfs to this directory
  TString fOutPDFName; // common pdf name, group pdf files together

  // Overplot variables

  Bool_t  fDrawOverPlot; // if kTrue, save overplot of all the th1ds
  TString fOverPlotTitle; // title specified by user for overplot, most likely the name of the trigger mask, channel, and logic
  Bool_t    fSaveOverPlotPDF; //if ktrue, save overplot to pdf file
  TCanvas   fOverCanvas; // Overplot canvas
  TLegend * fOverLegend; // Overplot legend
  Bool_t    fSaveOverPlotEPS; // if ktrue, save overplot eps file
  TString   fOutEPSName; // eps file name for overplot together

  // Overplot color array variables

  Bool_t    fUseColorArray; // if ktrue, use colorarray specified in a list file
  Color_t * fColors; // set in a list file, using enums of colors (just use AliRoot and type in the color, e.g. root [0] kBlue returns (const enum EColor)600)
  Int_t fNColors; // number of colors used

  // Overplot Scaling Variables

  Bool_t   fScaleAuto; // If kTRUE, set scaling automatically with GetMinAndMax()
  Double_t fMultMax; // Auto scaling factor, multiply the max
  Double_t fDivMin; // Auto scaling factor, divide the min

  Double_t fMinimum; // Maximum val in selected tgraphs, used for setting range for overplot
  Double_t fMaximum; // Minimum val in selected tgraphs, used for setting range for overplot

  // Plotting variables needed

  TCanvas * fCanvas; // new array in InitializeDrawSelected with the amount of fNDrawSelected
  TH1D *  fDrawPlot; // new array in InitializeDrawSelected with the amount of fNDrawSelected
  TLine ** fFillSeparationLine; // double array of Tlines - for each plot and for each change of fill number
    
  Int_t *fRunNumbers; // Runs to be used in macro
  Int_t *fFillNumbers; // Fill numbers to be associated to macro
  Int_t *fRawRunNumbers; // raw run numbers - taken from logbook
  Int_t *fRawFillNumbers;// raw fill numbers corresponding to runs - taken from logbook
  Int_t *fRawRunStats; //  raw stats for each run number - taken from logbook

  Int_t *fNDiffFills; // array of integers containing info
  Int_t *fBinArray; // array of histo bins corresponding to different RUNS

  Int_t fNRuns; // Number of Runs used in macro
  Int_t fNRawRuns; // Number of raw Runs in macro

  ClassDef(AliPSQAVisualization, 2)
};

#endif
