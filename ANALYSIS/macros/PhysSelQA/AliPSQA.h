/////////////////////////////////////
// Created by: Kevin McDermott     //
// email: kmcderm3@nd.edu          //
// CERN Summer Student 2012        //
// University of Notre Dame du Lac //
//                                 // 
// Revision: 1.0                   //
// Created on: August 6, 2012      //
/////////////////////////////////////

#ifndef ALIPSQA_H
#define ALIPSQA_H

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TMath.h"
#include "TList.h"
#include <Riostream.h>
#include "TSystem.h"
#include "TGaxis.h"
#include "TStyle.h"
#endif

class AliPSQA : public TObject
{
 public: 
  AliPSQA();
  ~AliPSQA();
  
  // Functions listed in order of their appearance in the macro/class

  // Initializers in constructor

  void InitializeRuns(const Char_t * listOfRuns);
  void InitializeTriggers(const Char_t * listOfTriggerNames, const Char_t * listOfTrigMaskChannel);
  void InitializePlots(const Char_t * listOfPlots);
  void SetQAColumnLabelForPlots(Int_t iplot);
  void SetQAColumnLabelForMatching(Int_t iplot);
  void InitializeRowTrigSet();
  void InitializeFilled();
  void InitializeValidRun();
  
  // In order of calls in macro

  void CacheFiles(); // Potential First call in macro
  void MakeDir(TString dir);

  Bool_t ComputeQAPerRun(); // Second call in macro
  TString GetRootFileName(Int_t irun, TString filename);
  void OpenTxtFile(Int_t irun, ofstream & outfile);
  Bool_t ProcessQA(Int_t itrigclass, Int_t itrigchannel, Int_t itriglogic, Int_t irun, TH2D* hStats, ofstream & outfile);
  Int_t SearchYLabels(Int_t itrigclass, Int_t itrigchannel, Int_t itriglogic, TH2D* hStats);
  void GetTrigMask(TString label, Int_t * trigMask);
  TString DecInttToBinTString(Int_t maskInt);
  void SearchXLabels(Int_t iplot, Int_t * ibinx, TH2D * hStats);
  void GetQAVals(TH2D * hstats, Int_t * ibinx, Int_t ibiny,  Double_t * qavals);
  void PlotQA(Int_t iplot, Int_t itrigclass, Int_t itrigchannel, Int_t itriglogic, Int_t irun, Double_t * qavals);  
  Double_t GetError(Double_t * qavals, Int_t iplot);
  Double_t UncNRatioCorrelated(Double_t * qavals);
  Double_t ErrorRatioNotSubsets(Double_t * qavals);
  void SaveQAValToTxtFile(Int_t iplot, Int_t itrigclass, Int_t itrigchannel, Int_t itriglogic,  Double_t * qavals, ofstream& outfile);
  void SaveQAPlots(); 
  void CleanUpGraphs(TGraphErrors * plot); 

  TList * GetBadFiles(){return fListBadFiles;}; // Potential third call in macro

  void DrawAllTGraphs(); // Potential fourth call in macro

  // Getters and Setters for macro conditions
  
  void    SetSaveCache(Bool_t save){fSaveCache = save;};
  Bool_t  GetSaveCache(){return fSaveCache;};

  void    SetCachePath(TString cpath){fCachePath = cpath;};
  TString GetCachePath(){return fCachePath;};

  void    SetLocalMode(Bool_t mode){fLocalMode = mode;};
  Bool_t  GetLocalMode(){return fLocalMode;};

  void    SetLocalPath(TString lpath){fLocalPath = lpath;};
  TString GetLocalPath(){return fLocalPath;};

  void    SetGridPath(TString gpath){fGridPath = gpath;};
  TString GetGridPath(){return fGridPath;}; 

  void    SetPSInput(TString input){fPSInput = input;};
  TString GetPSInput(){return fPSInput;};

  void    SetROOTInput(TString RI){fROOTInput = RI;};
  TString GetROOTInput(){return fROOTInput;};

  void    SetRootOutDirectory(TString out){fRootOutDirectory = out;};
  TString GetRootOutDirectory(){return fRootOutDirectory;};  

  void    SetTxtOutDirectory(TString out){fTxtOutDirectory = out;};
  TString GetTxtOutDirectory(){return fTxtOutDirectory;}; 

  void    SetOutRootName(TString name){fOutRootName = name;};
  TString GetOutRootName(){return fOutRootName;};

  void    SetSaveQAValToTxtFile(Bool_t savetxt){fSaveQAValToTxtFile = savetxt;};
  Bool_t  GetSaveQAValToTxtFile(){return fSaveQAValToTxtFile;};

  void    SetTxtFileName(TString txt){fTxtFileName = txt;};
  TString GetTxtFileName(){return fTxtFileName;};  

  void    SetSaveQAToROOTFile(Bool_t saverf){fSaveQAToROOTFile = saverf;};
  Bool_t  GetSaveQAToROOTFile(){return fSaveQAToROOTFile;};

  void    SetDrawAllTGraphs(Bool_t drawall){fDrawAllTGraphs = drawall;};
  Bool_t  GetDrawAllTGraphs(){return fDrawAllTGraphs;};

  // Legacy get/set

  void    SetQACycle(Int_t cycle){fQAcycle = cycle;};
  Int_t   GetQACycle(){return fQAcycle;};

 private: 
  TString fROOTInput; // name of root file to be analyzed
  Int_t *fRunNumbers; // Runs to be used in macro
  Int_t fNRuns; // Number of Runs used in macro
  TString * fPlots; // List of plots to be generated
  Int_t fNPlots; // Set number of plots to be used/qa numbers to be calculated per trig combo and run
  TString * fQAColumnsNumeratorP; // Sub strings from fPlots, used for naming plots and the key for the online interface, numerator variable
  TString * fQAColumnsDenominatorP; // Sub strings from fPlots, used for naming plots and the key for the online interface, denominator variable
  TString * fQAColumnsNumeratorM; // Sub strings from fPlots, used for searching of x-axis, numerator variable
  TString * fQAColumnsDenominatorM; // Sub strings from fPlots, used for searching of x-axis, denominator variable

  // Cache Variables

  Bool_t  fSaveCache; // Bool to save files to cache locally
  TString fCachePath; // cache subdirectory to be saved, subdirectory name specified in macro
  TString fCacheDirectory; // full path directory to save files, in the format of Form("%s/%09i/%s",fCachePath.Data(),fRunNumbers[irun],fPSInput.Data()) 
  
  // Grid location variables

  TString fGridPath; // Location on the Grid for files to be processed
  TString fPSInput; // Subdirectory for physics selection files to be used as input for macro
  Int_t fQAcycle; // ??

  // Local run variables

  Bool_t fLocalMode; // Run on the Grid or Locally
  TString fLocalPath; // Location of files to be processed locally
  
  // Variables for storing output

  Bool_t fSaveQAToROOTFile;
  TString fOutRootName;  // Name of macro output root file
  TString fRootOutDirectory; // Path for root output
  
  Bool_t fSaveQAValToTxtFile;
  TString fTxtFileName; // txt file name for each run
  TString fTxtOutDirectory; // Path for txt output files

  TList * fListBadFiles; // List of runs with no event data
  Int_t fNGoodRuns; // Number of Good Files

  Bool_t fDrawAllTGraphs;

  // Trigger Info

  TString * fTrigClasses; // Mapped trigger name from & trig mask
  Int_t fNTrigClasses; // Number of triggers used, given by number of trigger names
  

  /////////////////////////////////QUASI HARD-CODED ////////////////////////////////////
  ////////////////////// fTrigChannels[0] == Regular, fTrigChannels[1] == Fast, done this way in the case more partitions are needed, change the implementation, and otherwise just specify to loop over regular channels or fast channels, 


  TString * fTrigChannels; // Mapped names for fast/regular channels, room for more
  Int_t fNTrigChannels; // Number of trigger channels used
  
  ////////////////////////////////// !!!!!!!!!!!!!! Needs to be automated, hard code for now!!!!!!!!!!


  ///////////////////////////////  ///////////////////////////////  ///////////////////////////////
  static const Int_t fgkNTrigLogic = 18;
  ///////////////////////////////  ///////////////////////////////  ///////////////////////////////


  ////////////////////////////////// !!!!!!!!!!!!!! Needs to be automated, hard code for now!!!!!!!!!!

  TGraphErrors * fGraphs[15][30][2][fgkNTrigLogic]; // Plots to be generated
 
  Bool_t fRowTrigSet[30][2][fgkNTrigLogic]; // Check to make sure a specific combination of & and * is checked for a given run
  Bool_t fValidRun[30][2][fgkNTrigLogic];  // Check to see if the run produces any valid plots
  Bool_t fFilled[15][30][2][fgkNTrigLogic]; // Check to see if the run produces a specific plot

  ClassDef(AliPSQA,1)

};

#endif
