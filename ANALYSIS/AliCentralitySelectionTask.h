#ifndef ALICENTRALITYSELECTIONTASK_H
#define ALICENTRALITYSELECTIONTASK_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliCentralitySelectionTask
//   author: Alberica Toia
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TFile;
class TH1F;
class TList;
class TString;

class AliESDEvent;

class AliCentralitySelectionTask : public AliAnalysisTaskSE {

 public:

  AliCentralitySelectionTask();
  AliCentralitySelectionTask(const char *name);
  AliCentralitySelectionTask& operator= (const AliCentralitySelectionTask& ana);
  AliCentralitySelectionTask(const AliCentralitySelectionTask& c);
  virtual ~AliCentralitySelectionTask();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInput(const char* input)         {fAnalysisInput = input;}
  void SetMCInput()                        {fIsMCInput = kTRUE;}
  
  void SetPercentileFile(TString filename);
  void SetPercentileFile2(TString filename);
  void ReadCentralityHistos();
  void ReadCentralityHistos2();

  void AddPercentileFileToList(TString filename) { fFileList->Add(new TObjString(filename)); }
  void AddPercentileFile2ToList(TString filename) { fFileList2->Add(new TObjString(filename)); }

  Float_t GetCorrV0(const AliESDEvent* esd, float &v0CorrResc);
  Float_t GetCorrSPD2(Float_t spd2raw,Float_t zv);

 private:

  Int_t SetupRun(AliESDEvent* esd);

  Int_t    fDebug;	   	// Debug flag
  TString  fAnalysisInput; 	// "ESD", "AOD"
  Bool_t   fIsMCInput;          // true when input is MC
  TFile   *fFile;               // file that holds the centrality vs multiplicity 1d
  TFile   *fFile2;              // file that holds the centrality vs multiplicity 2d
  TString  fCentfilename;       // name of this file 1d
  TString  fCentfilename2;      // name of this file 2d
  
  TList*   fFileList;           //! list of input files names
  TList*   fFileList2;          //! list of input files 2 names
  Int_t    fCurrentRun;         // current run number
  Int_t    fRunNo;              // reference run number

  Float_t  fCentV0M;            // percentile centrality from V0
  Float_t  fCentFMD;            // percentile centrality from FMD
  Float_t  fCentTRK;            // percentile centrality from tracks
  Float_t  fCentTKL;            // percentile centrality from tracklets
  Float_t  fCentCL0;            // percentile centrality from clusters in layer 0
  Float_t  fCentCL1;            // percentile centrality from clusters in layer 0
  Float_t  fCentV0MvsFMD;       // percentile centrality from V0 vs FMD
  Float_t  fCentTKLvsV0M;       // percentile centrality from tracklets vs V0
  Float_t  fCentZEMvsZDC;       // percentile centrality from ZEM vs ZDC

  TH1F    *fHtempV0M;           // histogram with centrality vs multiplicity using V0
  TH1F    *fHtempFMD;           // histogram with centrality vs multiplicity using FMD
  TH1F    *fHtempTRK;           // histogram with centrality vs multiplicity using tracks
  TH1F    *fHtempTKL;           // histogram with centrality vs multiplicity using tracklets
  TH1F    *fHtempCL0;           // histogram with centrality vs multiplicity using clusters in layer 0
  TH1F    *fHtempCL1;           // histogram with centrality vs multiplicity using clusters in layer 0
  TH1F    *fHtempV0MvsFMD;           // histogram with centrality vs multiplicity using V0 vs FMD   
  TH1F    *fHtempTKLvsV0M;           // histogram with centrality vs multiplicity using tracklets vs V0
  TH1F    *fHtempZEMvsZDC;           // histogram with centrality vs multiplicity using ZEM vs ZDC 

  ClassDef(AliCentralitySelectionTask,1); 

};

#endif

