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
class TH2F;

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

 private:
  Int_t    fDebug;	   	// Debug flag
  TString  fAnalysisInput; 	// "ESD", "AOD"
  Bool_t   fIsMCInput;          // true when input is MC
  TFile   *fFile;               // file that holds the centrality vs multiplicity 1d
  TFile   *fFile2;               // file that holds the centrality vs multiplicity 2d
  TString  fCentfilename;       // name of this file 1d
  TString  fCentfilename2;       // name of this file 2d

  Float_t  fCentV0M;            // percentile centrality from V0
  Float_t  fCentFMD;            // percentile centrality from FMD
  Float_t  fCentTRK;            // percentile centrality from tracks
  Float_t  fCentTKL;            // percentile centrality from tracklets
  Float_t  fCentCL0;            // percentile centrality from clusters in layer 0
  Float_t  fCentV0MvsFMD;       // percentile centrality from V0 vs FMD
  Float_t  fCentTKLvsV0M;       // percentile centrality from tracklets vs V0
  Float_t  fCentZEMvsZDC;       // percentile centrality from ZEM vs ZDC

  TH1D    *fHtempV0M;           // histogram with centrality vs multiplicity using V0
  TH1D    *fHtempFMD;           // histogram with centrality vs multiplicity using FMD
  TH1D    *fHtempTRK;           // histogram with centrality vs multiplicity using tracks
  TH1D    *fHtempTKL;           // histogram with centrality vs multiplicity using tracklets
  TH1D    *fHtempCL0;           // histogram with centrality vs multiplicity using clusters in layer 0
  TH1D    *fHtempV0MvsFMD;           // histogram with centrality vs multiplicity using V0 vs FMD   
  TH1D    *fHtempTKLvsV0M;           // histogram with centrality vs multiplicity using tracklets vs V0
  TH1D    *fHtempZEMvsZDC;           // histogram with centrality vs multiplicity using ZEM vs ZDC 

  ClassDef(AliCentralitySelectionTask,1); 

};

#endif

