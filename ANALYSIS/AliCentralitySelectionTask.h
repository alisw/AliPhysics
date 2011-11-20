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
class TH2F;
class TList;
class TString;

class AliESDEvent;
class AliESDtrackCuts;

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

  void SetInput(const char* input)         {fAnalysisInput = input;}
  void SetMCInput()                        {fIsMCInput = kTRUE;}
  void SetPass(Int_t pass)                 {fPass = pass;}
  void DontUseScaling()                    {fUseScaling=kFALSE;}  
  void DontUseCleaning()                   {fUseCleaning=kFALSE;}

 private:

  Int_t SetupRun(AliESDEvent* const esd);
  Bool_t IsOutlierV0MSPD(Float_t spd, Float_t v0, Int_t cent) const;
  Bool_t IsOutlierV0MTPC(Int_t tracks, Float_t v0, Int_t cent) const;
  Bool_t IsOutlierV0MZDC(Float_t zdc, Float_t v0) const;
  Bool_t IsOutlierV0MZDCECal(Float_t zdc, Float_t v0) const;

  TString  fAnalysisInput; 	// "ESD", "AOD"
  Bool_t   fIsMCInput;          // true when input is MC
  Int_t    fPass;               // pass of reconstruction
  Int_t    fCurrentRun;         // current run number
  Bool_t   fUseScaling;         // flag to use scaling 
  Bool_t   fUseCleaning;        // flag to use cleaning  
  Float_t  fV0MScaleFactor;     // scale factor V0M
  Float_t  fSPDScaleFactor;     // scale factor SPD
  Float_t  fTPCScaleFactor;     // scale factor TPC
  Float_t  fV0MScaleFactorMC;   // scale factor V0M for MC
  Float_t  fV0MSPDOutlierPar0;  // outliers parameter
  Float_t  fV0MSPDOutlierPar1;  // outliers parameter
  Float_t  fV0MTPCOutlierPar0;  // outliers parameter
  Float_t  fV0MTPCOutlierPar1;  // outliers parameter
  Float_t  fV0MSPDSigmaOutlierPar0;  // outliers parameter
  Float_t  fV0MSPDSigmaOutlierPar1;  // outliers parameter
  Float_t  fV0MSPDSigmaOutlierPar2;  // outliers parameter
  Float_t  fV0MTPCSigmaOutlierPar0;  // outliers parameter
  Float_t  fV0MTPCSigmaOutlierPar1;  // outliers parameter
  Float_t  fV0MTPCSigmaOutlierPar2;  // outliers parameter  			   			   
  Float_t  fV0MZDCOutlierPar0;	     // outliers parameter
  Float_t  fV0MZDCOutlierPar1;	     // outliers parameter
  Float_t  fV0MZDCEcalOutlierPar0;   // outliers parameter
  Float_t  fV0MZDCEcalOutlierPar1;   // outliers parameter

  AliESDtrackCuts* fTrackCuts;  //! optional track cuts

  Float_t  fZVCut;              //! z-vertex cut (in cm)
  Float_t  fOutliersCut;        //! outliers cut (in n-sigma)
  Int_t    fQuality;            //! quality for centrality determination

  Bool_t   fCVHN;               //! if the event is central trigger
  Bool_t   fCVLN;               //! if the event is semicentral trigger

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
  TH2F    *fHtempZEMvsZDC;           // histogram with centrality vs multiplicity using ZEM vs ZDC 

  TList   *fOutputList; // output list
  
  TH1F *fHOutCentV0M     ;    //control histogram for centrality
  TH1F *fHOutCentV0M_CVHN;    //control histogram for centrality
  TH1F *fHOutCentV0M_CVLN;    //control histogram for centrality
  TH1F *fHOutCentFMD     ;    //control histogram for centrality
  TH1F *fHOutCentTRK     ;    //control histogram for centrality
  TH1F *fHOutCentTKL     ;    //control histogram for centrality
  TH1F *fHOutCentCL0     ;    //control histogram for centrality
  TH1F *fHOutCentCL1     ;    //control histogram for centrality
  TH1F *fHOutCentV0MvsFMD;    //control histogram for centrality
  TH1F *fHOutCentTKLvsV0M;    //control histogram for centrality
  TH1F *fHOutCentZEMvsZDC;    //control histogram for centrality
  TH2F *fHOutCentV0MvsCentCL1;    //control histogram for centrality
  TH2F *fHOutCentV0MvsCentTRK;    //control histogram for centrality
  TH2F *fHOutCentTRKvsCentCL1;    //control histogram for centrality
  TH2F *fHOutCentV0MvsCentZDC;    //control histogram for centrality

  TH1F *fHOutMultV0M ;        //control histogram for multiplicity
  TH1F *fHOutMultV0O ;        //control histogram for multiplicity
  TH1F *fHOutMultFMD ;        //control histogram for multiplicity
  TH1F *fHOutMultTRK ;        //control histogram for multiplicity
  TH1F *fHOutMultTKL ;        //control histogram for multiplicity
  TH1F *fHOutMultCL0 ;        //control histogram for multiplicity
  TH1F *fHOutMultCL1 ;        //control histogram for multiplicity

  TH2F *fHOutMultV0MvsZDN;    //control histogram for multiplicity
  TH2F *fHOutMultZEMvsZDN;    //control histogram for multiplicity
  TH2F *fHOutMultV0MvsZDC;    //control histogram for multiplicity
  TH2F *fHOutMultZEMvsZDC;    //control histogram for multiplicity
  TH2F *fHOutMultZEMvsZDCw;   //control histogram for multiplicity
  TH2F *fHOutMultV0MvsCL1;    //control histogram for multiplicity
  TH2F *fHOutMultV0MvsTRK;    //control histogram for multiplicity
  TH2F *fHOutMultTRKvsCL1;    //control histogram for multiplicity
  TH2F *fHOutMultV0MvsV0O;    //control histogram for multiplicity
  TH2F *fHOutMultV0OvsCL1;    //control histogram for multiplicity
  TH2F *fHOutMultV0OvsTRK;    //control histogram for multiplicity

  TH1F *fHOutCentV0Mqual1     ;    //control histogram for centrality quality 1
  TH1F *fHOutCentTRKqual1     ;    //control histogram for centrality quality 1
  TH1F *fHOutCentCL1qual1     ;    //control histogram for centrality quality 1
  TH2F *fHOutMultV0MvsCL1qual1;    //control histogram for multiplicity quality 1
  TH2F *fHOutMultV0MvsTRKqual1;    //control histogram for multiplicity quality 1
  TH2F *fHOutMultTRKvsCL1qual1;    //control histogram for multiplicity quality 1

  TH1F *fHOutCentV0Mqual2     ;    //control histogram for centrality quality 2
  TH1F *fHOutCentTRKqual2     ;    //control histogram for centrality quality 2
  TH1F *fHOutCentCL1qual2     ;    //control histogram for centrality quality 2
  TH2F *fHOutMultV0MvsCL1qual2;    //control histogram for multiplicity quality 2
  TH2F *fHOutMultV0MvsTRKqual2;    //control histogram for multiplicity quality 2
  TH2F *fHOutMultTRKvsCL1qual2;    //control histogram for multiplicity quality 2

  TH1F *fHOutQuality ;        //control histogram for quality
  TH1F *fHOutVertex ;         //control histogram for vertex

  ClassDef(AliCentralitySelectionTask, 12); 
};

#endif

