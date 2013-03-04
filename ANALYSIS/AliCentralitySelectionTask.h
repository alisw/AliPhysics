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
  void DontUseScaling()                    {fUseScaling=kFALSE;}  
  void DontUseCleaning()                   {fUseCleaning=kFALSE;}
  void SetFillHistos()                     {fFillHistos=kTRUE; DefineOutput(1, TList::Class());
}

 private:

  Int_t SetupRun(const AliVEvent* const esd);
  Bool_t IsOutlierV0MSPD(Float_t spd, Float_t v0, Int_t cent) const;
  Bool_t IsOutlierV0MTPC(Int_t tracks, Float_t v0, Int_t cent) const;
  Bool_t IsOutlierV0MZDC(Float_t zdc, Float_t v0) const;
  Bool_t IsOutlierV0MZDCECal(Float_t zdc, Float_t v0) const;

  TString  fAnalysisInput; 	// "ESD", "AOD"
  Bool_t   fIsMCInput;          // true when input is MC
  Int_t    fCurrentRun;         // current run number
  Bool_t   fUseScaling;         // flag to use scaling 
  Bool_t   fUseCleaning;        // flag to use cleaning  
  Bool_t   fFillHistos;         // flag to fill the QA histos
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
  AliESDtrackCuts* fEsdTrackCuts;  //! optional track cuts
  AliESDtrackCuts* fEsdTrackCutsExtra1;  //! optional track cuts
  AliESDtrackCuts* fEsdTrackCutsExtra2;  //! optional track cuts

  Float_t  fZVCut;              //! z-vertex cut (in cm)
  Float_t  fOutliersCut;        //! outliers cut (in n-sigma)
  Int_t    fQuality;            //! quality for centrality determination

  Bool_t   fIsSelected;         //! V0BG rejection

  Bool_t   fMSL;                //! 
  Bool_t   fMSH;                //! 
  Bool_t   fMUL;                //! 
  Bool_t   fMLL;                //! 

  Bool_t   fEJE;                //! 
  Bool_t   fEGA;                //! 
  Bool_t   fPHS;                //! 

  Bool_t   fMB;                 //! if the event is MB
  Bool_t   fCVHN;               //! if the event is central trigger
  Bool_t   fCVLN;               //! if the event is semicentral trigger
  Bool_t   fCVHNbit;            //! if the event is central trigger
  Bool_t   fCVLNbit;            //! if the event is semicentral trigger
  Bool_t   fCCENT;              //! if the event is central trigger
  Bool_t   fCSEMI;              //! if the event is semicentral trigger
  Bool_t   fCCENTbit;           //! if the event is central trigger
  Bool_t   fCSEMIbit;           //! if the event is semicentral trigger

  Float_t  fCentV0M;            // percentile centrality from V0
  Float_t  fCentV0A;            // percentile centrality from V0A
  Float_t  fCentV0C;            // percentile centrality from V0C
  Float_t  fCentV0MEq;          // percentile centrality from V0 equalized channel
  Float_t  fCentV0AEq;          // percentile centrality from V0A equalized channel
  Float_t  fCentV0CEq;          // percentile centrality from V0C equalized channel
  Float_t  fCentFMD;            // percentile centrality from FMD
  Float_t  fCentTRK;            // percentile centrality from tracks
  Float_t  fCentTKL;            // percentile centrality from tracklets
  Float_t  fCentCL0;            // percentile centrality from clusters in layer 0
  Float_t  fCentCL1;            // percentile centrality from clusters in layer 1
  Float_t  fCentCND;            // percentile centrality from candle 
  Float_t  fCentZNA;            // percentile centrality from ZNA 
  Float_t  fCentZNC;            // percentile centrality from ZNC 
  Float_t  fCentNPA;            // percentile centrality from Npart (MC) 
  Float_t  fCentV0MvsFMD;       // percentile centrality from V0 vs FMD
  Float_t  fCentTKLvsV0M;       // percentile centrality from tracklets vs V0
  Float_t  fCentZEMvsZDC;       // percentile centrality from ZEM vs ZDC

  Float_t  fCentV0Mtrue;        // percentile centrality from true (sim) V0A+V0C
  Float_t  fCentV0Atrue;        // percentile centrality from true (sim) V0A
  Float_t  fCentV0Ctrue;        // percentile centrality from true (sim) V0C
  Float_t  fCentV0MEqtrue;      // percentile centrality from true (sim) V0A+V0C equalized channel
  Float_t  fCentV0AEqtrue;      // percentile centrality from true (sim) V0A equalized channel
  Float_t  fCentV0CEqtrue;      // percentile centrality from true (sim) V0C equalized channel
  Float_t  fCentFMDtrue;        // percentile centrality from true (sim) FMD
  Float_t  fCentTRKtrue;        // percentile centrality from true (sim) tracks
  Float_t  fCentTKLtrue;        // percentile centrality from true (sim) tracklets
  Float_t  fCentCL0true;        // percentile centrality from true (sim) Clusters in layer 0
  Float_t  fCentCL1true;        // percentile centrality from true (sim) Clusters in layer 1
  Float_t  fCentCNDtrue;        // percentile centrality from true (sim) tracks (candle condition)
  Float_t  fCentZNAtrue;        // percentile centrality from true (sim) ZNA
  Float_t  fCentZNCtrue;        // percentile centrality from true (sim) ZNC


  TH1F    *fHtempV0M;           // histogram with centrality vs multiplicity using V0
  TH1F    *fHtempV0A;           // histogram with centrality vs multiplicity using V0A
  TH1F    *fHtempV0C;           // histogram with centrality vs multiplicity using V0C
  TH1F    *fHtempV0MEq;         // histogram with centrality vs multiplicity using V0  equalized channel
  TH1F    *fHtempV0AEq;         // histogram with centrality vs multiplicity using V0A equalized channel
  TH1F    *fHtempV0CEq;         // histogram with centrality vs multiplicity using V0C equalized channel
  TH1F    *fHtempFMD;           // histogram with centrality vs multiplicity using FMD
  TH1F    *fHtempTRK;           // histogram with centrality vs multiplicity using tracks
  TH1F    *fHtempTKL;           // histogram with centrality vs multiplicity using tracklets
  TH1F    *fHtempCL0;           // histogram with centrality vs multiplicity using clusters in layer 0
  TH1F    *fHtempCL1;           // histogram with centrality vs multiplicity using clusters in layer 1
  TH1F    *fHtempCND;           // histogram with centrality vs multiplicity using candle
  TH1F    *fHtempZNA;           // histogram with centrality vs multiplicity using ZNA
  TH1F    *fHtempZNC;           // histogram with centrality vs multiplicity using ZNC
  TH1F    *fHtempV0MvsFMD;      // histogram with centrality vs multiplicity using V0 vs FMD   
  TH1F    *fHtempTKLvsV0M;      // histogram with centrality vs multiplicity using tracklets vs V0
  TH2F    *fHtempZEMvsZDC;      // histogram with centrality vs multiplicity using ZEM vs ZDC 
  TH1F    *fHtempNPA;           // histogram with centrality vs multiplicity using Npart

  TH1F    *fHtempV0Mtrue;       // histogram with centrality true (sim) vs multiplicity using V0
  TH1F    *fHtempV0Atrue;       // histogram with centrality true (sim) vs multiplicity using V0A
  TH1F    *fHtempV0Ctrue;       // histogram with centrality true (sim) vs multiplicity using V0C
  TH1F    *fHtempV0MEqtrue;     // histogram with centrality true (sim) vs multiplicity using V0 equalized channel
  TH1F    *fHtempV0AEqtrue;     // histogram with centrality true (sim) vs multiplicity using V0A equalized channel
  TH1F    *fHtempV0CEqtrue;     // histogram with centrality true (sim) vs multiplicity using V0C equalized channel
  TH1F    *fHtempFMDtrue;       // histogram with centrality true (sim) vs multiplicity using FMD
  TH1F    *fHtempTRKtrue;       // histogram with centrality true (sim) vs multiplicity using tracks
  TH1F    *fHtempTKLtrue;       // histogram with centrality true (sim) vs multiplicity using tracklets
  TH1F    *fHtempCL0true;       // histogram with centrality true (sim) vs multiplicity using clusters in layer 0
  TH1F    *fHtempCL1true;       // histogram with centrality true (sim) vs multiplicity using clusters in layer 1
  TH1F    *fHtempCNDtrue;       // histogram with centrality true (sim) vs multiplicity using candle
  TH1F    *fHtempZNAtrue;       // histogram with centrality true (sim) vs multiplicity using ZNA
  TH1F    *fHtempZNCtrue;       // histogram with centrality true (sim) vs multiplicity using ZNC

  TList   *fOutputList; // output list
  

  TH1F *fHOutCentV0M    ;       //control histogram for centrality
  TH1F *fHOutCentV0A    ;       //control histogram for centrality
  TH1F *fHOutCentV0C    ;       //control histogram for centrality
  TH1F *fHOutCentV0MEq  ;       //control histogram for centrality
  TH1F *fHOutCentV0AEq  ;       //control histogram for centrality
  TH1F *fHOutCentV0CEq  ;       //control histogram for centrality
  TH1F *fHOutCentV0MCVHN;       //control histogram for centrality
  TH1F *fHOutCentV0MCVLN;       //control histogram for centrality
  TH1F *fHOutCentV0MCVHNinMB;   //control histogram for centrality
  TH1F *fHOutCentV0MCVLNinMB;   //control histogram for centrality
  TH1F *fHOutCentV0MCCENT;      //control histogram for centrality
  TH1F *fHOutCentV0MCSEMI;      //control histogram for centrality
  TH1F *fHOutCentV0MCCENTinMB;  //control histogram for centrality
  TH1F *fHOutCentV0MCSEMIinMB;  //control histogram for centrality
  TH1F *fHOutCentV0MMSL;        //control histogram for centrality
  TH1F *fHOutCentV0MMSH;        //control histogram for centrality
  TH1F *fHOutCentV0MMUL;        //control histogram for centrality
  TH1F *fHOutCentV0MMLL;        //control histogram for centrality
  TH1F *fHOutCentV0MEJE;        //control histogram for centrality
  TH1F *fHOutCentV0MEGA;        //control histogram for centrality
  TH1F *fHOutCentV0MPHS;        //control histogram for centrality
  TH1F *fHOutCentV0MMSLinMB;    //control histogram for centrality
  TH1F *fHOutCentV0MMSHinMB;    //control histogram for centrality
  TH1F *fHOutCentV0MMULinMB;    //control histogram for centrality
  TH1F *fHOutCentV0MMLLinMB;    //control histogram for centrality
  TH1F *fHOutCentV0MEJEinMB;    //control histogram for centrality
  TH1F *fHOutCentV0MEGAinMB;    //control histogram for centrality
  TH1F *fHOutCentV0MPHSinMB;    //control histogram for centrality
  TH1F *fHOutCentFMD     ;      //control histogram for centrality
  TH1F *fHOutCentTRK     ;      //control histogram for centrality
  TH1F *fHOutCentTKL     ;      //control histogram for centrality
  TH1F *fHOutCentCL0     ;      //control histogram for centrality
  TH1F *fHOutCentCL1     ;      //control histogram for centrality
  TH1F *fHOutCentCND     ;      //control histogram for centrality
  TH1F *fHOutCentNPA     ;      //control histogram for centrality
  TH1F *fHOutCentZNA     ;      //control histogram for centrality
  TH1F *fHOutCentZNC     ;      //control histogram for centrality
  TH1F *fHOutCentV0MvsFMD;      //control histogram for centrality
  TH1F *fHOutCentTKLvsV0M;      //control histogram for centrality
  TH1F *fHOutCentZEMvsZDC;      //control histogram for centrality
  TH2F *fHOutCentV0MvsCentCL1;  //control histogram for centrality
  TH2F *fHOutCentV0MvsCentTRK;  //control histogram for centrality
  TH2F *fHOutCentTRKvsCentCL1;  //control histogram for centrality
  TH2F *fHOutCentV0MvsCentZDC;  //control histogram for centrality
  TH2F *fHOutCentV0AvsCentV0C;  //control histogram for centrality
  TH2F *fHOutCentV0AvsCentTRK;  //control histogram for centrality
  TH2F *fHOutCentV0AvsCentCND;  //control histogram for centrality
  TH2F *fHOutCentV0AvsCentCL1;  //control histogram for centrality
  TH2F *fHOutCentV0CvsCentTRK;  //control histogram for centrality
  TH2F *fHOutCentV0CvsCentCND;  //control histogram for centrality
  TH2F *fHOutCentV0CvsCentCL1;  //control histogram for centrality
  TH2F *fHOutCentNPAvsCentV0A;  //control histogram for centrality
  TH2F *fHOutCentNPAvsCentV0C;  //control histogram for centrality
  TH2F *fHOutCentNPAvsCentTRK;  //control histogram for centrality
  TH2F *fHOutCentNPAvsCentCND;  //control histogram for centrality
  TH2F *fHOutCentNPAvsCentCL1;  //control histogram for centrality
  TH2F *fHOutCentZNAvsCentV0A;  //control histogram for centrality
  TH2F *fHOutCentZNAvsCentV0C;  //control histogram for centrality
  TH2F *fHOutCentZNAvsCentTRK;  //control histogram for centrality
  TH2F *fHOutCentZNAvsCentCND;  //control histogram for centrality
  TH2F *fHOutCentZNAvsCentCL1;  //control histogram for centrality

  TH2F *fHOutMultV0AC;          //control histogram for multiplicity
  TH1F *fHOutMultV0M ;          //control histogram for multiplicity
  TH1F *fHOutMultV0A ;          //control histogram for multiplicity
  TH1F *fHOutMultV0C ;          //control histogram for multiplicity
  TH1F *fHOutMultV0MEq ;        //control histogram for multiplicity
  TH1F *fHOutMultV0AEq ;        //control histogram for multiplicity
  TH1F *fHOutMultV0CEq ;        //control histogram for multiplicity
  TH1F *fHOutMultV0Mnc ;        //control histogram for multiplicity
  TH1F *fHOutMultV0Anc ;        //control histogram for multiplicity
  TH1F *fHOutMultV0Cnc ;        //control histogram for multiplicity
  TH1F *fHOutMultV0O ;          //control histogram for multiplicity
  TH2F *fHOutMultV0Cells ;      //control histogram for multiplicity
  TH1F *fHOutMultFMD ;          //control histogram for multiplicity
  TH1F *fHOutMultTRK ;          //control histogram for multiplicity
  TH1F *fHOutMultTKL ;          //control histogram for multiplicity
  TH1F *fHOutMultCL0 ;          //control histogram for multiplicity
  TH1F *fHOutMultCL1 ;          //control histogram for multiplicity
  TH1F *fHOutMultCND ;          //control histogram for multiplicity
  TH1F *fHOutMultNPA ;          //control histogram for multiplicity
  TH1F *fHOutMultZNA ;          //control histogram for multiplicity
  TH1F *fHOutMultZNC ;          //control histogram for multiplicity

  TH2F *fHOutMultV0MvsZDN;      //control histogram for multiplicity
  TH2F *fHOutMultZEMvsZDN;      //control histogram for multiplicity
  TH2F *fHOutMultV0MvsZDC;      //control histogram for multiplicity
  TH2F *fHOutMultZEMvsZDC;      //control histogram for multiplicity
  TH2F *fHOutMultZEMvsZDCw;     //control histogram for multiplicity
  TH2F *fHOutMultV0MvsCL1;      //control histogram for multiplicity
  TH2F *fHOutMultV0MvsTRK;      //control histogram for multiplicity
  TH2F *fHOutMultTRKvsCL1;      //control histogram for multiplicity
  TH2F *fHOutMultV0MvsV0O;      //control histogram for multiplicity
  TH2F *fHOutMultV0OvsCL1;      //control histogram for multiplicity
  TH2F *fHOutMultV0OvsTRK;      //control histogram for multiplicity
  TH2F *fHOutMultCL1vsTKL;      //control histogram for multiplicity

  TH1F *fHOutCentV0Mqual1     ; //control histogram for centrality quality 1
  TH1F *fHOutCentTRKqual1     ; //control histogram for centrality quality 1
  TH1F *fHOutCentCL1qual1     ; //control histogram for centrality quality 1
  TH2F *fHOutMultV0MvsCL1qual1; //control histogram for multiplicity quality 1
  TH2F *fHOutMultV0MvsTRKqual1; //control histogram for multiplicity quality 1
  TH2F *fHOutMultTRKvsCL1qual1; //control histogram for multiplicity quality 1

  TH1F *fHOutCentV0Mqual2     ; //control histogram for centrality quality 2
  TH1F *fHOutCentTRKqual2     ; //control histogram for centrality quality 2
  TH1F *fHOutCentCL1qual2     ; //control histogram for centrality quality 2
  TH2F *fHOutMultV0MvsCL1qual2; //control histogram for multiplicity quality 2
  TH2F *fHOutMultV0MvsTRKqual2; //control histogram for multiplicity quality 2
  TH2F *fHOutMultTRKvsCL1qual2; //control histogram for multiplicity quality 2

  TH1F *fHOutQuality ;          //control histogram for quality
  TH1F *fHOutVertex ;           //control histogram for vertex SPD
  TH1F *fHOutVertexT0 ;         //control histogram for vertex T0

  ClassDef(AliCentralitySelectionTask, 26); 
};

#endif

