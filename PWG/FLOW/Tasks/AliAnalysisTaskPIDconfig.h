#ifndef ALIANALYSISTASKPIDconfig_H
#define ALIANALYSISTASKPIDconfig_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskPIDconfig.h */
// Author: Naghmeh Mohammadi, 10/07/2014

//==============================================================================
//
//
//
//
//==============================================================================

#include <TVectorDfwd.h>

#include "AliAnalysisTaskSE.h"
#include "AliPID.h"  
#include "AliPIDResponse.h"
#include "AliCentrality.h"



class AliESDEvent;
class AliAODEvent;
class AliESDTrack;
class AliAODTrack;
class AliPIDResponse;
class TList;
class AliVEvent;
class TH1F;
class TH2F;
class TH3F; 
class TH3F;



class AliAnalysisTaskPIDconfig : public AliAnalysisTaskSE {
  
  
public:
  AliAnalysisTaskPIDconfig();
  AliAnalysisTaskPIDconfig(const char *name);
  virtual ~AliAnalysisTaskPIDconfig();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);

  void SetFilterBit(Double_t b){fFilterBit = b;}
  void SetCentralityPercentileMin(Int_t b){fCentralityPercentileMin = b;}
  void SetCentralityPercentileMax(Int_t b){fCentralityPercentileMax = b;}
  void SetCentralityEstimator(TString b){fCentralityEstimator = b;}
  void SetUseCentrality(Bool_t b=kTRUE){fUseCentrality = b;}
  void SetCentralityTrigger(Int_t b=AliVEvent::kMB){ fTriggerSelection = b;}
  void SetDCAxyCut(Int_t b){fDCAxyCut = b;}
  void SetDCAzCut(Int_t b){fDCAzCut = b;}
  void SetCutTPCmultiplicityOutliersAOD(Bool_t b){fCutTPCmultiplicityOutliersAOD = b;}
  void SetData2011(Bool_t b){fData2011 = b;}
  void CheckCentrality(AliVEvent *event,Bool_t &centralitypass); //to use only events with the correct centrality....
  void SetCuts(Bool_t b){fPIDcuts = b;}
  //void MultiplicityOutlierCut(AliVEvent *event,Bool_t &centralitypass,Int_t ntracks);
  void SetPIDcontoursList(TList* b){fContourCutList = b;}
  //TGraph* GetPIDcontours(TString specie, Double_t Plow, Double_t Phigh,Int_t centMin, Int_t centMax){}

protected:  


  
private:
    AliVEvent             *fVevent;
    AliESDEvent           *fESD;
    AliAODEvent           *fAOD;
    AliPIDResponse        *fPIDResponse;             //! PID response Handler
    Int_t                  fTriggerSelection;
    Int_t                  fCentralityPercentileMin;
    Int_t                  fCentralityPercentileMax;
    Double_t               fFilterBit;
    Double_t               fDCAxyCut;
    Double_t               fDCAzCut;
    Bool_t                 fData2011;
    Bool_t                 fTriggerMB;
    Bool_t                 fTriggerCentral;
    Bool_t                 fUseCentrality;
    Bool_t                 fCutTPCmultiplicityOutliersAOD;
    Bool_t                 fPIDcuts;
    TString                fCentralityEstimator;   //"V0M","TRK","TKL","ZDC","FMD"
    TList                 *fContourCutList;
    TList                 *fListQA;           // List of all lists
    TList                 *fListQAtpctof;     //! List with combined PID from TPC + TOF
    TList                 *fListQAInfo;
    TH1F                  *fhistCentralityPass;
    TH1F                  *fNoEvents;
    TH1F                  *fpVtxZ;
    TH2F                  *fhistDCABefore;
    TH2F                  *fhistDCAAfter;
    TH1F                  *fhistPhiDistBefore;
    TH1F                  *fhistPhiDistAfter;
    TH1F                  *fhistEtaDistBefore;
    TH1F                  *fhistEtaDistAfter;
    TH2F                  *fTPCvsGlobalMultBeforeOutliers;
    TH2F                  *fTPCvsGlobalMultAfterOutliers;
    TH2F                  *fTPCvsGlobalMultAfter;
    TH2F                  *fHistBetavsPTOFbeforePID;
    TH2F                  *fHistdEdxVsPTPCbeforePID;
    TH2F                  *fHistBetavsPTOFafterPID;
    TH2F                  *fHistdEdxVsPTPCafterPID;
    TH3F                  *fhistNsigmaP;
    TH3F                  *fhistNsigmaPt;



  //qa object initialisation
  void SetupTPCTOFqa();
  void SetupEventInfo();
  //

  AliAnalysisTaskPIDconfig(const AliAnalysisTaskPIDconfig &other);
  AliAnalysisTaskPIDconfig& operator=(const AliAnalysisTaskPIDconfig &other);
  
  ClassDef(AliAnalysisTaskPIDconfig,2)  // Task to properly set the PID response functions of all detectors
};

#endif
