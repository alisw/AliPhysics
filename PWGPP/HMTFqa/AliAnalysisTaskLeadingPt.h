#ifndef ALIANALYSISTASKLEADINGPT_H
#define ALIANALYSISTASKLEADINGPT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */
//Authors: Antonio Ortiz Velasquez, antonio.ortiz@nucleares.unam.mx
//modified by: Sergio Iga , sergio.arturo.iga.buitron@cern.ch
//ptleading distributions, vs nch analysis

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TProfile.h>
#include <TTreeStream.h>
#include <TRandom.h>
#include <TObject.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 
#include <AliESDtrackCuts.h>

#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"

class AliAnalysisUtils;
class AliPPVsMultUtils;

class AliAnalysisTaskLeadingPt : public AliAnalysisTaskSE {
 public:
 

  AliAnalysisTaskLeadingPt();
  AliAnalysisTaskLeadingPt(const char *name);
  virtual ~AliAnalysisTaskLeadingPt();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     
  //Double_t GetMinPt() { return fMinPt; }   
  //Int_t    GetTreeOption() { return fTreeOption; }  


  virtual void  SetTrigger(UInt_t ktriggerInt) {ftrigBit = ktriggerInt;}
  virtual void  SetTrackFilterGolden(AliAnalysisFilter* trackF) {fTrackFilterGolden = trackF;}
  virtual void  SetTrackFilterTPC(AliAnalysisFilter* trackF) {fTrackFilterTPC = trackF;}
  virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetPtLeadingCut(Double_t ptlCut){fPtCutLeading = ptlCut;}
  virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;} 
  virtual void  SetPileUpRejMV(Bool_t isrej) {fPileUpRejMV = isrej;} 
  virtual void  SetPileUpMvSettings(AliAnalysisUtils* fAliAnalysisUtils);
  virtual void  SetNcontributors(Int_t nContributors){fnContributors = nContributors;}


 private:
  virtual Float_t GetVertex(const AliVEvent* event) const;
  virtual void AnalyzeESD(AliESDEvent* esd); 
  virtual void AnalyzeAOD(AliAODEvent* aod); 

  ULong64_t GetEventIdAsLong(AliVHeader* header) const;

  AliESDEvent* fESD;                  //! ESD object
  AliAODEvent* fAOD;                  //! AOD object

  AliAnalysisFilter* fTrackFilterGolden;    //  Track Filter, set 2010 with golden cuts
  AliAnalysisFilter* fTrackFilterTPC; // track filter for TPC only tracks
  TString       fAnalysisType;        //  "ESD" or "AOD"
  UInt_t       ftrigBit;
  TRandom*      fRandom;              //! random number generator
  Bool_t        fPileUpRej;           // kTRUE is pile-up is rejected
  Bool_t        fPileUpRejMV;
  AliAnalysisUtils* fAliAnalysisUtils;
  Int_t fnContributors;



  //
  // Cuts and options
  //

  Double_t     fVtxCut;             // Vtx cut on z position in cm
  Double_t     fEtaCut;             // Eta cut used to select particles
  //
  // Help variables
  //
  Short_t      fTriggeredEventMB;   // 1 = triggered, 0 = not trigged (MC only)
  Short_t      fVtxStatus;          // -1 = no vtx, 0 = outside cut, 1 = inside cut
  Float_t      fZvtx;               // z vertex
  Int_t        fRun;                // run no
  ULong64_t    fEventId;            // unique event id
  Float_t      fPtCutLeading;

  //
  // Output objects
  //
  TList*        fListOfObjects;     //! Output list of objects
  TH1I*         fEvents;            //! No of accepted events
  TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
  TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts
  TH1F* fn1;

	//W/o PileUp Rejection
  TH1F* fnchL[7];
  TH1F* fptLL[7];
  TH1F* fptL[7];
  TH1F* fnchH[7];
  TH1F* fptLH[7];
  TH1F* fptH[7];
  TH2F* ftrackVsClusters;

	//SPD PileUp Rejection
  TH1F* fnchL_SPD[7];
  TH1F* fptLL_SPD[7];
  TH1F* fptL_SPD[7];
  TH1F* fnchH_SPD[7];
  TH1F* fptLH_SPD[7];
  TH1F* fptH_SPD[7];
  TH2F* ftrackVsClusters_SPD;

	//MV PileUp Rejection
  TH1F* fnchL_MV[7];
  TH1F* fptLL_MV[7];
  TH1F* fptL_MV[7];
  TH1F* fnchH_MV[7];
  TH1F* fptLH_MV[7];
  TH1F* fptH_MV[7];
  TH2F* ftrackVsClusters_MV;

	//SPD && MV PileUp Rejection
  TH1F* fnchL_SPDMV[7];
  TH1F* fptLL_SPDMV[7];
  TH1F* fptL_SPDMV[7];
  TH1F* fnchH_SPDMV[7];
  TH1F* fptLH_SPDMV[7];
  TH1F* fptH_SPDMV[7];
  TH2F* ftrackVsClusters_SPDMV;

  AliAnalysisTaskLeadingPt(const AliAnalysisTaskLeadingPt&);            // not implemented
  AliAnalysisTaskLeadingPt& operator=(const AliAnalysisTaskLeadingPt&); // not implemented

  //TTree*        fTree;              //! Debug tree 

  ClassDef(AliAnalysisTaskLeadingPt, 1);    //Analysis task for high pt analysis 
};

#endif
