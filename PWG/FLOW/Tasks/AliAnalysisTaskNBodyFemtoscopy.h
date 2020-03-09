/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************** 
* Femtoscopy task for N bodies * 
Laura Serksnyte
laura.serksnyte@cern.ch
**************************************/ 

#ifndef ALIANALYSISTASKNBODYFEMTOSCOPY_H
#define ALIANALYSISTASKNBODYFEMTOSCOPY_H

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TString.h"

//================================================================================================================

class AliAnalysisTaskNBodyFemtoscopy : public AliAnalysisTaskSE{
 public:
  
  AliAnalysisTaskNBodyFemtoscopy();
  AliAnalysisTaskNBodyFemtoscopy(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskNBodyFemtoscopy(); 
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  
  // 0.) Methods called in the constructor:
  virtual void InitializeArrays();
 
  // 1.) Methods called in UserCreateOutputObjects():
  virtual void BookAndNestAllLists();
  virtual void BookControlHistograms();
  virtual void BookFinalResultsHistograms();

  // 2.) Methods called in UserExec(Option_t *):
  
  Bool_t CommonEventCuts(AliVEvent *ave, Float_t centrality); // Common event and vertex cuts.
  Bool_t CommonTrackCuts(AliAODTrack *atrack); // Common track cuts.

  // 3.) Setters and getters:
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 
  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  void SetPtBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbinsPt = nbins;
   this->fMinBinPt = min;
   this->fMaxBinPt = max;
  };

  void SetCentralityEstimator(TString estimatorName)
  {
   this->fCentralityEstimator = estimatorName;
  };

  void SetCentralityBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbinsCentrality = nbins;
   this->fMinCentrality = min;
   this->fMaxCentrality = max;
  };

   void SetPtRange(Float_t ptMin, Float_t ptMax)
  {
   this->fPtRange[0] = ptMin;
   this->fPtRange[1] = ptMax;
  };
  void SetEtaRange(Float_t etaMin, Float_t etaMax)
  {
   this->fEtaRange[0] = etaMin;
   this->fEtaRange[1] = etaMax;
  };
  void SetPhiRange(Float_t phiMin, Float_t phiMax)
  {
   this->fPhiRange[0] = phiMin;
   this->fPhiRange[1] = phiMax;
  };
  void SetFilterBit(Int_t filterBitCheck)
  {
   this->fFilterBit = filterBitCheck;
  };
  void SetRejectEventsNoPrimaryVertex(Bool_t primvertex) {
  	this->fRejectEventsNoPrimaryVertex = primvertex;
  };
  void SetVertexZCut(Float_t VertexZMin, Float_t VertexZMax)
  {
   this->fVertexZ[0] = VertexZMin;
   this->fVertexZ[1] = VertexZMax;
  };


   private:
  AliAnalysisTaskNBodyFemtoscopy(const AliAnalysisTaskNBodyFemtoscopy& aatmpf);
  AliAnalysisTaskNBodyFemtoscopy& operator=(const AliAnalysisTaskNBodyFemtoscopy& aatmpf);
  
  // 0.) Base lists:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  // 1.) Control histograms:  
  TList *fControlHistogramsList; // list to hold all control histograms
  TH1F *fPtHist;                 // atrack->Pt()
  TH1F *fCentralityHist;         // ams->GetMultiplicityPercentile()
  TH1F *fNumberOfTracksHist;         // ams->GetMultiplicityPercentile()
  // After event cuts
  TH1F *fCentralityHistVzCut;         // ams->GetMultiplicityPercentile() after VzCut
  // After track cuts
  TH1F *fPtHistEtaCut;                 // atrack->Pt() after Eta cut
  TH1F *fPtHistEtaCutPTCut;            // atrack->Pt() after Eta and pT cut
  TH1F *fPtHistEtaCutPTCutPhiCut;      // atrack->Pt() after Eta, pT and Phi cut
  // Other variables
  Int_t fNbinsPt;                // number of bins
  Float_t fMinBinPt;             // min bin
  Float_t fMaxBinPt;             // min bin
  TString fCentralityEstimator;  // centrality estimator choice
  Float_t fNbinsCentrality;      // cenrality bin number
  Float_t fMinCentrality;        // min centrality
  Float_t fMaxCentrality;        // max centrality
  
  // 2.) Final results:
  TList *fFinalResultsList; // list to hold all histograms with final results

  // 3) Global cuts: track and event
  //Event 
  Bool_t fRejectEventsNoPrimaryVertex; //Reject events without primary vertex
  Bool_t fCutOnVertexZ;                //Apply vertexZ cut
  Float_t fVertexZ[2];                // VertexZMin = fVertexZ[0], VertexZMax = fVertexZ[1]
  //Track
  Bool_t fApplyCommonTrackCuts;       //Apply common track cuts
  Float_t fPtRange[2];                // ptMin = fPtRange[0], ptMax = fPtRange[1]
  Float_t fEtaRange[2];               // etaMin = etaRange[0], etaMax = etaRange[1]
  Float_t fPhiRange[2];               // phiMin = phiRange[0], phiMax = phiRange[1]
  Int_t fFilterBit;                   //  Filter bit choice


  // Increase this counter in each new version:
  ClassDef(AliAnalysisTaskNBodyFemtoscopy, 3);

};

//================================================================================================================

#endif




