
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
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"
#include "TExMap.h"
#include "TH1F.h"
#include "TString.h"
#include "TRandom3.h"

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
  virtual void BookEverything(); // use this to book all un-classified objects

  // 2.) Methods called in UserExec(Option_t *):
  
  Bool_t CommonEventCuts(AliVEvent *ave, Float_t centrality); // Common event and vertex cuts.
  Bool_t CommonTrackCuts(AliAODTrack *atrack); // Common track cuts.
  Bool_t GlobalTrackCuts(AliAODTrack *gtrack); // global track cuts, does not include the filterbit
  Bool_t Pion(AliAODTrack *atrack, Int_t charge = 1, Bool_t bPrimary = kTRUE); // check if particle is a pion
  Bool_t Kaon(AliAODTrack *atrack, Int_t charge = 1, Bool_t bPrimary = kTRUE);  // check if particle is a kaon
  Bool_t Proton(AliAODTrack *atrack, Int_t charge = 1, Bool_t bPrimary = kTRUE);  // check if particle is a proton
  virtual void GlobalTracksAODTEST(AliAODEvent *aAOD, Int_t index); // fill fGlobalTracksAODTEST - take only  normal global tracks
  virtual void CreateRandomIndices(AliAODEvent *aAOD);

  // 3.) Setters and getters:
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 
  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  // set binning ================================================================================================================
  void SetPtBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbinsPt = nbins;
   this->fMinBinPt = min;
   this->fMaxBinPt = max;
  };
  void SetCentralityBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbinsCentrality = nbins;
   this->fMinCentrality = min;
   this->fMaxCentrality = max;
  };
  void SetMultiplicityBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbinsMultiplicity = nbins;
   this->fMinBinMultiplicity = min;
   this->fMaxBinMultiplicity = max;
  };
  void SetTPCOnlyVsGlobalBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbinsTPCOnlyVsGlobal = nbins;
   this->fMinBinTPCOnlyVsGlobal = min;
   this->fMaxBinTPCOnlyVsGlobal = max;
  };

  // set ranges and limits ======================================================================================================
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
  void SetVertexZCut(Float_t VertexZMin, Float_t VertexZMax)
  {
   this->fVertexZ[0] = VertexZMin;
   this->fVertexZ[1] = VertexZMax;
  };
  void SetInclusiveSigmaCuts(Int_t pidFunction, Double_t sigmaValue)
  {
   // pidFunction: [0=Pion(...),1=Kaon(...),2=Proton(...)]
   fUseDefaultInclusiveSigmaCuts = kFALSE;
   this->fInclusiveSigmaCuts[pidFunction] = sigmaValue;
  };
  void SetExclusiveSigmaCuts(Int_t pidFunction, Int_t pidExclusive, Double_t sigmaValue)
  {
   // pidFunction: [0=Pion(...),1=Kaon(...),2=Proton(...)]
   fUseDefaultExclusiveSigmaCuts = kFALSE;
   this->fExclusiveSigmaCuts[pidFunction][pidExclusive] = sigmaValue;
  };


  // Set booleans ===============================================================================================================
  void SetRejectEventsNoPrimaryVertex(Bool_t primvertex) { this->fRejectEventsNoPrimaryVertex = primvertex;};
  void SetCutOnVertexZ(Bool_t vertexz) { this->fCutOnVertexZ = vertexz;};
  void SetApplyCommonTrackCuts(Bool_t commontrackcuts) { this->fApplyCommonTrackCuts = commontrackcuts;};
  void SetUseDefaultInclusiveSigmaCuts(Bool_t inclusivecuts) { this->fUseDefaultInclusiveSigmaCuts = inclusivecuts;};
  void SetUseDefaultExclusiveSigmaCuts(Bool_t exclusivecuts) { this->fUseDefaultExclusiveSigmaCuts = exclusivecuts;};
  void SetRejectFakeTracks(Bool_t rft) { this->fRejectFakeTracks = rft;};
  void SetProcessBothKineAndReco(Bool_t pbkar) { this->fProcessBothKineAndReco = pbkar;};

  // Set other things ===========================================================================================================
  void SetCentralityEstimator(TString estimatorName)
  {
   this->fCentralityEstimator = estimatorName;
  };
  void SetFilterBit(Int_t filterBitCheck)
  {
   this->fFilterBit = filterBitCheck;
  };
  


   private:
  AliAnalysisTaskNBodyFemtoscopy(const AliAnalysisTaskNBodyFemtoscopy& aatmpf);
  AliAnalysisTaskNBodyFemtoscopy& operator=(const AliAnalysisTaskNBodyFemtoscopy& aatmpf);
  
  AliPIDResponse *fPIDResponse; //! PID response object
  // 0.) Base lists, index array, ...:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)
  TArrayI *fRandomIndices;  // an array to keep randomized indices

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
  TH1F *fNumberOfTracksHistAfterAllCuts;         // count tracks that passed selection
  // Test PID 
  TH2F *fTestPIDTrueFalsePositive;         // check how many true positives 
  // Test global to TPC
  TH1F *fTestTPCOnlyVsGlobal;         // check how many false positives 

  // Binnning
  Int_t fNbinsPt;                // pt bin number
  Float_t fMinBinPt;             // min bin pt
  Float_t fMaxBinPt;             // min bin pt
  Int_t fNbinsCentrality;      // cenrality bin number
  Float_t fMinCentrality;        // min bin centrality
  Float_t fMaxCentrality;        // max bin centrality
  Int_t fNbinsMultiplicity;             // Multiplicity bin number
  Float_t fMinBinMultiplicity;             // min bin multiplicity
  Float_t fMaxBinMultiplicity;             // max bin multiplicity 
  Int_t fNbinsTPCOnlyVsGlobal;             // TPCOnlyVsGlobal bin number
  Float_t fMinBinTPCOnlyVsGlobal;             // min bin TPCOnlyVsGlobal 
  Float_t fMaxBinTPCOnlyVsGlobal;             // max bin TPCOnlyVsGlobal 

  // Limits and ranges
  Float_t fPtRange[2];                // ptMin = fPtRange[0], ptMax = fPtRange[1]
  Float_t fEtaRange[2];               // etaMin = etaRange[0], etaMax = etaRange[1]
  Float_t fPhiRange[2];               // phiMin = phiRange[0], phiMax = phiRange[1]
  Float_t fVertexZ[2];                // VertexZMin = fVertexZ[0], VertexZMax = fVertexZ[1]
  Double_t fInclusiveSigmaCuts[3];           // [PID function] 
  Double_t fExclusiveSigmaCuts[3][3];        // [PID function][PID exclusive]

  // Booleans
  Bool_t fRejectEventsNoPrimaryVertex; //Reject events without primary vertex
  Bool_t fCutOnVertexZ;                //Apply vertexZ cut
  Bool_t fApplyCommonTrackCuts;       //Apply common track cuts
  Bool_t fUseDefaultInclusiveSigmaCuts;		// use defaults 3 sigma cut for pid inclusive
  Bool_t fUseDefaultExclusiveSigmaCuts;		// use defaults 3 sigma cut for pid exclusive
  Bool_t fRejectFakeTracks; // if set to kFALSE, and if fMC is available, get the corresponding MC particle by taking absolute value of label
  Bool_t fProcessBothKineAndReco;

  // Other things
  TString fCentralityEstimator;  // centrality estimator choice
  Int_t fFilterBit;                   //  Filter bit choice
  
  // No setters, getters required
  TExMap *fGlobalTracksAODTEST[3];
  TClonesArray *fAllTracksTEST[3];   
  AliMCEvent *fMC; // placeholder for MC info
  

  // 2.) Final results:
  TList *fFinalResultsList; // list to hold all histograms with final results

  // Increase this counter in each new version:
  ClassDef(AliAnalysisTaskNBodyFemtoscopy, 5);

};

//================================================================================================================

#endif




