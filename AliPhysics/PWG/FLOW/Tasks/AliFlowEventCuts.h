/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowEventCuts:
// An event cut class
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#ifndef ALIFLOWEVENTCUTS_H
#define ALIFLOWEVENTCUTS_H

#include <float.h>
#include <limits.h>
#include "TNamed.h"

class AliVEvent;
class AliMCEvent;
class TBrowser;
class AliAnalysisUtils;
class AliMultSelection;
#include "TList.h"
#include "TH1.h"
#include "AliTriggerAnalysis.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimpleCuts.h"

class AliFlowEventCuts : public AliFlowEventSimpleCuts {

 public:
  enum refMultMethod { kTPConly, kSPDtracklets, kVZERO, kV0=kVZERO, kV0M=kVZERO, kSPD1clusters, kZDC, kV0C, kV0A};

  AliFlowEventCuts();
  AliFlowEventCuts(const char* name, const char* title = "AliFlowEventCuts");
  AliFlowEventCuts(const AliFlowEventCuts& someCuts);
  AliFlowEventCuts& operator=(const AliFlowEventCuts& someCuts);
  virtual  ~AliFlowEventCuts();
   
  virtual Bool_t IsSelected(TObject* obj, TObject *objmc);

  Bool_t PassesCuts(AliVEvent* event, AliMCEvent *mcevent);
  
  static AliFlowEventCuts* StandardCuts();
  
  void SetNumberOfTracksMax(Int_t value) {fNumberOfTracksMax=value;fCutNumberOfTracks=kTRUE;}
  void SetNumberOfTracksMin(Int_t value) {fNumberOfTracksMin=value;fCutNumberOfTracks=kTRUE;}
  void SetNumberOfTracksRange(Int_t min, Int_t max) {fNumberOfTracksMin=min;fNumberOfTracksMax=max;fCutNumberOfTracks=kTRUE;}
  void SetRefMultMax(Int_t value) {fRefMultMax=value;fCutRefMult=kTRUE;}
  void SetRefMultMin(Int_t value) {fRefMultMin=value;fCutRefMult=kTRUE;}
  void SetRefMultRange(Int_t min, Int_t max) {fRefMultMin=min;fRefMultMax=max;fCutRefMult=kTRUE;}
  void SetImpactParameterMax(Double_t value) {fImpactParameterMax=value;fCutImpactParameter=kTRUE;}
  void SetImpactParameterMin(Double_t value) {fImpactParameterMin=value;fCutImpactParameter=kTRUE;}
  void SetImpactParameterRange(Double_t min, Double_t max) {fImpactParameterMin=min;fImpactParameterMax=max;fCutImpactParameter=kTRUE;}
  void SetPrimaryVertexXrange(Double_t min, Double_t max)
       {fCutPrimaryVertexX=kTRUE; fPrimaryVertexXmin=min; fPrimaryVertexXmax=max;}
  void SetPrimaryVertexYrange(Double_t min, Double_t max)
       {fCutPrimaryVertexY=kTRUE; fPrimaryVertexYmin=min; fPrimaryVertexYmax=max;}
  void SetPrimaryVertexZrange(Double_t min, Double_t max)
       {fCutPrimaryVertexZ=kTRUE; fPrimaryVertexZmin=min; fPrimaryVertexZmax=max;}
  void SetNContributorsRange(Int_t min, Int_t max=INT_MAX) 
       {fCutNContributors=kTRUE; fNContributorsMin=min; fNContributorsMax=max;}
  void SetMeanPtRange(Double_t min, Double_t max) {fCutMeanPt=kTRUE; fMeanPtMax=max; fMeanPtMin=min;}
  void SetCutSPDvertexerAnomaly(Bool_t b=kTRUE) {fCutSPDvertexerAnomaly=b;}
  void SetCutZDCtiming(Bool_t c=kTRUE) {fCutZDCtiming=c;}
  void SetCutSPDTRKVtxZ(Bool_t b=kTRUE) {fCutSPDTRKVtxZ=b;}
  void SetCutTPCmultiplicityOutliers(Bool_t b=kTRUE) {fCutTPCmultiplicityOutliers=b;}  
  void SetCutTPCmultiplicityOutliersAOD(Bool_t b=kTRUE) {fCutTPCmultiplicityOutliersAOD=b;}

  //================================================//
   void SetSphericityCut(Double_t sphericityMin, Double_t sphericityMax) {
     fSphericityMin = sphericityMin;
     fSphericityMax = sphericityMax;
     fUseSphericityCut = kTRUE;
   }
   void SelectEventsWithHighPtParticles(Double_t gPt) {
     fMinPtOfTrigger = gPt;
     fUseMinPtOfTrigger = kTRUE;
   }
  //================================================//

  Int_t GetNumberOfTracksMax() const {return fNumberOfTracksMax;}
  Int_t GetNumberOfTracksMin() const {return fNumberOfTracksMin;}
  Int_t GetRefMultMax() const {return fRefMultMax;}
  Int_t GetRefMultMin() const {return fRefMultMin;}
  void SetRefMultMethod(refMultMethod m) {fRefMultMethod=m;}
  void SetRefMultMethod(AliESDtrackCuts::MultEstTrackType m) { fRefMultMethodAliESDtrackCuts=m; 
                                                               fUseAliESDtrackCutsRefMult=kTRUE; }
  refMultMethod GetRefMultMethod() const {return fRefMultMethod;}
  void SetRefMultCuts( AliFlowTrackCuts* cuts ) {fRefMultCuts=static_cast<AliFlowTrackCuts*>(cuts->Clone());}
  void SetMeanPtCuts( AliFlowTrackCuts* cuts ) {fMeanPtCuts=static_cast<AliFlowTrackCuts*>(cuts->Clone());}
  AliFlowTrackCuts* GetRefMultCuts() const {return fRefMultCuts;}
  void DefineHistograms();
  void SetQA(Bool_t b=kTRUE) {if (b) DefineHistograms();}
  TList* GetQA() const {return fQA;}
  TH1* QAbefore(Int_t i) {return static_cast<TH1*>(static_cast<TList*>(fQA->At(0))->At(i));}
  TH1* QAafter(Int_t i) {return static_cast<TH1*>(static_cast<TList*>(fQA->At(1))->At(i));}

  Int_t RefMult(AliVEvent* event, AliMCEvent *mcEvent = 0x0);
  //Int_t GetRefMult() {return fRefMult;}
  Int_t GetReferenceMultiplicity(AliVEvent* event, AliMCEvent *mcEvent) {return RefMult(event,mcEvent);}
  const char* CentrMethName(refMultMethod method) const;
  void SetCentralityPercentileMethod( refMultMethod m) {fCentralityPercentileMethod=m;}
  void SetUseCentralityUnchecked(Bool_t b=kTRUE) {fUseCentralityUnchecked=b;}

  Float_t GetCentrality(AliVEvent* event, AliMCEvent* mcEvent);
  void SetUsedDataset(Bool_t b=kTRUE) {fData2011=b;}    // confusing name, better use different interface
  void SetLHC10h(Bool_t b=kTRUE) {fData2011=(!b);}      // TODO let cut object determine runnumber and period
  void SetLHC11h(Bool_t b=kTRUE) {fData2011=b;}         // use this only as 'manual override'
  void SetCheckPileup(Bool_t b=kFALSE) {fCheckPileUp = b;}  // In case a pile-up rejection is required

  void Browse(TBrowser* b);
  Long64_t Merge(TCollection* list);  
  TH2F *GetCorrelationTPCvsGlobalMultiplicity() {return fhistTPCvsGlobalMult;}  


 private:
  TList* fQA; //QA
  Bool_t fCutNumberOfTracks;//cut on # of tracks
  Int_t fNumberOfTracksMax;  //limits
  Int_t fNumberOfTracksMin;  //limits
  Bool_t fCutRefMult; //cut on refmult
  refMultMethod fRefMultMethod; //how do we calculate refmult?
  Bool_t fUseAliESDtrackCutsRefMult; //use AliESDtrackCuts for refmult calculation
  AliESDtrackCuts::MultEstTrackType fRefMultMethodAliESDtrackCuts;
  Int_t fRefMultMax; //max refmult
  Int_t fRefMultMin; //min refmult
  AliFlowTrackCuts* fRefMultCuts; //cuts
  AliFlowTrackCuts* fMeanPtCuts; //mean pt cuts
  AliFlowTrackCuts* fStandardTPCcuts; //Standard TPC cuts
  AliFlowTrackCuts* fStandardGlobalCuts; //StandardGlobalCuts
  AliAnalysisUtils* fUtils;  //! analysis utils object
  AliMultSelection* fMultSelection; //! new centrality framework
  Bool_t fCutPrimaryVertexX; //cut on x of prim vtx
  Double_t fPrimaryVertexXmax; //max x prim vtx
  Double_t fPrimaryVertexXmin; //min x prim vtx
  Bool_t fCutPrimaryVertexY; //cut on y of prim vtx
  Double_t fPrimaryVertexYmax; //max y prim vtx
  Double_t fPrimaryVertexYmin; //min y prim vtx
  Bool_t fCutPrimaryVertexZ; //cut on z of prim vtx
  Double_t fPrimaryVertexZmax; //max z prim vtx
  Double_t fPrimaryVertexZmin; //min z prim vtx
  Bool_t fCutNContributors; //cut on number of contributors
  Int_t fNContributorsMax; //maximal number of contrib
  Int_t fNContributorsMin; //minimal number of contrib
  Bool_t fCutMeanPt; //cut on mean pt
  Double_t fMeanPtMax; //max mean pt
  Double_t fMeanPtMin; //min mean pt
  Bool_t fCutSPDvertexerAnomaly; //cut on the spd vertexer anomaly
  Bool_t fCutSPDTRKVtxZ; //require compatibility between SPDvertexz TRKvertexz
  Bool_t fCutTPCmultiplicityOutliers; //cut TPC multiplicity outliers
  Bool_t fCutTPCmultiplicityOutliersAOD; // cut TPC outliers in 10h or 11h aod
  Bool_t fUseCentralityUnchecked; //use the unchecked method
  refMultMethod fCentralityPercentileMethod; //where to get the percentile from
  Bool_t fCutZDCtiming;   //cut on ZDC timing
  AliTriggerAnalysis fTrigAna; //trigger analysis object
  Bool_t fCutImpactParameter; //cut on impact parameter (MC header)
  Double_t fImpactParameterMin; // min impact parameter
  Double_t fImpactParameterMax; // max impact parameter
  TH2F *fhistTPCvsGlobalMult; //!correlation between TPCMult and GlobalMult
  Bool_t fData2011; //2011 data is used
  Bool_t fCheckPileUp; //pile-up

  Double_t fSphericityMin;//sphericity min cut (currently only for AODs)
  Double_t fSphericityMax;//sphericity max cut (currently only for AODs)
  Bool_t fUseSphericityCut;//sphericity cut (currently only for AODs)

  Double_t fMinPtOfTrigger;//min pT of the trigger particle within the event
  Bool_t fUseMinPtOfTrigger;//usage of min pT cut for the trigger

  ClassDef(AliFlowEventCuts,8)
};

#endif


