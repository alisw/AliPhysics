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
class TBrowser;
#include "TList.h"
#include "TH1.h"
#include "AliTriggerAnalysis.h"
#include "AliFlowTrackCuts.h"

class AliFlowEventCuts : public TNamed {

 public:
  enum refMultMethod { kTPConly, kSPDtracklets, kV0, kSPD1clusters };

  AliFlowEventCuts();
  AliFlowEventCuts(const char* name, const char* title = "AliFlowEventCuts");
  AliFlowEventCuts(const AliFlowEventCuts& someCuts);
  AliFlowEventCuts& operator=(const AliFlowEventCuts& someCuts);
  virtual  ~AliFlowEventCuts();
  
  virtual Bool_t IsSelected(TObject* obj);

  Bool_t PassesCuts(AliVEvent* event);
  
  static AliFlowEventCuts* StandardCuts();
  
  void SetNumberOfTracksMax(Int_t value) {fNumberOfTracksMax=value;fCutNumberOfTracks=kTRUE;}
  void SetNumberOfTracksMin(Int_t value) {fNumberOfTracksMin=value;fCutNumberOfTracks=kTRUE;}
  void SetNumberOfTracksRange(Int_t min, Int_t max) {fNumberOfTracksMin=min;fNumberOfTracksMax=max;fCutNumberOfTracks=kTRUE;}
  void SetRefMultMax(Int_t value) {fRefMultMax=value;fCutRefMult=kTRUE;}
  void SetRefMultMin(Int_t value) {fRefMultMin=value;fCutRefMult=kTRUE;}
  void SetRefMultRange(Int_t min, Int_t max) {fRefMultMin=min;fRefMultMax=max;fCutRefMult=kTRUE;}
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
  void SetCutTPCmultiplicityOutliers(Bool_t b=kTRUE) {fCutTPCmultiplicityOutliers=b;}

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

  Int_t RefMult(AliVEvent* event);
  //Int_t GetRefMult() {return fRefMult;}
  Int_t GetReferenceMultiplicity(AliVEvent* event) {return RefMult(event);}
  const char* CentrMethName(refMultMethod method) const;
  void SetCentralityPercentileRange(Float_t min, Float_t max){ fCentralityPercentileMin=min;
                                                               fCentralityPercentileMax=max;
                                                               fCutCentralityPercentile=kTRUE; }
  void SetCentralityPercentileMethod( refMultMethod m) {fCentralityPercentileMethod=m;}
  void SetUseCentralityUnchecked(Bool_t b=kTRUE) {fUseCentralityUnchecked=b;}
  
  void Browse(TBrowser* b);
  Long64_t Merge(TCollection* list);  

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
  Bool_t fCutTPCmultiplicityOutliers; //cut TPC multiplicity outliers
  Bool_t fCutCentralityPercentile; //cut on centrality perc. from AliESDCentrality
  Bool_t fUseCentralityUnchecked; //use the unchecked method
  refMultMethod fCentralityPercentileMethod; //where to get the percentile from
  Float_t fCentralityPercentileMax; // max centr. perc
  Float_t fCentralityPercentileMin; // min centr. perc
  Bool_t fCutZDCtiming;   //cut on ZDC timing
  AliTriggerAnalysis fTrigAna; //trigger analysis object

  ClassDef(AliFlowEventCuts,4)
};

#endif


