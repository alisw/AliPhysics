//
// Author: Ionut-Cristian Arsene, 2015/08/07
// email: iarsene@cern.ch
//
// Event mixing handler
//
#ifndef ALIMIXINGHANDLER_H
#define ALIMIXINGHANDLER_H

#include <TNamed.h>
#include <TArrayF.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TString.h>

#include "AliHistogramManager.h"
#include "AliReducedVarManager.h"

class AliMixingHandler : public TNamed {

public:
  AliMixingHandler();
  AliMixingHandler(const Char_t* name, const Char_t* title);
  virtual ~AliMixingHandler();
  
  // setters
  void SetMixLikeSign(Bool_t flag) {fMixLikeSign = flag;}
  void SetPoolDepth(Int_t n) {fPoolDepth = n;}
  void SetMixingThreshold(Float_t fr) {fMixingThreshold = fr;}
  void SetDownscaleEvents(Float_t ds) {fDownscaleEvents = ds;}
  void SetDownscaleTracks(Float_t ds) {fDownscaleTracks = ds;}
  void SetNParallelCuts(Int_t n) {fNParallelCuts = n;}
  void SetCentralityLimits(Int_t n, const Float_t* arr)  {fCentralityLimits.Set(n,arr);}
  void SetEventVertexLimits(Int_t n, const Float_t* arr) {fEventVertexLimits.Set(n,arr);}
  void SetEventPlaneLimits(Int_t n, const Float_t* arr)  {fEventPlaneLimits.Set(n,arr);}
  void SetHistogramManager(AliHistogramManager* histos) {fHistos = histos;}
  void SetHistClassNames(const Char_t* names) {fHistClassNames = names;}
  void SetEventVariables(AliReducedVarManager::Variables centVar, AliReducedVarManager::Variables vtxVar, AliReducedVarManager::Variables epVar);
  
  // getters
  Int_t GetDepth() const {return fPoolDepth;}
  Float_t GetMixingThreshold() const {return fMixingThreshold;}
  Float_t GetDownscaleEvents() const {return fDownscaleEvents;}
  Float_t GetDownscaleTracks() const {return fDownscaleTracks;}
  Int_t GetNParallelCuts() const {return fNParallelCuts;}
  Int_t GetPoolSize(Int_t cut, Float_t centrality, Float_t vtxz, Float_t ep);
  Int_t GetPoolSize(Int_t cut, Int_t eventCategory);
  TString GetHistClassNames() const {return fHistClassNames;};
  
  void Init();
  Int_t FindEventCategory(Float_t centrality, Float_t vtxz, Float_t ep);
  Int_t GetEventPlaneBin(Int_t category);
  Int_t GetEventVertexBin(Int_t category);
  Int_t GetCentralityBin(Int_t category);
  void FillEvent(TList* leg1List, TList* leg2List, Float_t* values, Int_t type);
  Bool_t AcceptTrack();    // randomly accept/reject a track for mixing
  void RunLeftoverMixing(Int_t type);
  void PrintMixingLists(Int_t debug);  
  
private:
   AliMixingHandler(const AliMixingHandler& handler);             
   AliMixingHandler& operator=(const AliMixingHandler& handler);      
   
  // User options
  Int_t fPoolDepth;              // depth of the event mixing pool
  Float_t fMixingThreshold;      // within a (centrality,vtx,ep) mix all pools with entries > fMixingThreshold*fPoolDepth
  Float_t fDownscaleEvents;      // random downscale adding events to the pools
  Float_t fDownscaleTracks;      // random downscale adding tracks fo the pools
  
  TClonesArray fPoolsLeg1;         // array of pools
  TClonesArray fPoolsLeg2;         // array of pools
  Int_t fNParallelCuts;            // number of parallel cuts which are run
  TString fHistClassNames;         // name of the histogram classes for each cut, separated by a semicolon ";"
  TArrayI fPoolSize;               // counters for the pool sizes
  Bool_t fIsInitialized;           // check if the mixing handler is initialized
  Bool_t fMixLikeSign;             // mix or not like-sign tracks (default is true)
  
  TArrayF fCentralityLimits;
  TArrayF fEventVertexLimits;
  TArrayF fEventPlaneLimits;
  AliReducedVarManager::Variables fCentralityVariable;
  AliReducedVarManager::Variables fEventVertexVariable;
  AliReducedVarManager::Variables fEventPlaneVariable;
  
  AliHistogramManager* fHistos;    // histogram manager
  
  void RunEventMixing(TClonesArray* leg1Pool, TClonesArray* leg2Pool, ULong_t mixingMask, Int_t type, Float_t* values);
  ULong_t IncrementPoolSizes(TList* list1, TList* list2, Int_t eventCategory);
  void ResetPoolSizes(ULong_t mixingMask, Int_t category);  
  
  ClassDef(AliMixingHandler,1);
};

#endif
