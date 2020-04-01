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
#include "AliReducedInfoCut.h"

class AliMixingHandler : public TNamed {
   
public:
   enum Constants {
      kMixResonanceLegs=0,         // event mixing for resonance inv mass bkg
      kMixCorrelation,                 // event mixing for correlations
      kNMaxVariables = 10
   };

public:
  AliMixingHandler(Int_t mixingSetup=kMixResonanceLegs);
  AliMixingHandler(const Char_t* name, const Char_t* title, Int_t mixingSetup=kMixResonanceLegs);
  virtual ~AliMixingHandler();
  
  // setters
  void AddMixingVariable(AliReducedVarManager::Variables var, Int_t nBins, const Float_t* binLims);
  void AddMixingVariable(AliReducedVarManager::Variables var, Int_t nBins, const Double_t* binLims);
  void SetMixLikeSign(Bool_t flag) {fMixLikeSign = flag;}
  void SetMixLikePairs(Bool_t flag) {SetMixLikeSign(flag);}      // synonim function to SetMixLikeSign
  void SetPoolDepth(Int_t n) {fPoolDepth = n;}
  void SetMixingThreshold(Float_t fr) {fMixingThreshold = fr;}
  void SetDownscaleEvents(Float_t ds) {fDownscaleEvents = ds;}
  void SetDownscaleTracks(Float_t ds) {fDownscaleTracks = ds;}
  void SetNParallelCuts(Int_t n) {fNParallelCuts = n;}
  void SetNParallelPairCuts(Int_t n) {fNParallelPairCuts = n;}
  void SetHistogramManager(AliHistogramManager* histos) {fHistos = histos;}
  void SetHistClassNames(const Char_t* names) {fHistClassNames = names;}
  void AddCrossPairsCut(AliReducedInfoCut* cut) {fCrossPairsCuts.Add(cut);}
  void AddOppositeSignPairsCut(AliReducedInfoCut* cut) {fCrossPairsCuts.Add(cut);}    // synonim function to AddCrossPairsCut() used for charged legs
  void AddLikePairsLeg1Cut(AliReducedInfoCut* cut) {fLikePairsLeg1Cuts.Add(cut);}
  void AddLikeSignPairsPPCut(AliReducedInfoCut* cut) {fLikePairsLeg1Cuts.Add(cut);}  // synonim function to AddLikePairsLeg1Cut() used for charged legs
  void AddLikePairsLeg2Cut(AliReducedInfoCut* cut) {fLikePairsLeg2Cuts.Add(cut);}
  void AddLikeSignPairsMMCut(AliReducedInfoCut* cut) {fLikePairsLeg2Cuts.Add(cut);}  // synonim function to AddLikePairsLeg2Cut() used for charged legs
  void AddPairsCut(AliReducedInfoCut* cut) {
    fCrossPairsCuts.Add(cut);
    fLikePairsLeg1Cuts.Add(cut);
    fLikePairsLeg2Cuts.Add(cut);
  }
  
  // getters
  Int_t GetDepth() const {return fPoolDepth;}
  Float_t GetMixingThreshold() const {return fMixingThreshold;}
  Float_t GetDownscaleEvents() const {return fDownscaleEvents;}
  Float_t GetDownscaleTracks() const {return fDownscaleTracks;}
  Int_t GetNParallelCuts() const {return fNParallelCuts;}
  Int_t GetNParallelPairCuts() const {return fNParallelPairCuts;}
  Int_t GetPoolSize(Int_t cut, Float_t* values);
  Int_t GetPoolSize(Int_t cut, Int_t eventCategory) const;
  TString GetHistClassNames() const {return fHistClassNames;};
  Int_t GetNMixingVariables() const {return fNMixingVariables;}
  Int_t GetMixingSetup() const {return fMixingSetup;}
  
  void Init();
  Int_t FindEventCategory(Float_t* values);
  Int_t GetBinFromCategory(Int_t iVar, Int_t category) const;
  void FillEvent(TList* leg1List, TList* leg2List, Float_t* values, Int_t type=-1);
  Bool_t AcceptTrack();    // randomly accept/reject a track for mixing
  void RunLeftoverMixing(Int_t type=-1);
  void PrintMixingLists(Int_t debug);  
  ULong_t IsPairSelected(Float_t* values, Int_t pairType);
  
private:
   AliMixingHandler(const AliMixingHandler& handler);             
   AliMixingHandler& operator=(const AliMixingHandler& handler);      
   
  // User options
  Int_t    fMixingSetup;          //  see Constants for various options 
  Int_t fPoolDepth;              // depth of the event mixing pool
  // TODO: Add option to trigger the event mixing when a certain number of tracks in a given pool is reached
  // TODO: Add option for rolling buffer mixing
  Float_t fMixingThreshold;      // within a (centrality,vtx,ep) mix all pools with entries > fMixingThreshold*fPoolDepth
  Float_t fDownscaleEvents;      // random downscale adding events to the pools
  Float_t fDownscaleTracks;      // random downscale adding tracks fo the pools
  
  TClonesArray fPoolsLeg1;         // array of pools
  TClonesArray fPoolsLeg2;         // array of pools
  Int_t fNParallelCuts;            // number of parallel cuts which are run
  Int_t fNParallelPairCuts;        // number of parallel pair cuts which are run
  TString fHistClassNames;         // name of the histogram classes for each cut, separated by a semicolon ";"
  TArrayI fPoolSize;               // counters for the pool sizes
  Bool_t fIsInitialized;           // check if the mixing handler is initialized
  Bool_t fMixLikeSign;             // mix or not like-sign tracks (default is true)
  
  TArrayF fVariableLimits[kNMaxVariables];
  AliReducedVarManager::Variables fVariables[kNMaxVariables];
  Int_t  fNMixingVariables;
  
  AliHistogramManager* fHistos;    // histogram manager
  
  TList fCrossPairsCuts;         // cut object for cross pairs 
  TList fLikePairsLeg1Cuts;    // cut object for LEG1 like pairs
  TList fLikePairsLeg2Cuts;    // cut object for LEG2 like pairs
  
  void RunEventMixing(TClonesArray* leg1Pool, TClonesArray* leg2Pool, ULong_t mixingMask, Int_t type, Float_t* values);
  ULong_t IncrementPoolSizes(TList* list1, TList* list2, Int_t eventCategory);
  void ResetPoolSizes(ULong_t mixingMask, Int_t category);  
  
  ClassDef(AliMixingHandler,4);
};

#endif
