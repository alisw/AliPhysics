// Class for cutting on ALICE Var manager and other track specific information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
//   07/09/2016

#ifndef ALIREDUCEDTRACKCUT_H
#define ALIREDUCEDTRACKCUT_H

#include "AliReducedVarCut.h"

//_________________________________________________________________________
class AliReducedTrackCut : public AliReducedVarCut {

 public:
  AliReducedTrackCut();
  AliReducedTrackCut(const Char_t* name, const Char_t* title);
  virtual ~AliReducedTrackCut();

  void SetRejectKinks(Bool_t flag=kTRUE) {fRejectKinks = flag;}
  void SetRejectTaggedGamma(Bool_t flag = kTRUE) {fRejectTaggedGamma = flag;}
  void SetRejectTaggedPureGamma(Bool_t flag = kTRUE) {fRejectTaggedPureGamma = flag;}
  void SetRequestITSrefit(Bool_t flag = kTRUE) {fRequestITSrefit = flag;}
  void SetRequestTPCrefit(Bool_t flag = kTRUE) {fRequestTPCrefit = flag;}
  void SetITShitRequest(UChar_t hitMap, Bool_t useAND = kFALSE) {fCutOnITShitMap = hitMap; fUseANDonITShitMap = useAND; fRequestCutOnITShitMap = kTRUE;}
  void SetRequestITSanyCut(Int_t layer) {for(Int_t i=0;i<=layer;++i) fCutOnITShitMap |= (1<<i); fUseANDonITShitMap = kFALSE; fRequestCutOnITShitMap = kTRUE;}
  void SetRequestSPDfirst() {fCutOnITShitMap |= (1<<0); fUseANDonITShitMap = kFALSE; fRequestCutOnITShitMap = kTRUE;}
  void SetRequestSPDany() {fCutOnITShitMap |= (1<<0); fCutOnITShitMap |= (1<<1); fUseANDonITShitMap = kFALSE; fRequestCutOnITShitMap = kTRUE;}
  void SetRequestSPDboth() {fCutOnITShitMap |= (1<<0); fCutOnITShitMap |= (1<<1); fUseANDonITShitMap = kTRUE; fRequestCutOnITShitMap = kTRUE;}
  void SetITSsharedClsRejectMap(UChar_t sharedMap, Bool_t useAND=kFALSE) {fCutOnITSsharedClsMap = sharedMap; fUseANDonITSsharedClsMap = useAND; fRequestCutOnITSsharedClsMap = kTRUE;}
  void SetRejectSharedClsOnLayer(Int_t layer) {if(layer<0 || layer>5) return; fCutOnITSsharedClsMap |= (1<<layer); fUseANDonITSsharedClsMap = kFALSE; fRequestCutOnITSsharedClsMap = kTRUE;}
  void SetRejectSharedClsSPDany() {fCutOnITSsharedClsMap |= (1<<0); fCutOnITSsharedClsMap |= (1<<1); fUseANDonITSsharedClsMap = kFALSE; fRequestCutOnITSsharedClsMap = kTRUE;}
  void SetRejectSharedClsSPDboth() {fCutOnITSsharedClsMap |= (1<<0); fCutOnITSsharedClsMap |= (1<<1); fUseANDonITSsharedClsMap = kTRUE; fRequestCutOnITSsharedClsMap = kTRUE;}
  void SetRequestTOFout(Bool_t flag = kTRUE) {fRequestTOFout = flag;}  
  void SetRequestTRDmatch(Bool_t flag = kTRUE) {fRequestTRDonlineMatch = flag;}
  
  void SetTrackFilterMap(UInt_t map, Bool_t useAND=kTRUE) {fCutOnTrackFilterMap = map; fUseANDonTrackFilterMap=useAND;}
  void SetTrackFilterBit(Int_t bit) {if(bit>=0 && bit<32) fCutOnTrackFilterMap |= (UInt_t(1)<<bit);}
  void SetUseANDonTrackFilterMap(Bool_t useAND=kTRUE) {fUseANDonTrackFilterMap=useAND;}
  
  void SetTrackQualityFilterMap(UInt_t map, UInt_t exclusionMap=0, Bool_t useAND=kTRUE) {
     fCutOnTrackQualityMap = map; fCutOnTrackQualityMapExclude=exclusionMap; fUseANDonTrackQualityMap=useAND;}
  void SetTrackQualityFilterBit(Int_t bit, Bool_t exclude=kFALSE) {
     if(bit>=0 && bit<32) {
        fCutOnTrackQualityMap |= (UInt_t(1)<<bit);
        if(exclude) fCutOnTrackQualityMapExclude |= (UInt_t(1)<<bit);
     }
   }
  void SetUseANDonTrackQualityFilterMap(Bool_t useAND=kTRUE) {fUseANDonTrackQualityMap=useAND;}
  
  void SetRejectPureMC(Bool_t reject=kTRUE) {fRejectPureMC=reject;}
  void SetApplyReweightMCpt(Bool_t reweight=kTRUE) {fReweightMCpt=reweight;}
    
  void SetMCFilterMap(UInt_t map, Bool_t useAND=kTRUE) {fCutOnMCFilterMap = map; fUseANDonMCFilterMap=useAND;}
  void SetMCFilterBit(Int_t bit) {if(bit>=0 && bit<32) fCutOnMCFilterMap |= (UInt_t(1)<<bit);}
  void SetUseANDonMCFilterMap(Bool_t useAND=kTRUE) {fUseANDonMCFilterMap=useAND;}
  
  Bool_t GetRejectKinks() const {return fRejectKinks;}
  Bool_t GetRejectTaggedGamma() const {return fRejectTaggedGamma;}
  Bool_t GetRejectTaggedPureGamma() const {return fRejectTaggedPureGamma;}
  UInt_t GetTrackFilterMap() const {return fCutOnTrackFilterMap;}
  Bool_t GetUseANDonTrackFilterMap() const {return fUseANDonTrackFilterMap;}
  UInt_t GetTrackQualityFilterMap() const {return fCutOnTrackQualityMap;}
  Bool_t GetUseANDonTrackQualityFilterMap() const {return fUseANDonTrackQualityMap;}
  Bool_t GetRejectPureMC() const {return fRejectPureMC;}
  Bool_t GetApplyReweightMCpt() const {return fReweightMCpt;}
  UInt_t GetMCFilterMap() const {return fCutOnMCFilterMap;}
  Bool_t GetUseANDonMCFilterMap() const {return fUseANDonMCFilterMap;}
  Bool_t GetRequestITSrefit() const {return fRequestITSrefit;}
  Bool_t GetRequestTPCrefit() const {return fRequestTPCrefit;}
  UChar_t GetITShitMapRequest() const {return fCutOnITShitMap;}
  Bool_t GetUseANDonITShitMap() const {return fUseANDonITShitMap;}
  Bool_t GetUseCutOnITShitMap() const {return fRequestCutOnITShitMap;}
  Bool_t GetRequestTOFout() const {return fRequestTOFout;}
  Bool_t GetRequestTRDmatch() const {return fRequestTRDonlineMatch;}
  
  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(TObject* obj, Float_t* values);
  
 protected: 
      
  // Cuts on track specific quantities
  // global track quantities
  Bool_t    fRejectKinks;                       // if true, reject kinks
  Bool_t    fRejectTaggedGamma;       // if true, reject tagged gamma conversions
  Bool_t    fRejectTaggedPureGamma;  // if true, reject only the high purity tagged gamma conversions
  UInt_t    fCutOnTrackQualityMap;      // map encoding requests on the quality of the track (see bits in AliReducedBaseTrack::fQualityFlags)
  UInt_t    fCutOnTrackQualityMapExclude;    // if bits are enabled, corresponding requests on fCutOnTrackQualityMap are negated
  Bool_t    fUseANDonTrackQualityMap;   // if false, apply an OR on enabled positions; if true apply AND
  
  // selections on the track filter map 
  UInt_t    fCutOnTrackFilterMap;         // map encoding the various requests on track filters
  Bool_t    fUseANDonTrackFilterMap;  // if false, at least one of the enabled positions in the cut map should be on in the track filters map
   
   // Reject pure MC tracks
   Bool_t   fRejectPureMC;
   //Apply or not re-weighting on MC
   Bool_t fReweightMCpt;
   
   // selections on the MC signals bit map 
   UInt_t    fCutOnMCFilterMap;         // map encoding the various requests on MC signals
   Bool_t    fUseANDonMCFilterMap;  // if false, at least one of the enabled positions in the cut map should be on in the track MC signals map
   
   // ITS quantities
  Bool_t    fRequestITSrefit;                // if true, request kITSrefit flag to be on
  UChar_t   fCutOnITShitMap;               // hit map encoding various requests on cluster configurations
  Bool_t    fUseANDonITShitMap;        // if false, at least one of the enabled positions in the cut map should have corresponding clusters
                                                          // if true, all the enabled bits should have matching clusters
                                                          // The bits left to zero will not be used when evaluating the cut
                                                          // e.g.  100000 and FALSE: equivalent to SPDfirst
                                                          // e.g.  110000 and FALSE: equivalent to SPDany (at least one hit in the first two layers)
                                                          // e.g.  110000 and TRUE: equivalent to requesting both SPD layers to have a hit
   Bool_t   fRequestCutOnITShitMap;  // if true, apply the cut above
   UChar_t  fCutOnITSsharedClsMap;       // map encoding various requests on shared cluster configurations to reject
   Bool_t   fUseANDonITSsharedClsMap;      // if false, at least one of the enabled positions in the map should have corresponding shared clusters to reject
                                          // if true, all the enabled bits should have shared clusters to reject
                                          // The bits left to zero will not be used when evaluating the cut
                                          // e.g.  100000 and TRUE: reject tracks with a shared cls in the first layer
                                          // e.g.  110000 and TRUE: reject tracks with shared cls in both of the first 2 layers
                                          // e.g.  110000 and FALSE: reject tracks with shared cls in either of the first 2 layers
   Bool_t   fRequestCutOnITSsharedClsMap;  // if true, apply the cut above
   
   // TPC quantities
   Bool_t   fRequestTPCrefit;               // if true, request TPC refit      
   
   // TOF quantities
   Bool_t   fRequestTOFout;               // if true, request TOF out

   // TRD selections
   Bool_t   fRequestTRDonlineMatch;    // if true, request the track to be matched to a TRD online track
   
  AliReducedTrackCut(const AliReducedTrackCut &c);
  AliReducedTrackCut& operator= (const AliReducedTrackCut &c);
  
  ClassDef(AliReducedTrackCut,3);
};

#endif
