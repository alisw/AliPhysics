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
  void SetRequestTOFout(Bool_t flag = kTRUE) {fRequestTOFout = flag;}  
  void SetRequestTRDmatch(Bool_t flag = kTRUE) {fRequestTRDonlineMatch = flag;}
  
  Bool_t GetRejectKinks() const {return fRejectKinks;}
  Bool_t GetRejectTaggedGamma() const {return fRejectTaggedGamma;}
  Bool_t GetRejectTaggedPureGamma() const {return fRejectTaggedPureGamma;}
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
   
   // ITS quantities
  Bool_t    fRequestITSrefit;                // if true, request kITSrefit flag to be on
  UChar_t fCutOnITShitMap;               // hit map encoding various requests on cluster configurations
  Bool_t    fUseANDonITShitMap;        // if false, at least one of the enabled positions in the cut map should have corresponding clusters
                                                          // if true, all the enabled bits should have matching clusters
                                                          // The bits left to zero will not be used when evaluating the cut
                                                          // e.g.  100000 and FALSE: equivalent to SPDfirst
                                                          // e.g.  110000 and FALSE: equivalent to SPDany (at least one hit in the first two layers)
                                                          // e.g.  110000 and TRUE: equivalent to requesting both SPD layers to have a hit
   Bool_t   fRequestCutOnITShitMap;  // if true, apply the cut above
   
   // TPC quantities
   Bool_t   fRequestTPCrefit;               // if true, request TPC refit      
   
   // TOF quantities
   Bool_t   fRequestTOFout;               // if true, request TOF out

   // TRD selections
   Bool_t   fRequestTRDonlineMatch;    // if true, request the track to be matched to a TRD online track
   
  AliReducedTrackCut(const AliReducedTrackCut &c);
  AliReducedTrackCut& operator= (const AliReducedTrackCut &c);
  
  ClassDef(AliReducedTrackCut,2);
};

#endif
