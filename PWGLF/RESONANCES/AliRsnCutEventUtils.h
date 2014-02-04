//
// Class AliRsnCutEventUtils
//
// It works with ESD and AOD events.
//
// authors: F. Bellini (fbellini@cern.ch)

#ifndef ALIRSNCUTEVENTUTILS_H
#define ALIRSNCUTEVENTUTILS_H

#include "AliRsnCut.h"

class AliVVertex;
class AliAnalysisUtils;

class AliRsnCutEventUtils : public AliRsnCut {
 public:

  AliRsnCutEventUtils(const char *name = "cutEventUtils", Bool_t rmFirstEvInChunck = kFALSE, Bool_t checkPileUppA2013 = kTRUE);
  AliRsnCutEventUtils(const AliRsnCutEventUtils &copy);
  AliRsnCutEventUtils &operator=(const AliRsnCutEventUtils &copy);
  virtual ~AliRsnCutEventUtils() {;};
  
  void           SetRemovePileUppA2013(Bool_t doit = kTRUE) {fCheckPileUppA2013 = doit;}
  void           SetRemoveFirstEvtInChunk(Bool_t doit = kTRUE) {fIsRmFirstEvInChunck = doit;}
  void           SetUseMVPlpSelection(Bool_t useMVPlpSelection = kFALSE) { fUseMVPlpSelection = useMVPlpSelection;}
  void           SetUseVertexSelection2013pA(Bool_t zvtxpA2013 = kTRUE)   {fUseVertexSelection2013pA = zvtxpA2013;}
  Bool_t         IsSelected(TObject *object);
  AliAnalysisUtils* GetAnalysisUtils() { return fUtils; }
  void           SetAnalysisUtils(AliAnalysisUtils* utils){ fUtils = utils; }
  void           SetMinPlpContribMV(Int_t minPlpContribMV) { fMinPlpContribMV = minPlpContribMV;}
  void           SetMinPlpContribSPD(Int_t minPlpContribSPD) { fMinPlpContribSPD = minPlpContribSPD;}

 private:
  
  Bool_t              fIsRmFirstEvInChunck; // if kTRUE, remove the first event in the chunk (pA2013)
  Bool_t              fCheckPileUppA2013; // check and reject pileupped events (pA2013)
  Bool_t              fUseMVPlpSelection; // check for pile-up from multiple vtx 
  Int_t               fMinPlpContribMV; // min. n. of MV pile-up contributors
  Int_t               fMinPlpContribSPD; // min. n. of pile-up contributors from SPD
  Bool_t              fUseVertexSelection2013pA;// check and reject vertex of events for pA2013
 
  AliAnalysisUtils  * fUtils; //pointer to the AliAnalysisUtils object

  ClassDef(AliRsnCutEventUtils, 2)
    
    };

#endif
