#ifndef ALILIGHTNANALYSIS_H
#define ALILIGHTNANALYSIS_H

/*
 * AliLightNAnalysis.h
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliLightNEvent.h"
#include "AliLightNEventCuts.h"
#include "AliLightNTrack.h"
#include "AliLightNTrackCuts.h"
#include "Rtypes.h"
#include <vector>
#include "TClonesArray.h"
#include <iostream>
class AliLightNAnalysis {
 public:
  AliLightNAnalysis();
  void SetMVPileUp(bool mvPileUp){fMVPileUp=mvPileUp;};
  void SetEvtCutQA(bool setQA){fEvtCutQA=setQA;};
  void SetEventCutsProton(AliLightNEventCuts *cuts){fEvtCuts=cuts;};
  void SetEventCutsDeuteron(AliLightNEventCuts *cuts){fEvtCuts=cuts;};
  TList *GetEventCutHists(){return fEvtCuts->GetHistList();};
  void SetTrackCutsProton(AliLightNTrackCuts *cuts){fTrackCuts=cuts;};
  void SetTrackCutsDeuteron(AliLightNTrackCuts *cuts){fTrackCuts=cuts;};
  TList *GetTrackCutHists(){return fTrackCuts->GetQAHists();};
  TList *GetTrackCutHistsMC(){return fTrackCuts->GetMCQAHists();};
  void SetAntiTrackCutsProton(AliLightNTrackCuts *cuts){fAntiTrackCuts=cuts;};
  void SetAntiTrackCutsDeuteron(AliLightNTrackCuts *cuts){fAntiTrackCuts=cuts;};
  TList *GetAntitrackCutHists(){return fAntiTrackCuts->GetQAHists();};
  TList *GetAntitrackCutHistsMC(){return fAntiTrackCuts->GetMCQAHists();};
  TList *GetQAList() {return fQA;};
  void SetTrackBufferSize(int size){fTrackBufferSize=size;};
  void Init();
  TString ClassName() {return "AliLightNAnalysis";};
  void Make(AliAODEvent *evt);
  
  virtual ~AliLightNAnalysis();
 private:
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fMVPileUp;                           //!
  bool fEvtCutQA;                           //!
  TList *fQA;                               //!
  AliLightNTrack *fLightNTrack;          //!
  AliLightNEvent *fEvent;               //!
  AliLightNEventCuts *fEvtCuts;         //!
  AliLightNTrackCuts *fTrackCuts;       //!
  AliLightNTrackCuts *fAntiTrackCuts;   //!
  int fTrackBufferSize;
  AliAODTrack **fGTI;			//!
  ClassDef(AliLightNAnalysis,1)
};

#endif /* ALILIGHTNANALYSIS_H */
