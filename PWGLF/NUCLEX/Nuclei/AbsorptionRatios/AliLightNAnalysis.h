#ifndef ALILIGHTNANALYSIS_H
#define ALILIGHTNANALYSIS_H

/*
 * AliLightNAnalysis.h
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */

#include "AliLightNEvent.h"
#include "AliLightNEventCuts.h"
#include "AliLightNTrack.h"
#include "AliLightNTrackCuts.h"
#include "Rtypes.h"
#include <vector>
#include "TClonesArray.h"
#include <iostream>
#include "AliMultSelection.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
class AliLightNAnalysis {
public:
    AliLightNAnalysis();
    void SetMVPileUp(bool mvPileUp){fMVPileUp=mvPileUp;};
    void SetEvtCutQA(bool setQA){fEvtCutQA=setQA;};
    void SetEventCutsParticle(AliLightNEventCuts *cuts){fEvtCuts=cuts;};
    TList *GetEventCutHists(){return fEvtCuts->GetHistList();};
    void SetTrackCutsProton(AliLightNTrackCuts *cuts){fTrackCutsProton=cuts;};
    void SetTrackCutsDeuteron(AliLightNTrackCuts *cuts){fTrackCutsDeuteron=cuts;};
    TList *GetTrackCutHistsProton(){return fTrackCutsProton->GetQAHists();};
    TList *GetTrackCutHistsDeuteron(){return fTrackCutsDeuteron->GetQAHists();};
    void SetAntiTrackCutsProton(AliLightNTrackCuts *cuts){fAntiTrackCutsProton=cuts;};
    void SetAntiTrackCutsDeuteron(AliLightNTrackCuts *cuts){fAntiTrackCutsDeuteron=cuts;};
    TList *GetAntitrackCutHistsProton(){return fAntiTrackCutsProton->GetQAHists();};
    TList *GetAntitrackCutHistsDeuteron(){return fAntiTrackCutsDeuteron->GetQAHists();};
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
    AliLightNTrackCuts *fTrackCutsProton;       //!
    AliLightNTrackCuts *fTrackCutsDeuteron;       //!
    AliLightNTrackCuts *fAntiTrackCutsProton;   //!
    AliLightNTrackCuts *fAntiTrackCutsDeuteron;   //!
    int fTrackBufferSize;
    AliAODTrack **fGTI;			//!
    ClassDef(AliLightNAnalysis,1)
};

#endif /* ALILIGHTNANALYSIS_H */
