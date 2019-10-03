#ifndef ALIANALYSISTASKLIGHTN_H
#define ALIANALYSISTASKLIGHTN_H

/*
 * AliAnalysisTaskFemtoPlotration.h
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */


#include "AliAnalysisTaskSE.h"
#include "AliLightNAnalysis.h"
#include "AliLightNEventCuts.h"
#include "AliLightNTrackCuts.h"
#include "TList.h"
class AliAnalysisTaskLightN : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskLightN();
    AliAnalysisTaskLightN(const char *name,bool isMC);
    virtual ~AliAnalysisTaskLightN();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *){};
    void SetMVPileUp(bool mvPileUp){fMVPileUp=mvPileUp;};
    void SetEvtCutQA(bool setQA){fEvtCutQA=setQA;};
    void SetEventCutsParticle(AliLightNEventCuts *cuts){fEvtCutsParticle=cuts;};
    void SetTrackCutsProton(AliLightNTrackCuts *cuts){fTrackCutsProton=cuts;};
    void SetAntiTrackCutsProton(AliLightNTrackCuts *cuts){fAntiTrackCutsProton=cuts;};
    void SetTrackCutsDeuteron(AliLightNTrackCuts *cuts){fTrackCutsDeuteron=cuts;};
    void SetAntiTrackCutsDeuteron(AliLightNTrackCuts *cuts){fAntiTrackCutsDeuteron=cuts;};
    void SetTrackBufferSize(int size){fTrackBufferSize=size;};
    
private:
    int fTrackBufferSize;				//
    bool fMVPileUp;					//
    bool fEvtCutQA;					//
    bool fIsMC;					//
    TString fname;              //
    AliLightNAnalysis *fAnalysisParticle;		//!
    TList *fQA;					//!
    AliLightNEventCuts *fEvtCutsParticle;			//
    TList *fEvtHistListParticle;				//!
    AliLightNTrackCuts *fTrackCutsProton;			//
    AliLightNTrackCuts *fAntiTrackCutsProton;		//
    AliLightNTrackCuts *fTrackCutsDeuteron;			//
    AliLightNTrackCuts *fAntiTrackCutsDeuteron;		//
    TList *fTrackCutHistListProton;			//!
    TList *fTrackCutHistMCListProton;		//!
    TList *fAntiTrackCutHistListProton;		//!
    TList *fAntiTrackCutHistMCListProton;		//!
    TList *fTrackCutHistListDeuteron;		//!
    TList *fTrackCutHistMCListDeuteron;		//!
    TList *fAntiTrackCutHistListDeuteron;		//!
    TList *fAntiTrackCutHistMCListDeuteron;		//!
    ClassDef(AliAnalysisTaskLightN,1)
};

#endif /* ALIANALYSISTASKLIGHTN_H */

