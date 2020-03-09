/*
 * AliAnalysisTaskFemtoTutorial.h
 *
 *  Created on: 11 Dec 2019
 *      Author: Boris Bajtl
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_POMEGA_PENNE_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_POMEGA_PENNE_H_

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliAODTrack.h"
#include "TROOT.h"
#include "TSystem.h"
#include "AliPIDResponse.h"
#include "AliAODInputHandler.h"

#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"

class AliAnalysisTaskPOmegaPenne : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskPOmegaPenne();
    AliAnalysisTaskPOmegaPenne(const char *name, bool isMC);
    AliAnalysisTaskPOmegaPenne(const AliAnalysisTaskPOmegaPenne&);    // copy ctor
    AliAnalysisTaskPOmegaPenne& operator=(const AliAnalysisTaskPOmegaPenne&); // copy operator
    virtual ~AliAnalysisTaskPOmegaPenne();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *){};
    void SetEventCuts(            AliFemtoDreamEventCuts   *evtCuts )  { fEventCuts            =   evtCuts;  };
    void SetTrackCutsProton(      AliFemtoDreamTrackCuts   *trkCuts )  { fTrackCutsProton      =   trkCuts;  };
    void SetTrackCutsAntiProton(  AliFemtoDreamTrackCuts   *trkCuts )  { fTrackCutsAntiProton  =   trkCuts;  };
    void Setv0Cuts(               AliFemtoDreamv0Cuts      *v0Cuts  )  { fLambdaV0Cuts         =   v0Cuts;   }
    void SetAntiv0Cuts(           AliFemtoDreamv0Cuts      *v0Cuts  )  { fAntiLambdaV0Cuts     =   v0Cuts;   }
    void SetTrackCutsXion(        AliFemtoDreamCascadeCuts *cascCuts)  { fCascadeCutsXi        =   cascCuts; };
    void SetTrackCutsAntiXion(    AliFemtoDreamCascadeCuts *cascCuts)  { fCascadeCutsAntiXi    =   cascCuts; };
    void SetCollectionConfig(     AliFemtoDreamCollConfig  *config  )  { fConfig               =   config;   };
 private:
    void ResetGlobalTrackReference();
    void StoreGlobalTrackReference(AliAODTrack *track);
    bool                                fIsMC;                 //
    AliAODEvent                        *aaEvent;               //      UserExec:Current Event
    AliAODTrack                        *aaTrack;               //      UserExec:Current Track
    AliFemtoDreamEvent                 *fEvent;                //!
    AliFemtoDreamTrack                 *fTrack;                //!
    AliFemtoDreamCascade               *fCascade;              //!
    AliFemtoDreamEventCuts             *fEventCuts;            //
    AliFemtoDreamTrackCuts             *fTrackCutsProton;      //
    AliFemtoDreamTrackCuts             *fTrackCutsAntiProton;  //
    AliFemtoDreamv0                    *fv0;                   //!
    AliFemtoDreamv0Cuts                *fLambdaV0Cuts;         //
    AliFemtoDreamv0Cuts                *fAntiLambdaV0Cuts;     //
    AliFemtoDreamCascadeCuts           *fCascadeCutsXi;        //
    AliFemtoDreamCascadeCuts           *fCascadeCutsAntiXi;    //
    AliFemtoDreamCollConfig            *fConfig;               //
    AliFemtoDreamPairCleaner           *fPairCleaner;          //!
    AliFemtoDreamPartCollection        *fPartColl;             //!
    AliAODTrack                       **fGTI;                  //!
    int                                 fTrackBufferSize;      //
    // Output Container
    TList                              *tlEventCuts;           //!
    TList                              *tlTrackCutsProton;     //!
    TList                              *tlAntiTrackCutsProton; //!
    TList                              *tlLambdaList;          //!
    TList                              *tlAntiLambdaList;      //!
    TList                              *tlCascadeCutsXi;       //!
    TList                              *tlAntiCascadeCutsXi;   //!
    TList                              *tlResults;             //!
    TList                              *tlResultsQA;           //!
  
    ClassDef(AliAnalysisTaskPOmegaPenne,8)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_POMEGA_PENNE_H_ */
