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
    void SetEventCuts(            AliFemtoDreamEventCuts   *evtCuts         )  { fEventCuts             =   evtCuts;        };
    // void SetTrackCutsProton(      AliFemtoDreamTrackCuts   *trkCuts )  { fTrackCutsProton      =   trkCuts;  };
    // void SetTrackCutsAntiProton(  AliFemtoDreamTrackCuts   *trkCuts )  { fTrackCutsAntiProton  =   trkCuts;  };
    void Setv0Cuts(               AliFemtoDreamv0Cuts      *v0Cuts          )  { fLambdaV0Cuts          =   v0Cuts;         };
    void SetAntiv0Cuts(           AliFemtoDreamv0Cuts      *antiV0Cuts      )  { fAntiLambdaV0Cuts      =   antiV0Cuts;     };
    void SetTrackCutsXion(        AliFemtoDreamCascadeCuts *cascCuts        )  { fCascadeCutsXi         =   cascCuts;       };
    void SetTrackCutsAntiXion(    AliFemtoDreamCascadeCuts *antiCascCuts    )  { fCascadeCutsAntiXi     =   antiCascCuts;   };
    void SetCollectionConfig(     AliFemtoDreamCollConfig  *config          )  { fConfig                =   config;         };
    // #2
    void SetEventCuts2(            AliFemtoDreamEventCuts   *evtCuts2       )  { fEventCuts2            =   evtCuts2;       };
    void Setv0Cuts2(               AliFemtoDreamv0Cuts      *v0Cuts2        )  { fLambdaV0Cuts2         =   v0Cuts2;        };
    void SetAntiv0Cuts2(           AliFemtoDreamv0Cuts      *antiV0Cuts2    )  { fAntiLambdaV0Cuts2     =   antiV0Cuts2;    };
    void SetTrackCutsXion2(        AliFemtoDreamCascadeCuts *cascCuts2      )  { fCascadeCutsXi2        =   cascCuts2;      };
    void SetTrackCutsAntiXion2(    AliFemtoDreamCascadeCuts *antiCascCuts2  )  { fCascadeCutsAntiXi2    =   antiCascCuts2;  };
    // void SetCollectionConfig2(     AliFemtoDreamCollConfig  *config         )  { fConfig2               =   config;         };
    float CalculateInvMassHere(AliFemtoDreamv0 *v0, int PDGPosDaug, int PDGNegDaug);        // copied from AliFemtoDreamv0Cuts
    float CalculateInvMassLambda(TVector3 momPosDaughter, TVector3 momNegDaughter);
    float CalculateInvMassXi(TVector3 momBach, TVector3 momPosDaughter, TVector3 momNegDaughter);
    // void MixChildParticles(std::vector<AliFemtoDreamBasePart> XiVector, 
    //                        std::vector<AliFemtoDreamBasePart> LambdaVector, 
    //                        TList *outputLists,
    //                        bool checkSameParticleMixing = true);

    // void Setv0Cuts_rec(            AliFemtoDreamv0Cuts      *v0Cuts_rec        )  { fLambdaV0Cuts_rec        =   v0Cuts_rec;        };
    // void SetAntiv0Cuts_rec(        AliFemtoDreamv0Cuts      *v0AntiCuts_rec        )  { fAntiLambdaV0Cuts_rec        =   v0AntiCuts_rec;        };

 private:
    void ResetGlobalTrackReference();
    void StoreGlobalTrackReference(AliVTrack *track);
    bool                                fIsMC;                 //
    AliVEvent                          *VEvent;                //      UserExec:Current Event
    AliVTrack                          *VTrack;                //      UserExec:Current Track
    AliFemtoDreamEvent                 *fEvent;                //!
    AliFemtoDreamTrack                 *fTrack;                //!
    AliFemtoDreamEventCuts             *fEventCuts;            //
    AliFemtoDreamEventCuts             *fEventCuts2;           //
    // AliFemtoDreamTrackCuts             *fTrackCutsProton;      //
    // AliFemtoDreamTrackCuts             *fTrackCutsAntiProton;  //
    AliFemtoDreamv0                    *fv0;                   //!
    AliFemtoDreamCascade               *fCascade;              //!
    AliFemtoDreamv0Cuts                *fLambdaV0Cuts;         //
    AliFemtoDreamv0Cuts                *fAntiLambdaV0Cuts;     //
    AliFemtoDreamCascadeCuts           *fCascadeCutsXi;        //
    AliFemtoDreamCascadeCuts           *fCascadeCutsAntiXi;    //
    AliFemtoDreamv0                    *fv0_2;                 //!
    AliFemtoDreamCascade               *fCascade2;             //!
    AliFemtoDreamv0Cuts                *fLambdaV0Cuts2;        //
    AliFemtoDreamv0Cuts                *fAntiLambdaV0Cuts2;    //
    AliFemtoDreamCascadeCuts           *fCascadeCutsXi2;       //
    AliFemtoDreamCascadeCuts           *fCascadeCutsAntiXi2;   //
    AliFemtoDreamCollConfig            *fConfig;               //       Keep Lambda config
    // AliFemtoDreamCollConfig            *fConfig2;              //       Keep Xi config
    AliFemtoDreamPairCleaner           *fPairCleaner;          //!
    AliFemtoDreamPairCleaner           *fPairCleaner2;         //!
    AliFemtoDreamPartCollection        *fPartColl;             //!
    AliFemtoDreamPartCollection        *fPartColl2;            //!
    // AliFemtoDreamv0Cuts                *fLambdaV0Cuts_rec;     //!
    // AliFemtoDreamv0Cuts                *fAntiLambdaV0Cuts_rec; //!
    AliVTrack                          **fGTI;                 //!
    int                                 fTrackBufferSize;      //
    // ## Output Container
    TList                              *tlEventCuts;           //!
    // TList                              *tlTrackCutsProton;     //!
    // TList                              *tlAntiTrackCutsProton; //!
    TList                              *tlLambdaList;          //!
    TList                              *tlAntiLambdaList;      //!
    TList                              *tlCascadeCutsXi;       //!
    TList                              *tlAntiCascadeCutsXi;   //!
    // #2
    TList                              *tlEventCuts2;                       //!
    TList                              *tlLambdaList2;                      //!
    TList                              *tlAntiLambdaList2;                  //!
    TList                              *tlCascadeCutsXi2;                   //!
    TList                              *tlAntiCascadeCutsXi2;               //!
    TList                              *tlResults;                          //!      
    TList                              *tlResults2;                         //!
    TList                              *tlResultsQA;                        //!      PairCleaner - Keep Lambda
    TList                              *tlResultsQA2;                       //!      PairCleaner - Keep Xi
    // ## MC Container
    // TList                              *tlProtonMC;            //!
    // TList                              *tlAntiProtonMC;        //!
    TList                              *tlLambdaMC;                         //!
    TList                              *tlAntiLambdaMC;                     //!
    TList                              *tlRecombination;                    //!      Recombinations Lists and histos
    TH1F                               *hInvMassLambda_total;               //!
    TH1F                               *hInvMassLambda_shared_pion;         //!
    TH1F                               *hInvMassLambda_shared_proton;       //!
    TH1F                               *hInvMassXi_total;                   //!
    TH1F                               *hInvMassXi_shared_bach;             //!
    TH1F                               *hInvMassXi_shared_pi_daugh;         //!
    TH1F                               *hInvMassXi_shared_prot_daugh;       //!
    TH1F                               *hInvMassXi_shared_Lambda;           //!
    TH1F                               *hInvMassLambda_sanityCheck;         //!
    TH1                                *hInvMassXi_sanityCheck;             //!
    TH1F                               *fEvtCounter;                        //!
    AliFemtoDreamv0                    *fv0_recomb;                         //!
    TList                              *tlLambdaList_rec;                   //!
    TList                              *tlAntiLambdaList_rec;               //!

  
    ClassDef(AliAnalysisTaskPOmegaPenne,21)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_POMEGA_PENNE_H_ */
