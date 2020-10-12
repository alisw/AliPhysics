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
    // steering params from macro
    void SetMixBeforePC         ( bool                      choiceBefore    )  { fmixBeforePC           =   choiceBefore;   };
    void SetMixAfterPC          ( bool                      choiceAfter     )  { fmixAfterPC            =   choiceAfter;    };
    void SetFullBlastQA         ( bool                      choiceQA        )  { ffullBlastQA           =   choiceQA;       };
    void SetInvMassPairClean    ( bool                      choicePC        )  { fisInvMassPairClean    =   choicePC;      };
    void SetMultTrigger         ( TString                   choiceMult      )  { fmultTrigger           =   choiceMult;    };
    // Cuts #1
    void SetEventCuts(            AliFemtoDreamEventCuts   *evtCuts         )  { fEventCuts             =   evtCuts;        };
    void Setv0Cuts              ( AliFemtoDreamv0Cuts      *v0Cuts          )  { fLambdaV0Cuts          =   v0Cuts;         };
    void SetAntiv0Cuts          ( AliFemtoDreamv0Cuts      *antiV0Cuts      )  { fAntiLambdaV0Cuts      =   antiV0Cuts;     };
    void SetTrackCutsXion       ( AliFemtoDreamCascadeCuts *cascCuts        )  { fCascadeCutsXi         =   cascCuts;       };
    void SetTrackCutsAntiXion   ( AliFemtoDreamCascadeCuts *antiCascCuts    )  { fCascadeCutsAntiXi     =   antiCascCuts;   };
    // Cuts #2
    void SetEventCuts2          ( AliFemtoDreamEventCuts   *evtCuts2       )  { fEventCuts2            =   evtCuts2;       };
    void Setv0Cuts2             ( AliFemtoDreamv0Cuts      *v0Cuts2        )  { fLambdaV0Cuts2         =   v0Cuts2;        };
    void SetAntiv0Cuts2         ( AliFemtoDreamv0Cuts      *antiV0Cuts2    )  { fAntiLambdaV0Cuts2     =   antiV0Cuts2;    };
    void SetTrackCutsXion2      ( AliFemtoDreamCascadeCuts *cascCuts2      )  { fCascadeCutsXi2        =   cascCuts2;      };
    void SetTrackCutsAntiXion2  ( AliFemtoDreamCascadeCuts *antiCascCuts2  )  { fCascadeCutsAntiXi2    =   antiCascCuts2;  };
    // config 
    void SetCollectionConfig    ( AliFemtoDreamCollConfig  *config          )  { fConfig                =   config;         };

    // my analysis functions
    float CalculateInvMassLambda(TVector3 momNegDaughter, int PDGnegDaughter, TVector3 momPosDaughter, int PDGposDaughter);
    
    float CalculateInvMassLambda(AliFemtoDreamBasePart *lambdaParticle, bool isAntiParticle);
    
    float CalculateInvMassXi(TVector3 momBach, int PGGbach, TVector3 momPosDaughter, int PDGposDaughter, TVector3 momNegDaughter, int PDGnegDaughter);
    
    float CalculateInvMassXi(AliFemtoDreamBasePart *xiParticle, bool isAntiParticle);
    
    void CleanDecay(std::vector<AliFemtoDreamBasePart> *Decay, string particleSteering);

    void CleanDecayAndDecay(std::vector<AliFemtoDreamBasePart> *Decay1,
                            std::vector<AliFemtoDreamBasePart> *Decay2,
                            bool isAntiParticle);
    
    float WeightLambda(float pT);
    
    float WeightAntiLambda(float pT);
    
    float WeightXi(float pT);
    
    float WeightAntiXi(float pT);

float RelativePairMomentum(AliFemtoDreamBasePart *part1, const int pdg1, AliFemtoDreamBasePart *part2, const int pdg2);
    
 private:
    void ResetGlobalTrackReference();
    void StoreGlobalTrackReference(AliVTrack *track);
    // control vars
    bool                                fmixBeforePC;          //
    bool                                fmixAfterPC;           //
    bool                                fisInvMassPairClean;   //
    TString                             fmultTrigger;          //
    bool                                ffullBlastQA;          //
    AliVEvent                          *VEvent;                //      UserExec:Current Event
    AliVTrack                          *VTrack;                //      UserExec:Current Track
    AliFemtoDreamEvent                 *fEvent;                //!
    AliFemtoDreamTrack                 *fTrack;                //!
    AliFemtoDreamCollConfig            *fConfig;               //       Keep Lambda config
    AliVTrack                          **fGTI;                 //!
    int                                 fTrackBufferSize;      //
    //  #1
    AliFemtoDreamEventCuts             *fEventCuts;            //
    AliFemtoDreamv0                    *fv0;                   //!
    AliFemtoDreamCascade               *fCascade;              //!
    AliFemtoDreamv0Cuts                *fLambdaV0Cuts;         //
    AliFemtoDreamv0Cuts                *fAntiLambdaV0Cuts;     //
    AliFemtoDreamCascadeCuts           *fCascadeCutsXi;        //
    AliFemtoDreamCascadeCuts           *fCascadeCutsAntiXi;    //
    AliFemtoDreamPairCleaner           *fPairCleaner;          //!
    AliFemtoDreamPartCollection        *fPartColl;             //!
    // particle vectors
    std::vector<AliFemtoDreamBasePart>  vLambda;               //!
    std::vector<AliFemtoDreamBasePart>  vAntiLambda;           //!
    std::vector<AliFemtoDreamBasePart>  vXi;                   //!
    std::vector<AliFemtoDreamBasePart>  vAntiXi;               //!
    // #2
    AliFemtoDreamEventCuts             *fEventCuts2;           //
    AliFemtoDreamv0                    *fv0_2;                 //!
    AliFemtoDreamCascade               *fCascade2;             //!
    AliFemtoDreamv0Cuts                *fLambdaV0Cuts2;        //
    AliFemtoDreamv0Cuts                *fAntiLambdaV0Cuts2;    //
    AliFemtoDreamCascadeCuts           *fCascadeCutsXi2;       //
    AliFemtoDreamCascadeCuts           *fCascadeCutsAntiXi2;   //
    AliFemtoDreamPairCleaner           *fPairCleaner2;         //!
    AliFemtoDreamPartCollection        *fPartColl2;            //!
    // ## Output Container
    TList                              *tlEventCuts;           //!
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
    TList                              *tlLambdaMC;                         //!
    TList                              *tlAntiLambdaMC;                     //!
    TList                              *tlXiMC;                             //!
    TList                              *tlAntiXiMC;                         //!

    ////////////////////////////////////////////////////////////////////////////////
    // My recombination stuff ///////////////
    ///////////////////////////////////////////////////////////////////////////////
    std::vector<AliFemtoDreamBasePart>  vLambda_recomb;                     //!
    std::vector<AliFemtoDreamBasePart>  tmpLambda_recomb;                   //!
    std::vector<AliFemtoDreamBasePart>  tmpXi_recomb;                       //!
    std::vector<AliFemtoDreamBasePart>  tmpAntiLambda_recomb;               //!
    std::vector<AliFemtoDreamBasePart>  tmpAntiXi_recomb;                   //!
    TList                              *tlRecombination_before;             //!      Recombinations Lists and histos
    TList                              *tlRecombination_after;              //!      Recombinations Lists and histos
    //////////////////////
    // before Histos /////
    //////////////////////
    // particles ///
    TH1F                               *hInvMassLambda_sanityCheck_before;                  //!
    TH1F                               *hInvMassLambda_total_before;                        //!
    TH1F                               *hInvMassLambda_shared_pion_before;                  //!
    TH1F                               *hInvMassLambda_shared_proton_before;                //!
    TH1F                               *hInvMassLambda_shared_lambda_before;                //!
    TH1F                               *hInvMassXi_sanityCheck_before;                      //!
    TH1F                               *hInvMassXi_total_before;                            //!
    TH1F                               *hInvMassXi_shared_bach_before;                      //!
    TH1F                               *hInvMassXi_shared_pi_daugh_before;                  //!
    TH1F                               *hInvMassXi_shared_prot_daugh_before;                //!
    TH1F                               *hInvMassXi_shared_Lambda_before;                    //!
    TH1F                               *hInvMassXi_shared_pion_bach_prot_daugh_before;      //!
    TH1F                               *hInvMassXi_nothing_shared;                          //!
    // anti particles ///
    TH1F                               *hInvMassAntiLambda_sanityCheck_before;              //!
    TH1F                               *hInvMassAntiLambda_total_before;                    //!
    TH1F                               *hInvMassAntiLambda_shared_pion_before;              //!
    TH1F                               *hInvMassAntiLambda_shared_proton_before;            //!
    TH1F                               *hInvMassAntiLambda_shared_lambda_before;            //!
    TH1F                               *hInvMassAntiXi_sanityCheck_before;                  //!
    TH1F                               *hInvMassAntiXi_total_before;                        //!
    TH1F                               *hInvMassAntiXi_shared_bach_before;                  //!
    TH1F                               *hInvMassAntiXi_shared_pi_daugh_before;              //!
    TH1F                               *hInvMassAntiXi_shared_prot_daugh_before;            //!
    TH1F                               *hInvMassAntiXi_shared_Lambda_before;                //!
    TH1F                               *hInvMassAntiXi_shared_pion_bach_prot_daugh_before;  //!
    TH1F                               *hInvMassAntiXi_nothing_shared;                      //!
    
    TH1F                               *fEvtCounterBefore;                                  //!
    //////////////////////
    // after Histos /////
    //////////////////////
    TList                              *tlLambdaRecombination_after;                        //!
    TList                              *tlAntiLambdaRecombination_after;                    //!
    TList                              *tlXiRecombination_after;                            //!
    TList                              *tlAntiXiRecombination_after;                        //!
    // particles ///
    TH1F                               *hInvMassLambda_sanityCheck_after;                           //!
    TH1F                               *hInvMassLambda_pi_bach_Xi_after;                            //!
    TH1F                               *hInvMassLambda_pi_daugh_Xi_after;                           //!
    TH1F                               *hInvMassLambda_prot_Xi_after;                               //!
    TH1F                               *hInvMassLambda_full_lambda_from_Xi_after;                   //!
    TH1F                               *hInvMassXi_sanityCheck_after;                               //!
    TH1F                               *hInvMassXi_Lamda_pi_daugh_after;                            //!
    TH1F                               *hInvMassXi_Lamda_prot_daugh_after;                          //!
    TH1F                               *hInvMassXi_Lamda_pi_bach_after;                             //!
    TH1F                               *hInvMassXi_Lamda_full_after;                                //!
    TH1F                               *hInvMassXi_Lamda_pi_no_correctLambdaMass;                   //!
    TH1F                               *hInvMassXi_Lamda_prot_no_correctLambdaMass;                 //!
    // anti particles ///
    TH1F                               *hInvMassAntiLambda_sanityCheck_after;                       //!
    TH1F                               *hInvMassAntiLambda_pi_bach_Xi_after;                        //!
    TH1F                               *hInvMassAntiLambda_pi_daugh_Xi_after;                       //!
    TH1F                               *hInvMassAntiLambda_prot_Xi_after;                           //!
    TH1F                               *hInvMassAntiLambda_full_lambda_from_Xi_after;               //!
    TH1F                               *hInvMassAntiXi_sanityCheck_after;                           //!
    TH1F                               *hInvMassAntiXi_AntiLamda_antipi_daugh_after;                //!
    TH1F                               *hInvMassAntiXi_AntiLamda_antiprot_daugh_after;              //!
    TH1F                               *hInvMassAntiXi_AntiLamda_antipi_bach_after;                 //!
    TH1F                               *hInvMassAntiXi_AntiLamda_full_after;                        //!
    TH1F                               *hInvMassAntiXi_AntiLamda_antipi_no_correctAntiLambdaMass;   //!
    TH1F                               *hInvMassAntiXi_AntiLamda_antiprot_no_correctAntiLambdaMass; //!

    TH1F                               *fEvtCounterAfter;                                           //!
    /////////////////////
    // Inv Mass PC   /////
    //////////////////////
    // folder inv mass pairclean
    TList                              *tlInvMassPairClean;                                 //!<!
        //- folder decay
    TList                              *tlCleanDecay;                                       //!<!
            //-> Decay Diff To PDG Mass
    TH1F                               *hLambdaCleanedPartMassDiffToPDG_Decay;                    //!<!
    TH1F                               *hAntiLambdaCleanedPartMassDiffToPDG_Decay;                //!<!
    TH1F                               *hXiCleanedPartMassDiffToPDG_Decay;                        //!<!
    TH1F                               *hAntiXiCleanedPartMassDiffToPDG_Decay;                    //!<!
            //-> Decay Mass
    TH1F                               *hLambdaCleanedPartMass_Decay;                    //!<!
    TH1F                               *hAntiLambdaCleanedPartMass_Decay;                //!<!
    TH1F                               *hXiCleanedPartMass_Decay;                        //!<!
    TH1F                               *hAntiXiCleanedPartMass_Decay;                    //!<!
        //- folder decay and decay
    TList                              *tlCleanDecayAndDecay;                               //!<!
            //-> DecayAndDecay Diff To PDG Mass
    TH1F                               *hLambdaCleanedPartMassDiffToPDG_DecayDecay;                    //!<!
    TH1F                               *hAntiLambdaCleanedPartMassDiffToPDG_DecayDecay;                //!<!
    TH1F                               *hXiCleanedPartMassDiffToPDG_DecayDecay;                        //!<!
    TH1F                               *hAntiXiCleanedPartMassDiffToPDG_DecayDecay;                    //!<!
            //-> DecayAndDecay Mass                             
    TH1F                               *hLambdaCleanedPartMass_DecayDecay;                    //!<!
    TH1F                               *hAntiLambdaCleanedPartMass_DecayDecay;                //!<!
    TH1F                               *hXiCleanedPartMass_DecayDecay;                        //!<!
    TH1F                               *hAntiXiCleanedPartMass_DecayDecay;                    //!<!
    // CPA stuffs
    TH1F                              **hCPA_stuff;                                             //!<!
        //- folder CPA MC after Pairlcean
    TList                              *tlCPA_MC_afterPairClean;                                //!<!
    TList                              *tlLambda;                                               //!<!
    TList                              *tlAntiLambda;                                               //!<!
    TList                              *tlXi;                                               //!<!
    TList                              *tlAntiXi;                                               //!<!
    TH2F                               *CPAPtBinningPrim_lambda;                                 //!<!
    TH2F                               *CPAPtBinningMat_lambda;                                  //!<!
    TH2F                               *CPAPtBinningSec_lambda;                                  //!<!
    TH2F                               *CPAPtBinningCont_lambda;                                 //!<!
    TH2F                               *CPAPtBinningPrim_lambda_dump;                            //!<!
    TH2F                               *CPAPtBinningMat_lambda_dump;                             //!<!
    TH2F                               *CPAPtBinningSec_lambda_dump;                             //!<!
    TH2F                               *CPAPtBinningCont_lambda_dump;                            //!<!
    TH2F                               *CPAPtBinningPrim_antilambda;                             //!<!    
    TH2F                               *CPAPtBinningMat_antilambda;                              //!<!
    TH2F                               *CPAPtBinningSec_antilambda;                              //!<!
    TH2F                               *CPAPtBinningCont_antilambda;                             //!<!
    TH2F                               *CPAPtBinningPrim_xi;                                     //!<!
    TH2F                               *CPAPtBinningMat_xi;                                      //!<!
    TH2F                               *CPAPtBinningSec_xi;                                      //!<!
    TH2F                               *CPAPtBinningCont_xi;                                     //!<!    
    TH2F                               *CPAPtBinningPrim_xi_dump;                                //!<!
    TH2F                               *CPAPtBinningMat_xi_dump;                                 //!<!
    TH2F                               *CPAPtBinningSec_xi_dump;                                 //!<!
    TH2F                               *CPAPtBinningCont_xi_dump;                                //!<!
    TH2F                               *CPAPtBinningPrim_antixi;                                 //!<!
    TH2F                               *CPAPtBinningMat_antixi;                                  //!<!
    TH2F                               *CPAPtBinningSec_antixi;                                  //!<!
    TH2F                               *CPAPtBinningCont_antixi;                                 //!<!
    
    //////////////////////
    // weird stuff   /////
    //////////////////////
    TH1F                               *kStarXiLambda_unchanged;                                //!<!
    TH1F                               *kStarXiLambda_changed;                                  //!<!
    TH1F                               *kStarAntiXiAntiLambda_unchanged;                                //!<!
    TH1F                               *kStarAntiXiAntiLambda_changed;                                  //!<!
    
    ClassDef(AliAnalysisTaskPOmegaPenne,33)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_POMEGA_PENNE_H_ */
