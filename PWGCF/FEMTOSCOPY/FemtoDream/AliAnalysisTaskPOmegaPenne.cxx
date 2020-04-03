// /*
//  * AliAnalysisTaskPOmegaPenne.cxx
//  *
//  *  Created on: 11 Dec 2019
//  *      Author: Boris Bajtl
//  */

#include "AliAnalysisTaskPOmegaPenne.h"
#include <string.h>

ClassImp(AliAnalysisTaskPOmegaPenne)

    AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne() :  AliAnalysisTaskSE(),
                                                                fIsMC(false),
                                                                aaEvent(0),
                                                                aaTrack(0),
                                                                fEvent(0),
                                                                fTrack(0),
                                                                fCascade(0),
                                                                fEventCuts(0),
                                                                // fTrackCutsProton(0),
                                                                // fTrackCutsAntiProton(0),
                                                                fv0(0),
                                                                fLambdaV0Cuts(0),
                                                                fAntiLambdaV0Cuts(0),
                                                                fCascadeCutsXi(0),
                                                                fCascadeCutsAntiXi(0),
                                                                fConfig(0),
                                                                fPairCleaner(0),
                                                                fPartColl(0),
                                                                fGTI(0),
                                                                fTrackBufferSize(10000),
                                                                tlEventCuts(0),
                                                                // tlTrackCutsProton(0),
                                                                // tlAntiTrackCutsProton(0),
                                                                tlLambdaList(0),
                                                                tlAntiLambdaList(0),
                                                                tlCascadeCutsXi(0),
                                                                tlAntiCascadeCutsXi(0),
                                                                tlResults(0),
                                                                tlResultsQA(0),
                                                                // tlProtonMC(0),
                                                                // tlAntiProtonMC(0),
                                                                tlLambdaMC(0),
                                                                tlAntiLambdaMC(0)
{
}
AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne(const char *name, bool isMC) : AliAnalysisTaskSE(name),
                                                                                      fIsMC(isMC),
                                                                                      aaEvent(0),
                                                                                      aaTrack(0),
                                                                                      fEvent(0),
                                                                                      fTrack(0),
                                                                                      fCascade(0),
                                                                                      fEventCuts(0),
                                                                                    //   fTrackCutsProton(0),
                                                                                    //   fTrackCutsAntiProton(0),
                                                                                      fv0(0),
                                                                                      fLambdaV0Cuts(0),
                                                                                      fAntiLambdaV0Cuts(0),
                                                                                      fCascadeCutsXi(0),
                                                                                      fCascadeCutsAntiXi(0),
                                                                                      fConfig(0),
                                                                                      fPairCleaner(0),
                                                                                      fPartColl(0),
                                                                                      fGTI(0),
                                                                                      fTrackBufferSize(10000),
                                                                                      tlEventCuts(0),
                                                                                    //   tlTrackCutsProton(0),
                                                                                    //   tlAntiTrackCutsProton(0),
                                                                                      tlLambdaList(0),
                                                                                      tlAntiLambdaList(0),
                                                                                      tlCascadeCutsXi(0),
                                                                                      tlAntiCascadeCutsXi(0),
                                                                                      tlResults(0),
                                                                                      tlResultsQA(0),
                                                                                    //   tlProtonMC(0),
                                                                                    //   tlAntiProtonMC(0),
                                                                                      tlLambdaMC(0),
                                                                                      tlAntiLambdaMC(0)
{
    DefineOutput(1, TList::Class());    // Event Cuts
    // DefineOutput(2, TList::Class());    // Proton Track Cuts
    // DefineOutput(3, TList::Class());    // Anti Proton Track Cuts
    DefineOutput(2, TList::Class());    // Lambda Track Cuts
    DefineOutput(3, TList::Class());    // Anti Lambda Track Cuts
    DefineOutput(4, TList::Class());    // Xi Track Cuts
    DefineOutput(5, TList::Class());    // Anti Xi Track Cuts
    DefineOutput(6, TList::Class());    // Results
    DefineOutput(7, TList::Class());    // QA Results
    if (isMC)
    {
        // DefineOutput(10, TList::Class());    // MC Track Proton
        // DefineOutput(11, TList::Class());    // MC Track AntiProton
        DefineOutput(8, TList::Class());    // MC V0 - Lamba
        DefineOutput(9, TList::Class());    // MC AntiV0 - AntiLambda
    }
    
}
AliAnalysisTaskPOmegaPenne::~AliAnalysisTaskPOmegaPenne()       // Destructor
{
    delete aaEvent;
    delete aaTrack;
    delete fEvent;
    delete fTrack;
    delete fCascade;
    delete fEventCuts;
    // delete fTrackCutsProton;
    // delete fTrackCutsAntiProton;
    delete fv0;
    delete fLambdaV0Cuts;
    delete fAntiLambdaV0Cuts;
    delete fCascadeCutsXi;
    delete fCascadeCutsAntiXi;
    delete fConfig;
    delete fPairCleaner;
    delete fPartColl;
    delete *fGTI;
    delete tlEventCuts;
    // delete tlTrackCutsProton;
    // delete tlAntiTrackCutsProton;
    delete tlLambdaList;
    delete tlAntiLambdaList;
    delete tlCascadeCutsXi;
    delete tlAntiCascadeCutsXi;
    delete tlResults;
    delete tlResultsQA;
    // delete tlProtonMC;
    // delete tlAntiProtonMC;
    delete tlLambdaMC;
    delete tlAntiLambdaMC;
    }

// // Copy Constructor
AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne(const AliAnalysisTaskPOmegaPenne &obj) : AliAnalysisTaskSE(obj),
                                                                                                fIsMC(obj.fIsMC),
                                                                                                aaEvent(obj.aaEvent),
                                                                                                aaTrack(obj.aaTrack),
                                                                                                fEvent(obj.fEvent),
                                                                                                fTrack(obj.fTrack),
                                                                                                fCascade(obj.fCascade),
                                                                                                fEventCuts(obj.fEventCuts),
                                                                                                // fTrackCutsProton(obj.fTrackCutsProton),
                                                                                                // fTrackCutsAntiProton(obj.fTrackCutsAntiProton),
                                                                                                fv0(obj.fv0),
                                                                                                fLambdaV0Cuts(obj.fLambdaV0Cuts),
                                                                                                fAntiLambdaV0Cuts(obj.fAntiLambdaV0Cuts),
                                                                                                fCascadeCutsXi(obj.fCascadeCutsXi),
                                                                                                fCascadeCutsAntiXi(obj.fCascadeCutsAntiXi),
                                                                                                fConfig(obj.fConfig),
                                                                                                fPairCleaner(obj.fPairCleaner),
                                                                                                fPartColl(obj.fPartColl),
                                                                                                fGTI(obj.fGTI),
                                                                                                fTrackBufferSize(obj.fTrackBufferSize),
                                                                                                tlEventCuts(obj.tlEventCuts),
                                                                                                // tlTrackCutsProton(obj.tlTrackCutsProton),
                                                                                                // tlAntiTrackCutsProton(obj.tlAntiTrackCutsProton),
                                                                                                tlLambdaList(obj.tlLambdaList),
                                                                                                tlAntiLambdaList(obj.tlAntiLambdaList),
                                                                                                tlCascadeCutsXi(obj.tlCascadeCutsXi),
                                                                                                tlAntiCascadeCutsXi(obj.tlAntiCascadeCutsXi),
                                                                                                tlResults(obj.tlResults),
                                                                                                tlResultsQA(obj.tlResultsQA),
                                                                                                // tlProtonMC(obj.tlProtonMC),
                                                                                                // tlAntiProtonMC(obj.tlAntiProtonMC),
                                                                                                tlLambdaMC(obj.tlLambdaMC),
                                                                                                tlAntiLambdaMC(obj.tlAntiLambdaMC)

{
}

// AliAnalysisTaskPOmegaPenne& AliAnalysisTaskPOmegaPenne::operator=(const AliAnalysisTaskPOmegaPenne &other)
// {
//     AliAnalysisTaskSE::operator=(other);
//     this->fIsMC = other.fIsMC;
//     this->aaEvent = other.aaEvent;
//     this->aaTrack = other.aaTrack;
//     this->fOutput = other.fOutput;
//     this->fEvent = other.fEvent;
//     this->fTrack = other.fTrack;
//     this->fCascade = other.fCascade;
//     this->fEventCuts = other.fEventCuts;
//     this->fTrackCutsProton = other.fTrackCutsProton;
//     this->fTrackCutsAntiProton = other.fTrackCutsAntiProton;
//     this->fCascadeCutsXi = other.fCascadeCutsXi;
//     this->fCascadeCutsAntiXi = other.fCascadeCutsAntiXi;
//     this->fConfig = other.fConfig;
//     this->fPairCleaner = other.fPairCleaner;
//     this->fPartColl = other.fPartColl;
//     this->fGTI = other.fGTI;
//     this->fTrackBufferSize = other.fTrackBufferSize;

//     return *this;
// }

void AliAnalysisTaskPOmegaPenne::UserCreateOutputObjects()
{
   
    fEvent = new AliFemtoDreamEvent(false, true, GetCollisionCandidates());
    fTrack = new AliFemtoDreamTrack();
    fTrack->SetUseMCInfo(fIsMC);
    fGTI = new AliAODTrack *[fTrackBufferSize];
    
    fEventCuts->InitQA();

    // // Proton Cuts      ###########
    // if (!fTrackCutsProton){AliFatal("Track Cuts for Particle Proton not set!");}
    // fTrackCutsProton->Init();
    // fTrackCutsProton->SetName("Protons");
    // // ##

    // // AntiProton Cuts  ###########
    // if (!fTrackCutsAntiProton){AliFatal("Track Cuts for Particle AntiProton not set!");}
    // fTrackCutsAntiProton->Init();
    // fTrackCutsAntiProton->SetName("AntiProtons");
    // // ##
 
    // Lambda Cutys    ###########
    if (!fLambdaV0Cuts){AliFatal("Track Cuts for Particle Lambda not set!");}
    fLambdaV0Cuts->Init();
    fLambdaV0Cuts->SetName("Lambda");
    // ##

    // AntiLambda Cutys    ###########
    if (!fAntiLambdaV0Cuts){AliFatal("Track Cuts for Particle AntiLambda not set!");}
    fAntiLambdaV0Cuts->Init();
    fAntiLambdaV0Cuts->SetName("AntiLambda");
    // ##

    // V0 Candidates
    fv0 = new AliFemtoDreamv0();
    fv0->SetUseMCInfo(fIsMC);
    fv0->GetPosDaughter()->SetUseMCInfo(fIsMC); 
    fv0->GetNegDaughter()->SetUseMCInfo(fIsMC); 
    fv0->SetPDGCode(3122);
    fv0->SetPDGDaughterPos(2212);
    fv0->SetPDGDaughterNeg(211);
    // ##

    // Xi Cuts    ###########
    if (!fCascadeCutsXi){AliFatal("Track Cuts for Particle Xi not set!");}
    fCascadeCutsXi->Init();
    fCascadeCutsXi->SetName("Xi");
    // ##
    
    // AntiXi Cuts    ###########
    if (!fCascadeCutsAntiXi){AliFatal("Track Cuts for Particle AntiXi not set!");}
    fCascadeCutsAntiXi->Init();
    fCascadeCutsAntiXi->SetName("AntiXi");
    // ##

    // Cascade Cuts     #########
    fCascade = new AliFemtoDreamCascade();          // Initial Cascade Object
    fCascade->SetUseMCInfo(fIsMC);
    //PDG Codes should be set assuming Xi- to also work for Xi+
    fCascade->SetPDGCode(3312);
    fCascade->SetPDGDaugPos(2212);
    fCascade->GetPosDaug()->SetUseMCInfo(fIsMC);
    fCascade->SetPDGDaugNeg(211);
    fCascade->GetNegDaug()->SetUseMCInfo(fIsMC);
    fCascade->SetPDGDaugBach(211);
    fCascade->GetBach()->SetUseMCInfo(fIsMC);
    fCascade->Setv0PDGCode(3122);
    // ##

    fPairCleaner = new AliFemtoDreamPairCleaner(2, 2, false);
    fPartColl = new AliFemtoDreamPartCollection(fConfig, false);
    
    tlCascadeCutsXi = new TList();
    tlCascadeCutsXi->SetName("XiCascade");
    tlCascadeCutsXi->SetOwner();

    tlAntiCascadeCutsXi = new TList();
    tlAntiCascadeCutsXi->SetName("AntiXiCascade");
    tlAntiCascadeCutsXi->SetOwner();

    tlResultsQA = new TList();
    tlResultsQA->SetName("ResultsQA");
    tlResultsQA->SetOwner();

    // Connect Cuts to OutputContainers
    tlEventCuts             = fEventCuts->GetHistList();
    // tlTrackCutsProton       = fTrackCutsProton->GetQAHists();
    // tlAntiTrackCutsProton   = fTrackCutsAntiProton->GetQAHists();
    tlLambdaList            = fLambdaV0Cuts->GetQAHists();
    tlAntiLambdaList        = fAntiLambdaV0Cuts->GetQAHists();
    tlCascadeCutsXi         = fCascadeCutsXi->GetQAHists();
    tlAntiCascadeCutsXi     = fCascadeCutsAntiXi->GetQAHists();
    tlResults               = fPartColl->GetHistList();
    tlResultsQA->Add(         fPartColl->GetQAList());
    tlResultsQA->Add(         fPairCleaner->GetHistList());
    tlResultsQA->Add(         fEvent->GetEvtCutList());

    PostData(1, tlEventCuts);
    // PostData(2, tlTrackCutsProton);
    // PostData(3, tlAntiTrackCutsProton);
    PostData(2, tlLambdaList);
    PostData(3, tlAntiLambdaList);
    PostData(4, tlCascadeCutsXi);
    PostData(5, tlAntiCascadeCutsXi);
    PostData(6, tlResults);
    PostData(7, tlResultsQA);

    // if (fTrackCutsProton->GetIsMonteCarlo())
    // {
    //     tlProtonMC = fTrackCutsProton->GetMCQAHists();
    //     PostData(10, tlProtonMC);
    // }
    // if (fTrackCutsAntiProton->GetIsMonteCarlo())
    // {
    //     tlAntiProtonMC = fTrackCutsAntiProton->GetMCQAHists();
    //     PostData(11, tlAntiProtonMC);
    // }
    if (fLambdaV0Cuts->GetIsMonteCarlo())
    {
        tlLambdaMC = fLambdaV0Cuts->GetMCQAHists();
        PostData(8, tlLambdaMC);
    }
    if (fAntiLambdaV0Cuts->GetIsMonteCarlo())
    {
        tlAntiLambdaMC = fAntiLambdaV0Cuts->GetMCQAHists();
        PostData(9, tlAntiLambdaMC);
    }
}

// static std::vector<AliFemtoDreamBasePart> vProtons;         // Particle Vectors  
// static std::vector<AliFemtoDreamBasePart> vAntiProtons;     
static std::vector<AliFemtoDreamBasePart> vLambda;           
static std::vector<AliFemtoDreamBasePart> vAntiLambda;       
static std::vector<AliFemtoDreamBasePart> vXi;           
static std::vector<AliFemtoDreamBasePart> vAntiXi;       

void AliAnalysisTaskPOmegaPenne::UserExec(Option_t *)
{
    aaEvent = dynamic_cast<AliAODEvent *>(fInputEvent);
    
    if (!aaEvent)
    {
        AliWarning("No Input aaEvent");
    }
    else
    {
        fEvent->SetEvent(aaEvent);
        if (fEventCuts->isSelected(fEvent))
        {
            ResetGlobalTrackReference();
            for (int iTrack = 0; iTrack < aaEvent->GetNumberOfTracks(); ++iTrack)
            {
                aaTrack = dynamic_cast<AliAODTrack *>(aaEvent->GetTrack(iTrack));
                if (!aaTrack)
                {
                    AliFatal("No Standard AOD");
                    return;
                }
                StoreGlobalTrackReference(aaTrack);
            }
           
            fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

            // vProtons.clear();
            // vAntiProtons.clear();
            vXi.clear();
            vAntiXi.clear();
            vLambda.clear();
            vAntiLambda.clear();
        

            // for (int iTrack = 0; iTrack < aaEvent->GetNumberOfTracks(); ++iTrack)
            // {
            //     aaTrack = dynamic_cast<AliAODTrack *>(aaEvent->GetTrack(iTrack));
            //     if (!aaTrack)
            //     {
            //         AliFatal("No Standard AOD");
            //         return;
            //     }
            //     fTrack->SetTrack(aaTrack);

            //     // mark track (anti-)proton and/or (anti-)xi
            //     if (fTrackCutsProton->isSelected(fTrack))
            //     {
            //         vProtons.push_back(*fTrack);
            //     }
            //     if (fTrackCutsAntiProton->isSelected(fTrack))
            //     {
            //         vAntiProtons.push_back(*fTrack);
            //     }
            // }

            // ## Lambda Selection
            fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
            for (int iv0 = 0; iv0 < static_cast<TClonesArray *>(aaEvent->GetV0s())->GetEntriesFast(); ++iv0)
            {
                AliAODv0 *v0 = aaEvent->GetV0(iv0);
                fv0->Setv0(aaEvent, v0, fEvent->GetMultiplicity());
                fLambdaV0Cuts->isSelected(fv0);
                fAntiLambdaV0Cuts->isSelected(fv0);
            }

            // ## Xi selection
            for (int iCasc = 0; iCasc < static_cast<TClonesArray *>(aaEvent->GetCascades())->GetEntriesFast(); ++iCasc)
            {
                AliAODcascade *casc = aaEvent->GetCascade(iCasc);
                fCascade->SetCascade(aaEvent, casc);
                if (fCascadeCutsXi->isSelected(fCascade))
                {
                    vXi.push_back(*fCascade);
                }
                if (fCascadeCutsAntiXi->isSelected(fCascade))
                {
                    vAntiXi.push_back(*fCascade);
                }
            }                                                                         
            // remove double-matched tracks
            fPairCleaner->ResetArray();
            // fPairCleaner->CleanTrackAndDecay(&vProtons, &vXi, 0);
            // fPairCleaner->CleanTrackAndDecay(&vAntiProtons, &vAntiXi, 1);
            fPairCleaner->CleanTrackAndDecay(&vLambda, &vXi, 0);
            fPairCleaner->CleanTrackAndDecay(&vAntiLambda, &vAntiXi, 1);
            
            fPairCleaner->CleanDecay(&vXi, 0);
            fPairCleaner->CleanDecay(&vAntiXi, 1);
            
            // fPairCleaner->StoreParticle(vProtons);
            // fPairCleaner->StoreParticle(vAntiProtons);

            fPairCleaner->StoreParticle(vLambda);
            fPairCleaner->StoreParticle(vAntiLambda);

            fPairCleaner->StoreParticle(vXi);
            fPairCleaner->StoreParticle(vAntiXi);

            // lambdas nicht in storeparticle weil sonst mit setevent pairQA betrieben wird was wir nicht brauchen
            fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality());
            // soweit ich das richtig verstanden habe wird pairQA mit den teilchen gemacht die im pairCleaner 
            // sind und pdgCodes in der richtigen Reihenfolge vorhanden sind.


            PostData(1, tlEventCuts);
            // PostData(2, tlTrackCutsProton);
            // PostData(3, tlAntiTrackCutsProton);
            PostData(2, tlLambdaList);
            PostData(3, tlAntiLambdaList);
            PostData(4, tlCascadeCutsXi);
            PostData(5, tlAntiCascadeCutsXi);
            PostData(6, tlResults);
            PostData(7, tlResultsQA);
            if (fIsMC)
            {
                // PostData(10, tlProtonMC);
            
                // PostData(11, tlAntiProtonMC);
            
                PostData(8, tlLambdaMC);
            
                PostData(9, tlAntiLambdaMC);
            }
        }
    }
}

void AliAnalysisTaskPOmegaPenne::ResetGlobalTrackReference()
{
    //This method was inherited form H. Beck analysis
    for (UShort_t i = 0; i < fTrackBufferSize; i++)
    {
        fGTI[i] = nullptr;
        // std::fill(fGTI.begin(),fGTI.end(), nullptr);
    }
}

//  Stores TrackID in Global Track Reference Array 'fGTI' if ID > 0
//
void AliAnalysisTaskPOmegaPenne::StoreGlobalTrackReference(AliAODTrack *aaTrack)
{
    //This method was inherited form H. Beck analysis

    const int trackID = aaTrack->GetID();
    if (trackID < 0)
    {
        return;
    }
    if (trackID >= fTrackBufferSize)
    {
        printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID, fTrackBufferSize);
        return;
    }

    if (fGTI[trackID])
    {
        if ((!aaTrack->GetFilterMap()) && (!aaTrack->GetTPCNcls()))
        {
            return;
        }
        if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls())
        {
            printf("WARNING! global track info already there!");
            printf("    ###     TPCNcls track1 %u Track2 %u", (fGTI[trackID])->GetTPCNcls(), aaTrack->GetTPCNcls());
            printf("   ###     FilterMap Track1 %u track2 %u\n", (fGTI[trackID])->GetFilterMap(), aaTrack->GetFilterMap());
        }
    }
    fGTI[trackID] = aaTrack;
}
