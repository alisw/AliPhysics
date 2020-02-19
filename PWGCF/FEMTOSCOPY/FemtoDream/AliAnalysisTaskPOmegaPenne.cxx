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
                                                                fTrackCutsProton(0),
                                                                fTrackCutsAntiProton(0),
                                                                fv0(0),
                                                                fLambdaV0Cuts(0),
                                                                fAntiLambdaV0Cuts(0),
                                                                fCascadeCutsXion(0),
                                                                fCascadeCutsAntiXion(0),
                                                                fConfig(0),
                                                                fPairCleaner(0),
                                                                fPartColl(0),
                                                                fGTI(0),
                                                                fTrackBufferSize(10000),
                                                                tlEventCuts(0),
                                                                tlTrackCutsProton(0),
                                                                tlAntiTrackCutsProton(0),
                                                                tlLambdaList(0),
                                                                tlAntiLambdaList(0),
                                                                tlCascadeCutsXi(0),
                                                                tlAntiCascadeCutsXi(0),
                                                                tlResults(0),
                                                                tlResultsQA(0)
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
                                                                                      fTrackCutsProton(0),
                                                                                      fTrackCutsAntiProton(0),
                                                                                      fv0(0),
                                                                                      fLambdaV0Cuts(0),
                                                                                      fAntiLambdaV0Cuts(0),
                                                                                      fCascadeCutsXion(0),
                                                                                      fCascadeCutsAntiXion(0),
                                                                                      fConfig(0),
                                                                                      fPairCleaner(0),
                                                                                      fPartColl(0),
                                                                                      fGTI(0),
                                                                                      fTrackBufferSize(10000),
                                                                                      tlEventCuts(0),
                                                                                      tlTrackCutsProton(0),
                                                                                      tlAntiTrackCutsProton(0),
                                                                                      tlLambdaList(0),
                                                                                      tlAntiLambdaList(0),
                                                                                      tlCascadeCutsXi(0),
                                                                                      tlAntiCascadeCutsXi(0),
                                                                                      tlResults(0),
                                                                                      tlResultsQA(0)
{
    DefineOutput(1, TList::Class());    // Event Cuts
    DefineOutput(2, TList::Class());    // Proton Track Cuts
    DefineOutput(3, TList::Class());    // Anti Proton Track Cuts
    DefineOutput(4, TList::Class());    // Lambda Track Cuts
    DefineOutput(5, TList::Class());    // Anti Lambda Track Cuts
    DefineOutput(6, TList::Class());    // Xi Track Cuts
    DefineOutput(7, TList::Class());    // Anti Xi Track Cuts
    DefineOutput(8, TList::Class());    // Results
    DefineOutput(9, TList::Class());    // QA Results
}
AliAnalysisTaskPOmegaPenne::~AliAnalysisTaskPOmegaPenne()
{
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
                                                                                                fTrackCutsProton(obj.fTrackCutsProton),
                                                                                                fTrackCutsAntiProton(obj.fTrackCutsAntiProton),
                                                                                                fv0(obj.fv0),
                                                                                                fLambdaV0Cuts(obj.fLambdaV0Cuts),
                                                                                                fAntiLambdaV0Cuts(obj.fAntiLambdaV0Cuts),
                                                                                                fCascadeCutsXion(obj.fCascadeCutsXion),
                                                                                                fCascadeCutsAntiXion(obj.fCascadeCutsAntiXion),
                                                                                                fConfig(obj.fConfig),
                                                                                                fPairCleaner(obj.fPairCleaner),
                                                                                                fPartColl(obj.fPartColl),
                                                                                                fGTI(obj.fGTI),
                                                                                                fTrackBufferSize(obj.fTrackBufferSize),
                                                                                                tlEventCuts(obj.tlEventCuts),
                                                                                                tlTrackCutsProton(obj.tlTrackCutsProton),
                                                                                                tlAntiTrackCutsProton(obj.tlAntiTrackCutsProton),
                                                                                                tlLambdaList(obj.tlLambdaList),
                                                                                                tlAntiLambdaList(obj.tlAntiLambdaList),
                                                                                                tlCascadeCutsXi(obj.tlCascadeCutsXi),
                                                                                                tlAntiCascadeCutsXi(obj.tlAntiCascadeCutsXi),
                                                                                                tlResults(obj.tlResults),
                                                                                                tlResultsQA(obj.tlResultsQA)
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
//     this->fCascadeCutsXion = other.fCascadeCutsXion;
//     this->fCascadeCutsAntiXion = other.fCascadeCutsAntiXion;
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

    // Proton Cuts      ###########
    if (!fTrackCutsProton){AliFatal("Track Cuts for Particle Proton not set!");}
    fTrackCutsProton->Init();
    fTrackCutsProton->SetName("Protons");
    // ##

    // AntiProton Cuts  ###########
    if (!fTrackCutsAntiProton){AliFatal("Track Cuts for Particle AntiProton not set!");}
    fTrackCutsAntiProton->Init();
    fTrackCutsAntiProton->SetName("AntiProtons");
    // ##
 
    // Lambda Cutys    ###########
    fLambdaV0Cuts->Init();
    // ##

    // AntiLambda Cutys    ###########
    fAntiLambdaV0Cuts->Init();
    // ##

    // V0 Candidates
    fv0 = new AliFemtoDreamv0();
    fv0->SetPDGCode(3122);
    fv0->SetPDGDaughterPos(2212);
    fv0->SetPDGDaughterNeg(211);
    // ##

    // Xion Cuts    ###########
    if (!fCascadeCutsXion){AliFatal("Track Cuts for Particle Xi not set!");}
    fCascadeCutsXion->Init();
    fCascadeCutsXion->SetName("Xions");
    // ##
    
    // AntiXion Cuts    ###########
    if (!fCascadeCutsAntiXion){AliFatal("Track Cuts for Particle AntiXi not set!");}
    fCascadeCutsAntiXion->Init();
    fCascadeCutsAntiXion->SetName("AntiXions");
    // ##

    // Cascade Cuts     #########
    fCascade = new AliFemtoDreamCascade();          // Initial Cascade Object
    fCascade->SetUseMCInfo(fCascadeCutsXion->GetIsMonteCarlo() || fCascadeCutsAntiXion->GetIsMonteCarlo());
    //PDG Codes should be set assuming Xi- to also work for Xi+
    fCascade->SetPDGCode(3312);
    fCascade->SetPDGDaugPos(2212);
    fCascade->GetPosDaug()->SetUseMCInfo(fCascadeCutsXion->GetIsMonteCarlo() || fCascadeCutsAntiXion->GetIsMonteCarlo());
    fCascade->SetPDGDaugNeg(211);
    fCascade->GetNegDaug()->SetUseMCInfo(fCascadeCutsXion->GetIsMonteCarlo() || fCascadeCutsAntiXion->GetIsMonteCarlo());
    fCascade->SetPDGDaugBach(211);
    fCascade->GetBach()->SetUseMCInfo(fCascadeCutsXion->GetIsMonteCarlo() || fCascadeCutsAntiXion->GetIsMonteCarlo());
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
    tlTrackCutsProton       = fTrackCutsProton->GetQAHists();
    tlAntiTrackCutsProton   = fTrackCutsAntiProton->GetQAHists();
    tlLambdaList            = fLambdaV0Cuts->GetQAHists();
    tlAntiLambdaList        = fAntiLambdaV0Cuts->GetQAHists();
    tlCascadeCutsXi->Add(     fCascadeCutsXion->GetQAHists());
    tlAntiCascadeCutsXi->Add( fCascadeCutsAntiXion->GetQAHists());
    tlResults               = fPartColl->GetHistList();
    tlResultsQA->Add(         fPartColl->GetQAList());
    tlResultsQA->Add(         fPairCleaner->GetHistList());
    tlResultsQA->Add(         fEvent->GetEvtCutList());

    PostData(1, tlEventCuts);
    PostData(2, tlTrackCutsProton);
    PostData(3, tlAntiTrackCutsProton);
    PostData(4, tlLambdaList);
    PostData(5, tlAntiLambdaList);
    PostData(6, tlCascadeCutsXi);
    PostData(7, tlAntiCascadeCutsXi);
    PostData(8, tlResults);
    PostData(9, tlResultsQA);
}

static std::vector<AliFemtoDreamBasePart> vProtons;         // Particle Vectors  
static std::vector<AliFemtoDreamBasePart> vAntiProtons;     
static std::vector<AliFemtoDreamBasePart> vXions;           
static std::vector<AliFemtoDreamBasePart> vAntiXions;       

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

            vProtons.clear();
            vAntiProtons.clear();
            vXions.clear();
            vAntiXions.clear();
            
            for (int iTrack = 0; iTrack < aaEvent->GetNumberOfTracks(); ++iTrack)
            {
                aaTrack = dynamic_cast<AliAODTrack *>(aaEvent->GetTrack(iTrack));
                if (!aaTrack)
                {
                    AliFatal("No Standard AOD");
                    return;
                }
                fTrack->SetTrack(aaTrack);

                // mark track (anti-)proton and/or (anti-)xion
                if (fTrackCutsProton->isSelected(fTrack))
                {
                    vProtons.push_back(*fTrack);
                }
                if (fTrackCutsAntiProton->isSelected(fTrack))
                {
                    vAntiProtons.push_back(*fTrack);
                }
            }

            // ## Lambda Selection
            fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
            for (int iv0 = 0; iv0 < static_cast<TClonesArray *>(aaEvent->GetV0s())->GetEntriesFast(); ++iv0)
            {
                AliAODv0 *v0 = aaEvent->GetV0(iv0);
                fv0->Setv0(aaEvent, v0, fEvent->GetMultiplicity());
                fLambdaV0Cuts->isSelected(fv0);
                fAntiLambdaV0Cuts->isSelected(fv0);
            }

            for (int iCasc = 0; iCasc < static_cast<TClonesArray *>(aaEvent->GetCascades())->GetEntriesFast(); ++iCasc)
            {
                AliAODcascade *casc = aaEvent->GetCascade(iCasc);
                fCascade->SetCascade(aaEvent, casc);
                if (fCascadeCutsXion->isSelected(fCascade))
                {
                    vXions.push_back(*fCascade);
                }
                if (fCascadeCutsAntiXion->isSelected(fCascade))
                {
                    vAntiXions.push_back(*fCascade);
                }
            }                                                                         
            // remove double-matched tracks
            fPairCleaner->ResetArray();
            fPairCleaner->CleanTrackAndDecay(&vProtons, &vXions, 0);
            fPairCleaner->CleanTrackAndDecay(&vAntiProtons, &vAntiXions, 1);
            
            fPairCleaner->CleanDecay(&vXions, 0);
            fPairCleaner->CleanDecay(&vAntiXions, 1);
            
            fPairCleaner->StoreParticle(vProtons);
            fPairCleaner->StoreParticle(vAntiProtons);
            fPairCleaner->StoreParticle(vXions);
            fPairCleaner->StoreParticle(vAntiXions);

            // lambdas nicht in storeparticle weil sonst mit setevent pairQA betrieben wird was wir nicht brauchen
            fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

            PostData(1, tlEventCuts);
            PostData(2, tlTrackCutsProton);
            PostData(3, tlAntiTrackCutsProton);
            PostData(4, tlLambdaList);
            PostData(5, tlAntiLambdaList);
            PostData(6, tlCascadeCutsXi);
            PostData(7, tlAntiCascadeCutsXi);
            PostData(8, tlResults);
            PostData(9, tlResultsQA);
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
