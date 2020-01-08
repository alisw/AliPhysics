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
                                                                Event(nullptr),
                                                                track(nullptr),
                                                                fOutput(nullptr),
                                                                fEvent(nullptr),
                                                                fTrack(nullptr),
                                                                fEventCuts(nullptr),
                                                                fTrackCutsProton(nullptr),
                                                                fTrackCutsAntiProton(nullptr),
                                                                fTrackCutsKaon(nullptr),
                                                                fTrackCutsAntiKaon(nullptr),
                                                                fConfig(nullptr),
                                                                fPairCleaner(nullptr),
                                                                fPartColl(nullptr),
                                                                fGTI(nullptr),
                                                                fTrackBufferSize(2000)
{
}

AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne(const char *name, bool isMC) :   AliAnalysisTaskSE(name),
                                                                                        fIsMC(isMC),
                                                                                        Event(nullptr),
                                                                                        track(nullptr),
                                                                                        fOutput(nullptr),
                                                                                        fEvent(nullptr),
                                                                                        fTrack(nullptr),
                                                                                        fEventCuts(nullptr),
                                                                                        fTrackCutsProton(nullptr),
                                                                                        fTrackCutsAntiProton(nullptr),
                                                                                        fTrackCutsKaon(nullptr),
                                                                                        fTrackCutsAntiKaon(nullptr),
                                                                                        fConfig(nullptr),
                                                                                        fPairCleaner(nullptr),
                                                                                        fPartColl(nullptr),
                                                                                        fGTI(nullptr),
                                                                                        fTrackBufferSize(2000)
{
    DefineOutput(1, TList::Class());
}
AliAnalysisTaskPOmegaPenne::~AliAnalysisTaskPOmegaPenne()
{
    // TODO Auto-generated destructor stub
}

void AliAnalysisTaskPOmegaPenne::UserCreateOutputObjects()
{
    fOutput = new TList();
    fOutput->SetName("Output");
    fOutput->SetOwner();

    fEvent = new AliFemtoDreamEvent(false, true, GetCollisionCandidates());
    fOutput->Add(fEvent->GetEvtCutList());
    fTrack = new AliFemtoDreamTrack();
    fTrack->SetUseMCInfo(fIsMC);
    fGTI = new AliAODTrack *[fTrackBufferSize];
    
    fEventCuts->InitQA();
    fOutput->Add(fEventCuts->GetHistList());

    // Proton Cuts      ###########
    if (!fTrackCutsProton)
    {
        AliFatal("Track Cuts for Particle Proton not set!");
    }
    fTrackCutsProton->Init();
    fTrackCutsProton->SetName("Protonen");
    fOutput->Add(fTrackCutsProton->GetQAHists());
    if (fTrackCutsProton->GetIsMonteCarlo())
    {
        fTrackCutsProton->SetMCName("MCProtonen");
        fOutput->Add(fTrackCutsProton->GetMCQAHists());
    }
    // ##

    // Antiproton Cuts  ###########
    if (!fTrackCutsAntiProton)
    {
        AliFatal("Track Cuts for Particle AntiProton not set!");
    }
    fTrackCutsAntiProton->Init();
    fTrackCutsAntiProton->SetName("AntiProtonen");
    fOutput->Add(fTrackCutsAntiProton->GetQAHists());
    if (fTrackCutsAntiProton->GetIsMonteCarlo())
    {
        fTrackCutsAntiProton->SetMCName("MCAntiProtonen");
        fOutput->Add(fTrackCutsAntiProton->GetMCQAHists());
    }
    // ##

    // Kaon Cuts    ###########
    if (!fTrackCutsKaon)
    {
        AliFatal("Track Cuts for Particle Kaon not set!");
    }
    fTrackCutsKaon->Init();
    fTrackCutsKaon->SetName("Kaonen");
    fOutput->Add(fTrackCutsKaon->GetQAHists());
    if (fTrackCutsKaon->GetIsMonteCarlo())
    {
        fTrackCutsKaon->SetMCName("MCKaonen");
        fOutput->Add(fTrackCutsKaon->GetMCQAHists());
    }
    // ##
    
    // Antikaon Cuts    ###########
    if (!fTrackCutsAntiKaon)
    {
        AliFatal("Track Cuts for Particle AntiKaon not set!");
    }
    fTrackCutsAntiKaon->Init();
    fTrackCutsAntiKaon->SetName("AntiKaonen");
    fOutput->Add(fTrackCutsAntiKaon->GetQAHists());
    if (fTrackCutsAntiKaon->GetIsMonteCarlo())
    {
        fTrackCutsAntiKaon->SetMCName("MCAntiKaonen");
        fOutput->Add(fTrackCutsAntiKaon->GetMCQAHists());
    }
    // ##

    fPairCleaner = new AliFemtoDreamPairCleaner(0, 2, false);
    fOutput->Add(fPairCleaner->GetHistList());
    fPartColl = new AliFemtoDreamPartCollection(fConfig, false);
    fOutput->Add(fPartColl->GetHistList());
    fOutput->Add(fPartColl->GetQAList());
    PostData(1, fOutput);
}

void AliAnalysisTaskPOmegaPenne::UserExec(Option_t *)
{
    Event = dynamic_cast<AliAODEvent *>(fInputEvent);
    //AliAODEvent *Event = static_cast<AliAODEvent*>(InputEvent());
    
    if (!Event)
    {
        AliWarning("No Input Event");
    }
    else
    {
        fEvent->SetEvent(Event);
        if (fEventCuts->isSelected(fEvent))
        {
            ResetGlobalTrackReference();
            for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack)
            {
                track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));
                if (!track)
                {
                    AliFatal("No Standard AOD");
                    return;
                }
                StoreGlobalTrackReference(track);
            }
           
            fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

            std::vector<AliFemtoDreamBasePart> vProtons;         // Particle Vectors  
            std::vector<AliFemtoDreamBasePart> vAntiProtons;     
            std::vector<AliFemtoDreamBasePart> vKaons;           
            std::vector<AliFemtoDreamBasePart> vAntiKaons;       

            for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack)
            {
                track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));
                if (!track)
                {
                    AliFatal("No Standard AOD");
                    return;
                }
                fTrack->SetTrack(track);

                // mark track (anti-)proton and/or (anti-)kaon
                if (fTrackCutsProton->isSelected(fTrack))
                {
                    vProtons.push_back(*fTrack);
                }
                if (fTrackCutsAntiProton->isSelected(fTrack))
                {
                    vAntiProtons.push_back(*fTrack);
                }
                if (fTrackCutsKaon->isSelected(fTrack))
                {
                    vKaons.push_back(*fTrack);
                }
                if (fTrackCutsAntiKaon->isSelected(fTrack))
                {
                    vAntiKaons.push_back(*fTrack);
                }
            }
            
            // remove double-matched tracks
            fPairCleaner->CleanDecayAndDecay(&vProtons, &vKaons, 0);
            fPairCleaner->CleanDecayAndDecay(&vAntiProtons, &vAntiKaons, 0);
            fPairCleaner->ResetArray();
            fPairCleaner->StoreParticle(vProtons);
            fPairCleaner->StoreParticle(vAntiProtons);
            fPairCleaner->StoreParticle(vKaons);
            fPairCleaner->StoreParticle(vAntiKaons);

            fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

            PostData(1, fOutput);
        }
    }
}

void AliAnalysisTaskPOmegaPenne::ResetGlobalTrackReference()
{
    //This method was inherited form H. Beck analysis
    for (UShort_t i = 0; i < fTrackBufferSize; i++)
    {
        fGTI[i] = 0;
    }
}
void AliAnalysisTaskPOmegaPenne::StoreGlobalTrackReference(AliAODTrack *track)
{
    //This method was inherited form H. Beck analysis

    const int trackID = track->GetID();
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
        if ((!track->GetFilterMap()) && (!track->GetTPCNcls()))
        {
            return;
        }
        if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls())
        {
            printf("WARNING! global track info already there!");
            printf("    ###     TPCNcls track1 %u Track2 %u", (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
            printf("   ###     FilterMap Track1 %u track2 %u\n", (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
        }
    }
    fGTI[trackID] = track;
}
