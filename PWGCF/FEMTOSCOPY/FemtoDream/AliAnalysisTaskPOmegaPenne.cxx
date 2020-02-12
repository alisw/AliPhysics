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
                                                                Event(0),
                                                                track(0),
                                                                fOutput(0),
                                                                fEvent(0),
                                                                fTrack(0),
                                                                fCascade(0),
                                                                fEventCuts(0),
                                                                fTrackCutsProton(0),
                                                                fTrackCutsAntiProton(0),
                                                                fCascadeCutsXion(0),
                                                                fCascadeCutsAntiXion(0),
                                                                fConfig(0),
                                                                fPairCleaner(0),
                                                                fPartColl(0),
                                                                fGTI(0),
                                                                fTrackBufferSize(10000)
{
}
AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne(const char *name, bool isMC) :   AliAnalysisTaskSE(name),
                                                                                        fIsMC(isMC),
                                                                                        Event(0),
                                                                                        track(0),
                                                                                        fOutput(0),
                                                                                        fEvent(0),
                                                                                        fTrack(0),
                                                                                        fCascade(0),
                                                                                        fEventCuts(0),
                                                                                        fTrackCutsProton(0),
                                                                                        fTrackCutsAntiProton(0),
                                                                                        fCascadeCutsXion(0),
                                                                                        fCascadeCutsAntiXion(0),
                                                                                        fConfig(0),
                                                                                        fPairCleaner(0),
                                                                                        fPartColl(0),
                                                                                        fGTI(0),
                                                                                        fTrackBufferSize(10000)
{
    DefineOutput(1, TList::Class());
}
AliAnalysisTaskPOmegaPenne::~AliAnalysisTaskPOmegaPenne()
{
    // TODO Auto-generated destructor stub
    // |-> HÄ, ne gar nicht!... wenn der run zu ende ist hört das Objekt einfach auf zu evXionstieren. Scheiß auf den Destruktor!
    // pass lieber in der UserExec Methode auf!!!
}

// // Copy Constructor
// AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne(const AliAnalysisTaskPOmegaPenne& obj) : AliAnalysisTaskSE(obj),
//                                                                                                 fIsMC(obj.fIsMC),
//                                                                                                 Event(obj.Event),
//                                                                                                 track(obj.track),
//                                                                                                 fOutput(obj.fOutput),
//                                                                                                 fEvent(obj.fEvent),
//                                                                                                 fTrack(obj.fTrack),
//                                                                                                 fEventCuts(obj.fEventCuts),
//                                                                                                 fTrackCutsProton(obj.fTrackCutsProton),
//                                                                                                 fTrackCutsAntiProton(obj.fTrackCutsAntiProton),
//                                                                                                 fCascadeCutsXion(obj.fCascadeCutsXion),
//                                                                                                 fCascadeCutsAntiXion(obj.fCascadeCutsAntiXion),
//                                                                                                 fConfig(obj.fConfig),
//                                                                                                 fPairCleaner(obj.fPairCleaner),
//                                                                                                 fPartColl(obj.fPartColl),
//                                                                                                 fGTI(obj.fGTI),
//                                                                                                 fTrackBufferSize(10000)
// {
//     if (obj.vProtons)
//     {
//         vProtons = obj.vProtons;
//     }
//     if (obj.vAntiProtons)
//     {
//         vProtons = obj.vAntiProtons;
//     }
//     if (obj.vXions)
//     {
//         vProtons = obj.vXions;
//     }
//     if (obj.vAntiXions)
//     {
//         vProtons = obj.vAntiXions;
//     }
// }

// AliAnalysisTaskPOmegaPenne& AliAnalysisTaskPOmegaPenne::operator=(const AliAnalysisTaskPOmegaPenne& other)
// {
//     AliAnalysisTaskSE::operator=(other);
//     this.fIsMC = other.fIsMC;
//     this.Event = other.Event;
//     this.track = other.track;
//     this.fOutput = other.fOutput;
//     this.fEvent = other.fEvent;
//     this.fTrack = other.fTrack;
//     this.fEventCuts = other.fEventCuts;
//     this.fTrackCutsProton = other.fTrackCutsProton;
//     this.fTrackCutsAntiProton = other.fTrackCutsAntiProton;
//     this.fCascadeCutsAntiXion = other.fCascadeCutsAntiXion;
//     this.fConfig = other.fConfig;
//     this.fPairCleaner = other.fPairCleaner;
//     this.fPartColl = other.fPartColl;
//     this.fGTI = other.fGTI;
//     this.fTrackBufferSize = other.fTrackBufferSize;

//     return *this;
// }

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
    fTrackCutsProton->SetName("Protons");
    fOutput->Add(fTrackCutsProton->GetQAHists());
    if (fTrackCutsProton->GetIsMonteCarlo())
    {
        fTrackCutsProton->SetMCName("MCProtonen");
        fOutput->Add(fTrackCutsProton->GetMCQAHists());
    }
    // ##

    // AntiProton Cuts  ###########
    if (!fTrackCutsAntiProton)
    {
        AliFatal("Track Cuts for Particle AntiProton not set!");
    }
    fTrackCutsAntiProton->Init();
    fTrackCutsAntiProton->SetName("AntiProtons");
    fOutput->Add(fTrackCutsAntiProton->GetQAHists());
    if (fTrackCutsAntiProton->GetIsMonteCarlo())
    {
        fTrackCutsAntiProton->SetMCName("MCAntiProtonen");
        fOutput->Add(fTrackCutsAntiProton->GetMCQAHists());
    }
    // ##

    // Xion Cuts    ###########
    if (!fCascadeCutsXion)
    {
        AliFatal("Track Cuts for Particle Xion not set!");
    }
    fCascadeCutsXion->Init();
    fCascadeCutsXion->SetName("Xions");
    fOutput->Add(fCascadeCutsXion->GetQAHists());
    if (fCascadeCutsXion->GetIsMonteCarlo())
    {
        fCascadeCutsXion->SetMCName("MCXion");
        fOutput->Add(fCascadeCutsXion->GetMCQAHists());
    }
    // ##
    
    // AntiXion Cuts    ###########
    if (!fCascadeCutsAntiXion)
    {
        AliFatal("Track Cuts for Particle AntiXion not set!");
    }
    fCascadeCutsAntiXion->Init();
    fCascadeCutsAntiXion->SetName("AntiXions");
    fOutput->Add(fCascadeCutsAntiXion->GetQAHists());
    if (fCascadeCutsAntiXion->GetIsMonteCarlo())
    {
        fCascadeCutsAntiXion->SetMCName("MCAntiXions");
        fOutput->Add(fCascadeCutsAntiXion->GetMCQAHists());
    }
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
    fOutput->Add(fPairCleaner->GetHistList());
    fPartColl = new AliFemtoDreamPartCollection(fConfig, false);
    fOutput->Add(fPartColl->GetHistList());
    fOutput->Add(fPartColl->GetQAList());
    PostData(1, fOutput);
}

static std::vector<AliFemtoDreamBasePart> vProtons;         // Particle Vectors  
static std::vector<AliFemtoDreamBasePart> vAntiProtons;     
static std::vector<AliFemtoDreamBasePart> vXions;           
static std::vector<AliFemtoDreamBasePart> vAntiXions;       

void AliAnalysisTaskPOmegaPenne::UserExec(Option_t *)
{
    Event = dynamic_cast<AliAODEvent *>(fInputEvent);
    
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
                track = dynamic_cast<AliAODTrack *>(Event->GetTrack(iTrack));
                if (!track)
                {
                    AliFatal("No Standard AOD");
                    return;
                }
                StoreGlobalTrackReference(track);
            }
           
            fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

            vProtons.clear();
            vAntiProtons.clear();
            vXions.clear();
            vAntiXions.clear();
            
            for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack)
            {
                track = dynamic_cast<AliAODTrack *>(Event->GetTrack(iTrack));
                if (!track)
                {
                    AliFatal("No Standard AOD");
                    return;
                }
                fTrack->SetTrack(track);

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
            for (int iCasc = 0; iCasc < static_cast<TClonesArray *>(Event->GetCascades())->GetEntriesFast(); ++iCasc)
            {
                AliAODcascade *casc = Event->GetCascade(iCasc);
                fCascade->SetCascade(Event, casc);
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

//  Stores TrackID in Global Track Reference Array 'fGTI' if ID > 0
//
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
