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
                                                                fEventCuts(0),
                                                                fTrackCutsProton(0),
                                                                fTrackCutsAntiProton(0),
                                                                fv0(0),
                                                                fLambdaV0Cuts(0),
                                                                fAntiLambdaV0Cuts(0),
                                                                // fCascade(0),
                                                                // fCascadeCutsXion(0),
                                                                // fCascadeCutsAntiXion(0),
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
                                                                // tlCascadeCutsXi(0),
                                                                // tlAntiCascadeCutsXi(0),
                                                                tlPairCleaner(0),
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
                                                                                      fEventCuts(0),
                                                                                      fTrackCutsProton(0),
                                                                                      fTrackCutsAntiProton(0),
                                                                                      fv0(0),
                                                                                      fLambdaV0Cuts(0),
                                                                                      fAntiLambdaV0Cuts(0),
                                                                                      // fCascade(0),
                                                                                      // fCascadeCutsXion(0),
                                                                                      // fCascadeCutsAntiXion(0),
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
                                                                                      // tlCascadeCutsXi(0),
                                                                                      // tlAntiCascadeCutsXi(0),
                                                                                      tlPairCleaner(0),
                                                                                      tlResults(0),
                                                                                      tlResultsQA(0)
{
    DefineOutput(1, TList::Class());    // Event Cuts
    DefineOutput(2, TList::Class());    // Proton Track Cuts
    DefineOutput(3, TList::Class());    // Anti Proton Track Cuts
    DefineOutput(4, TList::Class());    // Lambda Track Cuts
    DefineOutput(5, TList::Class());    // Anti Lambda Track Cuts
    DefineOutput(6, TList::Class());    // Pair Cleaner
    DefineOutput(7, TList::Class());    // Results
    DefineOutput(8, TList::Class());    // QA Results
}
AliAnalysisTaskPOmegaPenne::~AliAnalysisTaskPOmegaPenne()
{
    // TODO Auto-generated destructor stub
    // |-> HÄ, ne gar nicht!... wenn der run zu ende ist hört das Objekt einfach auf zu eXistieren. Scheiß auf den Destruktor!
    // pass lieber in der UserExec Methode auf!!!
}

// // Copy Constructor
// AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne(const AliAnalysisTaskPOmegaPenne& obj) : AliAnalysisTaskSE(obj),
//                                                                                                 fIsMC(obj.fIsMC),
//                                                                                                 aaEvent(obj.aaEvent),
//                                                                                                 aaTrack(obj.aaTrack),
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
//     AliAnalysisTaskSE::operator=(other); // hier ist doof glaub ich und in der header datei muss das auch noch hinzu
//     this.fIsMC = other.fIsMC;
//     this.aaEvent = other.aaEvent;
//     this.aaTrack = other.aaTrack;
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
    // tlEventCuts = new TList();
    // tlEventCuts->SetName("EventCuts");
    // tlEventCuts->SetOwner();

    // tlTrackCutsProton = new TList();
    // tlTrackCutsProton->SetName("TrackCuts");
    // tlTrackCutsProton->SetOwner();

    // tlAntiTrackCutsProton = new TList();
    // tlAntiTrackCutsProton->SetName("AntiTrackCuts");
    // tlAntiTrackCutsProton->SetOwner();

    // tlCascadeCutsXi = new TList();
    // tlCascadeCutsXi->SetName("CascadeCuts");
    // tlCascadeCutsXi->SetOwner();

    // tlAntiCascadeCutsXi = new TList();
    // tlAntiCascadeCutsXi->SetName("AntiCascadeCuts");
    // tlAntiCascadeCutsXi->SetOwner();
                                                                                      
    tlPairCleaner = new TList();
    tlPairCleaner->SetName("PairCleaner");
    tlPairCleaner->SetOwner();

    // tlResults = new TList();
    // tlResults->SetName("Results");
    // tlResults->SetOwner();

    tlResultsQA = new TList();
    tlResultsQA->SetName("ResultsQA");
    tlResultsQA->SetOwner();


    fEvent = new AliFemtoDreamEvent(false, true, GetCollisionCandidates());
    fTrack = new AliFemtoDreamTrack();
    fTrack->SetUseMCInfo(fIsMC);
    fGTI = new AliAODTrack *[fTrackBufferSize];
    
    fEventCuts->InitQA();
    
    // Proton Cuts      ###########
    if (!fTrackCutsProton)
    {
        AliFatal("Track Cuts for Particle Proton not set!");
    }
    fTrackCutsProton->Init();
    fTrackCutsProton->SetName("Protons");
    // if (fTrackCutsProton->GetIsMonteCarlo())             // is bei mir eh kein monte carlo im moment
    // {
    //     fTrackCutsProton->SetMCName("MCProtonen");
    //     tlTrackCutsProton->Add(fTrackCutsProton->GetMCQAHists());
    // }
    // ##

    // AntiProton Cuts  ###########
    if (!fTrackCutsAntiProton)
    {
        AliFatal("Track Cuts for Particle AntiProton not set!");
    }
    fTrackCutsAntiProton->Init();
    // fTrackCutsAntiProton->SetName("AntiProtons");
    // if (fTrackCutsAntiProton->GetIsMonteCarlo())         // is bei mir eh kein monte carlo im moment
    // {
    //     fTrackCutsAntiProton->SetMCName("MCAntiProtonen");
    //     tlAntiTrackCutsProton->Add(fTrackCutsAntiProton->GetMCQAHists());
    // }
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
  
    // // Xion Cuts    ###########
    // if (!fCascadeCutsXion)
    // {
    //     AliFatal("Track Cuts for Particle Xion not set!");
    // }
    // fCascadeCutsXion->Init();
    // fCascadeCutsXion->SetName("Xis");
    // tlCascadeCutsXi->Add(fCascadeCutsXion->GetQAHists());
    // if (fCascadeCutsXion->GetIsMonteCarlo())
    // {
    //     fCascadeCutsXion->SetMCName("MCXion");
    //     tlCascadeCutsXi->Add(fCascadeCutsXion->GetMCQAHists());
    // }
    // // ##
    
    // // AntiXion Cuts    ###########
    // if (!fCascadeCutsAntiXion)
    // {
    //     AliFatal("Track Cuts for Particle AntiXion not set!");
    // }
    // fCascadeCutsAntiXion->Init();
    // fCascadeCutsAntiXion->SetName("AntiXis");
    // tlAntiCascadeCutsXi->Add(fCascadeCutsAntiXion->GetQAHists());
    // if (fCascadeCutsAntiXion->GetIsMonteCarlo())
    // {
    //     fCascadeCutsAntiXion->SetMCName("MCAntiXions");
    //     tlAntiCascadeCutsXi->Add(fCascadeCutsAntiXion->GetMCQAHists());
    // }
    // // ##

    // // Cascade Cuts     #########
    // fCascade = new AliFemtoDreamCascade();          // Initial Cascade Object
    // fCascade->SetUseMCInfo(fCascadeCutsXion->GetIsMonteCarlo() || fCascadeCutsAntiXion->GetIsMonteCarlo());
    // //PDG Codes should be set assuming Xi- to also work for Xi+
    // fCascade->SetPDGCode(3312);
    // fCascade->SetPDGDaugPos(2212);
    // fCascade->GetPosDaug()->SetUseMCInfo(fCascadeCutsXion->GetIsMonteCarlo() || fCascadeCutsAntiXion->GetIsMonteCarlo());
    // fCascade->SetPDGDaugNeg(211);
    // fCascade->GetNegDaug()->SetUseMCInfo(fCascadeCutsXion->GetIsMonteCarlo() || fCascadeCutsAntiXion->GetIsMonteCarlo());
    // fCascade->SetPDGDaugBach(211);
    // fCascade->GetBach()->SetUseMCInfo(fCascadeCutsXion->GetIsMonteCarlo() || fCascadeCutsAntiXion->GetIsMonteCarlo());
    // fCascade->Setv0PDGCode(3122);
    // // ##

    fPairCleaner = new AliFemtoDreamPairCleaner(2, 2, false);
    fPartColl = new AliFemtoDreamPartCollection(fConfig, false);

    tlEventCuts             = fEventCuts->GetHistList();
    tlTrackCutsProton       = fTrackCutsProton->GetQAHists();
    tlAntiTrackCutsProton   = fTrackCutsAntiProton->GetQAHists();
    tlLambdaList            = fLambdaV0Cuts->GetQAHists();
    tlAntiLambdaList        = fAntiLambdaV0Cuts->GetQAHists();
    tlResults               = fPartColl->GetHistList();
    tlResultsQA->Add(fPartColl->GetQAList());
    tlResultsQA->Add(fPairCleaner->GetHistList());
    tlResultsQA->Add(fEvent->GetEvtCutList());

    PostData(1, tlEventCuts);
    PostData(2, tlTrackCutsProton);
    PostData(3, tlAntiTrackCutsProton);
    PostData(4, tlLambdaList);
    PostData(5, tlAntiLambdaList);
    // PostData(4, tlCascadeCutsXi);
    // PostData(5, tlAntiCascadeCutsXi);
    PostData(6, tlPairCleaner);
    PostData(7, tlResults);
    PostData(8, tlResultsQA);
}

static std::vector<AliFemtoDreamBasePart> vProtons;         // Particle Vectors  
static std::vector<AliFemtoDreamBasePart> vAntiProtons;     
static std::vector<AliFemtoDreamBasePart> vLambdas;           
static std::vector<AliFemtoDreamBasePart> vAntiLambdas;       
// static std::vector<AliFemtoDreamBasePart> vXions;           
// static std::vector<AliFemtoDreamBasePart> vAntiXions;       

void AliAnalysisTaskPOmegaPenne::UserExec(Option_t *)
{
    aaEvent = dynamic_cast<AliAODEvent *>(fInputEvent);
    
    if (!aaEvent)
    {
        AliWarning("No Input Event");
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
            vLambdas.clear();
            vAntiLambdas.clear();            
            // vXions.clear();
            // vAntiXions.clear();
            

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
            fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
            for (int iv0 = 0; iv0 < static_cast<TClonesArray *>(aaEvent->GetV0s())->GetEntriesFast(); ++iv0)
            {
                AliAODv0 *v0 = aaEvent->GetV0(iv0);
                fv0->Setv0(aaEvent, v0, fEvent->GetMultiplicity());
                if (fLambdaV0Cuts->isSelected(fv0))
                {
                    vLambdas.push_back(*fv0);
                }
                if (fAntiLambdaV0Cuts->isSelected(fv0))
                {
                    vAntiLambdas.push_back(*fv0);
                }
            }
            // @@ Xi Cascade ##
            // for (int iCasc = 0; iCasc < static_cast<TClonesArray *>(aaEvent->GetCascades())->GetEntriesFast(); ++iCasc)
            // {
            //     AliAODcascade *casc = aaEvent->GetCascade(iCasc);
            //     fCascade->SetCascade(aaEvent, casc);
            //     if (fCascadeCutsXion->isSelected(fCascade))
            //     {
            //         vXions.push_back(*fCascade);
            //     }
            //     if (fCascadeCutsAntiXion->isSelected(fCascade))
            //     {
            //         vAntiXions.push_back(*fCascade);
            //     }
            // }                    
                                                                 
            // remove double-matched tracks
            fPairCleaner->ResetArray();
            fPairCleaner->CleanTrackAndDecay(&vProtons, &vLambdas, 0);
            fPairCleaner->CleanTrackAndDecay(&vAntiProtons, &vAntiLambdas, 1);
            // fPairCleaner->CleanTrackAndDecay(&vProtons, &vXions, 0);
            // fPairCleaner->CleanTrackAndDecay(&vAntiProtons, &vAntiXions, 1);

            fPairCleaner->CleanDecay(&vLambdas, 0);
            fPairCleaner->CleanDecay(&vAntiLambdas, 1);

            // fPairCleaner->CleanDecay(&vXions, 0);
            // fPairCleaner->CleanDecay(&vAntiXions, 1);
            
            fPairCleaner->StoreParticle(vProtons);
            fPairCleaner->StoreParticle(vAntiProtons);
            fPairCleaner->StoreParticle(vLambdas);
            fPairCleaner->StoreParticle(vAntiLambdas);
            if (fPairCleaner->GetCounter() > 0)
            {
                if (fConfig->GetUseEventMixing())
                {
                    fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());
                }
            }
            // fPairCleaner->StoreParticle(vXions);
            // fPairCleaner->StoreParticle(vAntiXions);

            // fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality());        // xi's mit getRefMult08() anstelle von GetMultiplicity()

            PostData(1, tlEventCuts);
            PostData(2, tlTrackCutsProton);
            PostData(3, tlAntiTrackCutsProton);
            PostData(4, tlLambdaList);
            PostData(5, tlAntiLambdaList);
            // PostData(4, tlCascadeCutsXi);
            // PostData(5, tlAntiCascadeCutsXi);
            PostData(6, tlPairCleaner);
            PostData(7, tlResults);
            PostData(8, tlResultsQA);
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
