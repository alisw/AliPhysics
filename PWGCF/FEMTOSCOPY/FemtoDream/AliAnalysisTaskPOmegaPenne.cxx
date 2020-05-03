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
                                                                VEvent(0),
                                                                VTrack(0),
                                                                fEvent(0),
                                                                fTrack(0),
                                                                fEventCuts(0),
                                                                fEventCuts2(0),
                                                                fv0(0),
                                                                fv0_2(0),
                                                                fCascade(0),
                                                                fCascade2(0),
                                                                fLambdaV0Cuts(0),
                                                                fLambdaV0Cuts2(0),
                                                                fAntiLambdaV0Cuts(0),
                                                                fAntiLambdaV0Cuts2(0),
                                                                fCascadeCutsXi(0),
                                                                fCascadeCutsXi2(0),
                                                                fCascadeCutsAntiXi(0),
                                                                fCascadeCutsAntiXi2(0),
                                                                fConfig(0),
                                                                fPairCleaner(0),
                                                                fPairCleaner2(0),
                                                                fPartColl(0),
                                                                fPartColl2(0),
                                                                fGTI(0),
                                                                fTrackBufferSize(10000),
                                                                tlEventCuts(0),
                                                                tlLambdaList(0),
                                                                tlAntiLambdaList(0),
                                                                tlCascadeCutsXi(0),
                                                                tlAntiCascadeCutsXi(0),
                                                                tlEventCuts2(0),
                                                                tlLambdaList2(0),
                                                                tlAntiLambdaList2(0),
                                                                tlCascadeCutsXi2(0),
                                                                tlAntiCascadeCutsXi2(0),
                                                                tlResults(0),
                                                                tlResults2(0),
                                                                tlResultsQA(0),
                                                                tlResultsQA2(0),
                                                                tlLambdaMC(0),
                                                                tlAntiLambdaMC(0)
{
}
AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne(const char *name, bool isMC) : AliAnalysisTaskSE(name),
                                                                                      fIsMC(isMC),
                                                                                      VEvent(0),
                                                                                      VTrack(0),
                                                                                      fEvent(0),
                                                                                      fTrack(0),
                                                                                      fEventCuts(0),
                                                                                      fEventCuts2(0),
                                                                                      fv0(0),
                                                                                      fv0_2(0),
                                                                                      fCascade(0),
                                                                                      fCascade2(0),
                                                                                      fLambdaV0Cuts(0),
                                                                                      fAntiLambdaV0Cuts(0),
                                                                                      fCascadeCutsXi(0),
                                                                                      fCascadeCutsAntiXi(0),
                                                                                      fLambdaV0Cuts2(0),
                                                                                      fAntiLambdaV0Cuts2(0),
                                                                                      fCascadeCutsXi2(0),
                                                                                      fCascadeCutsAntiXi2(0),
                                                                                      fConfig(0),
                                                                                      fPairCleaner(0),
                                                                                      fPairCleaner2(0),
                                                                                      fPartColl(0),
                                                                                      fPartColl2(0),
                                                                                      fGTI(0),
                                                                                      fTrackBufferSize(10000),
                                                                                      tlEventCuts(0),
                                                                                      tlLambdaList(0),
                                                                                      tlAntiLambdaList(0),
                                                                                      tlCascadeCutsXi(0),
                                                                                      tlAntiCascadeCutsXi(0),
                                                                                      tlEventCuts2(0),
                                                                                      tlLambdaList2(0),
                                                                                      tlAntiLambdaList2(0),
                                                                                      tlCascadeCutsXi2(0),
                                                                                      tlAntiCascadeCutsXi2(0),
                                                                                      tlResults(0),
                                                                                      tlResults2(0),
                                                                                      tlResultsQA(0),
                                                                                      tlResultsQA2(0),
                                                                                      tlLambdaMC(0),
                                                                                      tlAntiLambdaMC(0)
{
    DefineOutput(1, TList::Class());    // Event Cuts
    DefineOutput(2, TList::Class());    // Lambda Track Cuts
    DefineOutput(3, TList::Class());    // Anti Lambda Track Cuts
    DefineOutput(4, TList::Class());    // Xi Track Cuts
    DefineOutput(5, TList::Class());    // Anti Xi Track Cuts
    DefineOutput(6, TList::Class());    // Results - PairCleaner  - Keep Lambda
    DefineOutput(7, TList::Class());    // QA Results             - Keep Lambda

    DefineOutput(8, TList::Class());    // Event Cuts             - Keep Xi
    DefineOutput(9, TList::Class());    // Lambda Track Cuts      - Keep Xi
    DefineOutput(10, TList::Class());   // Anti Lambda Track Cuts - Keep Xi
    DefineOutput(11, TList::Class());   // Xi Track Cuts          - Keep Xi
    DefineOutput(12, TList::Class());   // Anti Xi Track Cuts     - Keep Xi
    DefineOutput(13, TList::Class());   // Results2 - PairCleaner - Keep Xi
    DefineOutput(14, TList::Class());   // QA Results2            - Keep Xi
    if (isMC)
    {
        DefineOutput(15, TList::Class());    // MC V0 - Lamba
        DefineOutput(16, TList::Class());    // MC AntiV0 - AntiLambda
    }
    
}
AliAnalysisTaskPOmegaPenne::~AliAnalysisTaskPOmegaPenne()       // Destructor
{
    delete VEvent;
    delete VTrack;
    delete fEvent;
    delete fTrack;
    delete fEventCuts;
    delete fEventCuts2;
    delete fv0;
    delete fv0_2;
    delete fCascade;
    delete fCascade2;
    delete fLambdaV0Cuts;
    delete fAntiLambdaV0Cuts;
    delete fCascadeCutsXi;
    delete fCascadeCutsAntiXi;
    delete fLambdaV0Cuts2;
    delete fAntiLambdaV0Cuts2;
    delete fCascadeCutsXi2;
    delete fCascadeCutsAntiXi2;
    delete fConfig;
    delete fPairCleaner;
    delete fPairCleaner2;
    delete fPartColl;
    delete fPairCleaner2;
    delete *fGTI;
    delete tlEventCuts;
    delete tlLambdaList;
    delete tlAntiLambdaList;
    delete tlCascadeCutsXi;
    delete tlAntiCascadeCutsXi;
    delete tlEventCuts2;
    delete tlLambdaList2;
    delete tlAntiLambdaList2;
    delete tlCascadeCutsXi2;
    delete tlAntiCascadeCutsXi2;
    delete tlResults;
    delete tlResults2;
    delete tlResultsQA;
    delete tlResultsQA2;
    delete tlLambdaMC;
    delete tlAntiLambdaMC;
    if (fGTI) delete fGTI;
}

// // Copy Constructor
AliAnalysisTaskPOmegaPenne::AliAnalysisTaskPOmegaPenne(const AliAnalysisTaskPOmegaPenne &obj) : AliAnalysisTaskSE(obj),
                                                                                                fIsMC(obj.fIsMC),
                                                                                                VEvent(obj.VEvent),
                                                                                                VTrack(obj.VTrack),
                                                                                                fEvent(obj.fEvent),
                                                                                                fTrack(obj.fTrack),
                                                                                                fEventCuts(obj.fEventCuts),
                                                                                                fEventCuts2(obj.fEventCuts2),
                                                                                                fv0(obj.fv0),
                                                                                                fv0_2(obj.fv0_2),
                                                                                                fCascade(obj.fCascade),
                                                                                                fCascade2(obj.fCascade2),
                                                                                                fLambdaV0Cuts(obj.fLambdaV0Cuts),
                                                                                                fAntiLambdaV0Cuts(obj.fAntiLambdaV0Cuts),
                                                                                                fCascadeCutsXi(obj.fCascadeCutsXi),
                                                                                                fCascadeCutsAntiXi(obj.fCascadeCutsAntiXi),
                                                                                                fLambdaV0Cuts2(obj.fLambdaV0Cuts2),
                                                                                                fAntiLambdaV0Cuts2(obj.fAntiLambdaV0Cuts2),
                                                                                                fCascadeCutsXi2(obj.fCascadeCutsXi2),
                                                                                                fCascadeCutsAntiXi2(obj.fCascadeCutsAntiXi2),
                                                                                                fConfig(obj.fConfig),
                                                                                                fPairCleaner(obj.fPairCleaner),
                                                                                                fPartColl(obj.fPartColl),
                                                                                                fPartColl2(obj.fPartColl2),
                                                                                                fPairCleaner2(obj.fPairCleaner2),
                                                                                                fGTI(obj.fGTI),
                                                                                                fTrackBufferSize(obj.fTrackBufferSize),
                                                                                                tlEventCuts(obj.tlEventCuts),
                                                                                                tlLambdaList(obj.tlLambdaList),
                                                                                                tlAntiLambdaList(obj.tlAntiLambdaList),
                                                                                                tlCascadeCutsXi(obj.tlCascadeCutsXi),
                                                                                                tlAntiCascadeCutsXi(obj.tlAntiCascadeCutsXi),
                                                                                                tlEventCuts2(obj.tlEventCuts2),
                                                                                                tlLambdaList2(obj.tlLambdaList2),
                                                                                                tlAntiLambdaList2(obj.tlAntiLambdaList2),
                                                                                                tlCascadeCutsXi2(obj.tlCascadeCutsXi2),
                                                                                                tlAntiCascadeCutsXi2(obj.tlAntiCascadeCutsXi2),
                                                                                                tlResults(obj.tlResults),
                                                                                                tlResults2(obj.tlResults2),
                                                                                                tlResultsQA(obj.tlResultsQA),
                                                                                                tlResultsQA2(obj.tlResultsQA2),
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
   
    fEvent = new AliFemtoDreamEvent(true, true, GetCollisionCandidates(), false);
    fTrack = new AliFemtoDreamTrack();
    fTrack->SetUseMCInfo(fIsMC);
    fGTI = new AliVTrack *[fTrackBufferSize];
    
    fEventCuts->InitQA();
 
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

    fEventCuts2->InitQA();
    // ############################################# NUMMER 2 - only Xi left alive ############################
    // Lambda Cutys    ###########
    if (!fLambdaV0Cuts2){AliFatal("Track Cuts for Particle Lambda not set!");}
    fLambdaV0Cuts2->Init();
    fLambdaV0Cuts2->SetName("Lambda");
    // ##

    // AntiLambda Cutys    ###########
    if (!fAntiLambdaV0Cuts2){AliFatal("Track Cuts for Particle AntiLambda not set!");}
    fAntiLambdaV0Cuts2->Init();
    fAntiLambdaV0Cuts2->SetName("AntiLambda");
    // ##

    // V0 Candidates
    fv0_2 = new AliFemtoDreamv0();
    fv0_2->SetUseMCInfo(fIsMC);
    fv0_2->GetPosDaughter()->SetUseMCInfo(fIsMC); 
    fv0_2->GetNegDaughter()->SetUseMCInfo(fIsMC); 
    fv0_2->SetPDGCode(3122);
    fv0_2->SetPDGDaughterPos(2212);
    fv0_2->SetPDGDaughterNeg(211);
    // ##

    // Xi Cuts    ###########
    if (!fCascadeCutsXi2){AliFatal("Track Cuts for Particle Xi not set!");}
    fCascadeCutsXi2->Init();
    fCascadeCutsXi2->SetName("Xi");
    // ##
    
    // AntiXi Cuts    ###########
    if (!fCascadeCutsAntiXi2){AliFatal("Track Cuts for Particle AntiXi not set!");}
    fCascadeCutsAntiXi2->Init();
    fCascadeCutsAntiXi2->SetName("AntiXi");
    // ##

    // Cascade Cuts     #########
    fCascade2 = new AliFemtoDreamCascade();          // Initial Cascade Object
    fCascade2->SetUseMCInfo(fIsMC);
    //PDG Codes should be set assuming Xi- to also work for Xi+
    fCascade2->SetPDGCode(3312);
    fCascade2->SetPDGDaugPos(2212);
    fCascade2->GetPosDaug()->SetUseMCInfo(fIsMC);
    fCascade2->SetPDGDaugNeg(211);
    fCascade2->GetNegDaug()->SetUseMCInfo(fIsMC);
    fCascade2->SetPDGDaugBach(211);
    fCascade2->GetBach()->SetUseMCInfo(fIsMC);
    fCascade2->Setv0PDGCode(3122);
    // ##
    // ############################################# ENDE - NUMMER 2 - only Xi left alive ######################


    fPairCleaner = new AliFemtoDreamPairCleaner(0, 4, false);       // keep Lambdas
    fPairCleaner2 = new AliFemtoDreamPairCleaner(0, 4, false);      // keep Xi
    fPartColl = new AliFemtoDreamPartCollection(fConfig, false);
    fPartColl2 = new AliFemtoDreamPartCollection(fConfig, false);
    
    tlCascadeCutsXi = new TList();
    tlCascadeCutsXi->SetName("XiCascade");
    tlCascadeCutsXi->SetOwner();

    tlAntiCascadeCutsXi = new TList();
    tlAntiCascadeCutsXi->SetName("AntiXiCascade");
    tlAntiCascadeCutsXi->SetOwner();

    tlResultsQA = new TList();
    tlResultsQA->SetName("ResultsQA");
    tlResultsQA->SetOwner();

    tlResultsQA2 = new TList();
    tlResultsQA2->SetName("ResultsQA");
    tlResultsQA2->SetOwner();

    // Connect Cuts to OutputContainers
    tlEventCuts             = fEventCuts->GetHistList();
    tlLambdaList            = fLambdaV0Cuts->GetQAHists();
    tlAntiLambdaList        = fAntiLambdaV0Cuts->GetQAHists();
    tlCascadeCutsXi         = fCascadeCutsXi->GetQAHists();
    tlAntiCascadeCutsXi     = fCascadeCutsAntiXi->GetQAHists();
    tlResults               = fPartColl->GetHistList();
    tlResultsQA->Add(         fPartColl->GetQAList());
    tlResultsQA->Add(         fPairCleaner->GetHistList());
    tlResultsQA->Add(         fEvent->GetEvtCutList());
    
    tlEventCuts2             = fEventCuts2->GetHistList();
    tlLambdaList2            = fLambdaV0Cuts2->GetQAHists();
    tlAntiLambdaList2        = fAntiLambdaV0Cuts2->GetQAHists();
    tlCascadeCutsXi2         = fCascadeCutsXi2->GetQAHists();
    tlAntiCascadeCutsXi2     = fCascadeCutsAntiXi2->GetQAHists();
    tlResults2              = fPartColl2->GetHistList();
    tlResultsQA2->Add(        fPartColl2->GetQAList());
    tlResultsQA2->Add(        fPairCleaner2->GetHistList());
    tlResultsQA2->Add(        fEvent->GetEvtCutList());

    PostData(1, tlEventCuts);
    PostData(2, tlLambdaList);
    PostData(3, tlAntiLambdaList);
    PostData(4, tlCascadeCutsXi);
    PostData(5, tlAntiCascadeCutsXi);
    PostData(6, tlResults);
    PostData(7, tlResultsQA);

    PostData(8, tlEventCuts2);
    PostData(9, tlLambdaList2);
    PostData(10, tlAntiLambdaList2);
    PostData(11, tlCascadeCutsXi2);
    PostData(12, tlAntiCascadeCutsXi2);
    PostData(13, tlResults2);
    PostData(14, tlResultsQA2);

    if (fLambdaV0Cuts->GetIsMonteCarlo())
    {
        tlLambdaMC = fLambdaV0Cuts->GetMCQAHists();
        PostData(15, tlLambdaMC);
    }
    if (fAntiLambdaV0Cuts->GetIsMonteCarlo())
    {
        tlAntiLambdaMC = fAntiLambdaV0Cuts->GetMCQAHists();
        PostData(16, tlAntiLambdaMC);
    }
}

static std::vector<AliFemtoDreamBasePart> vLambda;           // keep Lambda after PairCleaner
static std::vector<AliFemtoDreamBasePart> vAntiLambda;       
static std::vector<AliFemtoDreamBasePart> vXi;           
static std::vector<AliFemtoDreamBasePart> vAntiXi;       

static std::vector<AliFemtoDreamBasePart> vLambda2;             // keep Xi after PairCleaner
static std::vector<AliFemtoDreamBasePart> vAntiLambda2;         
static std::vector<AliFemtoDreamBasePart> vXi2;                 
static std::vector<AliFemtoDreamBasePart> vAntiXi2;             

void AliAnalysisTaskPOmegaPenne::UserExec(Option_t *)
{
    VEvent = static_cast<AliVEvent *>(fInputEvent);
    
    if (!VEvent)
    {
        AliWarning("No Input VEvent");
    }
    else
    {
        fEvent->SetEvent(VEvent);
        if (fEventCuts->isSelected(fEvent))
        {
            ResetGlobalTrackReference();
            for (int iTrack = 0; iTrack < VEvent->GetNumberOfTracks(); ++iTrack)
            {
                VTrack = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
                if (!VTrack)
                {
                    AliFatal("No Standard AOD");
                    return;
                }
                StoreGlobalTrackReference(VTrack);
            }
           
            // fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);      // für protonen

            vXi.clear();
            vAntiXi.clear();
            vLambda.clear();
            vAntiLambda.clear();
        
            vXi2.clear();
            vAntiXi2.clear();
            vLambda2.clear();
            vAntiLambda2.clear();


            // for (int iTrack = 0; iTrack < VEvent->GetNumberOfTracks(); ++iTrack)
            // {
            //     aaTrack = dynamic_cast<AliAODTrack *>(VEvent->GetTrack(iTrack));
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

            // irgendwie benötigt um GetV0s() und GetCascade() zu holen
            AliAODEvent* aodEvent = static_cast<AliAODEvent*>(VEvent); // caste input event auf ein AODEvent 
            
            // ## Lambda Selection ## keep Lambdas
            fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
            fv0_2->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
            
            for (int iv0 = 0; iv0 < static_cast<TClonesArray *>(aodEvent->GetV0s())->GetEntriesFast(); ++iv0)
            {
                AliAODv0 *v0 = aodEvent->GetV0(iv0);
                fv0->Setv0(VEvent, v0, fEvent->GetMultiplicity());

                // ## Lambda Selection 1 ## keep Lambda
                if (fLambdaV0Cuts->isSelected(fv0)) 
                {
                    vLambda.push_back(*fv0);
                    vLambda[vLambda.size() - 1].SetCPA(1.0);
                }
                if (fAntiLambdaV0Cuts->isSelected(fv0)) 
                {
                    vAntiLambda.push_back(*fv0);
                    vAntiLambda[vAntiLambda.size() - 1].SetCPA(1.0);
                }

                fv0_2->Setv0(VEvent, v0, fEvent->GetMultiplicity());
                // ## Lambda Selection 2 ## keep Xi            
                if (fLambdaV0Cuts2->isSelected(fv0_2)) 
                {
                    vLambda2.push_back(*fv0_2);
                    vLambda2[vLambda2.size() - 1].SetCPA(0.5);
                }
                if (fAntiLambdaV0Cuts2->isSelected(fv0_2)) 
                {
                    vAntiLambda2.push_back(*fv0_2);
                    vAntiLambda2[vAntiLambda2.size() - 1].SetCPA(0.5);
                }
            }

            // ## Xi selection
            for (int iCasc = 0; iCasc < static_cast<TClonesArray *>(aodEvent->GetCascades())->GetEntriesFast(); ++iCasc)
            {
                AliAODcascade *casc = aodEvent->GetCascade(iCasc);
                fCascade->SetCascade(VEvent, casc);

                // ## Xi selection 1 ### keep Lambda
                if (fCascadeCutsXi->isSelected(fCascade))
                {
                    vXi.push_back(*fCascade);
                    vXi[vXi.size() - 1].SetCPA(0.5);
                }
                if (fCascadeCutsAntiXi->isSelected(fCascade))
                {
                    vAntiXi.push_back(*fCascade);
                    vAntiXi[vAntiXi.size() -1].SetCPA(0.5);
                }

                fCascade2->SetCascade(VEvent, casc);
                // ## Xi selection 2 ### keep Xi
                if (fCascadeCutsXi2->isSelected(fCascade2))
                {
                    vXi2.push_back(*fCascade2);
                    vXi2[vXi2.size() - 1].SetCPA(1.0);
                }
                if (fCascadeCutsAntiXi2->isSelected(fCascade2))
                {
                    vAntiXi2.push_back(*fCascade2);
                    vAntiXi2[vAntiXi2.size() -1].SetCPA(1.0);
                }
            }

            // remove double-matched tracks
            fPairCleaner->ResetArray();
            fPairCleaner2->ResetArray();
            
            fPairCleaner->CleanDecayAndDecay(&vXi, &vLambda,  0);
            fPairCleaner->CleanDecayAndDecay(&vAntiXi, &vAntiLambda, 1);
            fPairCleaner->CleanDecay(&vLambda, 2);
            fPairCleaner->CleanDecay(&vAntiLambda, 3);
            
            fPairCleaner2->CleanDecayAndDecay(&vXi2, &vLambda2,  0);
            fPairCleaner2->CleanDecayAndDecay(&vAntiXi2, &vAntiLambda2, 1);
            fPairCleaner2->CleanDecay(&vLambda2, 2);
            fPairCleaner2->CleanDecay(&vAntiLambda2, 3);

            // fPairCleaner->CleanDecay(&vXi, 0);
            // fPairCleaner->CleanDecay(&vAntiXi, 1);

            fPairCleaner->StoreParticle(vLambda); 
            fPairCleaner->StoreParticle(vAntiLambda);

            fPairCleaner->StoreParticle(vXi);
            fPairCleaner->StoreParticle(vAntiXi);

            fPairCleaner2->StoreParticle(vLambda2); 
            fPairCleaner2->StoreParticle(vAntiLambda2);

            fPairCleaner2->StoreParticle(vXi2);
            fPairCleaner2->StoreParticle(vAntiXi2);

            fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality()); // proton xi and lambda analysis
            fPartColl2->SetEvent(fPairCleaner2->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality()); // proton xi and lambda analysis
            // soweit ich das richtig verstanden habe wird pairQA mit den teilchen gemacht die im pairCleaner 
            // sind und pdgCodes in der richtigen Reihenfolge vorhanden sind.


            PostData(1, tlEventCuts);
            PostData(2, tlLambdaList);
            PostData(3, tlAntiLambdaList);
            PostData(4, tlCascadeCutsXi);
            PostData(5, tlAntiCascadeCutsXi);
            PostData(6, tlResults);
            PostData(7, tlResultsQA);

            PostData(8, tlEventCuts2);
            PostData(9, tlLambdaList2);
            PostData(10, tlAntiLambdaList2);
            PostData(11, tlCascadeCutsXi2);
            PostData(12, tlAntiCascadeCutsXi2);
            PostData(13, tlResults2);
            PostData(14, tlResultsQA2);
            if (fIsMC)
            {
                PostData(15, tlLambdaMC);
            
                PostData(16, tlAntiLambdaMC);
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
void AliAnalysisTaskPOmegaPenne::StoreGlobalTrackReference(AliVTrack *vTrack)
{
    //This method was inherited form H. Beck analysis
    AliNanoAODTrack *nanoTrack = static_cast<AliNanoAODTrack*>(vTrack);
    const int trackID = vTrack->GetID();
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
        if ((!nanoTrack->GetFilterMap()) && (!vTrack->GetTPCNcls()))
        {
            return;
        }
        if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap() || fGTI[trackID]->GetTPCNcls())
        {
            printf("WARNING! global track info already there!");
            printf("    ###     TPCNcls track1 %u Track2 %u", (fGTI[trackID])->GetTPCNcls(), vTrack->GetTPCNcls());
            printf("   ###     FilterMap Track1 %u track2 %u\n", dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap(), nanoTrack->GetFilterMap());
        }
    }
    fGTI[trackID] = vTrack;
}
