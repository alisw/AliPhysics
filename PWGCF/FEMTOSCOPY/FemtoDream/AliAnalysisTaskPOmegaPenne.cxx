// /*
//  * AliAnalysisTaskPOmegaPenne.cxx
//  *
//  *  Created on: 11 Dec 2019
//  *      Author: Boris Bajtl
//  */



#include "AliAnalysisTaskPOmegaPenne.h"
#include <string.h>
#include "AliNanoAODTrack.h"
#include "TDatabasePDG.h"

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
                                                                tlAntiLambdaMC(0),
                                                                tlRecombination(0),
                                                                hInvMassLambda_total(0),
                                                                hInvMassLambda_shared_pion(0),
                                                                hInvMassLambda_shared_proton(0),
                                                                hInvMassXi_total(0),
                                                                hInvMassXi_shared_bach(0),
                                                                hInvMassXi_shared_pi_daugh(0),
                                                                hInvMassXi_shared_prot_daugh(0),
                                                                hInvMassXi_shared_Lambda(0),
                                                                hInvMassLambda_sanityCheck(0),
                                                                hInvMassXi_sanityCheck(0),
                                                                fEvtCounter(0)
                                                                // fLambdaV0Cuts_rec(0),
                                                                // fAntiLambdaV0Cuts_rec(0)
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
                                                                                      tlRecombination(0),
                                                                                      tlLambdaMC(0),
                                                                                      tlAntiLambdaMC(0),
                                                                                      hInvMassLambda_total(0),
                                                                                      hInvMassLambda_shared_pion(0),
                                                                                      hInvMassLambda_shared_proton(0),
                                                                                      hInvMassXi_total(0),
                                                                                      hInvMassXi_shared_bach(0),
                                                                                      hInvMassXi_shared_pi_daugh(0),
                                                                                      hInvMassXi_shared_prot_daugh(0),
                                                                                      hInvMassXi_shared_Lambda(0),
                                                                                      hInvMassLambda_sanityCheck(0),
                                                                                      hInvMassXi_sanityCheck(0),
                                                                                      fEvtCounter(0)
                                                                                    //   fLambdaV0Cuts_rec(0),
                                                                                    //   fAntiLambdaV0Cuts_rec(0)
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
    DefineOutput(15, TList::Class());       // reconstruction from daugthers histograms

    if (isMC)
    {
        DefineOutput(16, TList::Class());    // MC V0 - Lamba
        DefineOutput(17, TList::Class());    // MC AntiV0 - AntiLambda
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
    delete tlRecombination;
    delete hInvMassLambda_total;
    delete hInvMassLambda_shared_pion;
    delete hInvMassLambda_shared_proton;
    delete hInvMassXi_total;
    delete hInvMassXi_shared_bach;
    delete hInvMassXi_shared_pi_daugh;
    delete hInvMassXi_shared_prot_daugh;
    delete hInvMassXi_shared_Lambda;
    delete hInvMassLambda_sanityCheck;
    delete hInvMassXi_sanityCheck;
                                                                
    delete fEvtCounter;
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
                                                                                                tlAntiLambdaMC(obj.tlAntiLambdaMC),
                                                                                                tlRecombination(obj.tlRecombination),
                                                                                                hInvMassLambda_total(obj.hInvMassLambda_total),
                                                                                                hInvMassLambda_shared_pion(obj.hInvMassLambda_shared_pion),
                                                                                                hInvMassLambda_shared_proton(obj.hInvMassLambda_shared_proton),
                                                                                                hInvMassXi_total(obj.hInvMassXi_total),
                                                                                                hInvMassXi_shared_bach(obj.hInvMassXi_shared_bach),
                                                                                                hInvMassXi_shared_pi_daugh(obj.hInvMassXi_shared_pi_daugh),
                                                                                                hInvMassXi_shared_prot_daugh(obj.hInvMassXi_shared_prot_daugh),
                                                                                                hInvMassXi_shared_Lambda(obj.hInvMassXi_shared_Lambda),
                                                                                                hInvMassLambda_sanityCheck(obj.hInvMassLambda_sanityCheck),
                                                                                                hInvMassXi_sanityCheck(obj.hInvMassXi_sanityCheck),
                                                                                                fEvtCounter(obj.fEvtCounter)
                                                                                                // fLambdaV0Cuts_rec(obj.fLambdaV0Cuts_rec),
                                                                                                // fAntiLambdaV0Cuts_rec(obj.fAntiLambdaV0Cuts_rec)

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
    fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());


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


    // ############################################# Recombination Cuts ######################
    // Lambda Cutys    ###########
    // if (!fLambdaV0Cuts_rec){AliFatal("Track Cuts for Particle Lambda_recombination not set!");}
    // fLambdaV0Cuts_rec->Init();
    // fLambdaV0Cuts_rec->SetName("Lambda_rec");
    // ##
    // AntiLambda Cutys    ###########
    // if (!fAntiLambdaV0Cuts_rec){AliFatal("Track Cuts for Particle AntiLambda_recombination not set!");}
    // fAntiLambdaV0Cuts_rec->Init();
    // fAntiLambdaV0Cuts_rec->SetName("AntiLambda_rec");
    // ##
    // ############################################# ENDE - Recombination Cuts ######################


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


    tlRecombination = new TList();        // Lambda and Xi recombination statistic histogramms for interchanged daughters
    tlRecombination->SetName("Recombination");
    tlRecombination->SetOwner();

    // hInvMassLambda_total = new TH1F("InvariantMassLambda", "Invariant Mass LAMBDA", 400, 1.112, 1.120); 
    hInvMassLambda_sanityCheck = new TH1F("InvariantMassLambdaSanityCheck", "Invariant Mass LAMBDA Sanity Check", 500, 1.04, 2.60);           // mit meiner funktion ausgerechnete invariante masse aus den selektierten Teilchen
    hInvMassLambda_total = new TH1F("InvariantMassLambdatotal", "Invariant Mass LAMBDA total", 500, 1.0, 3.0); 
    hInvMassLambda_shared_pion = new TH1F("InvariantMassLambdaSharedPion", "Invariant Mass LAMBDA shared Pion", 500, 1.0, 3.0);
    hInvMassLambda_shared_proton = new TH1F("InvariantMassLambdaSharedProton", "Invariant Mass LAMBDA shared Proton", 500, 1.0, 3.0);
    hInvMassXi_sanityCheck = new TH1F("InvariantMassXiSanityCheck", "Invariant Mass XI Sanity Check", 600, 1.200, 1.600);                   // mit meiner funktion ausgerechnete invariante masse aus den selektierten Teilchen
    hInvMassXi_total = new TH1F("InvariantMassXiTotal", "Invariant Mass XI total", 600, 0.700, 2.500);                                      
    hInvMassXi_shared_bach = new TH1F("InvariantMassXiSharedBach", "Invariant Mass XI shared Bachelor Pi", 600, 0.700, 2.500); 
    hInvMassXi_shared_pi_daugh = new TH1F("InvariantMassXiSharedPiDaugh", "Invariant Mass XI shared Pi Daugh", 600, 0.700, 2.500); 
    hInvMassXi_shared_prot_daugh = new TH1F("InvariantMassXiSharedProtDaugh", "Invariant Mass XI shared Prot Daugh", 600, 0.700, 2.500); 
    hInvMassXi_shared_Lambda = new TH1F("InvariantMassXiSharedLambda", "Invariant Mass XI shared Lambda", 600, 0.700, 2.500); 

    fEvtCounter = new TH1F("EventCounter", "Event Counter", 7, 0, 7);
    fEvtCounter->GetXaxis()->SetBinLabel(1, "Prot_Lambda + pi_Xi1");        // reconstruct Lambda
    fEvtCounter->GetXaxis()->SetBinLabel(2, "Prot_Lambda + pi_Xi2");        //
    fEvtCounter->GetXaxis()->SetBinLabel(3, "Prot_Xi + pi_Lambda");         //
    fEvtCounter->GetXaxis()->SetBinLabel(4, "Lambda + pi_Xi1");             // reconstruct Xi
    fEvtCounter->GetXaxis()->SetBinLabel(5, "Lambda + pi_Xi2");             //
    fEvtCounter->GetXaxis()->SetBinLabel(6, "Lambda + pi_Lambda");          //
    fEvtCounter->GetXaxis()->SetBinLabel(7, "prot_Lambda + pi_Lambda");     // reconstruct Lambda from other Lambda

    tlRecombination->Add(hInvMassLambda_sanityCheck);
    tlRecombination->Add(hInvMassLambda_total);
    tlRecombination->Add(hInvMassLambda_shared_pion);
    tlRecombination->Add(hInvMassLambda_shared_proton);
    tlRecombination->Add(hInvMassXi_sanityCheck);
    tlRecombination->Add(hInvMassXi_total);
    tlRecombination->Add(hInvMassXi_shared_bach);
    tlRecombination->Add(hInvMassXi_shared_pi_daugh);
    tlRecombination->Add(hInvMassXi_shared_prot_daugh);
    tlRecombination->Add(hInvMassXi_shared_Lambda);
    tlRecombination->Add(fEvtCounter);

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
    tlResults2               = fPartColl2->GetHistList();
    tlResultsQA2->Add(        fPartColl2->GetQAList());
    tlResultsQA2->Add(        fPairCleaner2->GetHistList());
    tlResultsQA2->Add(        fEvent->GetEvtCutList());

    PostData(1, tlEventCuts);           // cuts keeping Lambda
    PostData(2, tlLambdaList);
    PostData(3, tlAntiLambdaList);
    PostData(4, tlCascadeCutsXi);
    PostData(5, tlAntiCascadeCutsXi);
    PostData(6, tlResults);
    PostData(7, tlResultsQA);

    PostData(8, tlEventCuts2);          //  cuts keeping Xi
    PostData(9, tlLambdaList2);
    PostData(10, tlAntiLambdaList2);
    PostData(11, tlCascadeCutsXi2);
    PostData(12, tlAntiCascadeCutsXi2);
    PostData(13, tlResults2);
    PostData(14, tlResultsQA2);

    PostData(15, tlRecombination);         // reconstruction from daugthers histograms


    if (fLambdaV0Cuts->GetIsMonteCarlo())
    {
        tlLambdaMC = fLambdaV0Cuts->GetMCQAHists();
        PostData(16, tlLambdaMC);
    }
    if (fAntiLambdaV0Cuts->GetIsMonteCarlo())
    {
        tlAntiLambdaMC = fAntiLambdaV0Cuts->GetMCQAHists();
        PostData(17, tlAntiLambdaMC);
    }
}         


void AliAnalysisTaskPOmegaPenne::UserExec(Option_t *)
{
    // VEvent = dynamic_cast<AliVEvent *>(fInputEvent);
    VEvent = fInputEvent;
    
    if (!fInputEvent)
    {
        AliWarning("No Input VEvent");
        return;
    }

    fEvent->SetEvent(fInputEvent);
    if (fEventCuts->isSelected(fEvent))
    {
        ResetGlobalTrackReference();
        for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)
        {
            VTrack = dynamic_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
            if (!VTrack)
            {
                AliFatal("No Standard AOD");
                return;
            }
            StoreGlobalTrackReference(VTrack);
        }
        
        int counter = 0;

        std::vector<AliFemtoDreamBasePart> vLambda;         // keep Xi after OPairCleaner
        std::vector<AliFemtoDreamBasePart> vAntiLambda;
        std::vector<AliFemtoDreamBasePart> vXi;
        std::vector<AliFemtoDreamBasePart> vAntiXi;
        std::vector<AliFemtoDreamBasePart> vLambda2;        // keep Lambda after OPairCleaner
        std::vector<AliFemtoDreamBasePart> vAntiLambda2;
        std::vector<AliFemtoDreamBasePart> vXi2;
        std::vector<AliFemtoDreamBasePart> vAntiXi2;

        // irgendwie benötigt um GetV0s() und GetCascade() zu holen
        AliAODEvent *aodEvent = dynamic_cast<AliAODEvent *>(fInputEvent); // caste input event auf ein AODEvent

        // ## Lambda Selection ## keep Lambdas
        fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
        fv0_2->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

        for (int iv0 = 0; iv0 < dynamic_cast<TClonesArray *>(aodEvent->GetV0s())->GetEntriesFast(); ++iv0)
        {
            AliAODv0 *v0 = aodEvent->GetV0(iv0);
            fv0->Setv0(fInputEvent, v0, fEvent->GetMultiplicity());

            // ## Lambda Selection 1 ## keep Xi
            if (fLambdaV0Cuts->isSelected(fv0))
            {
                vLambda.push_back(*fv0);
                vLambda[vLambda.size() - 1].SetCPA(0.5);
            }
            if (fAntiLambdaV0Cuts->isSelected(fv0))
            {
                vAntiLambda.push_back(*fv0);
                vAntiLambda[vAntiLambda.size() - 1].SetCPA(0.5);
            }

            fv0_2->Setv0(fInputEvent, v0, fEvent->GetMultiplicity());
            // ## Lambda Selection 2 ## keep Lambda
            if (fLambdaV0Cuts2->isSelected(fv0_2))
            {
                vLambda2.push_back(*fv0_2);
                vLambda2[vLambda2.size() - 1].SetCPA(1.0);
            }
            if (fAntiLambdaV0Cuts2->isSelected(fv0_2))
            {
                vAntiLambda2.push_back(*fv0_2);
                vAntiLambda2[vAntiLambda2.size() - 1].SetCPA(1.0);
            }
        }
        // ## Xi selection
        for (int iCasc = 0; iCasc < dynamic_cast<TClonesArray *>(aodEvent->GetCascades())->GetEntriesFast(); ++iCasc)
        {
            AliAODcascade *casc = aodEvent->GetCascade(iCasc);
            fCascade->SetCascade(fInputEvent, casc);

            // ## Xi selection 1 ### keep Xi
            if (fCascadeCutsXi->isSelected(fCascade))
            {
                vXi.push_back(*fCascade);
                vXi[vXi.size() - 1].SetCPA(1.0);
            }
            if (fCascadeCutsAntiXi->isSelected(fCascade))
            {
                vAntiXi.push_back(*fCascade);
                vAntiXi[vAntiXi.size() - 1].SetCPA(1.0);
            }

            fCascade2->SetCascade(fInputEvent, casc);
            // ## Xi selection 2 ### keep Lambda
            if (fCascadeCutsXi2->isSelected(fCascade2))
            {
                vXi2.push_back(*fCascade2);
                vXi2[vXi2.size() - 1].SetCPA(0.5);
            }
            if (fCascadeCutsAntiXi2->isSelected(fCascade2))
            {
                vAntiXi2.push_back(*fCascade2);
                vAntiXi2[vAntiXi2.size() - 1].SetCPA(0.5);
            }
        }
        

        //###########################################
        // Lambda - Lambda recombinations
        //##########################################
        std::vector<AliFemtoDreamBasePart> vLambda_recomb;
        std::vector<AliFemtoDreamBasePart> tmpLambda_recomb(0); // recombination Vector for the loop
        
        // ein lambda mit allen höheren kombinieren (siehe zweite schleife)
        // schleife läuft nur bis zum vorletzten lambda
        for (size_t iterLamb = 0; iterLamb + 1 < vLambda.size(); iterLamb++)       
        {
            if(vLambda.size() == 1) break;        // abbrechen wenn Lambda nur ein Teilchen enthält oder 

            // recombiniere lambda[iterLamb] mit den darauf folgenden Lambdas 
            // - dadurch werden nicht doppelt Lambdas aber im moment noch doppelt Tracks wenn sie sich zwei Lambdas teilen  
            // tausche nur den Impuls der für die invariante Masse benötigt wird
            // 
            // GetMomentum(0) - Lambda
            // GetMomentum(1) - Pion
            // GetMomentum(2) - Proton

            for (size_t iterUpwards = iterLamb + 1; iterUpwards < vLambda.size(); iterUpwards++) 
            {
                tmpLambda_recomb.clear();
                // check for shared tracks
                if(vLambda[iterLamb].GetIDTracks().size() < 2 || vLambda[iterUpwards].GetIDTracks().size() < 2 ) continue;    // failsafe if the Lambda has no 2 tracks

                if ( vLambda[iterLamb].GetIDTracks()[0] == vLambda[iterUpwards].GetIDTracks()[0])
                {
                    tmpLambda_recomb.push_back(vLambda[iterLamb]);
                    tmpLambda_recomb[ 0 ].SetMomentum(2, vLambda[iterUpwards].GetMomentum(2));
                    vLambda_recomb.push_back(tmpLambda_recomb[0]);
                    fEvtCounter->Fill(6);
                    for (size_t iterLamb_recomb = 0; iterLamb_recomb < vLambda_recomb.size(); iterLamb_recomb++)
                    {
                        hInvMassLambda_shared_pion->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1), 
                                                                                vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                        hInvMassLambda_total->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1), 
                                                                          vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                    }
                }
                if(vLambda[iterLamb].GetIDTracks()[1] == vLambda[iterUpwards].GetIDTracks()[1])
                {
                    tmpLambda_recomb.push_back(vLambda[iterLamb]);
                    tmpLambda_recomb[ 0 ].SetMomentum(1, vLambda[iterUpwards].GetMomentum(1));
                    vLambda_recomb.push_back(tmpLambda_recomb[0]);
                    fEvtCounter->Fill(6);
                    for (size_t iterLamb_recomb = 0; iterLamb_recomb < vLambda_recomb.size(); iterLamb_recomb++)
                    {
                        hInvMassLambda_shared_proton->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1), 
                                                                                  vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                        hInvMassLambda_total->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1), 
                                                                          vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                    }
                }
                else 
                {
                    // save recombination lambda twice for each for manipulation of each track
                    tmpLambda_recomb.push_back(vLambda[iterLamb]);                    
                    tmpLambda_recomb.push_back(vLambda[iterLamb]);
                    // take next lambdas (iterUpwards) and manipulate the two lambdas before
                    tmpLambda_recomb[0].SetMomentum(1, vLambda[iterUpwards].GetMomentum(1));
                    tmpLambda_recomb[1].SetMomentum(2, vLambda[iterUpwards].GetMomentum(2));
                    vLambda_recomb.push_back(tmpLambda_recomb[0]);
                    vLambda_recomb.push_back(tmpLambda_recomb[1]);
                    fEvtCounter->Fill(6);
                    fEvtCounter->Fill(6);
                    for (size_t iterLamb_recomb = 0; iterLamb_recomb < vLambda_recomb.size(); iterLamb_recomb++)
                    {
                        hInvMassLambda_shared_proton->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1), vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                        hInvMassLambda_total->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1), vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                    }
                }
            }
        }
        vLambda_recomb.clear();

        //###########################################
        // Lambda - Xi recombinations
        //##########################################
        std::vector<AliFemtoDreamBasePart> vXi_recomb; 
        std::vector<AliFemtoDreamBasePart> tmpXi_recomb(0);     // temporary recombination vector to calculate new invMasses

        for (size_t iterLamb = 0; iterLamb < vLambda.size(); iterLamb++)   // ein lambda mit allen Xi's kombinieren (siehe zweite schleife)
        {
            if( !vLambda.size() || !vXi.size() ) break;        // abbrechen wenn lambda oder Xi leer ist/sind

            // GetIDTracks()
            // [0] - negativeDaughter
            // [1] - positiveDaughter
            if(vLambda[iterLamb].GetIDTracks().size() < 2) 
            {
                continue;    // failsafe if the Lambda has no 2 tracks
            }

            // recombiniere vLambda[iterLamb] mit jeder Tochter der Xi's
            // - nur Impuls manipulation damit invariante Masse ausgerechnet werden kann
            // ## XI
            // GetMomentum(0) - Xi
            // GetMomentum(1) - Pi-Daughter
            // GetMomentum(2) - Proton-Daughter
            // GetMomentum(3) - Pi-Bachelor
            // Hinweis>>Cascade initialisiert AliFemtoBasePart.fP mit 4. d.h. es sollte sich beim Impulsvektor um alle Zerfallsprodukte handeln
            // GetIDTracks()
            // [0] - negativeDaughter
            // [1] - positiveDaughter
            // [2] - Bachelor
            for (size_t iterXi = 0; iterXi < vXi.size(); iterXi++) 
            {
                if(vXi[iterXi].GetMomenta().size() < 4) 
                {
                    // std::cout << " Momentum vXi zu kleine dimension" << std::endl;
                    continue;   // failsafe, falls gespeichertes Xi keine 4 Momenta besitzt
                }
                // reset temporary recombination vectors
                tmpLambda_recomb.clear();
                tmpXi_recomb.clear();
                
                // safe recombination lambda three times for each following lambda 
                // - for all combinations - Xi_1pi-Lambda_prot ; Xi_2pi-Lambda_prot ; Xi_prot-Lambda_pi
                tmpLambda_recomb.push_back(vLambda[iterLamb]);
                tmpLambda_recomb.push_back(vLambda[iterLamb]);
                tmpLambda_recomb.push_back(vLambda[iterLamb]);
                
                if(tmpLambda_recomb.size() >= 3 && vXi[iterXi].GetMomenta().size() >=3)
                {
                    // take Xi's constituents and manipulate the three lambdas before
                    tmpLambda_recomb[0].SetMomentum(1, vXi[iterXi].GetMomentum(0)); // Bachelor Xi-Pion mit Lambda-Proton
                    tmpLambda_recomb[1].SetMomentum(1, vXi[iterXi].GetMomentum(2)); // Daughter Xi-Pion mit Lambda-Proton
                    tmpLambda_recomb[2].SetMomentum(2, vXi[iterXi].GetMomentum(3)); // Daughter Xi-Proton mit Lambda-Pion
                    vLambda_recomb.push_back(tmpLambda_recomb[0]);
                    vLambda_recomb.push_back(tmpLambda_recomb[1]);
                    vLambda_recomb.push_back(tmpLambda_recomb[2]);
                }
                // ## Xi pairing
                if (vXi[iterXi].GetIDTracks()[2] == vLambda[iterLamb].GetIDTracks()[0])     // ## ## Bachelor shared ## ##
                {
                    tmpXi_recomb.push_back(vXi[iterXi]);
                    tmpXi_recomb.push_back(vXi[iterXi]);
                    tmpXi_recomb.push_back(vXi[iterXi]);
                    tmpXi_recomb[0].SetMomentum(1, vLambda[iterLamb].GetMomentum(1));   // set Pi-Daughter
                    tmpXi_recomb[1].SetMomentum(2, vLambda[iterLamb].GetMomentum(2));   // set Proton-Daughter
                    tmpXi_recomb[2].SetMomentum(1, vLambda[iterLamb].GetMomentum(1));   // set full Lambda
                    tmpXi_recomb[2].SetMomentum(2, vLambda[iterLamb].GetMomentum(2));   // set full Lambda
                    for (int i = 1; i < tmpXi_recomb.size(); i++)
                    {
                        vXi_recomb.push_back(tmpXi_recomb[i]);
                    }
                    for (size_t j = 0; j < tmpXi_recomb.size(); j++)
                    {
                        float invMassToStore = CalculateInvMassXi( tmpXi_recomb[j].GetMomentum(3),  // Bach
                                                                   tmpXi_recomb[j].GetMomentum(2),  // posDaught
                                                                   tmpXi_recomb[j].GetMomentum(1) ); // negDaught)
                        hInvMassXi_shared_bach->Fill(invMassToStore);
                        hInvMassXi_total->Fill(invMassToStore);
                    }
                }
                if (vXi[iterXi].GetIDTracks()[0] == vLambda[iterLamb].GetIDTracks()[0])    // ## ## pion daughter shared ## ##
                {
                    if (vXi[iterXi].GetIDTracks()[1] == vLambda[iterLamb].GetIDTracks()[1]) // ## ## and daughter proton shared -> full lambda shared ## ##
                    {
                    tmpXi_recomb.push_back(vXi[iterXi]);
                    tmpXi_recomb[0].SetMomentum(3, vLambda[iterLamb].GetMomentum(1));   // set only Bachelor
                    vXi_recomb.push_back(tmpXi_recomb[0]);

                    hInvMassXi_shared_Lambda->Fill(CalculateInvMassXi( tmpXi_recomb[0].GetMomentum(3),
                                                                       tmpXi_recomb[0].GetMomentum(2),
                                                                       tmpXi_recomb[0].GetMomentum(1) ));
                    }
                    else // ## ## only daughter pion shared ## ##
                    {
                        tmpXi_recomb.push_back(vXi[iterXi]);
                        tmpXi_recomb.push_back(vXi[iterXi]);
                        tmpXi_recomb[0].SetMomentum(3, vLambda[iterLamb].GetMomentum(1)); // set Bachelor
                        tmpXi_recomb[1].SetMomentum(2, vLambda[iterLamb].GetMomentum(2)); // set Proton-Daughter
                        for (int i = 1; i < tmpXi_recomb.size(); i++)
                        {
                            vXi_recomb.push_back(tmpXi_recomb[i]);
                        }
                        for (size_t j = 0; j < tmpXi_recomb.size(); j++)
                        {
                            float invMassToStore = CalculateInvMassXi( tmpXi_recomb[j].GetMomentum(3),  // Bach
                                                                       tmpXi_recomb[j].GetMomentum(2),  // posDaught
                                                                       tmpXi_recomb[j].GetMomentum(1) ); // negDaught)

                            hInvMassXi_shared_pi_daugh->Fill(invMassToStore);
                            hInvMassXi_total->Fill(invMassToStore);
                        }
                    }
                }
                if (vXi[iterXi].GetIDTracks()[1] == vLambda[iterLamb].GetIDTracks()[1] && vXi[iterXi].GetIDTracks()[0] != vLambda[iterLamb].GetIDTracks()[0])     // ## ## only daughter proton shared ## ##
                {
                    tmpXi_recomb.push_back(vXi[iterXi]);
                    tmpXi_recomb.push_back(vXi[iterXi]);
                    tmpXi_recomb[0].SetMomentum(3, vLambda[iterLamb].GetMomentum(1)); // set Bachelor
                    tmpXi_recomb[1].SetMomentum(1, vLambda[iterLamb].GetMomentum(1)); // set Pi-Daughter
                    for (int i = 1; i < tmpXi_recomb.size(); i++)
                    {
                        vXi_recomb.push_back(tmpXi_recomb[i]);
                    }
                    for (size_t j = 0; j < tmpXi_recomb.size(); j++)
                    {
                        hInvMassXi_shared_prot_daugh->Fill(CalculateInvMassXi(tmpXi_recomb[j].GetMomentum(3),   // Bach
                                                                              tmpXi_recomb[j].GetMomentum(2),   // posDaught
                                                                              tmpXi_recomb[j].GetMomentum(1))); // negDaught)
                        hInvMassXi_total->Fill(CalculateInvMassXi(tmpXi_recomb[j].GetMomentum(3),               // Bach
                                                                  tmpXi_recomb[j].GetMomentum(2),               // posDaught
                                                                  tmpXi_recomb[j].GetMomentum(1)));             // negDaught
                    }
                }
                else   // ## ## nothing shared ## ##
                {
                    // get the Xi and manipulate the Bachelor and Daughters
                    for (int j = 0; j < 4; j++)
                    {
                        tmpXi_recomb.push_back(vXi[iterXi]);
                    }

                    tmpXi_recomb[0].SetMomentum(3, vLambda[iterLamb].GetMomentum(1)); // set Bachelor
                    tmpXi_recomb[1].SetMomentum(1, vLambda[iterLamb].GetMomentum(1)); // set Pi-Daughter
                    tmpXi_recomb[2].SetMomentum(2, vLambda[iterLamb].GetMomentum(2)); // set Proton-Daughter
                    tmpXi_recomb[3].SetMomentum(1, vLambda[iterLamb].GetMomentum(1)); // set full Lambda
                    tmpXi_recomb[3].SetMomentum(2, vLambda[iterLamb].GetMomentum(2)); // set full Lambda
                    for (int i = 1; i < tmpXi_recomb.size(); i++)
                    {
                        vXi_recomb.push_back(tmpXi_recomb[i]);
                    }
                    for (size_t j = 0; j < tmpXi_recomb.size(); j++)
                    {
                        hInvMassXi_total->Fill(CalculateInvMassXi( tmpXi_recomb[j].GetMomentum(3),               // Bach
                                                                   tmpXi_recomb[j].GetMomentum(2),               // posDaught
                                                                   tmpXi_recomb[j].GetMomentum(1)) );             // negDaught
                    }
                }
            }
        }
        vXi_recomb.clear();
        vLambda_recomb.clear();
        
        //###########################################
        // Cleanup and Postdata
        //##########################################
        
        // remove double-matched tracks
        fPairCleaner->ResetArray();
        fPairCleaner2->ResetArray();

        // #1
        fPairCleaner->CleanDecayAndDecay(&vXi, &vLambda, 0);
        fPairCleaner->CleanDecayAndDecay(&vAntiXi, &vAntiLambda, 1);
        fPairCleaner->CleanDecay(&vLambda, 2);
        fPairCleaner->CleanDecay(&vAntiLambda, 3);

        fPairCleaner->CleanDecay(&vXi, 0);
        fPairCleaner->CleanDecay(&vAntiXi, 1);

        fPairCleaner->StoreParticle(vLambda);
        fPairCleaner->StoreParticle(vAntiLambda);
        fPairCleaner->StoreParticle(vXi);
        fPairCleaner->StoreParticle(vAntiXi);

        // #2
        fPairCleaner2->CleanDecayAndDecay(&vXi2, &vLambda2, 0);
        fPairCleaner2->CleanDecayAndDecay(&vAntiXi2, &vAntiLambda2, 1);
        fPairCleaner2->CleanDecay(&vLambda2, 2);
        fPairCleaner2->CleanDecay(&vAntiLambda2, 3);

        fPairCleaner2->StoreParticle(vLambda2);
        fPairCleaner2->StoreParticle(vAntiLambda2);
        fPairCleaner2->StoreParticle(vXi2);
        fPairCleaner2->StoreParticle(vAntiXi2);

        fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality()); // proton xi and lambda analysis
        // fPartColl2->SetEvent(fPairCleaner2->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality()); // proton xi and lambda analysis
        // soweit ich das richtig verstanden habe wird pairQA mit den teilchen gemacht die im pairCleaner
        // sind und pdgCodes in der richtigen Reihenfolge vorhanden sind.

        for (size_t j = 0; j < vXi.size(); j++)
        {
            if(!vXi[j].UseParticle()) continue;
            hInvMassXi_sanityCheck->Fill(CalculateInvMassXi( vXi[j].GetMomentum(3),               // Bach
                                                             vXi[j].GetMomentum(2),               // posDaught
                                                             vXi[j].GetMomentum(1)) );            // negDaught
        }
        for (size_t k = 0; k < vLambda.size(); k++)
        {
            if(!vLambda[k].UseParticle()) continue;
            hInvMassLambda_sanityCheck->Fill(CalculateInvMassLambda( vLambda[k].GetMomentum(1),               // posDaught
                                                                     vLambda[k].GetMomentum(2) ));             // negDaught
        }
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

        PostData(15, tlRecombination); // reconstruction from daugthers histograms

        if (fIsMC)
        {
            PostData(16, tlLambdaMC);

            PostData(17, tlAntiLambdaMC);
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
    AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack*>(vTrack);
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
float AliAnalysisTaskPOmegaPenne::CalculateInvMassLambda(TVector3 momPosDaughter, TVector3 momNegDaughter)
{
    float invMass = 0;
    
    float massDP = TDatabasePDG::Instance()->GetParticle(2212)->Mass(); // Proton
    float massDN = TDatabasePDG::Instance()->GetParticle(211)->Mass();  // Pion
    float EDaugP = TMath::Sqrt(
        massDP * massDP + 
        momPosDaughter.X() * momPosDaughter.X() + 
        momPosDaughter.Y() * momPosDaughter.Y() + 
        momPosDaughter.Z() * momPosDaughter.Z()
    );
    float EDaugN = TMath::Sqrt(
        massDN * massDN + 
        momNegDaughter.X() * momNegDaughter.X() + 
        momNegDaughter.Y() * momNegDaughter.Y() + 
        momNegDaughter.Z() * momNegDaughter.Z()
    );
    float energysum = EDaugP + EDaugN;
    float pSum2 = 
        ( momNegDaughter.X() + momPosDaughter.X() ) * 
        ( momNegDaughter.X() + momPosDaughter.X() ) 
        +
        ( momNegDaughter.Y() + momPosDaughter.Y() ) * 
        ( momNegDaughter.Y() + momPosDaughter.Y() ) 
        + 
        ( momNegDaughter.Z() + momPosDaughter.Z() ) * 
        ( momNegDaughter.Z() + momPosDaughter.Z() )
    ;
    invMass = TMath::Sqrt(energysum * energysum - pSum2);
    return invMass;
}

float AliAnalysisTaskPOmegaPenne::CalculateInvMassXi(TVector3 momBach, TVector3 momPosDaughter, TVector3 momNegDaughter)
{
    float massPosDaugh = TDatabasePDG::Instance()->GetParticle(2212)->Mass();  // Proton
    float massNegDaugh = TDatabasePDG::Instance()->GetParticle(211)->Mass();   // Pion
    float massBach = massNegDaugh;
    float massV0 = TDatabasePDG::Instance()->GetParticle(3122)->Mass();  // Lambda
    
    TVector3 PtotV0 = (momPosDaughter + momNegDaughter);
    float Ev0 = ::sqrt(massV0 * massV0 + PtotV0.Mag2());

    float EBach = ::sqrt(massBach + momBach.Mag2());

    float Ptot2Casc = (PtotV0 + momBach).Mag2();

    return ::sqrt(pow(Ev0 + EBach,2) - Ptot2Casc);
}

float AliAnalysisTaskPOmegaPenne::CalculateInvMassHere(AliFemtoDreamv0 *v0, int PDGPosDaug, int PDGNegDaug)     // copied from AliFemtoDreamv0Cuts
{
    Double_t invMass = 0;
    
    float massDP = TDatabasePDG::Instance()->GetParticle(PDGPosDaug)->Mass();
    float massDN = TDatabasePDG::Instance()->GetParticle(PDGNegDaug)->Mass();
    float EDaugP = TMath::Sqrt(
        massDP * massDP + 
        v0->GetPosDaughter()->GetMomentum().X() * v0->GetPosDaughter()->GetMomentum().X() + 
        v0->GetPosDaughter()->GetMomentum().Y() * v0->GetPosDaughter()->GetMomentum().Y() + 
        v0->GetPosDaughter()->GetMomentum().Z() * v0->GetPosDaughter()->GetMomentum().Z()
    );
    float EDaugN = TMath::Sqrt(
        massDN * massDN + 
        v0->GetNegDaughter()->GetMomentum().X() * v0->GetNegDaughter()->GetMomentum().X() + 
        v0->GetNegDaughter()->GetMomentum().Y() * v0->GetNegDaughter()->GetMomentum().Y() + 
        v0->GetNegDaughter()->GetMomentum().Z() * v0->GetNegDaughter()->GetMomentum().Z())
    ;
    float energysum = EDaugP + EDaugN;
    float pSum2 = 
        (v0->GetNegDaughter()->GetMomentum().X() + v0->GetPosDaughter()->GetMomentum().X()) * 
        (v0->GetNegDaughter()->GetMomentum().X() + v0->GetPosDaughter()->GetMomentum().X()) +

        (v0->GetNegDaughter()->GetMomentum().Y() + v0->GetPosDaughter()->GetMomentum().Y()) * 
        (v0->GetNegDaughter()->GetMomentum().Y() + v0->GetPosDaughter()->GetMomentum().Y()) + 

        (v0->GetNegDaughter()->GetMomentum().Z() + v0->GetPosDaughter()->GetMomentum().Z()) * 
        (v0->GetNegDaughter()->GetMomentum().Z() + v0->GetPosDaughter()->GetMomentum().Z())
        ;
    invMass = TMath::Sqrt(energysum * energysum - pSum2);
    return invMass;
}

void AliAnalysisTaskPOmegaPenne::MixChildParticles(std::vector<AliFemtoDreamBasePart> XiVector, std::vector<AliFemtoDreamBasePart> LambdaVector, TList *outputLists, bool checkSameParticleMixing)
{
    std::vector<AliFemtoDreamBasePart> vLambda_recomb(0);
    std::vector<AliFemtoDreamBasePart> vXi_recomb(0);
    std::vector<AliFemtoDreamBasePart> tmpLambda_recomb(0);

    //###########################################
    // Lambda - Lambda recombinations
    //###########################################

    // ein lambda mit allen indexen aufwärts kombinieren (siehe zweite schleife)
    // schleife läuft nur bis zum vorletzten lambda
    for (size_t iterLamb = 0; iterLamb + 1 < LambdaVector.size(); iterLamb++)
    {
        if (LambdaVector.size() == 1)
            break; // abbrechen wenn Lambda nur ein Teilchen enthält

        // // recombiniere lambda[iterLamb] mit den darauf folgenden Lambdas
        // // - dadurch werden nicht doppelt Lambdas aber im moment noch doppelt Tracks wenn sie sich zwei Lambdas teilen
        // // tausche nur den Impuls der für die invariante Masse benötigt wird
        // //
        // // GetMomentum(0) - Lambda
        // // GetMomentum(1) - Pion
        // // GetMomentum(2) - Proton
        for (size_t iterUpwards = iterLamb + 1; iterUpwards < LambdaVector.size(); iterUpwards++)
        {
            tmpLambda_recomb.clear();
            // check for shared tracks
            if (LambdaVector[iterLamb].GetIDTracks().size() < 2 || LambdaVector[iterUpwards].GetIDTracks().size() < 2)
                continue; // failsafe if the Lambda has no 2 tracks

            if (LambdaVector[iterLamb].GetIDTracks()[0] == LambdaVector[iterUpwards].GetIDTracks()[0])
            {
                tmpLambda_recomb.push_back(LambdaVector[iterLamb]);
                tmpLambda_recomb[0].SetMomentum(2, LambdaVector[iterUpwards].GetMomentum(2));
                vLambda_recomb.push_back(tmpLambda_recomb[0]);
                fEvtCounter->Fill(6);
                for (size_t iterLamb_recomb = 0; iterLamb_recomb < vLambda_recomb.size(); iterLamb_recomb++)
                {
                    hInvMassLambda_shared_pion->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1),
                                                                            vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                    hInvMassLambda_total->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1),
                                                                      vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                }
            }
            if (LambdaVector[iterLamb].GetIDTracks()[1] == LambdaVector[iterUpwards].GetIDTracks()[1])
            {
                tmpLambda_recomb.push_back(LambdaVector[iterLamb]);
                tmpLambda_recomb[0].SetMomentum(1, LambdaVector[iterUpwards].GetMomentum(1));
                vLambda_recomb.push_back(tmpLambda_recomb[0]);
                fEvtCounter->Fill(6);
                for (size_t iterLamb_recomb = 0; iterLamb_recomb < vLambda_recomb.size(); iterLamb_recomb++)
                {
                    hInvMassLambda_shared_proton->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1),
                                                                              vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                    hInvMassLambda_total->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1),
                                                                      vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                }
            }
            else
            {
                // save recombination lambda twice for each for manipulation of each track
                tmpLambda_recomb.push_back(LambdaVector[iterLamb]);
                tmpLambda_recomb.push_back(LambdaVector[iterLamb]);
                // take next lambdas (iterUpwards) and manipulate the two lambdas before
                tmpLambda_recomb[0].SetMomentum(1, LambdaVector[iterUpwards].GetMomentum(1));
                tmpLambda_recomb[1].SetMomentum(2, LambdaVector[iterUpwards].GetMomentum(2));
                vLambda_recomb.push_back(tmpLambda_recomb[0]);
                vLambda_recomb.push_back(tmpLambda_recomb[1]);
                fEvtCounter->Fill(6);
                fEvtCounter->Fill(6);
                for (size_t iterLamb_recomb = 0; iterLamb_recomb < vLambda_recomb.size(); iterLamb_recomb++)
                {
                    hInvMassLambda_shared_proton->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1), vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                    hInvMassLambda_total->Fill(CalculateInvMassLambda(vLambda_recomb[iterLamb_recomb].GetMomentum(1), vLambda_recomb[iterLamb_recomb].GetMomentum(2)));
                }
            }
        }
    }
    vLambda_recomb.clear();

    //###########################################
    // Lambda - Xi recombinations
    //##########################################
    std::vector<AliFemtoDreamBasePart> tmpXi_recomb(0); // temporary recombination vector to calculate new invMasses

    for (size_t iterLamb = 0; iterLamb < LambdaVector.size(); iterLamb++) // ein lambda mit allen Xi's kombinieren (siehe zweite schleife)
    {
        if (!LambdaVector.size() || !XiVector.size())
            break; // abbrechen wenn lambda oder Xi leer ist/sind

        // recombiniere LambdaVector[iterLamb] mit jeder Tochter der Xi's
        // - nur Impuls manipulation damit invariante Masse ausgerechnet werden kann
        // ## XI
        // GetMomentum(0) - Xi
        // GetMomentum(1) - Pi-Daughter
        // GetMomentum(2) - Proton-Daughter
        // GetMomentum(3) - Pi-Bachelor
        // Hinweis>>Cascade initialisiert AliFemtoBasePart.fP mit 4. d.h. es sollte sich beim Impulsvektor um alle Zerfallsprodukte handeln
        for (size_t iterXi = 0; iterXi < XiVector.size(); iterXi++)
        {
            // reset temporary recombination vectors
            tmpLambda_recomb.clear();
            tmpXi_recomb.clear();

            // safe recombination lambda three times for each following lambda
            // - for all combinations - Xi_1pi-Lambda_prot ; Xi_2pi-Lambda_prot ; Xi_prot-Lambda_pi
            tmpLambda_recomb.push_back(LambdaVector[iterLamb]);
            tmpLambda_recomb.push_back(LambdaVector[iterLamb]);
            tmpLambda_recomb.push_back(LambdaVector[iterLamb]);

            if (tmpLambda_recomb.size() >= 3 && XiVector[iterXi].GetMomenta().size() >= 3)
            {
                // take Xi's constituents and manipulate the three lambdas before
                tmpLambda_recomb[0].SetMomentum(1, XiVector[iterXi].GetMomentum(0)); // Bachelor Xi-Pion mit Lambda-Proton
                tmpLambda_recomb[1].SetMomentum(1, XiVector[iterXi].GetMomentum(2)); // Daughter Xi-Pion mit Lambda-Proton
                tmpLambda_recomb[2].SetMomentum(2, XiVector[iterXi].GetMomentum(3)); // Daughter Xi-Proton mit Lambda-Pion
                vLambda_recomb.push_back(tmpLambda_recomb[0]);
                vLambda_recomb.push_back(tmpLambda_recomb[1]);
                vLambda_recomb.push_back(tmpLambda_recomb[2]);
            }
            // ## Xi pairing
            if (XiVector[iterXi].GetIDTracks()[0] == LambdaVector[iterLamb].GetIDTracks()[0]) // ## ## Bachelor shared ## ##
            {
                tmpXi_recomb.push_back(XiVector[iterXi]);
                tmpXi_recomb.push_back(XiVector[iterXi]);
                tmpXi_recomb.push_back(XiVector[iterXi]);
                tmpXi_recomb[0].SetMomentum(1, LambdaVector[iterLamb].GetMomentum(1)); // set Pi-Daughter
                tmpXi_recomb[1].SetMomentum(2, LambdaVector[iterLamb].GetMomentum(2)); // set Proton-Daughter
                tmpXi_recomb[2].SetMomentum(1, LambdaVector[iterLamb].GetMomentum(1)); // set full Lambda
                tmpXi_recomb[2].SetMomentum(2, LambdaVector[iterLamb].GetMomentum(2)); // set full Lambda
                for (int i = 1; i < tmpXi_recomb.size(); i++)
                {
                    vXi_recomb.push_back(tmpXi_recomb[i]);
                }
                for (size_t j = 0; j < tmpXi_recomb.size(); j++)
                {
                    float invMassToStore = CalculateInvMassXi(tmpXi_recomb[j].GetMomentum(3),  // Bach
                                                              tmpXi_recomb[j].GetMomentum(2),  // posDaught
                                                              tmpXi_recomb[j].GetMomentum(1)); // negDaught)
                    hInvMassXi_shared_bach->Fill(invMassToStore);
                    hInvMassXi_total->Fill(invMassToStore);
                }
            }
            if (XiVector[iterXi].GetIDTracks()[1] == LambdaVector[iterLamb].GetIDTracks()[1]) // ## ## pion daughter shared ## ##
            {
                if (XiVector[iterXi].GetIDTracks()[2] == LambdaVector[iterLamb].GetIDTracks()[2]) // ## ## and daughter proton shared -> full lambda shared ## ##
                {
                    tmpXi_recomb.push_back(XiVector[iterXi]);
                    tmpXi_recomb[0].SetMomentum(3, LambdaVector[iterLamb].GetMomentum(1)); // set only Bachelor
                    vXi_recomb.push_back(tmpXi_recomb[0]);

                    hInvMassXi_shared_Lambda->Fill(CalculateInvMassXi(tmpXi_recomb[0].GetMomentum(3),
                                                                      tmpXi_recomb[0].GetMomentum(2),
                                                                      tmpXi_recomb[0].GetMomentum(1)));
                }
                else // ## ## only daughter pion shared ## ##
                {
                    tmpXi_recomb.push_back(XiVector[iterXi]);
                    tmpXi_recomb.push_back(XiVector[iterXi]);
                    tmpXi_recomb[0].SetMomentum(3, LambdaVector[iterLamb].GetMomentum(1)); // set Bachelor
                    tmpXi_recomb[1].SetMomentum(2, LambdaVector[iterLamb].GetMomentum(2)); // set Proton-Daughter
                    for (int i = 1; i < tmpXi_recomb.size(); i++)
                    {
                        vXi_recomb.push_back(tmpXi_recomb[i]);
                    }
                    for (size_t j = 0; j < tmpXi_recomb.size(); j++)
                    {
                        float invMassToStore = CalculateInvMassXi(tmpXi_recomb[j].GetMomentum(3),  // Bach
                                                                  tmpXi_recomb[j].GetMomentum(2),  // posDaught
                                                                  tmpXi_recomb[j].GetMomentum(1)); // negDaught)

                        hInvMassXi_shared_pi_daugh->Fill(invMassToStore);
                        hInvMassXi_total->Fill(invMassToStore);
                    }
                }
            }
            if (XiVector[iterXi].GetIDTracks()[2] == LambdaVector[iterLamb].GetIDTracks()[2] && XiVector[iterXi].GetIDTracks()[1] != LambdaVector[iterLamb].GetIDTracks()[1]) // ## ## only daughter proton shared ## ##
            {
                tmpXi_recomb.push_back(XiVector[iterXi]);
                tmpXi_recomb.push_back(XiVector[iterXi]);
                tmpXi_recomb[0].SetMomentum(3, LambdaVector[iterLamb].GetMomentum(1)); // set Bachelor
                tmpXi_recomb[1].SetMomentum(1, LambdaVector[iterLamb].GetMomentum(1)); // set Pi-Daughter
                for (int i = 1; i < tmpXi_recomb.size(); i++)
                {
                    vXi_recomb.push_back(tmpXi_recomb[i]);
                }
                for (size_t j = 0; j < tmpXi_recomb.size(); j++)
                {
                    hInvMassXi_shared_prot_daugh->Fill(CalculateInvMassXi(tmpXi_recomb[j].GetMomentum(3),   // Bach
                                                                          tmpXi_recomb[j].GetMomentum(2),   // posDaught
                                                                          tmpXi_recomb[j].GetMomentum(1))); // negDaught)
                    hInvMassXi_total->Fill(CalculateInvMassXi(tmpXi_recomb[j].GetMomentum(3),               // Bach
                                                              tmpXi_recomb[j].GetMomentum(2),               // posDaught
                                                              tmpXi_recomb[j].GetMomentum(1)));             // negDaught
                }
            }
            else // ## ## nothing shared ## ##
            {
                // get the Xi and manipulate the Bachelor and Daughters
                for (int j = 0; j < 4; j++)
                {
                    tmpXi_recomb.push_back(XiVector[iterXi]);
                }
                tmpXi_recomb[0].SetMomentum(3, LambdaVector[iterLamb].GetMomentum(1)); // set Bachelor
                tmpXi_recomb[1].SetMomentum(1, LambdaVector[iterLamb].GetMomentum(1)); // set Pi-Daughter
                tmpXi_recomb[2].SetMomentum(2, LambdaVector[iterLamb].GetMomentum(2)); // set Proton-Daughter
                tmpXi_recomb[3].SetMomentum(1, LambdaVector[iterLamb].GetMomentum(1)); // set full Lambda
                tmpXi_recomb[3].SetMomentum(2, LambdaVector[iterLamb].GetMomentum(2)); // set full Lambda
                for (int i = 1; i < tmpXi_recomb.size(); i++)
                {
                    vXi_recomb.push_back(tmpXi_recomb[i]);
                }
                for (size_t j = 0; j < tmpXi_recomb.size(); j++)
                {
                    hInvMassXi_total->Fill(CalculateInvMassXi(tmpXi_recomb[j].GetMomentum(3),   // Bach
                                                              tmpXi_recomb[j].GetMomentum(2),   // posDaught
                                                              tmpXi_recomb[j].GetMomentum(1))); // negDaught
                }
            }
        }
    }
    vXi_recomb.clear();
    vLambda_recomb.clear();
}